#pragma once

#include <TBranch.h>
#include <TChain.h>
#include <concepts>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <type_traits>

template <typename... branch_type> class root_chain {
private:
    template <class T> struct to_array { typedef T type; };
    // template <class T, std::size_t N> struct to_array<T[N]> { typedef std::array<typename to_array<T>::type, N> type; };
    template <class T> using to_array_t = typename to_array<T>::type;
    std::unique_ptr<TChain> chain{};
    std::tuple<to_array_t<branch_type> *...> data;
    std::array<TBranch *, sizeof...(branch_type)> b_add;

    template <class T, T I, T... J> inline void do_delete(std::integer_sequence<T, I, J...>) {
        do_delete(std::integer_sequence<std::size_t, I>{});
        do_delete(std::integer_sequence<std::size_t, J...>{});
    }

    template <class T, T I> inline void do_delete(std::integer_sequence<T, I>) {
        typedef typename std::tuple_element<I, std::tuple<to_array_t<branch_type>...>>::type thistype;
        if (std::is_trivially_copyable<thistype>::value)
            std::default_delete<thistype>{}(std::get<I>(data));
        // simliarI to std::is_array<thistype>::value?delete [] std::get<>(data):delete std::get<I>(data);
    }

    template <class T, T I, T... J> inline void do_setaddress(std::integer_sequence<T, I, J...>, const std::array<const char *, sizeof...(branch_type)> names) {
        do_setaddress(std::integer_sequence<std::size_t, I>{}, names);
        do_setaddress(std::integer_sequence<std::size_t, J...>{}, names);
    }

    template <class T, T I> inline void do_setaddress(std::integer_sequence<T, I>, const std::array<const char *, sizeof...(branch_type)> names) {
        b_add[I] = chain->GetBranch(std::get<I>(names));
        assert(b_add[I]);
        typedef typename std::tuple_element<I, std::tuple<to_array_t<branch_type>...>>::type thistype;
        if (std::is_trivially_copyable<thistype>::value) {
            // normal case, assign space and SetBranchAddress, the ugly expression
            // is for properly dealing with arrays, the assigned space is then freed
            // by me (not ROOT)
            std::get<I>(data) = reinterpret_cast<thistype *>(new thistype);
            chain->SetBranchAddress(std::get<I>(names),
                                    reinterpret_cast<decltype(new typename std::tuple_element<I, std::tuple<branch_type...>>::type)>(std::get<I>(data)), &b_add[I]);
            // following method is broken, don't know why:
            // b_add[I]->SetAddress(std::get<I>(data));
        } else {
            // ok, seems for non-trivially copyable containers, ROOT will only accept
            // T** pointer, and assign space on its own, free on its own
            std::get<I>(data) = nullptr;
            chain->SetBranchAddress(std::get<I>(names), &std::get<I>(data), &b_add[I]);
        }
    }

    template <typename T, T... IDs> std::tuple<to_array_t<branch_type> &...> get_elements(std::integer_sequence<T, IDs...>, std::size_t id) {
        auto curr = chain->LoadTree(id);
        if (curr < 0) {
            throw;
        }
        for (auto &i : b_add)
            i->GetEntry(curr);
        // chain->GetEntry(curr);
        return {*std::get<IDs>(data)...};
        // return std::forward_as_tuple(*std::get<IDs>(data)...);
    }

public:
    template <typename T>
    // requires std::is_nothrow_convertible_v<typename std::remove_cvref_t<T>::value_type, std::string>
    root_chain(T &&file_list, const char *tree_name, const std::array<const char *, sizeof...(branch_type)> names) {
        chain = std::make_unique<TChain>(tree_name);
        for (const auto &i : file_list) {
            chain->Add(i.c_str());
        }
        do_setaddress(std::make_integer_sequence<std::size_t, sizeof...(branch_type)>{}, names);
    }

    root_chain(const root_chain &) = delete;
    root_chain(root_chain &&) = default;

    ~root_chain() {
        chain->ResetBranchAddresses();
        do_delete(std::make_integer_sequence<std::size_t, sizeof...(branch_type)>{});
        // cleaning of chain is done by unique_ptr
        // the chain will then clean up TBranches
    }

    decltype(auto) get_elements(std::size_t id) { return get_elements(std::make_integer_sequence<std::size_t, sizeof...(branch_type)>{}, id); }

    decltype(auto) get_entries() { return chain->GetEntries(); }

    // used by range-based-for
    class iter {
    private:
        std::size_t num;
        root_chain &se;

    public:
        iter(int m_num, root_chain &root_chain_self) : num(m_num), se(root_chain_self) {}
        iter &operator++() {
            num++;
            return *this;
        }
        bool operator!=(const iter &other) const { return num != other.num; }
        decltype(auto) operator*() { return se.get_elements(num); }
    };

    decltype(auto) operator[](std::size_t id) { return get_elements(id); }

    decltype(auto) begin() { return iter(0, *this); }
    decltype(auto) end() { return iter(get_entries(), *this); }
};

template <typename... branch_type> class chain_runner {
private:
    std::vector<root_chain<branch_type...>> chains{};
    size_t nthread{};

public:
    template <typename T>
    // requires std::is_nothrow_convertible_v<typename std::remove_cvref_t<T>::value_type, std::string>
    chain_runner(T &&file_list, const char *tree_name, const std::array<const char *, sizeof...(branch_type)> names, size_t nthread_ = 16) {
        nthread = std::min(nthread_, file_list.size());
        auto split_list = std::make_unique<std::vector<std::string>[]>(nthread);
        chains.reserve(nthread);
        for (size_t i = 0; i < file_list.size(); i++) {
            split_list[i % nthread].push_back(file_list[i]);
        }
        for (size_t i = 0; i < nthread; i++) {
            chains.emplace_back(split_list[i], tree_name, names);
        }
    }
    size_t get_entries() {
        size_t entries = 0;
        for (auto &i : chains) {
            entries += i.get_entries();
        }
        return entries;
    }

    template <typename RunManager, bool finalize_unpar = true, typename... Args> [[nodiscard]] RunManager run(Args &&...args) {
        RunManager run(std::forward<Args>(args)...);
        std::vector<std::thread> execute_threads;
        execute_threads.reserve(nthread);
        std::mutex finalize_lock;
        for (auto &chain : chains) {
            execute_threads.emplace_back([&] {
                auto thread_obj = run.get_thread_object();
                for (const auto &s : chain) {
                    std::apply([&thread_obj](auto &&...m_args) { thread_obj.run(std::forward<decltype(m_args)>(m_args)...); }, s);
                }
                if constexpr (finalize_unpar) {
                    std::lock_guard<std::mutex> lock(finalize_lock);
                    thread_obj.finalize();
                } else {
                    thread_obj.finalize();
                }
            });
        }
        for (auto &i : execute_threads) {
            i.join();
        }
        run.finalize();
        return run;
    }
};
