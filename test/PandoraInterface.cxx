/**
 *  @file   LArRecoMP/test/PandoraInterface.cc
 *
 *  @brief  Implementation of the lar reco mp application
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Helpers/XmlHelper.h"
#include "Managers/PluginManager.h"
#include "TLorentzVector.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArPersistency/EventReadingAlgorithm.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"
#include "Pandora/ObjectCreation.h"
#include "Objects/MCParticle.h"
#include <cstddef>

#ifdef LIBTORCH_DL
#include "larpandoradlcontent/LArDLContent.h"
#endif

#include "PandoraInterface.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <getopt.h>
#include <iostream>
#include <string>

#include "chain_helper.h"
#include "TGraph2D.h"

using namespace pandora;
using namespace lar_reco;

int main(int argc, char *argv[])
{
    int errorNo(0);
    const Pandora *pPrimaryPandora(nullptr);

    try
    {
        Parameters parameters;

        if (!ParseCommandLine(argc, argv, parameters))
            return 1;

#ifdef MONITORING
        TApplication *pTApplication = new TApplication("LArReco", &argc, argv);
        pTApplication->SetReturnFromRun(kTRUE);
#endif
        CreatePandoraInstances(parameters, pPrimaryPandora);

        if (!pPrimaryPandora)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ProcessEvents(parameters, pPrimaryPandora);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
        errorNo = 1;
    }
    catch (const StopProcessingException &)
    {
        // Exit gracefully
        errorNo = 0;
    }
    catch (...)
    {
        std::cerr << "Unknown exception: " << std::endl;
        errorNo = 1;
    }

    MultiPandoraApi::DeletePandoraInstances(pPrimaryPandora);
    return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_reco
{

void CreatePandoraInstances(const Parameters &parameters, const Pandora *&pPrimaryPandora)
{
    pPrimaryPandora = new Pandora();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPrimaryPandora));
#ifdef LIBTORCH_DL
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*pPrimaryPandora));
#endif
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPrimaryPandora));

    if (!pPrimaryPandora)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    MultiPandoraApi::AddPrimaryPandoraInstance(pPrimaryPandora);

    ProcessExternalParameters(parameters, pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraApi::SetLArTransformationPlugin(*pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPrimaryPandora, parameters.m_settingsFile));
}

//------------------------------------------------------------------------------------------------------------------------------------------

class JUNOEvent {
  public:
    static JUNOEvent &instance() {
        static JUNOEvent m_instance{
            "/var/home/yan/code/rdf_framework_juno/build/process.root"};
        return m_instance;
    }
    decltype(auto) operator[](std::size_t id) {
        std::cerr << "id: " << id << std::endl;
        return junoEDMfile.get_elements(id - 1);
    }
    size_t size() { return junoEDMfile.get_entries(); }

  private:
    JUNOEvent(std::string path)
        : junoEDMfile(
              std::vector<std::string>{path}, "treeout",
              {"plotnpe", "plothittime", "initpoint", "exitpoint", "initP"}) {}
    root_chain<TGraph2D, TGraph2D, TLorentzVector, TLorentzVector,
               TLorentzVector>
        junoEDMfile;
};

// class junoTrackP : public object_creation::MCParticle::Parameters{
// public:
//     int PDG{};
//     double start_theta{};
//     double start_phi{};
//     double start_t{};
//     double end_theta{};
//     double end_phi{};
//     double end_t{};
// };

class JUNOparticle : public pandora::MCParticle{
    public:
    JUNOparticle (const object_creation::MCParticle::Parameters &parameters) : pandora::MCParticle(parameters) {}
};

class junoTrackF
    : public pandora::ObjectFactory<object_creation::MCParticle::Parameters,
                                    object_creation::MCParticle::Object> {
  public:
    virtual Parameters *NewParameters() const { return new Parameters; }
    virtual StatusCode Read(Parameters &, FileReader &) const {
        return pandora::STATUS_CODE_SUCCESS;
    }
    virtual StatusCode Write(const Object *const, FileWriter &) const {
        return pandora::STATUS_CODE_SUCCESS;
    }

  private:
    virtual StatusCode Create(const Parameters &parameters,
                              const Object *&pObject) const {
        pObject = new JUNOparticle(parameters);
        return pandora::STATUS_CODE_SUCCESS;
    }
};

void InitializeTrack(const TLorentzVector &start, const TLorentzVector &end,
                     const TLorentzVector &p,
                     const Pandora *const pPrimaryPandora) {
    object_creation::MCParticle::Parameters parameters;
    parameters.m_energy = p.E();
    parameters.m_momentum.Set({(float)p.Px(), (float)p.Py(), (float)p.Pz()});
    parameters.m_vertex.Set(
        {(float)start.T(), (float)start.Theta(), (float)start.Phi()});
    parameters.m_endpoint.Set(
        {(float)end.T(), (float)end.Theta(), (float)end.Phi()});
    parameters.m_particleId = 13;
    parameters.m_mcParticleType = pandora::MC_VIEW_W;
    parameters.m_pParentAddress = nullptr;
    try {
        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=,
            PandoraApi::MCParticle::Create(*pPrimaryPandora, parameters,
                                           junoTrackF{}));
    } catch (const pandora::StatusCodeException &) {
        std::cerr
            << "CreatePandoraMCParticles - unable to create mc "
               "neutrino, insufficient or invalid information supplied "
            << std::endl;
    }
}

class JUNOHit : public pandora::CaloHit {
  public:
    JUNOHit(const object_creation::CaloHit::Parameters &parameters)
        : pandora::CaloHit(parameters) {}
};

class JUNOCarloHitF
    : public pandora::ObjectFactory<object_creation::CaloHit::Parameters,
                                    object_creation::CaloHit::Object> {
  public:
    virtual Parameters *NewParameters() const { return new Parameters; }
    virtual StatusCode Read(Parameters &, FileReader &) const {
        return pandora::STATUS_CODE_SUCCESS;
    }
    virtual StatusCode Write(const Object *const, FileWriter &) const {
        return pandora::STATUS_CODE_SUCCESS;
    }

  protected:
    virtual StatusCode Create(const Parameters &parameters,
                              const Object *&pObject) const {
        pObject = new JUNOHit(parameters);
        return pandora::STATUS_CODE_SUCCESS;
    }
};

void InitializeCarlo(double t, double charge, double theta, double phi, const Pandora *const pPandora){
    object_creation::CaloHit::Parameters caloHitParameters;
    caloHitParameters.m_cellSize0 = 1;
    caloHitParameters.m_cellSize1 = 1;
    caloHitParameters.m_cellThickness = 0.5;
    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
    caloHitParameters.m_time = 0.;
    caloHitParameters.m_nCellRadiationLengths = 1;
    caloHitParameters.m_nCellInteractionLengths = 1;
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
    caloHitParameters.m_layer = 0;
    caloHitParameters.m_isInOuterSamplingLayer = false;
    caloHitParameters.m_inputEnergy = charge;
    caloHitParameters.m_mipEquivalentEnergy = 1;
    caloHitParameters.m_electromagneticEnergy = 1;
    caloHitParameters.m_hadronicEnergy = 1;
    caloHitParameters.m_hitType = pandora::TPC_VIEW_W;
    const double wpos_cm(
            pPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoW(theta, phi));
    caloHitParameters.m_positionVector.Set(pandora::CartesianVector(t, 0, wpos_cm));
    caloHitParameters.m_expectedDirection.Set(pandora::CartesianVector(0, 0, 1));
    caloHitParameters.m_cellNormalVector.Set(pandora::CartesianVector(0, 0, 1));
    caloHitParameters.m_pParentAddress = nullptr;
    
    // caloHitParameters.m_positionVector = pandora::CartesianVector(t, 0, theta - phi);
    try {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                PandoraApi::CaloHit::Create(*pPandora,
                                                            caloHitParameters,
                                                            JUNOCarloHitF{}));
    } catch (const pandora::StatusCodeException &) {
        return;
    }
}

void ProcessEvents(const Parameters &parameters,
                   const Pandora *const pPrimaryPandora) {
    int nEvents(0);

    while ((nEvents++ < parameters.m_nEventsToProcess) ||
           (0 > parameters.m_nEventsToProcess)) {
        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl
                      << "   PROCESSING EVENT: " << (nEvents - 1) << std::endl
                      << std::endl;
        auto &&[plotnpe, plothittime, initpoint, exitpoint, initP] =
            JUNOEvent::instance()[nEvents];
        InitializeTrack(initpoint, exitpoint, initP, pPrimaryPandora);
        auto nhit = plotnpe.GetN();
        for (int i = 0; i < nhit; ++i) {
            InitializeCarlo(plothittime.GetZ()[i], plotnpe.GetZ()[i],
                            plotnpe.GetX()[i], plotnpe.GetY()[i],
                            pPrimaryPandora);
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                PandoraApi::ProcessEvent(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                PandoraApi::Reset(*pPrimaryPandora));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int c(0);
    std::string recoOption;

    while ((c = getopt(argc, argv, "r:i:e:g:n:s:pNh")) != -1)
    {
        switch (c)
        {
            case 'r':
                recoOption = optarg;
                break;
            case 'i':
                parameters.m_settingsFile = optarg;
                break;
            case 'e':
                parameters.m_eventFileNameList = optarg;
                break;
            case 'g':
                parameters.m_geometryFileName = optarg;
                break;
            case 'n':
                parameters.m_nEventsToProcess = atoi(optarg);
                break;
            case 's':
                parameters.m_nEventsToSkip = atoi(optarg);
                break;
            case 'p':
                parameters.m_printOverallRecoStatus = true;
                break;
            case 'N':
                parameters.m_shouldDisplayEventNumber = true;
                break;
            case 'h':
            default:
                return PrintOptions();
        }
    }
    parameters.m_nEventsToProcess = JUNOEvent::instance().size();
    return ProcessRecoOption(recoOption, parameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl
              << "./bin/PandoraInterface " << std::endl
              << "    -r RecoOption          (required) [Full, AllHitsCR, AllHitsNu, CRRemHitsSliceCR, CRRemHitsSliceNu, AllHitsSliceCR, AllHitsSliceNu]"
              << std::endl
              << "    -i Settings            (required) [algorithm description: xml]" << std::endl
              << "    -e EventFileList       (optional) [colon-separated list of files: xml/pndr]" << std::endl
              << "    -g GeometryFile        (optional) [detector geometry description: xml/pndr]" << std::endl
              << "    -n NEventsToProcess    (optional) [no. of events to process]" << std::endl
              << "    -s NEventsToSkip       (optional) [no. of events to skip in first file]" << std::endl
              << "    -p                     (optional) [print status]" << std::endl
              << "    -N                     (optional) [print event numbers]" << std::endl
              << std::endl;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProcessRecoOption(const std::string &recoOption, Parameters &parameters)
{
    std::string chosenRecoOption(recoOption);
    std::transform(chosenRecoOption.begin(), chosenRecoOption.end(), chosenRecoOption.begin(), ::tolower);

    if ("full" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = true;
    }
    else if ("allhitscr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("nostitchingcr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsnu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else
    {
        std::cout << "LArReco, Unrecognized reconstruction option: " << recoOption << std::endl << std::endl;
        return PrintOptions();
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessExternalParameters(const Parameters &parameters, const Pandora *const pPandora)
{
    auto *const pEventReadingParameters = new lar_content::EventReadingAlgorithm::ExternalEventReadingParameters;
    pEventReadingParameters->m_geometryFileName = parameters.m_geometryFileName;
    pEventReadingParameters->m_eventFileNameList = parameters.m_eventFileNameList;
    if (parameters.m_nEventsToSkip.IsInitialized())
        pEventReadingParameters->m_skipToEvent = parameters.m_nEventsToSkip.Get();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArEventReading", pEventReadingParameters));

    auto *const pEventSteeringParameters = new lar_content::MasterAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = parameters.m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunStitching = parameters.m_shouldRunStitching;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = parameters.m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = parameters.m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = parameters.m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = parameters.m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldPerformSliceId = parameters.m_shouldPerformSliceId;
    pEventSteeringParameters->m_printOverallRecoStatus = parameters.m_printOverallRecoStatus;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));

#ifdef LIBTORCH_DL
    auto *const pEventSettingsParametersCopy = new lar_content::MasterAlgorithm::ExternalSteeringParameters(*pEventSteeringParameters);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
        pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "LArDLMaster", pEventSettingsParametersCopy));
#endif
}

} // namespace lar_reco
