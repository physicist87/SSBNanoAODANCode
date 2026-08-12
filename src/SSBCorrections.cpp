#include "../interface/SSBCorrections.h"
#include "../TextReader/TextReader.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "correction.h"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <vector>
#include "TLorentzVector.h"

namespace {
// Runs a correctionlib load; on failure, logs a Warning and returns nullptr
// instead of propagating. Only for non-essential corrections - required
// ones (jec_, jer_, lepton SFs) should crash loudly if missing, not fall back.
template <typename T, typename LoaderFn>
std::shared_ptr<const T> LoadOptionalCorrection(LoaderFn&& loader, const std::string& context, Logger& logger) {
    try {
        return loader();
    } catch (const std::exception& e) {
        logger.Warning() << context << ": " << e.what() << std::endl;
        return nullptr;
    }
}
}  // namespace

BTagAlgo ParseBTagAlgo(const std::string& jetBtagConfig) {
    if (jetBtagConfig.find("deepCSV") != std::string::npos) return BTagAlgo::DeepCSV;
    if (jetBtagConfig.find("deepJet") != std::string::npos) return BTagAlgo::DeepJet;
    if (jetBtagConfig.find("UParT") != std::string::npos)   return BTagAlgo::UParTAK4;
    if (jetBtagConfig.find("pfCSVV2") != std::string::npos) return BTagAlgo::CSVv2;
    return BTagAlgo::Unknown;
}

std::string BTagAlgoToString(BTagAlgo algo) {
    switch (algo) {
        case BTagAlgo::DeepCSV:  return "DeepCSV";
        case BTagAlgo::DeepJet:  return "DeepJet";
        // Must be "UParTAK4", not "UParT" - matches btagEff_UParTAK4.root's
        // filename and its "eff_UParTAK4_<flav>_<wp>" histogram names.
        case BTagAlgo::UParTAK4: return "UParTAK4";
        case BTagAlgo::CSVv2:    return "CSVv2";
        default:                 return "Unknown";
    }
}

namespace {
// JES full-uncertainty-set NP name -> correctionlib key, from the official
// JERC tutorial's JecConfigAK4.json (not derived by substitution - per-source
// naming doesn't follow the L1FastJet/L2Relative pattern). Full set only.
const std::map<std::string, std::string> kJesFullSetUnc_2016Pre = {
    {"CMS_scale_j_AbsoluteMPFBias", "Summer20UL16APVNanoV15_V1_MC_AbsoluteMPFBias_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteScale", "Summer20UL16APVNanoV15_V1_MC_AbsoluteScale_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteStat_2016", "Summer20UL16APVNanoV15_V1_MC_AbsoluteStat_AK4PFPuppi"},
    {"CMS_scale_j_FlavorQCD", "Summer20UL16APVNanoV15_V1_MC_FlavorQCD_AK4PFPuppi"},
    {"CMS_scale_j_Fragmentation", "Summer20UL16APVNanoV15_V1_MC_Fragmentation_AK4PFPuppi"},
    {"CMS_scale_j_PileUpDataMC", "Summer20UL16APVNanoV15_V1_MC_PileUpDataMC_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtBB", "Summer20UL16APVNanoV15_V1_MC_PileUpPtBB_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC1", "Summer20UL16APVNanoV15_V1_MC_PileUpPtEC1_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC2", "Summer20UL16APVNanoV15_V1_MC_PileUpPtEC2_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtHF", "Summer20UL16APVNanoV15_V1_MC_PileUpPtHF_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtRef", "Summer20UL16APVNanoV15_V1_MC_PileUpPtRef_AK4PFPuppi"},
    {"CMS_scale_j_RelativeFSR", "Summer20UL16APVNanoV15_V1_MC_RelativeFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC1_2016", "Summer20UL16APVNanoV15_V1_MC_RelativeJEREC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC2_2016", "Summer20UL16APVNanoV15_V1_MC_RelativeJEREC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJERHF", "Summer20UL16APVNanoV15_V1_MC_RelativeJERHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtBB", "Summer20UL16APVNanoV15_V1_MC_RelativePtBB_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC1_2016", "Summer20UL16APVNanoV15_V1_MC_RelativePtEC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC2_2016", "Summer20UL16APVNanoV15_V1_MC_RelativePtEC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtHF", "Summer20UL16APVNanoV15_V1_MC_RelativePtHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativeBal", "Summer20UL16APVNanoV15_V1_MC_RelativeBal_AK4PFPuppi"},
    {"CMS_scale_j_RelativeSample_2016", "Summer20UL16APVNanoV15_V1_MC_RelativeSample_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatEC_2016", "Summer20UL16APVNanoV15_V1_MC_RelativeStatEC_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatFSR_2016", "Summer20UL16APVNanoV15_V1_MC_RelativeStatFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatHF_2016", "Summer20UL16APVNanoV15_V1_MC_RelativeStatHF_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionECAL", "Summer20UL16APVNanoV15_V1_MC_SinglePionECAL_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionHCAL", "Summer20UL16APVNanoV15_V1_MC_SinglePionHCAL_AK4PFPuppi"},
    {"CMS_scale_j_TimePtEta_2016", "Summer20UL16APVNanoV15_V1_MC_TimePtEta_AK4PFPuppi"},
};

const std::map<std::string, std::string> kJesFullSetUnc_2016Post = {
    {"CMS_scale_j_AbsoluteMPFBias", "Summer20UL16NanoV15_V1_MC_AbsoluteMPFBias_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteScale", "Summer20UL16NanoV15_V1_MC_AbsoluteScale_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteStat_2016", "Summer20UL16NanoV15_V1_MC_AbsoluteStat_AK4PFPuppi"},
    {"CMS_scale_j_FlavorQCD", "Summer20UL16NanoV15_V1_MC_FlavorQCD_AK4PFPuppi"},
    {"CMS_scale_j_Fragmentation", "Summer20UL16NanoV15_V1_MC_Fragmentation_AK4PFPuppi"},
    {"CMS_scale_j_PileUpDataMC", "Summer20UL16NanoV15_V1_MC_PileUpDataMC_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtBB", "Summer20UL16NanoV15_V1_MC_PileUpPtBB_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC1", "Summer20UL16NanoV15_V1_MC_PileUpPtEC1_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC2", "Summer20UL16NanoV15_V1_MC_PileUpPtEC2_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtHF", "Summer20UL16NanoV15_V1_MC_PileUpPtHF_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtRef", "Summer20UL16NanoV15_V1_MC_PileUpPtRef_AK4PFPuppi"},
    {"CMS_scale_j_RelativeFSR", "Summer20UL16NanoV15_V1_MC_RelativeFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC1_2016", "Summer20UL16NanoV15_V1_MC_RelativeJEREC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC2_2016", "Summer20UL16NanoV15_V1_MC_RelativeJEREC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJERHF", "Summer20UL16NanoV15_V1_MC_RelativeJERHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtBB", "Summer20UL16NanoV15_V1_MC_RelativePtBB_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC1_2016", "Summer20UL16NanoV15_V1_MC_RelativePtEC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC2_2016", "Summer20UL16NanoV15_V1_MC_RelativePtEC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtHF", "Summer20UL16NanoV15_V1_MC_RelativePtHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativeBal", "Summer20UL16NanoV15_V1_MC_RelativeBal_AK4PFPuppi"},
    {"CMS_scale_j_RelativeSample_2016", "Summer20UL16NanoV15_V1_MC_RelativeSample_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatEC_2016", "Summer20UL16NanoV15_V1_MC_RelativeStatEC_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatFSR_2016", "Summer20UL16NanoV15_V1_MC_RelativeStatFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatHF_2016", "Summer20UL16NanoV15_V1_MC_RelativeStatHF_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionECAL", "Summer20UL16NanoV15_V1_MC_SinglePionECAL_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionHCAL", "Summer20UL16NanoV15_V1_MC_SinglePionHCAL_AK4PFPuppi"},
    {"CMS_scale_j_TimePtEta_2016", "Summer20UL16NanoV15_V1_MC_TimePtEta_AK4PFPuppi"},
};

const std::map<std::string, std::string> kJesFullSetUnc_2017 = {
    {"CMS_scale_j_AbsoluteMPFBias", "Summer20UL17NanoV15_V1_MC_AbsoluteMPFBias_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteScale", "Summer20UL17NanoV15_V1_MC_AbsoluteScale_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteStat_2017", "Summer20UL17NanoV15_V1_MC_AbsoluteStat_AK4PFPuppi"},
    {"CMS_scale_j_FlavorQCD", "Summer20UL17NanoV15_V1_MC_FlavorQCD_AK4PFPuppi"},
    {"CMS_scale_j_Fragmentation", "Summer20UL17NanoV15_V1_MC_Fragmentation_AK4PFPuppi"},
    {"CMS_scale_j_PileUpDataMC", "Summer20UL17NanoV15_V1_MC_PileUpDataMC_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtBB", "Summer20UL17NanoV15_V1_MC_PileUpPtBB_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC1", "Summer20UL17NanoV15_V1_MC_PileUpPtEC1_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC2", "Summer20UL17NanoV15_V1_MC_PileUpPtEC2_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtHF", "Summer20UL17NanoV15_V1_MC_PileUpPtHF_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtRef", "Summer20UL17NanoV15_V1_MC_PileUpPtRef_AK4PFPuppi"},
    {"CMS_scale_j_RelativeFSR", "Summer20UL17NanoV15_V1_MC_RelativeFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC1_2017", "Summer20UL17NanoV15_V1_MC_RelativeJEREC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC2_2017", "Summer20UL17NanoV15_V1_MC_RelativeJEREC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJERHF", "Summer20UL17NanoV15_V1_MC_RelativeJERHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtBB", "Summer20UL17NanoV15_V1_MC_RelativePtBB_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC1_2017", "Summer20UL17NanoV15_V1_MC_RelativePtEC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC2_2017", "Summer20UL17NanoV15_V1_MC_RelativePtEC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtHF", "Summer20UL17NanoV15_V1_MC_RelativePtHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativeBal", "Summer20UL17NanoV15_V1_MC_RelativeBal_AK4PFPuppi"},
    {"CMS_scale_j_RelativeSample_2017", "Summer20UL17NanoV15_V1_MC_RelativeSample_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatEC_2017", "Summer20UL17NanoV15_V1_MC_RelativeStatEC_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatFSR_2017", "Summer20UL17NanoV15_V1_MC_RelativeStatFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatHF_2017", "Summer20UL17NanoV15_V1_MC_RelativeStatHF_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionECAL", "Summer20UL17NanoV15_V1_MC_SinglePionECAL_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionHCAL", "Summer20UL17NanoV15_V1_MC_SinglePionHCAL_AK4PFPuppi"},
    {"CMS_scale_j_TimePtEta_2017", "Summer20UL17NanoV15_V1_MC_TimePtEta_AK4PFPuppi"},
};

const std::map<std::string, std::string> kJesFullSetUnc_2018 = {
    {"CMS_scale_j_AbsoluteMPFBias", "Summer20UL18NanoV15_V1_MC_AbsoluteMPFBias_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteScale", "Summer20UL18NanoV15_V1_MC_AbsoluteScale_AK4PFPuppi"},
    {"CMS_scale_j_AbsoluteStat_2018", "Summer20UL18NanoV15_V1_MC_AbsoluteStat_AK4PFPuppi"},
    {"CMS_scale_j_FlavorQCD", "Summer20UL18NanoV15_V1_MC_FlavorQCD_AK4PFPuppi"},
    {"CMS_scale_j_Fragmentation", "Summer20UL18NanoV15_V1_MC_Fragmentation_AK4PFPuppi"},
    {"CMS_scale_j_PileUpDataMC", "Summer20UL18NanoV15_V1_MC_PileUpDataMC_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtBB", "Summer20UL18NanoV15_V1_MC_PileUpPtBB_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC1", "Summer20UL18NanoV15_V1_MC_PileUpPtEC1_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtEC2", "Summer20UL18NanoV15_V1_MC_PileUpPtEC2_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtHF", "Summer20UL18NanoV15_V1_MC_PileUpPtHF_AK4PFPuppi"},
    {"CMS_scale_j_PileUpPtRef", "Summer20UL18NanoV15_V1_MC_PileUpPtRef_AK4PFPuppi"},
    {"CMS_scale_j_RelativeFSR", "Summer20UL18NanoV15_V1_MC_RelativeFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC1_2018", "Summer20UL18NanoV15_V1_MC_RelativeJEREC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJEREC2_2018", "Summer20UL18NanoV15_V1_MC_RelativeJEREC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativeJERHF", "Summer20UL18NanoV15_V1_MC_RelativeJERHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtBB", "Summer20UL18NanoV15_V1_MC_RelativePtBB_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC1_2018", "Summer20UL18NanoV15_V1_MC_RelativePtEC1_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtEC2_2018", "Summer20UL18NanoV15_V1_MC_RelativePtEC2_AK4PFPuppi"},
    {"CMS_scale_j_RelativePtHF", "Summer20UL18NanoV15_V1_MC_RelativePtHF_AK4PFPuppi"},
    {"CMS_scale_j_RelativeBal", "Summer20UL18NanoV15_V1_MC_RelativeBal_AK4PFPuppi"},
    {"CMS_scale_j_RelativeSample_2018", "Summer20UL18NanoV15_V1_MC_RelativeSample_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatEC_2018", "Summer20UL18NanoV15_V1_MC_RelativeStatEC_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatFSR_2018", "Summer20UL18NanoV15_V1_MC_RelativeStatFSR_AK4PFPuppi"},
    {"CMS_scale_j_RelativeStatHF_2018", "Summer20UL18NanoV15_V1_MC_RelativeStatHF_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionECAL", "Summer20UL18NanoV15_V1_MC_SinglePionECAL_AK4PFPuppi"},
    {"CMS_scale_j_SinglePionHCAL", "Summer20UL18NanoV15_V1_MC_SinglePionHCAL_AK4PFPuppi"},
    {"CMS_scale_j_TimePtEta_2018", "Summer20UL18NanoV15_V1_MC_TimePtEta_AK4PFPuppi"},
};

// Resolves a CMS NP name to its correctionlib key for the given RunPeriod.
// Empty return = not found (caller warns and skips the systematic).
std::string LookupJesFullSetUncertaintyKey(const std::string& runPeriod, const std::string& npName) {
    const std::map<std::string, std::string>* table = nullptr;
    if (runPeriod.find("2016PreVFP") != std::string::npos) {
        table = &kJesFullSetUnc_2016Pre;
    } else if (runPeriod.find("2016PostVFP") != std::string::npos) {
        table = &kJesFullSetUnc_2016Post;
    } else if (runPeriod.find("2017") != std::string::npos) {
        table = &kJesFullSetUnc_2017;
    } else if (runPeriod.find("2018") != std::string::npos) {
        table = &kJesFullSetUnc_2018;
    } else {
        return "";
    }
    auto it = table->find(npName);
    return (it != table->end()) ? it->second : "";
}
}  // namespace

SSBCorrections::SSBCorrections(TextReader* reader, const std::string inputfileName) {
    logger_.Debug() << "TextReader in SSBCorrections" << std::endl;
    logger_.Debug() << "Current directory: " << std::filesystem::current_path() << std::endl;
    reader->PrintoutVariables();


    // JSON root: new CAT cvmfs > legacy jsonpog-integration cvmfs > local copy.
    // Config *Path values are relative to whichever is picked.
    std::string jsonDir;
    const std::string catCvmfsPath  = "/cvmfs/cms-griddata.cern.ch/cat/metadata/";
    const std::string oldCvmfsPath  = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";
    if (std::filesystem::exists(catCvmfsPath)) {
        logger_.Info() << "Using new CAT CVMFS path for JSONs: " << catCvmfsPath << std::endl;
        jsonDir = catCvmfsPath;
    } else if (std::filesystem::exists(oldCvmfsPath)) {
        logger_.Info() << "New CAT CVMFS path not found - using legacy jsonpog-integration CVMFS path: "
                  << oldCvmfsPath << std::endl;
        jsonDir = oldCvmfsPath;
    } else {
        std::string localPath = std::filesystem::current_path().string() + "/jsonpog-integration/POG/";
        logger_.Warning() << "No CVMFS JSON path found. Falling back to local path: " << localPath << std::endl;
        jsonDir = localPath;
    }

    std::string puw_path     = reader->GetText("PUWeightPath");
    std::string jec_path     = reader->GetText("JECPath");
    std::string jer_path     = reader->GetText("JERPath");
    std::string jmar_path     = reader->GetText("JMARPath");
    std::string jveto_path    = reader->GetText("JetVetoPath");
    //std::string jveto_key    = reader->GetText("JetVetoName");
    //std::string jveto_map_key = reader->GetText("JetVetoKey");
    std::string jveto_name    = reader->GetText("JetVetoName");
    std::string jveto_map_key = reader->GetText("JetVetoKey");
    std::string jveto_type     = reader->GetText("JetVetoType");
    std::string jer_sf_path  = reader->GetText("JERSFPath");
    std::string jec_name     = reader->GetText("JECName");
    std::string jer_name     = reader->GetText("JERName");
    std::string jer_res_name = reader->GetText("JERResName");
    std::string jet_btag_conf = reader->GetText("Jet_btag");  // e.g., "deepCSVL", "deepJetM"
    std::string muon_path    = reader->GetText("MuonSFPath");
    std::string elec_path    = reader->GetText("ElecSFPath");
    std::string RunPeriod    = reader->GetText("RunRange");
    std::string puJson  =  "";
    btag_sf_type_            = reader->GetText("BTagSFType");

    // Muon Infor. ID  ISO // 
    std::string muon_id_corName  = reader->GetText( "MuonIDSFName"  );
    std::string muon_iso_corName = reader->GetText( "MuonIsoSFName" );
    
    // Electron ID ISO // 
    std::string ele_sf_name_      = reader->GetText( "ElecIDSFName"  );  
    //std::string ele_reco_sf_name_    = reader->GetText( "ElecRecoSFName" );

    // Trigger SF // 
    std::string Trig_sf_name_      = reader->GetText( "TrigSFFile"  );  
    std::string Trig_sf_histname_  = reader->GetText( "TrigSFHist"  );  



    //if inputfileName    

    if (RunPeriod.find("2016PreVFP") != std::string::npos) {
        year_ = "2016preVFP";
        puJson = "Collisions16_UltraLegacy_goldenJSON";
        
    } else if (RunPeriod.find("2016PostVFP") != std::string::npos) {
        year_ = "2016postVFP";
        puJson = "Collisions16_UltraLegacy_goldenJSON";
    } else if (RunPeriod.find("2017") != std::string::npos) {
        year_ = "2017";
        puJson = "Collisions17_UltraLegacy_goldenJSON";
    } else if (RunPeriod.find("2018") != std::string::npos) {
        year_ = "2018";
        puJson = "Collisions18_UltraLegacy_goldenJSON";
    } 
    else {
        logger_.Error() << "Unknown RunPeriod: " << RunPeriod << std::endl;
        year_ = "2018"; // fallback or throw error
    }


    bool is_data = inputfileName.find("Data") != std::string::npos;

    std::string era = "Unknown";

    if (is_data) {
   	 // For data: extract year and era from the input file name, e.g. "Run2016F-HIPM"
    	size_t run_pos = inputfileName.find("Run");
    	if (run_pos != std::string::npos && run_pos + 7 <= inputfileName.size()) {
        	// Extract year: 4 digits following "Run"
        	std::string year_digits = inputfileName.substr(run_pos + 3, 4); // e.g., "2016", "2017"

        	// Extract era: characters immediately after "RunYYYY"
        	era = inputfileName.substr(run_pos + 7); // e.g., "B", "F-HIPM"
        	size_t delim = era.find_first_of("._/");
        	if (delim != std::string::npos) {
            		era = era.substr(0, delim); // Trim suffix after era
        	}
    	}
    } else {
        // For MC: use the RunPeriod directly (e.g., "2016PreVFP", "2016PostVFP", "2017", "2018")
       era = "MC";  // Dummy placeholder; era is not needed for MC in most JEC logic
    }

    logger_.Debug() << "era : " << era <<  std::endl;
    logger_.Debug() << "jec_name: " << jec_name << std::endl;
    // Optional overrides; default to NanoAODv15 Puppi ("AK4PFPuppi"/"L1L2L3Res").
    std::string jetAlgoTag = reader->Check("JetAlgoTag") ? reader->GetText("JetAlgoTag") : "AK4PFPuppi";
    std::string jecLevelTag = reader->Check("JECLevel") ? reader->GetText("JECLevel") : "L1L2L3Res";
    // Keep the un-expanded base tag to also derive the L1FastJet-only name
    // below - jec_name gets overwritten with the full name next line.
    std::string jecBaseNameForL1 = jec_name;
    jec_name = ExpandJECName(jec_name, RunPeriod, era, is_data, jetAlgoTag, jecLevelTag);
    logger_.Debug() << "after ExpandJECName jec_name: " << jec_name << std::endl;
    std::string jecL1Name = ExpandJECName(jecBaseNameForL1, RunPeriod, era, is_data, jetAlgoTag, "L1FastJet");
    


    logger_.Debug() << "jsonDir + puw_path : " << jsonDir + puw_path << std::endl;
    auto puset = correction::CorrectionSet::from_file(jsonDir + puw_path);
    
    pu_weight_ = puset->at(puJson);


    auto jec_set = correction::CorrectionSet::from_file(jsonDir + jec_path);
    jec_ = jec_set->compound().at(jec_name);

    // Determine whether this compound correction needs a 5th "run" input (v15
    // Puppi DATA does; MC and old v9/PFchs bake era into the name instead).
    // Read this from the correction's own schema (inputs()) instead of trial-
    // evaluating with dummy values: a dummy run like "1" is always outside a
    // run-binned correction's valid range and throws a Binning bounds error
    // that looks just like a "wrong signature" failure, so probing can
    // misdetect a genuine 5-input (DATA) correction as 4-input (MC-only).
    jec_needs_run_ = false;
    {
        const auto& jec_inputs = jec_->inputs();
        std::ostringstream oss;
        for (const auto& var : jec_inputs) {
            oss << var.name() << " ";
            if (var.name() == "run") jec_needs_run_ = true;
        }
        logger_.Info() << "JEC compound correction '" << jec_name << "' declares "
                  << jec_inputs.size() << " input(s): " << oss.str()
                  << (jec_needs_run_ ? "- 'run' present, will pass the event run number."
                                      : "- no 'run' input.") << std::endl;
    }

    // L1FastJet-only correction - Type-1 MET baseline, same file as jec_.
    jec_l1_ = LoadOptionalCorrection<correction::Correction>(
        [&]() { return jec_set->at(jecL1Name); },
        "Could not load L1FastJet-only correction '" + jecL1Name + "' from " + jec_path +
        ". Type-1 MET recomputation will fall back to a (fullyCorrected-raw) delta instead "
        "of the CMS-recommended (fullyCorrected-L1only) one - this is a known-less-accurate "
        "fallback, not the standard recipe.", logger_);

    // JESSys/JESSysDir are new config keys; Check() first so an old config
    // without them just gets "no JES systematic" instead of a GetText/GetBool error.
    {
        std::string jesSys = reader->Check("JESSys") ? reader->GetText("JESSys") : "nominal";
        jes_sys_dir_ = reader->Check("JESSysDir") ? reader->GetText("JESSysDir") : "up";
        if (jesSys != "nominal" && !jesSys.empty()) {
            std::string uncKey = LookupJesFullSetUncertaintyKey(RunPeriod, jesSys);
            if (uncKey.empty()) {
                logger_.Warning() << "JESSys=\"" << jesSys << "\" not found in the Full-set "
                          << "NP-name table for RunPeriod=\"" << RunPeriod << "\" (see "
                          << "LookupJesFullSetUncertaintyKey in this file) - no JES systematic "
                          << "will be applied. Check spelling (e.g. \"CMS_scale_j_AbsoluteScale\") "
                          << "and that this RunPeriod has a Full-set table." << std::endl;
            } else if (jes_sys_dir_ != "up" && jes_sys_dir_ != "down") {
                logger_.Warning() << "JESSysDir=\"" << jes_sys_dir_ << "\" is not \"up\" or "
                          << "\"down\" - no JES systematic will be applied." << std::endl;
            } else {
                jes_unc_source_ = LoadOptionalCorrection<correction::Correction>(
                    [&]() { return jec_set->at(uncKey); },
                    "Could not load JES uncertainty-source correction '" + uncKey +
                    "' (JESSys=" + jesSys + ") from " + jec_path + ". No JES systematic will be applied.", logger_);
                if (jes_unc_source_) {
                    logger_.Info() << "Loaded JES uncertainty source '" << uncKey
                              << "' for JESSys=" << jesSys << " JESSysDir=" << jes_sys_dir_
                              << " - applying to GetCorrectedJetPt() (physics jets + Type-1 MET)."
                              << std::endl;
                }
            }
        }
    }

    auto jer_set = correction::CorrectionSet::from_file(jsonDir + jer_path);
    jer_ = jer_set->at(jer_res_name);

    auto jer_sf_set = correction::CorrectionSet::from_file(jsonDir + jer_sf_path);
    jer_sf_ = jer_sf_set->at(jer_name);

    // JER SF uncertainty is a separate correction from jer_sf_ (not a tag
    // input to it); combined as sf*(1+-sfUnc). Optional - some campaigns don't publish it.
    std::string jerSfUncName = jer_name;
    {
        size_t pos = jerSfUncName.find("ScaleFactor");
        if (pos != std::string::npos) {
            jerSfUncName.replace(pos, std::string("ScaleFactor").size(), "SFUncertainty");
        }
    }
    jer_sfunc_ = LoadOptionalCorrection<correction::Correction>(
        [&]() { return jer_sf_set->at(jerSfUncName); },
        "Could not load JER SF uncertainty correction (derived from JERName by replacing "
        "'ScaleFactor' with 'SFUncertainty'). JER SF up/down variations will not be available.", logger_);

    // CMS's official JERSmear tool, usually a separate file from
    // jet_jerc.json.gz; SmearJER() falls back in-house if unavailable.
    std::string jerSmearPath = reader->Check("JERSmearPath") ? reader->GetText("JERSmearPath") : "";
    std::string jerSmearName = reader->Check("JERSmearName") ? reader->GetText("JERSmearName") : "JERSmear";
    if (!jerSmearPath.empty()) {
        jer_smear_ = LoadOptionalCorrection<correction::Correction>(
            [&]() {
                auto jerSmearSet = correction::CorrectionSet::from_file(jsonDir + jerSmearPath);
                return jerSmearSet->at(jerSmearName);
            },
            "Could not load JERSmearPath=" + jerSmearPath +
            ". Falling back to an in-house JER smearing implementation.", logger_);
        if (jer_smear_) {
            logger_.Info() << "Loaded JER smearing tool '" << jerSmearName << "' from "
                      << jerSmearPath << " - using it for JER smearing (CMS-recommended)." << std::endl;
        }
    } else {
        logger_.Warning() << "JERSmearPath not set in config - falling back to an in-house "
                  << "JER smearing implementation instead of CMS's official JERSmear tool. "
                  << "Set JERSmearPath (e.g. 'JME/Run2-2018-UL-NanoAODv15/latest/jer_smear.json.gz') "
                  << "to use the recommended approach." << std::endl;
        jer_smear_ = nullptr;
    }

    // JMAR (PU jet ID SF) - optional; GetPUJetIDSFAndEff() returns 1.0 if unloaded.
    pujetid_sf_ = LoadOptionalCorrection<correction::Correction>(
        [&]() {
            auto jmar_sf_set = correction::CorrectionSet::from_file(jsonDir + jmar_path);
            return jmar_sf_set->at("PUJetID_eff");
        },
        "Could not load JMARPath=" + jmar_path +
        ". PU jet ID SF will be unavailable (GetPUJetIDSFAndEff() returns 1.0 - no PU-ID "
        "reweighting applied).", logger_);
    if (pujetid_sf_) {
        logger_.Info() << "Loaded PU jet ID SF from " << jmar_path << std::endl;
    }

    jveto_name_ = jveto_name;
    jveto_key_ = jveto_map_key;
    jveto_type_ = jveto_type;	

    if (!jveto_path.empty() && !jveto_name.empty()) {
        jetvetomap_ = LoadOptionalCorrection<correction::Correction>(
            [&]() {
                auto jetveto_set = correction::CorrectionSet::from_file(jsonDir + jveto_path);
                return jetveto_set->at(jveto_name);
            },
            "Failed to load jet veto map", logger_);
        if (jetvetomap_) {
            logger_.Info() << "Loaded jet veto map: " << jveto_name
                      << " with key: " << jveto_map_key
                      << " type: " << jveto_type << std::endl;
        }
    } else {
       logger_.Info() << "Jet veto map not configured, skipping..." << std::endl;
       jetvetomap_ = nullptr;
    }

    if (btag_sf_type_ != "comb" && btag_sf_type_ != "mujets") {
        logger_.Warning() << "Invalid BTagSFType: " << btag_sf_type_
                  << ". Using default 'comb'." << std::endl;
        btag_sf_type_ = "comb";
    }

    // Analysis type specific recommendation
    if (btag_sf_type_ == "comb") {
        logger_.Info() << "Using 'comb' SF (QCD + ttbar enriched regions)" << std::endl;
    } else {
        logger_.Info() << "Using 'mujets' SF (QCD enriched regions, bias avoidance)" << std::endl;
    }


    // Single canonical parse (ParseBTagAlgo, in SSBCorrections.h) shared with
    // Analysis.cpp's own btag_algo_ parsing, so the two can't drift out of
    // sync the way they did once before (see NOTES.md).
    BTagAlgo btag_algo_kind = ParseBTagAlgo(jet_btag_conf);
    if (btag_algo_kind == BTagAlgo::Unknown) {
        logger_.Warning() << "Unknown b-tag algorithm in Jet_btag: " << jet_btag_conf << std::endl;
    }
    std::string btag_algo = BTagAlgoToString(btag_algo_kind);
    std::string btag_wp = "";

    char last = jet_btag_conf.back();
    if (last == 'L' || last == 'l') btag_wp = "Loose";
    else if (last == 'M' || last == 'm') btag_wp = "Medium";
    else if (last == 'T' || last == 't') btag_wp = "Tight";
    else {
        logger_.Warning() << "Unknown WP in Jet_btag: " << jet_btag_conf << std::endl;
    }

    // BTagSFJsonPath config override takes priority; else guess the CAT or
    // legacy jsonpog-integration path scheme from jsonDir.
    std::string btag_sf_json;
    if (reader->Check("BTagSFJsonPath")) {
        btag_sf_json = reader->GetText("BTagSFJsonPath");
    } else if (jsonDir == catCvmfsPath) {
        btag_sf_json = "BTV/Run2-" + year_ + "-UL-NanoAODv15/latest/btagging.json.gz";
    } else {
        btag_sf_json = "BTV/" + year_ + "_UL/btagging.json.gz";
    }

    // UParT's tagger name isn't verified against btagging.json.gz; override via
    // config key BTagTaggerName if the guessed name is wrong.
    std::string btag_tagger;
    if (reader->Check("BTagTaggerName")) {
        btag_tagger = reader->GetText("BTagTaggerName");
    } else if (btag_algo_kind == BTagAlgo::DeepJet) {
        btag_tagger = "deepJet_" + btag_sf_type_;
    } else if (btag_algo_kind == BTagAlgo::UParTAK4) {
        // UParT only ships a single "comb" heavy-flavor SF ("UParTAK4_comb",
        // no "_mujets" variant) - BTagSFType has no effect on the tagger name here.
        if (btag_sf_type_ != "comb") {
            logger_.Info() << "btag_algo=UParT: BTagSFType='" << btag_sf_type_
                      << "' has no effect - your btagging.json.gz only provides a 'comb' "
                      << "b/c-jet SF for UParTAK4 (no 'mujets' variant exists)." << std::endl;
        }
        btag_tagger = "UParTAK4_comb";
    } else {
        btag_tagger = "deepCSV_" + btag_sf_type_;
    }
    std::string btag_eff_path = "";
    if (!is_data) {
        std::string process_subdir = GetProcessSubDir(inputfileName);
        btag_eff_path = "CorrectionFiles/BTag/UL" + RunPeriod
                      + "/" + process_subdir
                      + "/btagEff_" + btag_algo + ".root";
        logger_.Info() << "btag_eff_path: " << btag_eff_path << std::endl;
    }

    InitBtagSFCorrection(jsonDir + btag_sf_json, btag_tagger);
    if (!is_data) {
        //LoadMCBtagEfficiencies(btag_eff_path, btag_algo);
        LoadMCBtagEfficiencies(btag_eff_path, btag_algo, btag_wp);

    }


    // Load muon SF
    //auto muon_set = CorrectionSet::from_file(jsonDir+muon_path);
    logger_.Debug() << "jsonDir+muon_path " << jsonDir+muon_path << std::endl;
    auto muon_set = correction::CorrectionSet::from_file(jsonDir + muon_path);

    logger_.Debug() << "muon_id_corName : " << muon_id_corName << std::endl;
    logger_.Debug() << "muon_iso_corName : " << muon_iso_corName << std::endl;

    muon_id_   = muon_set->at(muon_id_corName);
    muon_iso_  = muon_set->at(muon_iso_corName);

    auto cset = correction::CorrectionSet::from_file(jsonDir + elec_path);
    ele_sf_ = cset->at(ele_sf_name_);
    std::string TrigSFPath = std::filesystem::current_path().string() + "/CorrectionFiles/Trig/";
    {
        std::unique_ptr<TFile> f_trg(TFile::Open((TrigSFPath + Trig_sf_name_).c_str()));
        if (!f_trg || f_trg->IsZombie()) {
            throw std::runtime_error("SSBCorrections: could not open trigger SF file: " + TrigSFPath + Trig_sf_name_);
        }
        auto* h = dynamic_cast<TH2*>(f_trg->Get(Trig_sf_histname_.c_str()));
        if (!h) {
            throw std::runtime_error("SSBCorrections: trigger SF histogram '" + Trig_sf_histname_ +
                                      "' not found (or not a TH2-derived type) in " + TrigSFPath + Trig_sf_name_);
        }
        H_trig.reset(static_cast<TH2*>(h->Clone()));
        H_trig->SetDirectory(nullptr);  // detach from f_trg's TDirectory before f_trg closes below
    }

    std::string rochesterCorrFile;
    logger_.Debug() << "RunPeriod : " << RunPeriod << std::endl;
    if (RunPeriod.find("2016Pre")  != std::string::npos) {    rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2016aUL.txt";}
    else if (RunPeriod.find("2016Post") != std::string::npos) { rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2016bUL.txt";}
    else if (RunPeriod.find("2017")     != std::string::npos) { rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2017UL.txt";}
    else if (RunPeriod.find("2018")     != std::string::npos) { rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2018UL.txt";}
    else {
        logger_.Error() << "Unknown RunPeriod in rochesterCorrFile : " << RunPeriod << std::endl;
    }

    logger_.Debug() << "rochesterCorrFile : " << rochesterCorrFile << std::endl;
    std::ifstream file(rochesterCorrFile);
    if(!file.good()){
            logger_.Error() << "Rochester file is not found: " << rochesterCorrFile << std::endl;
            throw std::runtime_error("SSBCorrections: Rochester correction file not found: " + rochesterCorrFile);
    }

    rc.init(rochesterCorrFile);
}

SSBCorrections::~SSBCorrections() {
    // H_trig and eff_histograms_ are unique_ptr-owned - cleaned up
    // automatically, nothing to do here manually.
    logger_.Info() << "Destructor called - cleaning up " << eff_histograms_.size()
              << " efficiency histograms..." << std::endl;
}

double SSBCorrections::GetCorrectedJetPt(double raw_pt, double eta, double area, double rho, unsigned int run_number) const {
    double sf = jec_needs_run_
        ? jec_->evaluate({area, eta, raw_pt, rho, static_cast<double>(run_number)})
        : jec_->evaluate({area, eta, raw_pt, rho});
    double corrected_pt = raw_pt * sf;

    // JES full-uncertainty-set shift, applied after nominal JEC so every
    // caller of GetCorrectedJetPt() picks it up consistently. No-op if unset.
    if (jes_unc_source_) {
        double scale = jes_unc_source_->evaluate({eta, corrected_pt});
        double factor = (jes_sys_dir_ == "up") ? (1.0 + scale) : (1.0 - scale);
        corrected_pt *= factor;
    }

    return corrected_pt;
}

double SSBCorrections::GetL1CorrectedJetPt(double raw_pt, double eta, double area, double rho) const {
    if (!jec_l1_) {
        // Not loaded - fall back to the less-accurate (corrected-raw) delta.
        return raw_pt;
    }
    double c1 = jec_l1_->evaluate({area, eta, raw_pt, rho});
    return raw_pt * c1;
}

double SSBCorrections::GetCorrectedJetMass(double raw_mass, double raw_pt, double eta, double area, double rho, unsigned int run_number) const {
    //double sf = jec_->evaluate({eta, raw_pt, area});
    double sf = jec_needs_run_
        ? jec_->evaluate({area, eta, raw_pt, rho, static_cast<double>(run_number)})
        : jec_->evaluate({area, eta, raw_pt, rho});
    return raw_mass * sf;
}

double SSBCorrections::GetJER(double eta, double pt, double rho) const {
    return jer_->evaluate({eta, pt, rho});
}

namespace {
// Deterministic per-(event, jet) seed for JER smearing - reproducible across
// re-runs, and reusable for up/down passes to isolate the systematic from
// jet-to-jet statistical fluctuation. Not bit-identical to any official recipe.
UInt_t ComputeJerSeed(ULong64_t event, double eta, double phi) {
    ULong64_t etaBits = static_cast<ULong64_t>(std::llround(eta * 1.0e4));
    ULong64_t phiBits = static_cast<ULong64_t>(std::llround(phi * 1.0e4));
    ULong64_t mixed = (event * 2654435761ULL)
                     ^ (etaBits * 40503ULL)
                     ^ (phiBits * 2246822519ULL);
    UInt_t seed = static_cast<UInt_t>(mixed & 0xFFFFFFFFULL);
    return seed == 0 ? 1u : seed;  // 0 means "seed from entropy" to TRandom3
}
}

double SSBCorrections::SmearJER(double reco_pt, double gen_pt, double gen_eta, double gen_phi,
                                 double eta, double phi, double rho,
                                 ULong64_t event, const std::string& jer_tag) const {
    // jer_sf_ is a plain 2-input (eta, pt) schema; up/down comes from the
    // separate jer_sfunc_ correction, combined as sf*(1+-unc).
    double sf = jer_sf_->evaluate({eta, reco_pt});
    if (jer_sfunc_ && jer_tag != "nominal") {
        double sfUnc = jer_sfunc_->evaluate({eta, reco_pt});
        if (jer_tag == "up")        sf = sf * (1.0 + sfUnc);
        else if (jer_tag == "down") sf = sf * (1.0 - sfUnc);
    }
    double resolution = jer_->evaluate({eta, reco_pt, rho});

    // Re-validate the gen match: dR<0.2 AND |reco_pt-gen_pt|<3*resolution*reco_pt,
    // checked here since JERSmear's evaluate() schema has no room for gen eta/phi.
    bool matched = false;
    if (gen_pt >= 0.0) {
        TLorentzVector jetDir, genDir;
        jetDir.SetPtEtaPhiM(1.0, eta, phi, 0.0);
        genDir.SetPtEtaPhiM(1.0, gen_eta, gen_phi, 0.0);
        double dR = jetDir.DeltaR(genDir);
        matched = (dR < 0.2) && (std::abs(reco_pt - gen_pt) < 3.0 * resolution * reco_pt);
    }
    double genPtForSmear = matched ? gen_pt : -1.0;

    if (jer_smear_) {
        // CMS's official JERSmear tool - deterministic hashprng-based smearing.
        double smear = jer_smear_->evaluate({reco_pt, eta, genPtForSmear, rho,
                                              static_cast<int>(event), resolution, sf});
        double corr = (std::isfinite(smear) && smear > 0.0) ? smear : 1.0;
        return std::max(0.0, reco_pt * corr);
    }

    // Fallback if jer_smear_ isn't loaded: in-house hybrid method, deterministic
    // via a per-(event,jet) seeded TRandom3.
    if (matched) {
        double delta_pt = reco_pt - genPtForSmear;
        double smeared_pt = genPtForSmear + sf * delta_pt;
        return std::max(0.0, smeared_pt);
    }

    // sf can be < 1 in some eta/pt bins (MC resolution needs to be improved,
    // i.e. scaled down to match data) - clamp to 0 rather than letting
    // sf*sf-1 go negative and produce NaN from sqrt().
    double stochasticTerm = std::sqrt(std::max(sf * sf - 1.0, 0.0));
    TRandom3 rnd(ComputeJerSeed(event, eta, phi));
    double smear_factor = 1.0 + stochasticTerm * rnd.Gaus(0, 1) * resolution;
    return std::max(0.0, reco_pt * smear_factor);
}

float SSBCorrections::GetPUJetIDSFAndEff(float pt, float eta, bool passPU, bool genMatched, const std::string& wp, const std::string& syst, bool getEff) const {
    if (!pujetid_sf_) {
        logger_.Error() << "[GetPUJetIDSFAndEff] PUJetID correction not loaded." << std::endl;
        return 1.0;
    }
    
    if (pt >= 50.0) return 1.0;
    
    try {
        std::string eval_type = getEff ? "MCEff" : syst;  // Efficiency or Scale Factor
        std::variant<double, std::vector<double>> val = pujetid_sf_->evaluate({eta, pt, eval_type, wp});
        return std::get<double>(val);
    } catch (const std::exception& e) {
        logger_.Error() << "[GetPUJetIDSFAndEff] Evaluation failed: " << e.what() << std::endl;
        return 1.0;
    }
}


double SSBCorrections::GetMuonRecoSF(double pt, double eta) const {
    std::variant<double, std::vector<double>> val = muon_reco_->evaluate({pt, eta});
    return std::get<double>(val);
}

double SSBCorrections::GetMuonIDSF(double pt, double eta, const std::string& syst_tag) const {
    std::variant<double, std::vector<double>> val = muon_id_->evaluate({eta, pt, syst_tag});
    return std::get<double>(val);
}

double SSBCorrections::GetMuonIsoSF(double pt, double eta, const std::string& syst_tag) const {
    std::variant<double, std::vector<double>> val = muon_iso_->evaluate({eta, pt, syst_tag});
    return std::get<double>(val);
}


double SSBCorrections::DoubleMuon_IDIsoEff(TLorentzVector lep1, TLorentzVector lep2,
                                           TString muidsys, TString muisosys, TString tracksys) const {
    float pt1 = std::min(lep1.Pt(), 119.999);
    float pt2 = std::min(lep2.Pt(), 119.999);
    float abseta1 = std::abs(lep1.Eta());
    float abseta2 = std::abs(lep2.Eta());

    std::string IDSyst = "nominal", IsoSyst = "nominal";

    if (muidsys.Contains("up", TString::kIgnoreCase)) IDSyst = "up";
    else if (muidsys.Contains("down", TString::kIgnoreCase)) IDSyst = "down";

    if (muisosys.Contains("up", TString::kIgnoreCase)) IsoSyst = "up";
    else if (muisosys.Contains("down", TString::kIgnoreCase)) IsoSyst = "down";
    double mu1id  = GetMuonIDSF(pt1, abseta1, IDSyst);
    double mu2id  = GetMuonIDSF(pt2, abseta2, IDSyst);
    double mu1iso = GetMuonIsoSF(pt1, abseta1, IsoSyst);
    double mu2iso = GetMuonIsoSF(pt2, abseta2, IsoSyst);

    // track SF is 1.0, so you don't need to apply them... 
    double mu1trk =1.0;// TrackSF(lep1->Eta());
    double mu2trk =1.0;// TrackSF(lep2->Eta());
/*
    if (tracksys.Contains("up", TString::kIgnoreCase)) {
        mu1trk += TrackSFErr(lep1->Eta(), tracksys);
        mu2trk += TrackSFErr(lep2->Eta(), tracksys);
    } else if (tracksys.Contains("down", TString::kIgnoreCase)) {
        mu1trk -= TrackSFErr(lep1->Eta(), tracksys);
        mu2trk -= TrackSFErr(lep2->Eta(), tracksys);
    }
*/
    return mu1id * mu2id * mu1iso * mu2iso * mu1trk * mu2trk;
}

double SSBCorrections::DoubleElec_Eff(
    const TLorentzVector& lep1, const TLorentzVector& lep2,
    double ele1sueta, double ele2sueta, 
    const std::string& id_wp,        // "Tight", "Medium", "Loose"
    const std::string& id_syst,      // "nominal", "up", "down"
    const std::string& reco_syst     // "nominal", "up", "down"
) const {
    
    // Clamp to the JSON's valid pt/eta ranges.
    float lep1pt = std::clamp(static_cast<float>(lep1.Pt()), 10.0f, 999.0f);
    float lep2pt = std::clamp(static_cast<float>(lep2.Pt()), 10.0f, 999.0f);
    float lep1sueta_clamped = std::clamp(static_cast<float>(ele1sueta), -3.0f, 3.0f);
    float lep2sueta_clamped = std::clamp(static_cast<float>(ele2sueta), -3.0f, 3.0f);

    std::string actual_id_wp = id_wp;
    if (actual_id_wp.empty()) {
        actual_id_wp = "Tight";  // Set default working point
        logger_.Warning() << "Empty ID working point, using default: " << actual_id_wp << std::endl;
    }
    
    // Use working point directly (not with "ID" prefix)
    float ele1id = GetElectronSF(actual_id_wp, lep1sueta_clamped, lep1pt, id_syst);
    float ele2id = GetElectronSF(actual_id_wp, lep2sueta_clamped, lep2pt, id_syst);
    
    // Step 3: Calculate Reco SF for each electron using GetElectronSF
    // "Reco" type will automatically choose RecoAbove20/RecoBelow20 based on pt
    float ele1reco = GetElectronSF("Reco", lep1sueta_clamped, lep1pt, reco_syst);
    float ele2reco = GetElectronSF("Reco", lep2sueta_clamped, lep2pt, reco_syst);
    
    // Step 4: Isolation SF is 1.0 (typically included in ID for electrons)
    float ele1iso = 1.0f;
    float ele2iso = 1.0f;
    
    // Step 5: Compute total double electron efficiency
    double doubleEleff = static_cast<double>(ele1id) * static_cast<double>(ele2id) * 
                        static_cast<double>(ele1iso) * static_cast<double>(ele2iso) * 
                        static_cast<double>(ele1reco) * static_cast<double>(ele2reco);
    
    return doubleEleff;
}

double SSBCorrections::MuonElec_Eff(const TLorentzVector& muon, const TLorentzVector& electron,
                                    double muon_eta, double electron_sueta,
                                    const std::string& mu_id_syst, 
                                    const std::string& mu_iso_syst,
                                    const std::string& ele_id_wp,
                                    const std::string& ele_id_syst, 
                                    const std::string& ele_reco_syst) const {
    
    // Step 1: Apply pt/eta limits based on JSON ranges
    // Muon: pt range typically up to ~120 GeV in measurements
    float muon_pt = std::clamp(static_cast<float>(muon.Pt()), 15.0f, 119.999f);
    float muon_abseta = std::clamp(std::abs(static_cast<float>(muon_eta)), 0.0f, 2.4f);
    
    // Electron: pt range 10.0 to Infinity in JSON, eta range within detector acceptance
    float electron_pt = std::clamp(static_cast<float>(electron.Pt()), 10.0f, 999.0f);
    float electron_sueta_clamped = std::clamp(static_cast<float>(electron_sueta), -3.0f, 3.0f);

    // Step 2: Calculate Muon ID SF
    float mu_id = GetMuonIDSF(muon_pt, muon_abseta, mu_id_syst);
    
    // Step 3: Calculate Muon Iso SF  
    float mu_iso = GetMuonIsoSF(muon_pt, muon_abseta, mu_iso_syst);
    
    // Step 4: Muon tracking SF (typically 1.0 for current analyses)
    float mu_trk = 1.0f;

    // Step 5: Calculate Electron ID SF
    std::string actual_ele_id_wp = ele_id_wp;
    if (actual_ele_id_wp.empty()) {
        actual_ele_id_wp = "Tight";  // Set default working point
        logger_.Warning() << "Empty electron ID working point, using default: " << actual_ele_id_wp << std::endl;
    }
    
    float ele_id = GetElectronSF(actual_ele_id_wp, electron_sueta_clamped, electron_pt, ele_id_syst);
    
    // Step 6: Calculate Electron Reco SF
    // "Reco" type will automatically choose RecoAbove20/RecoBelow20 based on pt
    float ele_reco = GetElectronSF("Reco", electron_sueta_clamped, electron_pt, ele_reco_syst);
    
    // Step 7: Electron isolation SF is 1.0 (typically included in ID for electrons)
    float ele_iso = 1.0f;

    // Step 8: Compute total muon-electron efficiency
    double muonelec_eff = static_cast<double>(mu_id) * static_cast<double>(mu_iso) * 
                         static_cast<double>(mu_trk) * static_cast<double>(ele_id) * 
                         static_cast<double>(ele_iso) * static_cast<double>(ele_reco);

    return muonelec_eff;
}




float SSBCorrections::GetElectronSF(const std::string& sf_type, float eta, float pt, const std::string& syst) const {

    // Map systematic names: "nominal" -> "sf", "up" -> "sfup", "down" -> "sfdown"
    std::string valtype = "sf";  // default

    if (syst == "up") {
        valtype = "sfup";
    } else if (syst == "down") {
        valtype = "sfdown";
    } else if (syst == "nominal" || syst.empty()) {
        valtype = "sf";
    } else {
        // Handle unknown systematic values
        logger_.Warning() << "Unknown systematic: '" << syst << "', using nominal (sf)" << std::endl;
        valtype = "sf";
    }

    // Map working point names from config to JSON format
    std::string working_point = sf_type;
    if (sf_type == "Reco") {
        working_point = (pt >= 20.0f) ? "RecoAbove20" : "RecoBelow20";
    } else if (sf_type.find("SCB") == 0) {
        // Convert SCB (Scale factor Cut-Based) format to JSON format
        if (sf_type == "SCBVeto") working_point = "Veto";
        else if (sf_type == "SCBLoose") working_point = "Loose";
        else if (sf_type == "SCBMedium") working_point = "Medium";
        else if (sf_type == "SCBTight") working_point = "Tight";
        else {
            logger_.Warning() << "Unknown SCB working point: " << sf_type << ", using Tight" << std::endl;
            working_point = "Tight";
        }
    } else if (sf_type.find("MVA") == 0) {
        // Convert MVA format to JSON format if needed
        if (sf_type == "MVALoose") working_point = "wp90noiso";
        else if (sf_type == "MVAMedium") working_point = "wp80noiso";
        else if (sf_type == "MVATight") working_point = "wp90iso";
        else {
            logger_.Warning() << "Unknown MVA working point: " << sf_type << ", using wp90iso" << std::endl;
            working_point = "wp90iso";
        }
    }
    // For other cases (direct JSON names like "Tight", "Medium"), use as-is

    try {
        if (!ele_sf_) {
            logger_.Error() << "[GetElectronSF] ele_sf_ pointer is null!" << std::endl;
            return 1.0;
        }

        // Validate parameters before evaluation
        if (year_.empty()) {
            logger_.Error() << "[GetElectronSF] Year is empty!" << std::endl;
            return 1.0;
        }

        if (working_point.empty()) {
            logger_.Error() << "[GetElectronSF] Working point is empty!" << std::endl;
            return 1.0;
        }

        // Use the correct order from JSON: {year, ValType, WorkingPoint, eta, pt}
        float result = ele_sf_->evaluate({year_, valtype, working_point, eta, pt});

        // Sanity check on result
        if (result <= 0.0 || result > 10.0) {
            logger_.Warning() << "[GetElectronSF] Unusual SF value: " << result
                      << " for parameters: " << working_point << ", eta=" << eta << ", pt=" << pt << std::endl;
        }

        return result;

    } catch (const std::exception& e) {
        logger_.Error() << "[GetElectronSF] Evaluation failed: " << e.what() << std::endl;
        logger_.Error() << "[GetElectronSF] Parameters: year='" << year_ << "', valtype='" << valtype
                  << "', working_point='" << working_point << "', eta=" << eta << ", pt=" << pt << std::endl;
        return 1.0;
    }
}

float SSBCorrections::GetPUWeight(float nTrueInt, const std::string& systTag) const {
    std::string variation = systTag;
    
    if (!pu_weight_) {
        logger_.Error() << "[GetPUWeight] PU correction not loaded." << std::endl;
        return 1.0;
    }

    // If systTag is "Central" or empty, treat as "nominal"
    if (systTag == "Central" || systTag == "") {
        variation = "nominal";
    }

    // If the variation is not valid, fallback to "nominal"
    if (variation != "nominal" && variation != "up" && variation != "down") {
        logger_.Warning() << "[GetPUWeight] Unrecognized systematic '" << variation
                  << "', defaulting to nominal." << std::endl;
        variation = "nominal";
    }

    try {
        std::variant<double, std::vector<double>> val = pu_weight_->evaluate({ nTrueInt, variation });
        double weight = std::get<double>(val);
        return weight;
    } catch (const std::exception& e) {
        logger_.Error() << "[GetPUWeight] Evaluation failed: " << e.what() << std::endl;
        return 1.0;
    }
}


float SSBCorrections::MatchGenPt(const TLorentzVector& reco_jet,
                                  const std::vector<TLorentzVector>& gen_jets,
                                  float maxDR,
                                  float* out_eta,
                                  float* out_phi) const {
    float minDR = maxDR;
    float matched_genpt = -1.0;

    for (const auto& gen_jet : gen_jets) {
        float dR = reco_jet.DeltaR(gen_jet);
        if (dR < minDR) {
            minDR = dR;
            matched_genpt = gen_jet.Pt();
            if (out_eta) *out_eta = gen_jet.Eta();
            if (out_phi) *out_phi = gen_jet.Phi();
        }
    }
    return matched_genpt;
}

std::vector<TLorentzVector> SSBCorrections::ApplyJetCorrections(
    const std::vector<TLorentzVector>& rawJets,
    const std::vector<float>& rawFactors,
    const std::vector<float>& areas,
    float rho,
    bool isData,
    bool applyJES,
    bool applyJER,
    const std::vector<TLorentzVector>& genJets,
    const std::vector<int>& genJetIndices,
    unsigned int run_number,
    ULong64_t event_number,
    const std::string& jerSysTag
) const {
    std::vector<TLorentzVector> correctedJets;
    correctedJets.reserve(rawJets.size());

    for (size_t i = 0; i < rawJets.size(); ++i) {
        double eta  = rawJets[i].Eta();
        double phi  = rawJets[i].Phi();

        double raw_pt   = rawJets[i].Pt() * (1.0 - rawFactors[i]);
        double raw_mass = rawJets[i].M()  * (1.0 - rawFactors[i]);

        double corrected_pt = raw_pt;
        double corrected_mass = raw_mass;

        if (applyJES) {
            corrected_pt   = GetCorrectedJetPt(raw_pt, eta, areas[i], rho, run_number);
            corrected_mass = GetCorrectedJetMass(raw_mass, raw_pt, eta, areas[i], rho, run_number);
        }

        // Snapshot pt before JER so the mass rescale below only folds in the
        // JER ratio, not the JEC factor a second time.
        double pt_before_jer = corrected_pt;

        if (!isData && applyJER) {
            // -1.0 = no gen match; SmearJER() falls back to stochastic smearing.
            float matched_genpt = -1.0;
            float matched_geneta = 0.0;
            float matched_genphi = 0.0;
            if (genJetIndices.size() > i && genJetIndices[i] >= 0 &&
                static_cast<size_t>(genJetIndices[i]) < genJets.size()) {
                matched_genpt  = genJets[genJetIndices[i]].Pt();
                matched_geneta = genJets[genJetIndices[i]].Eta();
                matched_genphi = genJets[genJetIndices[i]].Phi();
            }
            corrected_pt = SmearJER(corrected_pt, matched_genpt, matched_geneta, matched_genphi,
                                     eta, phi, rho, event_number, jerSysTag);
        }

        if (pt_before_jer > 0) {
            corrected_mass *= (corrected_pt / pt_before_jer);
        }

        TLorentzVector corr_jet;
        corr_jet.SetPtEtaPhiM(corrected_pt, eta, phi, corrected_mass);
        correctedJets.push_back(corr_jet);
    }

    return correctedJets;
}

TLorentzVector SSBCorrections::ApplyType1METWithCorrT1(
    double raw_met_pt,
    double raw_met_phi,
    const std::vector<TLorentzVector>& jetsAsStored,
    const std::vector<float>& jetRawFactors,
    const std::vector<float>& jetAreas,
    const std::vector<float>& jetMuonSubtrFactors,
    const std::vector<float>& jetChEmEF,
    const std::vector<float>& jetNeEmEF,
    const std::vector<int>& jetGenJetIndices,
    const std::vector<float>& corrT1RawPt,
    const std::vector<float>& corrT1Eta,
    const std::vector<float>& corrT1Phi,
    const std::vector<float>& corrT1Area,
    const std::vector<float>& corrT1MuonSubtrFactor,
    float rho,
    bool isData,
    bool applyJES,
    bool applyJER,
    const std::vector<TLorentzVector>& genJets,
    unsigned int run_number,
    ULong64_t event_number,
    const std::string& jerSysTag
) const {
    double met_px = raw_met_pt * std::cos(raw_met_phi);
    double met_py = raw_met_pt * std::sin(raw_met_phi);

    // Per-jet: muon-subtract raw pt, apply JES/JER, require Type-1 selection
    // (pt>15, |eta|<5.2, chEmEF+neEmEF<0.90), accumulate delta into MET.
    auto accumulate = [&](double eta, double phi, double rawPtRaw, double area,
                          double muonSubtrFactor, double chEmEF, double neEmEF,
                          float matchedGenPt, float matchedGenEta, float matchedGenPhi) {
        double rawPtNoMu = rawPtRaw * (1.0 - muonSubtrFactor);
        if (rawPtNoMu <= 0.0) return;

        double l1PtNoMu = rawPtNoMu;
        double corrPtNoMu = rawPtNoMu;
        if (applyJES) {
            l1PtNoMu   = GetL1CorrectedJetPt(rawPtNoMu, eta, area, rho);
            corrPtNoMu = GetCorrectedJetPt(rawPtNoMu, eta, area, rho, run_number);
        }
        if (!isData && applyJER) {
            corrPtNoMu = SmearJER(corrPtNoMu, matchedGenPt, matchedGenEta, matchedGenPhi,
                                   eta, phi, rho, event_number, jerSysTag);
        }

        // Type-1 selection - only jets passing this contribute to MET.
        bool passSel = (corrPtNoMu > 15.0) && (std::abs(eta) < 5.2) && ((chEmEF + neEmEF) < 0.90);
        if (!passSel) return;

        met_px += (l1PtNoMu - corrPtNoMu) * std::cos(phi);
        met_py += (l1PtNoMu - corrPtNoMu) * std::sin(phi);
    };

    // jetsAsStored holds NanoAOD's already-corrected pt, so undo rawFactor
    // first. Gen match reuses Jet_genJetIdx (same as ApplyJetCorrections) so
    // this path can't disagree with the physics-jet path's match.
    for (size_t i = 0; i < jetsAsStored.size(); ++i) {
        double rawFactor = (i < jetRawFactors.size()) ? jetRawFactors[i] : 0.0f;
        double area       = (i < jetAreas.size())       ? jetAreas[i]       : 0.5f;
        double muonSubtr  = (i < jetMuonSubtrFactors.size()) ? jetMuonSubtrFactors[i] : 0.0f;
        double chEmEF     = (i < jetChEmEF.size())      ? jetChEmEF[i]      : 0.0f;
        double neEmEF     = (i < jetNeEmEF.size())      ? jetNeEmEF[i]      : 0.0f;
        double rawPtRaw   = jetsAsStored[i].Pt() * (1.0 - rawFactor);

        float matchedGenPt = -1.0f;
        float matchedGenEta = 0.0f;
        float matchedGenPhi = 0.0f;
        if (!isData && applyJER && i < jetGenJetIndices.size() &&
            jetGenJetIndices[i] >= 0 &&
            static_cast<size_t>(jetGenJetIndices[i]) < genJets.size()) {
            matchedGenPt  = genJets[jetGenJetIndices[i]].Pt();
            matchedGenEta = genJets[jetGenJetIndices[i]].Eta();
            matchedGenPhi = genJets[jetGenJetIndices[i]].Phi();
        }
        accumulate(jetsAsStored[i].Eta(), jetsAsStored[i].Phi(), rawPtRaw, area, muonSubtr, chEmEF, neEmEF,
                   matchedGenPt, matchedGenEta, matchedGenPhi);
    }

    // CorrT1METJet_ jets: rawPt is already raw, no EM-fraction branches exist
    // (pass 0/0, per the tutorial), and there's no genJetIdx equivalent so
    // this falls back to the DeltaR-based MatchGenPt() rematch.
    for (size_t i = 0; i < corrT1RawPt.size(); ++i) {
        double area      = (i < corrT1Area.size())               ? corrT1Area[i]               : 0.5f;
        double muonSubtr = (i < corrT1MuonSubtrFactor.size())    ? corrT1MuonSubtrFactor[i]    : 0.0f;

        float matchedGenPt = -1.0f;
        float matchedGenEta = 0.0f;
        float matchedGenPhi = 0.0f;
        if (!isData && applyJER) {
            TLorentzVector recoJet4v;
            recoJet4v.SetPtEtaPhiM(corrT1RawPt[i], corrT1Eta[i], corrT1Phi[i], 0.0);
            matchedGenPt = MatchGenPt(recoJet4v, genJets, 0.2f, &matchedGenEta, &matchedGenPhi);
        }
        accumulate(corrT1Eta[i], corrT1Phi[i], corrT1RawPt[i], area, muonSubtr, 0.0, 0.0,
                   matchedGenPt, matchedGenEta, matchedGenPhi);
    }

    double corrected_met_pt  = std::sqrt(met_px * met_px + met_py * met_py);
    double corrected_met_phi = std::atan2(met_py, met_px);

    TLorentzVector correctedMET;
    correctedMET.SetPtEtaPhiM(corrected_met_pt, 0, corrected_met_phi, 0);
    return correctedMET;
}


bool SSBCorrections::ShouldVetoJet(const TLorentzVector& jet, double chEmEF, double neEmEF) const {
    // Only apply for 2018 data/MC
    if (year_ != "2018") {
        return false;
    }

    // Check if jetvetomap is loaded
    if (!jetvetomap_) {
        logger_.Warning() << "Jet veto map not loaded, skipping HEM veto" << std::endl;
        return false;
    }

    // pt/jetId pre-selection already applied by the only caller (JetSelector).
    if ((chEmEF + neEmEF) >= 0.90) {
        return false;
    }

    float eta = jet.Eta();
    float phi = jet.Phi();

    try {
        // Use configured key
        std::variant<double, std::vector<double>> val = jetvetomap_->evaluate({
            jveto_key_, eta, phi
        });
        
        double veto_flag = std::get<double>(val);

        // Return true if jet should be vetoed (non-zero value)
        bool should_veto = (veto_flag > 0.0);

        return should_veto;

    } catch (const std::exception& e) {
        logger_.Warning() << "Jet veto map evaluation failed for jet (eta="
                  << eta << ", phi=" << phi << ") with key=" << jveto_key_
                  << ": " << e.what() << std::endl;
        return false;
    }
}


void SSBCorrections::InitBtagSFCorrection(const std::string& json_path, 
                                          const std::string& tagger_name) {
    logger_.Info() << "Loading b-tagging SF from JSON: " << json_path << std::endl;
    
    auto cset = correction::CorrectionSet::from_file(json_path);

    // UParT ships only one heavy-flavor SF ("UParTAK4_comb", no "_mujets"
    // variant) and names light-flavor "UParTAK4_light" instead of "_incl".
    bool isUParT = (tagger_name.rfind("UParTAK4", 0) == 0);

    std::string heavy_flavor_name = tagger_name;  // e.g., "deepJet_comb" or "UParTAK4_comb"
    std::string light_flavor_name = tagger_name;

    if (isUParT) {
        // heavy_flavor_name is already the exact correction name as-is
        // ("UParTAK4_comb"); only the light-flavor suffix needs swapping.
        size_t pos = light_flavor_name.find("comb");
        if (pos != std::string::npos) {
            light_flavor_name.replace(pos, 4, "light");
        } else {
            logger_.Error() << "Expected 'comb' in tagger_name: " << tagger_name << std::endl;
            return;
        }
    } else {
        // Heavy flavor correction name (config-based selection)
        size_t pos = heavy_flavor_name.find("comb");
        if (pos != std::string::npos) {
            heavy_flavor_name.replace(pos, 4, btag_sf_type_);
        } else {
            logger_.Error() << "Expected 'comb' in tagger_name: " << tagger_name << std::endl;
            return;
        }

        // Light flavor correction name (always "incl")
        pos = light_flavor_name.find("comb");
        if (pos != std::string::npos) {
            light_flavor_name.replace(pos, 4, "incl");
        }
    }

    // Load corrections with generic keys
    btag_corrections_["heavy"] = cset->at(heavy_flavor_name);
    btag_corrections_["light"] = cset->at(light_flavor_name);

    logger_.Info() << "Loaded corrections: "
              << heavy_flavor_name << " and " << light_flavor_name << std::endl;

    // Lets Jet_btag="UParTM" alone determine the cut via GetBtagWPCut(),
    // instead of also requiring an explicit BTagDiscCut. Optional - falls
    // back to the config value if this fails to load.
    if (isUParT) {
        std::string wpValuesName = tagger_name.substr(0, tagger_name.find('_')) + "_wp_values";
        btag_wp_values_ = LoadOptionalCorrection<correction::Correction>(
            [&]() { return cset->at(wpValuesName); },
            "Could not load '" + wpValuesName +
            "'. BTagDiscCut must be set explicitly in the config for this tagger.", logger_);
        if (btag_wp_values_) {
            logger_.Info() << "Loaded WP-cut lookup correction ('" << wpValuesName
                      << "') - Jet_btag's working point letter alone can now determine "
                      << "BTagDiscCut if it's not set explicitly in the config." << std::endl;
        }
    }
}

double SSBCorrections::GetBtagWPCut(const std::string& wp) const {
    if (!btag_wp_values_) {
        return -1.0;
    }
    try {
        // Single input: the working point string, e.g. "L"/"M"/"T" - the
        // same convention already confirmed working for GetBtagSF()'s own
        // evaluate() call against this same btagging.json.gz.
        return btag_wp_values_->evaluate({wp});
    } catch (const std::exception& e) {
        logger_.Warning() << "GetBtagWPCut('" << wp << "') failed: " << e.what()
                  << ". If this is an input-count/type mismatch, the wp_values correction's "
                  << "actual schema differs from the single-string-input assumption here - "
                  << "check with `python3 -c \"import correctionlib; "
                  << "c=correctionlib.CorrectionSet.from_file('btagging.json.gz'); "
                  << "print(c['UParTAK4_wp_values'].inputs)\"`." << std::endl;
        return -1.0;
    }
}

std::string SSBCorrections::getBtagCorrectionName(int flavor) const {
    // 0 = light flavor (u,d,s,g), others = heavy flavor (b,c)
    return (flavor == 0) ? "light" : "heavy";
}

float SSBCorrections::GetBtagSF(float pt, float eta, int flav, 
                                const std::string& wp, const std::string& syst) const {
    float pt_clamped = std::clamp(pt, 20.0f, 1000.0f);
    float eta_abs = std::fabs(eta);

    std::string corr_name = getBtagCorrectionName(flav);  // "heavy" or "light"
    auto it = btag_corrections_.find(corr_name);
    
    if (it == btag_corrections_.end()) {
        logger_.Error() << "B-tag correction not found: " << corr_name << std::endl;
        return 1.0;
    }

    try {
        return it->second->evaluate({syst, wp, flav, eta_abs, pt_clamped});
    } catch (const std::exception& e) {
        logger_.Warning() << "GetBtagSF failed: " << e.what() << std::endl;
        return 1.0;
    }
}

float SSBCorrections::ComputeBTagEventWeight(const std::vector<float>& pts,
                                             const std::vector<float>& etas,
                                             const std::vector<int>& flavs,
                                             const std::vector<bool>& isTagged,
                                             const std::string& algo,
                                             const std::string& wp,
                                             const std::string& syst) const {
    
    // Check input vector sizes
    if (pts.size() != etas.size() || 
        pts.size() != flavs.size() || 
        pts.size() != isTagged.size()) {
        logger_.Warning() << "ComputeBTagEventWeight: input vectors have different sizes" << std::endl;
        return 1.0;
    }
    // Initialize probabilities
    float p_mc = 1.0;
    float p_data = 1.0;
    
    // Loop over all jets
    for (size_t i = 0; i < pts.size(); ++i) {
        float pt   = pts[i];
        float eta  = etas[i];
        int flav   = flavs[i];
        bool tagged = isTagged[i];
         
        // Get MC efficiency
        float eff = GetMCBtagEfficiency(pt, eta, flav, algo, wp);
        
        // Get scale factor
        float sf = GetBtagSF(pt, eta, flav, wp, syst);
        
        // Avoid division by zero
        if (eff < 1e-5) eff = 1e-5;
        if (eff > 1.0 - 1e-5) eff = 1.0 - 1e-5;
        
        // Calculate probabilities
        if (tagged) {
            // Tagged jet contribution
            p_mc *= eff;
            p_data *= (eff * sf);
        } else {
            // Not tagged jet contribution
            p_mc *= (1.0f - eff);
            p_data *= (1.0f - eff * sf);
        }
    }
    
    // Calculate final event weight
    float btag_evt_weight = (p_mc > 1e-10) ? p_data / p_mc : 1.0;
    
    return btag_evt_weight;
}


void SSBCorrections::LoadMCBtagEfficiencies(const std::string& filepath,
                                             const std::string& algo,
                                             const std::string& wp) {
    std::unique_ptr<TFile> f(TFile::Open(filepath.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        logger_.Error() << "Failed to open efficiency file: " << filepath << std::endl;
        return;
    }

    for (const std::string& flav : {"b", "c", "l"}) {
        std::string name = "eff_" + algo + "_" + flav + "_" + wp;
        auto* hist = dynamic_cast<TH2*>(f->Get(name.c_str()));
        if (hist) {
            std::unique_ptr<TH2> hist_copy(static_cast<TH2*>(hist->Clone((name + "_copy").c_str())));
            hist_copy->SetDirectory(0);
            // Assigning to the map slot destroys any previous unique_ptr there automatically.
            eff_histograms_[algo + "_" + flav + "_" + wp] = std::move(hist_copy);
            logger_.Debug() << "Loaded and copied hist: " << name << std::endl;
        } else {
            logger_.Warning() << "Histogram not found: " << name << std::endl;
        }
    }
    f->Close();
}

float SSBCorrections::GetMCBtagEfficiency(float pt, float eta, int flav, const std::string& algo, const std::string& wp) const {
    std::string flav_str = "l";
    if (flav == 5) flav_str = "b";
    else if (flav == 4) flav_str = "c";

    // Convert single character WP to full name for efficiency lookup
    std::string wp_full = wp;
    if (wp == "L") wp_full = "Loose";
    else if (wp == "M") wp_full = "Medium";
    else if (wp == "T") wp_full = "Tight";

    std::string key = algo + "_" + flav_str + "_" + wp_full;
    auto it = eff_histograms_.find(key);
    if (it == eff_histograms_.end()) {
        logger_.Warning() << "Efficiency hist not found: " << key << std::endl;
        return 1.0;
    }

    TH2* hist = it->second.get();
    int bin_x = hist->GetXaxis()->FindBin(pt);
    int bin_y = hist->GetYaxis()->FindBin(eta);
    float eff = hist->GetBinContent(bin_x, bin_y);

    return std::clamp(eff, 0.0f, 1.0f);
}

double SSBCorrections::RochesterCorrectionData(TString year, int Q, double pt, double eta, double phi, int s,int m) const{
        double correction;
        correction = rc.kScaleDT(Q,pt,eta,phi,s,m); return correction;
}

double SSBCorrections::RochesterCorrectionMC(TString year, int Q, double pt, double eta,double phi,int genID,double genPt,int nl, int s,int m) const{

    double correction = 1.0; double u =1.0;

    bool genMatch = genID != -1;

    if(genMatch) correction = rc.kSpreadMC(Q,pt,eta,phi,genPt,s,m);
    else if(!genMatch){u = gRandom->Rndm(); correction = rc.kSmearMC(Q,pt,eta,phi,nl,u,s,m);} //Random number is needed when gen-mathcing is failed
    return correction;
}


TLorentzVector SSBCorrections::METXYCorrection(const TLorentzVector& type1_met,
                                               int runnb, TString year, bool isMC, int npv, bool isUL, bool ispuppi
                                               ) const {


    std::pair<double, double> correctedMET = METXYCorr_Met_MetPhi(type1_met.Pt(),type1_met.Phi(),runnb,year,isMC,npv,isUL,ispuppi);

    double met_pt  = correctedMET.first;
    double met_phi = correctedMET.second;

    TLorentzVector corrected_met;
    corrected_met.SetPtEtaPhiM(met_pt,0,met_phi,0);
    return corrected_met;
}


TLorentzVector SSBCorrections::METXYCorrection_corrlib(const TLorentzVector& type1_met,
                                               const std::string& era,
                                               bool isData,
                                               int npv) const {
    if (!metphi_corr_) {
        logger_.Error() << "[METXYCorrection] MET correction object not loaded." << std::endl;
        return type1_met;
    }

    std::string data_tag = isData ? "data" : "mc";

    double corr_x = 0.0;
    double corr_y = 0.0;

    try {
        std::variant<double, std::vector<double>> val_x = metphi_corr_->evaluate({era, data_tag, npv, "x"});
        std::variant<double, std::vector<double>> val_y = metphi_corr_->evaluate({era, data_tag, npv, "y"});

        if (std::holds_alternative<double>(val_x)) {
            corr_x = std::get<double>(val_x);
        } else {
            logger_.Warning() << "[METXYCorrection] unexpected type for x correction" << std::endl;
        }

        if (std::holds_alternative<double>(val_y)) {
            corr_y = std::get<double>(val_y);
        } else {
            logger_.Warning() << "[METXYCorrection] unexpected type for y correction" << std::endl;
        }
    } catch (const std::exception& e) {
        logger_.Error() << "[METXYCorrection] Evaluation failed: " << e.what() << std::endl;
        return type1_met;
    }

    double met_x = type1_met.Px() - corr_x;
    double met_y = type1_met.Py() - corr_y;

    TLorentzVector corrected_met;
    corrected_met.SetPxPyPzE(met_x, met_y, 0, std::sqrt(met_x * met_x + met_y * met_y));
    return corrected_met;
}
/// Trigger SF 
double SSBCorrections::TrigDiMuon_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_) {
    double pt1 = lep1.Pt();
    double pt2 = lep2.Pt();

    double leading_pt = std::max(pt1, pt2);
    double subleading_pt = std::min(pt1, pt2);

    return GetTrgEff(leading_pt, subleading_pt, Sys_);
}

double SSBCorrections::TrigDiElec_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_) {
    double pt1 = lep1.Pt();
    double pt2 = lep2.Pt();

    double leading_pt = std::max(pt1, pt2);
    double subleading_pt = std::min(pt1, pt2);

    return GetTrgEff(leading_pt, subleading_pt, Sys_);
}

double SSBCorrections::TrigMuElec_Eff(TLorentzVector muon, TLorentzVector elec, TString Sys_) {
    double pt_mu = muon.Pt();   // Assuming lep1 is muon
    double pt_ele = elec.Pt();  // Assuming lep2 is electron

    return GetTrgEff(pt_ele, pt_mu, Sys_);
}

double SSBCorrections::GetTrgEff(double pt1, double pt2, TString Sys_) {
    // Clamp values below the histogram max range (assumed 500)
    double max_pt = 499.999;

    pt1 = std::min(pt1, max_pt);
    pt2 = std::min(pt2, max_pt);

    int xbin = H_trig->GetXaxis()->FindBin(pt1);
    int ybin = H_trig->GetYaxis()->FindBin(pt2);

    double trgsf = H_trig->GetBinContent(xbin, ybin);
    double trgsferr = 0.0;

    if (Sys_ == "nominal" || Sys_ == "Central") {trgsferr = 0.0;}
    else if (Sys_ == "Up" || Sys_ == "up")
        trgsferr = H_trig->GetBinError(xbin, ybin);
    else if (Sys_ == "Down" || Sys_ == "down")
        trgsferr = -H_trig->GetBinError(xbin, ybin);

    return trgsf + trgsferr;
}

std::string SSBCorrections::ExpandJECName(const std::string& base_jec_name, const std::string runPeriod, const std::string& era, bool is_data,
                                           const std::string& jetAlgo, const std::string& jecLevel) {
    //  "Summer19UL16_v7" -> prefix = "Summer19UL16", version = "7"
    size_t pos = base_jec_name.find("_V");
    if (pos == std::string::npos) {
        logger_.Error() << "Invalid jec_name format: " << base_jec_name << std::endl;
        return base_jec_name; // fallback
    }

    std::string prefix = base_jec_name.substr(0, pos);             // "Summer19UL16"
    std::string version = base_jec_name.substr(pos + 2);           // "7"
    // Config-driven (JetAlgoTag/JECLevel), defaulting to Puppi - v15's Jet_
    // collection is AK4 Puppi, not the old PFchs this was hardcoded to.
    std::string suffix = "_" + jecLevel + "_" + jetAlgo;

    std::string expanded = prefix;  // prefix

    // v15 Puppi compounds use a single name per year (era resolved inside the
    // correction via the "run" input - see jec_needs_run_), so skip the
    // per-era name-building chain below; only PFchs (v9) still needs it.
    bool useEraInName = (jetAlgo.find("Puppi") == std::string::npos);

    if (useEraInName) {
    if (runPeriod.find("16Pre") != std::string::npos) {
        if (is_data) {
            if (era == "B" || era == "C" || era == "D") {
                expanded += "_RunBCD";
            } else if (era == "E") {
                expanded += "_RunEF";
            } else if (era.find("F") == 0) {
                expanded += "_RunEF";
            } else if (era == "G" || era == "H") {
                expanded += "_RunGH";
            } else {
                expanded += "_RunFGH";  // fallback
            }
        } else {
            expanded += (era == "B" || era == "C" || era == "D" || era == "E") ? "APV" : "";
        }
    } else if (runPeriod.find("16Post") != std::string::npos){
        if (is_data) {
	    if (era.find("F") == 0) {
                expanded += "_RunFGH";
            } else if (era == "G" || era == "H") {
                expanded += "_RunFGH";
            } else {
                expanded += "_RunFGH";  // fallback
            }
        } else {
            expanded += (era == "F" || era == "G" || era == "H") ? "APV" : "";
        }

    } else if (runPeriod.find("2017") != std::string::npos) {
        if (is_data) {
            if (era == "B" ) {
                expanded += "_RunB";
            } else if (era == "C") {
                expanded += "_RunC";
            } else if (era == "D" ) {
                expanded += "_RunD";
            } else if (era == "E") {
                expanded += "_RunE";
            } else if (era == "F") {
                expanded += "_RunF";
            } else {
                expanded += "_RunB";
            }

        }
    } else if (runPeriod.find("2018") != std::string::npos) {
        if (is_data) {
            if (era == "A" ) {
                expanded += "_RunA";
            } else if (era == "B") {
                expanded += "_RunB";
            } else if (era == "C" ) {
                expanded += "_RunC";
            } else if (era == "D") {
                expanded += "_RunD";
            } else {
                expanded += "_RunB";
            }
        }
    }
    } // end useEraInName

    // Add version and data/MC tag
    expanded += "_V" + version + "_";
    expanded += is_data ? "DATA" : "MC";
    expanded += suffix;

    return expanded;
}
std::string SSBCorrections::GetJetVetoType() const {
	return jveto_type_; 
}

std::string SSBCorrections::GetProcessSubDir(const std::string& inputfileName) const {
    // DY madgraph (must check before plain DY)
    if (inputfileName.find("DYJetsToLL") != std::string::npos) {
        if (inputfileName.find("madgraph") != std::string::npos) {
            if (year_ == "2018") return "DY_madgraph";
            else                 return "DY";
        }
        return "DY";
    }
    // WJets -> use DY efficiency
    if (inputfileName.find("WJetsToLNu") != std::string::npos) return "DY";
    // TTbar
    if (inputfileName.find("TTbar_Signal") != std::string::npos) return "TTbar_Signal";
    if (inputfileName.find("TTbar_")       != std::string::npos) return "TTbarOther";
    if (inputfileName.find("TTJets_")      != std::string::npos) return "TTbarOther";
    // SingleTop
    if (inputfileName.find("ST_") != std::string::npos) return "SingleTop";
    // Diboson
    if (inputfileName.find("WW") != std::string::npos ||
        inputfileName.find("WZ") != std::string::npos ||
        inputfileName.find("ZZ") != std::string::npos) return "Diboson";
    // TTV
    if (inputfileName.find("TTW") != std::string::npos ||
        inputfileName.find("TTZ") != std::string::npos) return "TTV";

    logger_.Warning() << "GetProcessSubDir: no match for '" << inputfileName
              << "', falling back to TTbar_Signal" << std::endl;
    return "TTbar_Signal";
}
