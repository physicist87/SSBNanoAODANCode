#include "../interface/Analysis.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
// Reads a True/False config key with a validated default. A typo (e.g.
// "Ture") throws instead of silently being treated as "off" by a later
// `== "True"` check - unlike the old "== \"DUMMY\"" sentinel pattern, which
// only caught a fully-missing key, not a misspelled value.
TString ReadValidatedBoolConfig(TextReader* reader, const std::string& key,
                                 const char* defaultText, Logger& logger) {
    std::string raw = reader->Check(key) ? reader->GetText(key) : "DUMMY";
    if (raw == "DUMMY") {
        logger.Warning() << key << " not set in config - defaulting to " << defaultText << "." << std::endl;
        return TString(defaultText);
    }
    if (raw == "True" || raw == "False") return TString(raw);
    throw std::runtime_error("Invalid config value for '" + key + "': \"" + raw +
                              "\" (expected \"True\" or \"False\").");
}

// Same idea as ReadValidatedBoolConfig, for a fixed set of string choices
// (e.g. METtype: "PF"/"Puppi") instead of True/False.
TString ReadValidatedChoiceConfig(TextReader* reader, const std::string& key,
                                   const std::vector<std::string>& allowed,
                                   const char* defaultText, Logger& logger) {
    std::string raw = reader->Check(key) ? reader->GetText(key) : "DUMMY";
    if (raw == "DUMMY") {
        logger.Warning() << key << " not set in config - defaulting to " << defaultText << "." << std::endl;
        return TString(defaultText);
    }
    for (const auto& a : allowed) if (raw == a) return TString(raw);
    std::string allowedList;
    for (size_t i = 0; i < allowed.size(); ++i) {
        allowedList += "\"" + allowed[i] + "\"";
        if (i + 1 < allowed.size()) allowedList += "/";
    }
    throw std::runtime_error("Invalid config value for '" + key + "': \"" + raw +
                              "\" (expected one of " + allowedList + ").");
}
}  // namespace

// Constructor: initialize TTreeReader with TChain and branch list file
Analysis::Analysis(TChain *inputChain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt= -1)
    : chain(inputChain), fReader(inputChain), branchReader_(inputChain, fReader), NumEvt(NumEvt), outdir(seDirName), outfile(outputName){
    if (!chain) {
        throw std::runtime_error("Error: Invalid TChain pointer!");
    }
    FileName_ = SetInputFileName(inputName);
    isData = TString(FileName_).Contains("Data");
    logger_.Info() << "FileName_ : " << FileName_ << std::endl;
    logger_.Info() << "configFile : " << configFile << std::endl;
    std::string confDir = "./configs/";
    std::string confpath = "";
    confpath = confDir+configFile;
    SSBConfReader = std::make_unique<TextReader>();
    SSBConfReader->ReadFile(confpath);
    SSBConfReader->ReadVariables();
    SSBConfReader->PrintoutVariables();
    SSBCorr = std::make_unique<SSBCorrections>(SSBConfReader.get(), FileName_.Data());
    SSBCPVCal = std::make_unique<SSBCPVCalc>();
    // Initialize branches based on branch list file (NanoAODBranchReader owns
    // the actual maps/type-detection now - see interface/NanoAODBranchReader.h)
    branchReader_.InitBranches(branchListFile, isData);

    // Fail fast: check the branches this job actually needs (channel/METtype/
    // data-vs-MC dependent) right now, instead of discovering a NanoAOD
    // version/branch-list mismatch mid-event-loop via a MissingBranchError
    // or a wrong-looking default value.
    {
        std::string channel = SSBConfReader->GetText("Channel");
        std::string metType = SSBConfReader->Check("METtype") ? SSBConfReader->GetText("METtype") : "Puppi";
        std::vector<std::string> required = {
            "run", "luminosityBlock", "event", "PV_npvsGood",
            "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_rawFactor", "Jet_area",
        };
        if (metType == "PF") {
            required.push_back("RawMET_pt");
            required.push_back("RawMET_phi");
        } else {
            required.push_back("RawPuppiMET_pt");
            required.push_back("RawPuppiMET_phi");
        }
        if (channel.find("muon") != std::string::npos || channel.find("muel") != std::string::npos) {
            required.insert(required.end(), {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge"});
        }
        if (channel.find("elec") != std::string::npos || channel.find("muel") != std::string::npos) {
            required.insert(required.end(), {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_mass", "Electron_charge"});
        }
        if (!isData) {
            required.push_back("Generator_weight");
            required.push_back("Pileup_nTrueInt");
        }
        branchReader_.ValidateRequiredBranches(required);
    }

    cutflowName[selectionIndex(SelectionStage::AfterTriggerAndPV)] = "Step_0";
    cutflowName[selectionIndex(SelectionStage::AfterDileptonCuts)] = "Step_1" ;
    cutflowName[selectionIndex(SelectionStage::AfterZVeto)] = "Step_2";
    cutflowName[selectionIndex(SelectionStage::AfterJetMultiplicity)] = "Step_3";
    cutflowName[selectionIndex(SelectionStage::AfterMETCut)] = "Step_4";
    cutflowName[selectionIndex(SelectionStage::AfterBTagMultiplicity)] = "bTagged Jet >= 1";
    cutflowName[selectionIndex(SelectionStage::Reserved6)] = "bTagged Jet >= 2";
    cutflowName[selectionIndex(SelectionStage::Reserved7)] = "bTagged Jet == 2";
    cutflowName[selectionIndex(SelectionStage::AfterTopReconstruction)] = "Top-Recon.";
    cutflowName[selectionIndex(SelectionStage::Reserved9)] = "Top-Pt-Rewight";
 
    pi = TMath::Pi();
    Start();
}

Analysis::~Analysis() {
    // fout needs Write()/Close() before release; the others would free
    // themselves anyway, but log lines here confirm teardown order.
    if (fout) {
        fout->Write();
        logger_.Info() << "output file Name : " <<  fout->GetName() << std::endl;
        fout->Close(); // Close the file before releasing
        fout.reset();
        logger_.Info() << "fout successfully deleted." << std::endl;
    }

    if (SSBConfReader) {
        SSBConfReader.reset();
        logger_.Info() << "SSBConfReader successfully deleted." << std::endl;
    }

    if (SSBCorr) {
        SSBCorr.reset();
        logger_.Info() << "SSBCorr successfully deleted." << std::endl;
    }
    if (SSBCPVCal) {
        SSBCPVCal.reset();
        logger_.Info() << "SSBCPVCal successfully deleted." << std::endl;
    }
    logger_.Info() << "Analysis destructor completed." << std::endl;
}

// Version-agnostic NanoAOD branch reading lives in NanoAODBranchReader -
// Analysis calls through branchReader_ instead of touching the maps directly.

// NanoAODv15 has no Jet_jetId; PUPPI jet ID is reconstructed from energy fractions.
bool Analysis::PassConfiguredJetId(int idx) const {
    JetID jetId(branchReader_, RunPeriod);
    return jetId.PassConfigured(idx, JetId);
}



void Analysis::SetVariables() {
    logger_.Info() << "Set varibles " << std::endl;

    Lumi = SSBConfReader->GetNumber( "Luminosity" );
    RunPeriod = SSBConfReader->GetText( "RunRange" );
    Decaymode = SSBConfReader->GetText( "Channel" ); // Channel
    XsecTable_ = SSBConfReader->GetText( "XSecTablesName" );

    std::string blindStr = SSBConfReader->GetText( "isBlind" );
    isBlind = (blindStr == "True" || blindStr == "true");
    logger_.Info() << "[Blind] isData=" << isData << " isBlind=" << isBlind << std::endl;

    // Debug event tracer - optional, off by default; Check()-guarded since
    // DebugMode/DebugRun/DebugLumi/DebugEvent aren't in any existing config.
    {
        bool debugMode = SSBConfReader->Check("DebugMode") && SSBConfReader->GetBool("DebugMode");
        unsigned int debugRun = SSBConfReader->Check("DebugRun") ? SSBConfReader->GetNumberUint("DebugRun") : 0;
        unsigned int debugLumi = SSBConfReader->Check("DebugLumi") ? SSBConfReader->GetNumberUint("DebugLumi") : 0;
        unsigned long long debugEvent = SSBConfReader->Check("DebugEvent")
            ? static_cast<unsigned long long>(SSBConfReader->GetNumberDouble("DebugEvent")) : 0;
        debugFilter_.Configure(debugMode, debugRun, debugLumi, debugEvent);
        if (debugMode) {
            logger_.Info() << "[DebugEventFilter] enabled - tracing run=" << debugRun
                      << " lumi=" << debugLumi << " event=" << debugEvent
                      << " (0 = match any value for that field)" << std::endl;
        }
    }


    num_dleptrig = SSBConfReader->Size( "dileptrigger" );
    num_sleptrig = SSBConfReader->Size( "singleleptrigger" );

    for(int i =0; i < num_dleptrig; ++i)
    {
        logger_.Info() << SSBConfReader->GetText("dileptrigger",i+1) << std::endl;
        std::string tmptrg = SSBConfReader->GetText("dileptrigger",i+1);
        DLtrigName.push_back( removeSubstring( tmptrg, "_v") );
        trigName.push_back( removeSubstring( tmptrg, "_v")  );
    }

    for(int i =0; i < num_sleptrig; ++i)
    {
        logger_.Info() << SSBConfReader->GetText("singleleptrigger",i+1) << std::endl;
        std::string tmptrg = SSBConfReader->GetText("singleleptrigger",i+1);
        SLtrigName.push_back( removeSubstring(tmptrg,"_v") );
        trigName.push_back( removeSubstring(tmptrg,"_v") );
    }

    for (int i = 0; i < trigName.size(); ++i){
        logger_.Info() << "test trigName[i] " << trigName[i] << std::endl;
        triggerList[trigName[i]] =DeepCopy<bool>( branchReader_.boolSingles[trigName[i]]);
    }

    for (int i = 0; i < SSBConfReader->Size("METFilters"); ++i){
       std::string tmpnoisefl = SSBConfReader->GetText("METFilters",i+1);
       if(branchReader_.boolSingles[tmpnoisefl] ==NULL) {logger_.Error() << "Error!!!" << tmpnoisefl << std::endl;}
       noiseFilters[tmpnoisefl] = DeepCopy<bool>( branchReader_.boolSingles[tmpnoisefl]);
    }
    // Kinematic cut variables for Object 
    muon_pt     = SSBConfReader->GetNumber( "MuonPt_cut"     );// Muon pT Cut 
    muon_eta    = SSBConfReader->GetNumber( "MuonEta_cut"    );// Muon Eta Cut 
    muon_isocut  = SSBConfReader->GetNumber( "MuonIso_cut" );// 

    elec_pt     = SSBConfReader->GetNumber( "ElecPt_cut"     );
    elec_eta    = SSBConfReader->GetNumber( "ElecEta_cut"    );
    elec_isocut  = SSBConfReader->GetNumber( "ElecIso_cut" );


    MuonIsoType = SSBConfReader->GetText( "MuonIso_type" );
    MuonId      = SSBConfReader->GetText( "Muon_ID"      );

    ElecIsoType = SSBConfReader->GetText( "ElecIso_type" );
    ElecId      = SSBConfReader->GetText( "Elec_ID"      );

    JetId      = SSBConfReader->GetText("Jet_ID");
    JetbTag    = SSBConfReader->GetText("Jet_btag");
 
    
    // Kinematic variables for Object 
    jet_pt     = SSBConfReader->GetNumber( "Jet_pt"  );
    jet_eta    = SSBConfReader->GetNumber( "Jet_eta" );
    met_cut    = SSBConfReader->GetNumber( "MET_cut" );
 
    veto_muoniso_type = SSBConfReader->GetText( "VetoMuonIso" );
    veto_muoniso_cut  = SSBConfReader->GetNumber( "VetoMuonIsocut" );
    veto_muonid       = SSBConfReader->GetText( "VetoMuonId"  );
    veto_eleciso_type = SSBConfReader->GetText( "VetoElecIso" );
    veto_elecid       = SSBConfReader->GetText( "VetoElecId"  );

    // DoJES is a new key; default true (JES's old always-on behavior) so an
    // unmigrated config doesn't silently get JES=off from GetBool()'s default.
    dojes = SSBConfReader->Check("DoJES") ? SSBConfReader->GetBool("DoJES") : true;
    if (!SSBConfReader->Check("DoJES")) {
        logger_.Warning() << "DoJES not set in config - defaulting to true (JES applied), "
                              "matching this code's previous always-on behavior. "
                              "Add DoJES : \"true\"/\"false\" explicitly to silence this." << std::endl;
    }
    dojer = SSBConfReader->Check("DoJER") ? SSBConfReader->GetBool("DoJER") : true;
    if (!SSBConfReader->Check("DoJER")) {
        logger_.Warning() << "DoJER not set in config - defaulting to true (JER applied), "
                              "matching this code's previous always-on behavior. "
                              "Add DoJER : \"true\"/\"false\" explicitly to silence this." << std::endl;
    }
    logger_.Info() << "DoJES: " << (dojes ? "true" : "false")
                    << "  DoJER: " << (dojer ? "true" : "false") << std::endl;

    // New key; default "nominal" (previous behavior), same Check()-guard pattern as DoJES.
    JERSys = SSBConfReader->Check("JERSys") ? SSBConfReader->GetText("JERSys") : "nominal";
    if (!SSBConfReader->Check("JERSys")) {
        logger_.Warning() << "JERSys not set in config - defaulting to \"nominal\". "
                              "Add JERSys : \"nominal\"/\"up\"/\"down\" explicitly to silence this." << std::endl;
    }
    logger_.Info() << "JERSys: " << JERSys << std::endl;

    PileUpSys    = SSBConfReader->GetText("PileupSys");
    L1PreFireSys = SSBConfReader->GetText("L1PreFireSys");
    TrigSFSys    = SSBConfReader->GetText("TrigSFSys");    
    LepIdSFSys   = SSBConfReader->GetText("LepIDSFSys");
    LepIsoSFSys  = SSBConfReader->GetText("LepIsoSFSys");
    LepRecoSFSys  = SSBConfReader->GetText("LepRecoSFSys");
    LepTrackSFSys  = SSBConfReader->GetText("LepTrackSFSys");


    puid_wp_ = SSBConfReader->GetText("PUIDWorkingPoint");
    PUIDSFSys = SSBConfReader->GetText("PUIDSFSys"); 
    apply_puid_ = SSBConfReader->GetBool("ApplyPUID");
    puid_pt_threshold_ = SSBConfReader->GetNumber("PUIDPtThreshold");

    // Jet_puId doesn't exist in v15; force PUID off instead of rejecting every jet.
    if (apply_puid_ && !branchReader_.JetPuIdAvailable()) {
        logger_.Warning() << "ApplyPUID=True in config but Jet_puId branch is unavailable "
                  << "in this input file - disabling PU-jet-ID for this job." << std::endl;
        apply_puid_ = false;
    }

    // Validate PUID settings
    if (puid_wp_ != "L" && puid_wp_ != "M" && puid_wp_ != "T") {
        logger_.Warning() << "Invalid PUID working point: " << puid_wp_
                  << ". Using default 'L'." << std::endl;
        puid_wp_ = "L";
    }

    if (puid_pt_threshold_ <= 0 || puid_pt_threshold_ > 200) {
        logger_.Warning() << "Invalid PUID pT threshold: " << puid_pt_threshold_
                  << ". Using default 50.0 GeV." << std::endl;
        puid_pt_threshold_ = 50.0;
    }

    // Print PUID configuration
    logger_.Info() << "PUID Configuration:" << std::endl;
    logger_.Info() << "  Working Point: " << puid_wp_ << std::endl;
    logger_.Info() << "  Systematic: " << PUIDSFSys << std::endl;
    logger_.Info() << "  Apply PUID: " << (apply_puid_ ? "true" : "false") << std::endl;
    logger_.Info() << "  pT Threshold: " << puid_pt_threshold_ << " GeV" << std::endl;

    // NanoAODv15 default: PUPPI MET. Set METtype=PF in the config explicitly if you want PFMET.
    METtype = ReadValidatedChoiceConfig(SSBConfReader.get(), "METtype", {"PF", "Puppi"}, "Puppi", logger_);
    logger_.Info() << "  MET type: " << METtype << std::endl;

    applyMETXY = ReadValidatedBoolConfig(SSBConfReader.get(), "applyMETXY", "False", logger_);
    logger_.Info() << "  apply MET XY correction: " << applyMETXY << std::endl;

    applyRochester = ReadValidatedBoolConfig(SSBConfReader.get(), "applyRochester", "False", logger_);
    logger_.Info() << "  apply Rochester correction: " << applyRochester << std::endl;
}

void Analysis::SetObjectVariable() {
    muons_pt  = branchReader_.RequiredFloatArray("Muon_pt");
    muons_eta = branchReader_.RequiredFloatArray("Muon_eta"); 
    muons_phi = branchReader_.RequiredFloatArray("Muon_phi");
    muons_M   = branchReader_.RequiredFloatArray("Muon_mass");
    
    muons_Id  = branchReader_.RequiredBoolArray("Muon_looseId");
    muons_iso = branchReader_.RequiredFloatArray("Muon_pfRelIso03_all");

    /// Muon ID
    if      (TString(MuonId).Contains( "Loose"  ) )  { muons_Id = branchReader_.RequiredBoolArray("Muon_looseId"); }
    else if (TString(MuonId).Contains( "Medium"  ) ) { muons_Id = branchReader_.RequiredBoolArray("Muon_mediumId"); }
    else if (TString(MuonId).Contains( "Tight"  ) )  { muons_Id = branchReader_.RequiredBoolArray("Muon_tightId"); }
    else { logger_.Error() << "Muon ID Error" << std::endl; }

    if (TString(MuonIsoType).Contains("PFIsodbeta03")) {
        if (branchReader_.OptionalFloatArray("Muon_pfRelIso03_all") == nullptr) {
            logger_.Error() << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
            return;
        }
        muons_iso = branchReader_.RequiredFloatArray("Muon_pfRelIso03_all");
    
    } else if (TString(MuonIsoType).Contains("PFIsodbeta04")) {
        if (branchReader_.OptionalFloatArray("Muon_pfRelIso04_all") == nullptr) {
            logger_.Error() << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
            return;
        }
        muons_iso = branchReader_.RequiredFloatArray("Muon_pfRelIso04_all");
    } else {
            logger_.Error() << "Muon Iso type Error" << std::endl;
            return;
    }

    elecs_pt  = branchReader_.RequiredFloatArray("Electron_pt");
    elecs_eta = branchReader_.RequiredFloatArray("Electron_eta");
    elecs_phi = branchReader_.RequiredFloatArray("Electron_phi");
    elecs_M   = branchReader_.RequiredFloatArray("Electron_mass");

    if (TString(ElecIsoType).Contains("PFIsoRho03")) {
        elecs_iso = branchReader_.RequiredFloatArray("Electron_pfRelIso03_all");
    } else if (TString(ElecIsoType).Contains("PFIsoRho04")) {
        elecs_iso = branchReader_.RequiredFloatArray("Electron_pfRelIso03_all");
        logger_.Warning() << "No PFIsoRho04 in NanoAOD..." << std::endl;
    } else {
        logger_.Error() << "Electron Iso type Error" << std::endl;
    }

    // Electron_cutBased is UChar_t in v15, not in intVectors; use find(), not .at() (crashes on v15).
    {
        auto it = branchReader_.intVectors.find("Electron_cutBased");
        elecs_scbId = (it != branchReader_.intVectors.end()) ? it->second.get() : nullptr;
    }
    if (TString(ElecId).Contains("SCBLoose")) {
        //elecIdVariant = branchReader_.RequiredIntArray("Electron_cutBased");
        eleid_scbcut = 2;
    } else if (TString(ElecId).Contains("SCBMedium")) {
        eleid_scbcut = 3;
    } else if (TString(ElecId).Contains("SCBTight")) {
        eleid_scbcut = 4;
    } else if (TString(ElecId).Contains("SCBVeto")) {
        eleid_scbcut = 1;
    } else if (TString(ElecId).Contains("MVALoose")) {
        elecs_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WPL");
        eleid_scbcut = 2;
    } else if (TString(ElecId).Contains("MVAMedium")) {
        elecs_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WP90");
        eleid_scbcut = 3;
    } else if (TString(ElecId).Contains("MVATight")) {
        elecs_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WP80");
        eleid_scbcut = 4;
    } else if (TString(ElecId).Contains("MVAVeto")) {
        elecs_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WPL");
        eleid_scbcut = 1;
    } else {
        logger_.Error() << "Electron ID Error" << std::endl;
    }


    // Veto-lepton collections (used for jet cleaning & third-lepton veto).
    if (TString(veto_muoniso_type).Contains("PFIsodbeta03")) {
        if (branchReader_.OptionalFloatArray("Muon_pfRelIso03_all") == nullptr) {
            logger_.Error() << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
            return;
        }
        muonsveto_iso = branchReader_.RequiredFloatArray("Muon_pfRelIso03_all");
    } else if (TString(veto_muoniso_type).Contains("PFIsodbeta04")) {
        if (branchReader_.OptionalFloatArray("Muon_pfRelIso04_all") == nullptr) {
            logger_.Error() << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
            return;
        }
        muonsveto_iso = branchReader_.RequiredFloatArray("Muon_pfRelIso04_all");
    } else {
        logger_.Error() << "Muon Iso type Error" << std::endl;
        return;
    }

    
    if (muonsveto_iso != nullptr) {
    } else {
        logger_.Error() << "Error: muonsveto_iso is null after assignment!" << std::endl;
    }


    /// Muon ID
    if      (TString(veto_muonid).Contains( "Loose"  ) )  { muonsveto_Id = branchReader_.RequiredBoolArray("Muon_looseId"); }
    else if (TString(veto_muonid).Contains( "Medium"  ) ) { muonsveto_Id = branchReader_.RequiredBoolArray("Muon_mediumId");}
    else if (TString(veto_muonid).Contains( "Tight"  ) )  { muonsveto_Id = branchReader_.RequiredBoolArray("Muon_tightId"); }
    else { logger_.Error() << "Muon ID Error" << std::endl; }
    
    if (muonsveto_Id == nullptr) {
        logger_.Error() << "Error: muonsveto_Id is null after assignment!" << std::endl;
    }


    // Same UChar_t-in-v15 issue as elecs_scbId above; find()-guarded, not .at().
    {
        auto it = branchReader_.intVectors.find("Electron_cutBased");
        elecsveto_scbId = (it != branchReader_.intVectors.end()) ? it->second.get() : nullptr;
    }
    if (!branchReader_.BranchIsAvailable("Electron_cutBased")) {logger_.Error() <<"Error: Electron_cutBased branch not found at all!" <<std::endl;}

    if (TString(veto_elecid).Contains("SCBLoose")) {
        //elecIdVariant = branchReader_.RequiredIntArray("Electron_cutBased");
        elevetoid_scbcut = 2;
    } else if (TString(veto_elecid).Contains("SCBMedium")) {
        elevetoid_scbcut = 3;
    } else if (TString(veto_elecid).Contains("SCVTight")) {
        elevetoid_scbcut = 4;
    } else if (TString(veto_elecid).Contains("SCBVeto")) {
        elevetoid_scbcut = 1;
    } else if (TString(veto_elecid).Contains("MVALoose")) {
        elecsveto_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WPL");
        elevetoid_scbcut = 2;
    } else if (TString(veto_elecid).Contains("MVAMedium")) {
        elecsveto_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WP90");
        elevetoid_scbcut = 3;
    } else if (TString(veto_elecid).Contains("MVATight")) {
        elecsveto_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WP80");
        elevetoid_scbcut = 4;
    } else if (TString(veto_elecid).Contains("MVAVeto")) {
        elecsveto_mvaId = branchReader_.RequiredBoolArray("Electron_mvaFall17V2Iso_WPL");
        elevetoid_scbcut = 1;
    } else {
        logger_.Error() << "Electron ID in veto selection Error " << veto_elecid << std::endl;
    }

    if (TString(veto_eleciso_type).Contains("PFIsoRho03")) {
        if (branchReader_.floatVectors.find("Electron_pfRelIso03_all") != branchReader_.floatVectors.end()) {
            auto* ptr = branchReader_.RequiredFloatArray("Electron_pfRelIso03_all");
            if (ptr) {
                elecsveto_iso = ptr;
            } else {
                logger_.Error() << "Error: 'Electron_pfRelIso03_all' is a null unique_ptr." << std::endl;
            }
        } else {
            logger_.Error() << "Error: 'Electron_pfRelIso03_all' not found in branchReader_.floatVectors." << std::endl;
        }
    } else if (TString(veto_eleciso_type).Contains("PFIsoRho04")) {
        logger_.Warning() << "No PFIsoRho04 in NanoAOD..." << std::endl;
    } else {
        logger_.Error() << "Electron Iso type Error" << std::endl;
    }

    if (elecsveto_iso == nullptr) {
        logger_.Error() << "Error: elecsveto_iso is null after assignment!" << std::endl;
    }


    jets_pt  = branchReader_.RequiredFloatArray("Jet_pt");
    jets_eta = branchReader_.RequiredFloatArray("Jet_eta");
    jets_phi = branchReader_.RequiredFloatArray("Jet_phi");
    jets_M   = branchReader_.RequiredFloatArray("Jet_mass");
    // v15 has no Jet_jetId; jet ID is now evaluated per-jet via
    // PassConfiguredJetId() (interface/JetID.h) instead.
    jets_Id = nullptr;
    jets_puId = branchReader_.OptionalIntArray("Jet_puId"); // may be null (removed in v15)

    if (!isData){
        gen_jets_pt  = branchReader_.RequiredFloatArray("GenJet_pt");
        gen_jets_eta = branchReader_.RequiredFloatArray("GenJet_eta");
        gen_jets_phi = branchReader_.RequiredFloatArray("GenJet_phi");
        gen_jets_M   = branchReader_.RequiredFloatArray("GenJet_mass");
    }

    // "PFLoose"/"PFLooseLepVeto" were valid pre-v15 but have no Puppi formula.
    if (JetId != "PFTight" && JetId != "PFTightLepVeto") {
        logger_.Warning() << "Jet_ID='" << JetId << "' is not supported for NanoAODv15 Puppi jets. "
                  << "Use PFTight or PFTightLepVeto in the config." << std::endl;
    }


    // Set b-tag discriminator and threshold based on algorithm, WP and RunPeriod
    
    // First set the appropriate b-tag discriminator variable
    if (TString(JetbTag).Contains("deepCSV")) {
        jets_btag = branchReader_.RequiredFloatArray("Jet_btagDeepB");
    }
    else if (TString(JetbTag).Contains("deepJet")) {
        jets_btag = branchReader_.RequiredFloatArray("Jet_btagDeepFlavB");
    }
    else if (TString(JetbTag).Contains("UParT")) {
        // UParT (Unified Particle Transformer) AK4 b-tag score, new in NanoAODv15.
        jets_btag = branchReader_.RequiredFloatArray("Jet_btagUParTAK4B");
    }
    else if (TString(JetbTag).Contains("pfCSVV2")) {
        jets_btag = branchReader_.RequiredFloatArray("Jet_btagCSVV2");
    }
    else if (TString(JetbTag).Contains("CSV") && !TString(JetbTag).Contains("deep")) {
        jets_btag = branchReader_.RequiredFloatArray("Jet_btagCSV"); // Modify with appropriate variable name if needed
    }
    else if (TString(JetbTag).Contains("CISV")) {
        jets_btag = branchReader_.RequiredFloatArray("Jet_btagCISV"); // Modify with appropriate variable name if needed
    }
    else {
        logger_.Error() << "Error: Unknown b-tagging algorithm in " << JetbTag << std::endl;
        jets_btag = nullptr;
        bdisccut = -1.0;
        return;
    }

    // v15 removed DeepCSV branches; catch a resulting null here instead of
    // segfaulting later at jets_btag->At(...).
    if (jets_btag == nullptr) {
        logger_.Error() << "jets_btag is null: the branch for Jet_btag='" << JetbTag
                  << "' was not found in this input file. NanoAODv15 dropped the DeepCSV "
                  << "taggers - set Jet_btag to a deepJet (or PNet/UParT) working point in "
                  << "your config." << std::endl;
        bdisccut = -1.0;
        return;
    }


    // Shared with SSBCorrections.cpp's constructor so the two can't drift apart.
    btag_algo_ = ParseBTagAlgo(JetbTag.Data());  // JetbTag is a TString - ParseBTagAlgo takes std::string
    if (btag_algo_ == BTagAlgo::Unknown) {
        logger_.Warning() << "Unknown b-tagging algorithm in JetbTag: " << JetbTag << std::endl;
    }
    
    // Parse working point from JetbTag
    btag_wp_ = "M";  // Default to Medium
    if (TString(JetbTag).Contains("L")) {
        btag_wp_ = "L";
    } else if (TString(JetbTag).Contains("M")) {
        btag_wp_ = "M";
    } else if (TString(JetbTag).Contains("T")) {
        btag_wp_ = "T";
    } else {
        logger_.Warning() << "Unknown working point in JetbTag: " << JetbTag << std::endl;
    }
    
 
    // New taggers (e.g. UParT) have no hardcoded cut below; read "BTagDiscCut"
    // from config if present, which also lets it override the hardcoded cuts.
    double configBTagDiscCut = SSBConfReader->GetNumber("BTagDiscCut", /*PrintError=*/false);
    bool haveConfigBTagDiscCut = SSBConfReader->Check("BTagDiscCut");
    if (haveConfigBTagDiscCut) {
        bdisccut = configBTagDiscCut;
        // Print once, not per-event (SetObjectVariable runs every event).
        static bool printedBTagDiscCutInfo = false;
        if (!printedBTagDiscCutInfo) {
            logger_.Info() << "Using BTagDiscCut from config: " << bdisccut << std::endl;
            printedBTagDiscCutInfo = true;
        }
    } else if (btag_algo_ == BTagAlgo::UParTAK4) {
        // Look up the cut from btagging.json.gz's UParTAK4_wp_values instead of
        // requiring it in config; falls back to hard failure if unavailable.
        double lookedUpCut = SSBCorr->GetBtagWPCut(btag_wp_);
        static bool printedBTagWPCutInfo = false;
        if (lookedUpCut > 0.0) {
            bdisccut = lookedUpCut;
            if (!printedBTagWPCutInfo) {
                logger_.Info() << "BTagDiscCut not set in config - looked up " << bdisccut
                          << " from btagging.json.gz's UParTAK4_wp_values for WP='" << btag_wp_
                          << "'." << std::endl;
                printedBTagWPCutInfo = true;
            }
        } else {
            logger_.Error() << "Jet_btag='" << JetbTag << "' (UParT) has no hardcoded working-point "
                      << "cut value, 'BTagDiscCut' was not set in the config, AND the "
                      << "UParTAK4_wp_values lookup didn't return a usable cut (see any "
                      << "[WARNING] above from GetBtagWPCut/InitBtagSFCorrection for why). Add e.g. "
                      << "'BTagDiscCut : 0.xxxx' to the config with the CMS BTV-recommended UParT "
                      << "AK4 " << btag_wp_ << " cut for this campaign as a manual fallback."
                      << std::endl;
            bdisccut = -1.0; // fails NumbJetCut/BTaggingSFApply safely rather than tagging everything
        }
    }

    if (!haveConfigBTagDiscCut)
    if (TString(RunPeriod).Contains("2016PreVFP")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.2027;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.6001;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.8819;
            else logger_.Warning() << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0508;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.2598;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.6502;
            else logger_.Warning() << "Unknown deepJet working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("pfCSVV2")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.5426;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.8484;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.9535;
            else logger_.Warning() << "Unknown pfCSVV2 working point!" << std::endl;
        }
    }
    else if (TString(RunPeriod).Contains("2016PostVFP")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.1918;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.5847;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.8767;
            else logger_.Warning() << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0480;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.2489;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.6377;
            else logger_.Warning() << "Unknown deepJet working point!" << std::endl;
        }
    }
    else if (TString(RunPeriod).Contains("2017")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.1355;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.4506;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7738;
            else logger_.Warning() << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0532;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.3040;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7476;
            else logger_.Warning() << "Unknown deepJet working point!" << std::endl;
        }
    }
    else if (TString(RunPeriod).Contains("2018")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.1208;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.4168;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7665;
            else logger_.Warning() << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0490;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.2783;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7100;
            else logger_.Warning() << "Unknown deepJet working point!" << std::endl;
        }
    }
    else {
        logger_.Error() << "Error: Unsupported run period: " << RunPeriod << std::endl;
        jets_btag = nullptr;
        bdisccut = -1.0;
    }
    
    // MET_pt/MET_phi were renamed to PFMET_pt/PFMET_phi in NanoAODv15 -
    // try both names so the same "PF" config value works on v9 and v15 files.
    if (METtype == "PF"){
        met_pt  = branchReader_.GetFloatSinglePtrByAlias({"MET_pt", "PFMET_pt"});
        met_phi = branchReader_.GetFloatSinglePtrByAlias({"MET_phi", "PFMET_phi"});}
    else if (METtype == "Puppi"){
	met_pt  = branchReader_.RequiredFloatSingle("PuppiMET_pt");
        met_phi  = branchReader_.RequiredFloatSingle("PuppiMET_phi");}
    else {
        // Unrecognized/unset METtype: default to PuppiMET (NanoAODv15 default).
        logger_.Warning() << "METtype='" << METtype << "' not recognized - defaulting to PuppiMET." << std::endl;
        met_pt  = branchReader_.RequiredFloatSingle("PuppiMET_pt");
        met_phi = branchReader_.RequiredFloatSingle("PuppiMET_phi");
    }

    object_variables_set_ = true; 

}



// Shared per-stage control-plot filler; see interface/Analysis.h for details.
void Analysis::FillControlPlots(SelectionStage stage, bool includeJets, bool includeBJets) {
    const std::size_t i = selectionIndex(stage);

    FillHisto( h_DiLepMass[i], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
    FillHisto( h_Num_PV[i],    num_pv,                evt_weight_ );
    FillHisto( h_Lep1pt[i] ,   (Lep1).Pt()  , evt_weight_ );
    FillHisto( h_Lep1eta[i],   (Lep1).Eta() , evt_weight_ );
    FillHisto( h_Lep1phi[i],   (Lep1).Phi() , evt_weight_ );
    FillHisto( h_Lep2pt[i] ,   (Lep2).Pt()  , evt_weight_ );
    FillHisto( h_Lep2eta[i],   (Lep2).Eta() , evt_weight_ );
    FillHisto( h_Lep2phi[i],   (Lep2).Phi() , evt_weight_ );
    FillHisto( h_METpt[i]  ,   Met.Pt()  , evt_weight_ );
    FillHisto( h_METphi[i] ,   Met.Phi() , evt_weight_ );
    FillHisto( h_Num_Jets[i] , v_jet_idx.size(), evt_weight_ );

    if (includeJets) {
        FillHisto( h_Jet1pt[i] ,   (Jet1).Pt()  , evt_weight_ );
        FillHisto( h_Jet1eta[i],   (Jet1).Eta() , evt_weight_ );
        FillHisto( h_Jet1phi[i],   (Jet1).Phi() , evt_weight_ );
        FillHisto( h_Jet2pt[i] ,   (Jet2).Pt()  , evt_weight_ );
        FillHisto( h_Jet2eta[i],   (Jet2).Eta() , evt_weight_ );
        FillHisto( h_Jet2phi[i],   (Jet2).Phi() , evt_weight_ );
    }

    if (includeBJets) {
        FillHisto( h_Num_bJets[i], v_bjet_idx.size(), evt_weight_ );
    }
}

// Event loop function
void Analysis::Loop() {
    if (!chain) {
        logger_.Error() << "Error: TChain is null!" << std::endl;
        return;
    }

    Long64_t nEntries = chain->GetEntries();
    fReader.Restart(); // Reset the reader to the beginning

    MCSF();

    for (Long64_t ientry = 0; ientry < ((NumEvt == -1 || NumEvt > nEntries) ? nEntries : NumEvt); ++ientry) {


	current_entry_ = ientry;

        if (!fReader.Next()) { // Move to the next entry
            logger_.Error() << "Error: Failed to read entry " << ientry << std::endl;
            break;
        }

        // Debug event tracer: only prints when debugFilter_ matches this event.
        // Intentionally bypasses Logger's level gate - tracing is an explicit
        // ask, not a verbosity setting.
        unsigned int       traceRun  = branchReader_.RequiredUIntValue("run");
        unsigned int       traceLumi = branchReader_.RequiredUIntValue("luminosityBlock");
        unsigned long long traceEvt  = branchReader_.RequiredULongValue("event");
        bool isTracedEvent = debugFilter_.Matches(traceRun, traceLumi, traceEvt);
        if (isTracedEvent) {
            std::cout << "[TRACE " << traceRun << ":" << traceLumi << ":" << traceEvt
                      << "] entry=" << ientry << " isData=" << isData << std::endl;
        }
        // Prints which cut a traced event fails, right before its `continue`.
        auto traceCutFail = [&](const char* cutName) {
            if (isTracedEvent) {
                std::cout << "[TRACE " << traceRun << ":" << traceLumi << ":" << traceEvt
                          << "] failed: " << cutName << std::endl;
            }
        };

        evt_weight_ = 1.;
        MCSFApply();

        GenWeightApply();
        PUWeightApply();
        L1PreFireApply();

        SetObjectVariable(); //
        LeptonSelector(); 

        LeptonOrder();
        JetSelector();
	if (isjetveto_event_) {
            if (isTracedEvent) {
                std::cout << "[TRACE " << traceRun << ":" << traceLumi << ":" << traceEvt
                          << "] skipped: isjetveto_event_" << std::endl;
            }
            continue; // Skip entire event having jetveto map jet
        }
        PUIDSFApply();  // Apply PUID event weight using collected information

        JetOrder();
        bJetSelector();
        //METDefiner();

        if (isTracedEvent) {
            std::cout << "[TRACE " << traceRun << ":" << traceLumi << ":" << traceEvt
                      << "] after JetSelector/PUIDSFApply/bJetSelector: "
                      << "nLepton=" << v_lepton_idx.size()
                      << " nJet=" << v_jet_idx.size()
                      << " nBJet=" << v_bjet_idx.size()
                      << " evt_weight=" << evt_weight_ << std::endl;
        }

        //if (i > 100) break; //%lld supports Long64_t
        if (ientry % 10000 == 0) {
            logger_.Info() << "Processed " << ientry << " events" << std::endl;
        }


        if ( METFilterAPP() == false ) {traceCutFail("METFilterAPP"); continue;}
        if ( Trigger() == false ) {traceCutFail("Trigger"); continue;}

        // PV_npvsGood is UChar_t in NanoAODv15 (was Int_t in v9) - read it
        // through the version-agnostic single-value accessor.
        if (branchReader_.GetIntSingleValue("PV_npvsGood") < 1) {
            traceCutFail("PV_npvsGood>=1");
            continue;
        }

        FillHisto( h_Num_PV[selectionIndex(SelectionStage::AfterTriggerAndPV)]    , branchReader_.GetIntSingleValue("PV_npvsGood") , evt_weight_ );
        if ( NumIsoLeptons(2) == false ) {traceCutFail("NumIsoLeptons(2)"); continue;}

        if (ThirdLeptonVeto() == false ) {traceCutFail("ThirdLeptonVeto"); continue;}
        if (LeptonsPtAddtional() == false ) {traceCutFail("LeptonsPtAddtional"); continue;}
        if (DiLeptonMassCut() == false) {traceCutFail("DiLeptonMassCut"); continue;}

        LeptonSFApply();
        TriggerSFApply();

        num_pv = static_cast<int>(branchReader_.GetIntSingleValue("PV_npvsGood"));
        // Stage 1: no jets, no b-jets (h_Num_bJets was already commented
        // out here before this consolidation - left out on purpose).
        FillControlPlots(SelectionStage::AfterDileptonCuts, /*includeJets=*/false, /*includeBJets=*/false);

        if (ZVetoCut() == false) {traceCutFail("ZVetoCut"); continue;}

        FillControlPlots(SelectionStage::AfterZVeto, /*includeJets=*/false, /*includeBJets=*/false);

        if (NumJetCut(v_jet_idx) == false) {traceCutFail("NumJetCut"); continue;}

        FillControlPlots(SelectionStage::AfterJetMultiplicity, /*includeJets=*/true, /*includeBJets=*/false);

        if (METCut(Met) == false) {traceCutFail("METCut"); continue;}

        FillControlPlots(SelectionStage::AfterMETCut, /*includeJets=*/true, /*includeBJets=*/false);

        BTaggingSFApply();

        if (NumbJetCut(v_bjet_idx) == false) {traceCutFail("NumbJetCut"); continue;}

        FillControlPlots(SelectionStage::AfterBTagMultiplicity, /*includeJets=*/true, /*includeBJets=*/false);
        SetUpKINObs();
        if (!isKinSol) {
            traceCutFail("KinSolv (no solution)");
        }
        if (isKinSol)
        {
            if (isTracedEvent) {
                std::cout << "[TRACE " << traceRun << ":" << traceLumi << ":" << traceEvt
                          << "] survived to Top-Recon. stage: evt_weight=" << evt_weight_
                          << " Top.M()=" << Top.M() << " AnTop.M()=" << AnTop.M() << std::endl;
            }
            //FillHisto( h_HT[selectionIndex(SelectionStage::AfterTopReconstruction)]     , AllJetpt   , evt_weight_);
            FillControlPlots(SelectionStage::AfterTopReconstruction, /*includeJets=*/true, /*includeBJets=*/true);
            if ( Top.Pt() > AnTop.Pt() ) { Top1 = Top; Top2 = AnTop; }
            else { Top1 = AnTop; Top2 = Top; }
            
            FillHisto( h_TopMass      , Top.M()         , evt_weight_ );
            FillHisto( h_Toppt        , Top.Pt()        , evt_weight_ );
            FillHisto( h_Topphi       , Top.Phi()       , evt_weight_ );
            FillHisto( h_TopRapidity  , Top.Rapidity()  , evt_weight_ );
            FillHisto( h_TopEnergy    , Top.Energy()    , evt_weight_ );
            FillHisto( h_AnTopMass    , AnTop.M()       , evt_weight_ );
            FillHisto( h_AnToppt      , AnTop.Pt()      , evt_weight_ );
            FillHisto( h_AnTopphi     , AnTop.Phi()     , evt_weight_ );
            FillHisto( h_AnTopRapidity, AnTop.Rapidity(), evt_weight_ );
            FillHisto( h_AnTopEnergy  , AnTop.Energy()  , evt_weight_ );
            
            FillHisto( h_W1Mass , W1.M()  , evt_weight_ );
            FillHisto( h_W2Mass , W2.M()  , evt_weight_ );
            
            FillHisto( h_W1Mt , W1.Mt()  , evt_weight_ );
            FillHisto( h_W2Mt , W2.Mt()  , evt_weight_ );
            
            //FillHisto( h_bJet1Energy , bJet1.Energy()  , evt_weight_ );
            //FillHisto( h_bJet2Energy , bJet2.Energy()  , evt_weight_ );
            
            FillHisto( h_bJetEnergy   , bJet.Energy()   , evt_weight_ );
            FillHisto( h_AnbJetEnergy , AnbJet.Energy() , evt_weight_ );
            FillHisto( h_bJetPt       , bJet.Pt()   , evt_weight_ );
            FillHisto( h_AnbJetPt     , AnbJet.Pt() , evt_weight_ );
            FillHisto( h_LepEnergy    , Lep.Energy()    , evt_weight_ );
            FillHisto( h_AnLepEnergy  , AnLep.Energy()  , evt_weight_ );
            FillHisto( h_NuEnergy     , Nu.Energy()     , evt_weight_ );
            FillHisto( h_AnNuEnergy   , AnNu.Energy()   , evt_weight_ );

            std::vector<double> v_recocp_O;
            v_recocp_O.push_back( SSBCPVCal->getO1Vari( Top, AnTop, AnLep, Lep )  );
            v_recocp_O.push_back( SSBCPVCal->getO2Vari( Top, AnTop, bJet, AnbJet ) );
            v_recocp_O.push_back( SSBCPVCal->getO3Vari( bJet, AnbJet, AnLep, Lep ) );
            v_recocp_O.push_back( SSBCPVCal->getO4Vari( AnbJet, bJet, AnLep, Lep ) );
            v_recocp_O.push_back( SSBCPVCal->getO5Vari( bJet , AnbJet, AnLep, Lep ) );
            v_recocp_O.push_back( SSBCPVCal->getO6Vari( bJet , AnbJet, AnLep, Lep ) );
            v_recocp_O.push_back( SSBCPVCal->getO7Vari( Top , AnTop, AnLep, Lep ) );
            v_recocp_O.push_back( SSBCPVCal->getO8Vari( Top, AnTop, bJet , AnbJet, AnLep, Lep ) );
            v_recocp_O.push_back( SSBCPVCal->getO9Vari( bJet , AnbJet, AnLep, Lep )  );
            v_recocp_O.push_back( SSBCPVCal->getO10Vari( bJet , AnbJet, AnLep, Lep )  );
            v_recocp_O.push_back( SSBCPVCal->getO11Vari( bJet , AnbJet, AnLep, Lep )  );
            v_recocp_O.push_back( SSBCPVCal->getO12Vari( bJet , AnbJet, AnLep, Lep )  );
            v_recocp_O.push_back( SSBCPVCal->getO13Vari( bJet , AnbJet, AnLep, Lep )  );

            // Blind CP observables for data: sign is randomized per-observable,
            // deterministically seeded from run/lumi/event.
            if (isData && isBlind) {
                unsigned int      run_num  = branchReader_.RequiredUIntValue("run");
                unsigned int      lumi_num = branchReader_.RequiredUIntValue("luminosityBlock");
                unsigned long long evt_num = branchReader_.RequiredULongValue("event");

                uint64_t seed = (uint64_t)run_num  * 100000000ULL
                              + (uint64_t)lumi_num * 1000000ULL
                              + (evt_num & 0xFFFFFFFFULL);
                seed ^= (seed >> 33);
                seed *= 0xff51afd7ed558ccdULL;
                seed ^= (seed >> 33);
                seed *= 0xc4ceb9fe1a85ec53ULL;
                seed ^= (seed >> 33);

                for (int i = 0; i < (int)v_recocp_O.size(); ++i) {
                    if ((seed >> i) & 1ULL) v_recocp_O[i] = -v_recocp_O[i];
                }
            }

            for (int i = 0; i < (int)v_recocp_O.size(); ++i)
            {
               FillHisto( h_Reco_CPO_[i],         v_recocp_O[i], evt_weight_ );
               FillHisto( h_Reco_CPO_ReRange_[i],  v_recocp_O[i], evt_weight_ );
            }

                  
 
        }


    }// end of event iteration //
    logger_.Info() << "End Loop !!" << std::endl;
    fallbackCounter_.PrintSummary();
}


template <typename T>
std::unique_ptr<TTreeReaderValue<T>> Analysis::DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src) {
    if (src) {
        // If source pointer is valid, make a deep copy
        return std::make_unique<TTreeReaderValue<T>>(*src);
    } else {
	logger_.Warning() << "no source !!! "  << std::endl;
        // If source is nullptr, return nullptr
        return nullptr;
    }
}



std::string Analysis::removeSubstring(std::string &str, const std::string &keyword) {
    size_t pos = str.find(keyword);  // find specific keyword ! 
    if (pos != std::string::npos) {
        str.erase(pos);  // 
    }
    logger_.Debug() << "str : " << str << std::endl;
    return str; 
}

bool Analysis::METFilterAPP()
{
   bool metfilt_ = true;
   for ( const auto &pair : noiseFilters )
   { 
      TString METFiltName = pair.first;
//      if ( !(TString(FileName_).Contains("Data")) && 
//            TString(METFiltName).Contains("Flag_eeBadScFilter")  ) {continue;}
            if ( !(**pair.second) ){ metfilt_ = false; }
   }
   return metfilt_;
   //return true;
}

// Trigger Requirement Function
bool Analysis::SelTrigger(std::vector<std::string> v_sel)
{
    std::string trgName = "";
    int ptrigindex = 0;
    bool passtrig_ = false;

    for (int j = 0; j < v_sel.size(); j++)
    {
        trgName = v_sel[j];

        // trgName cout in triggerList 
        auto it = triggerList.find(trgName);
        if (it != triggerList.end() && it->second) {
            if (**(it->second)) {ptrigindex++;}
        } else {
            logger_.Error() << "Error: Trigger " << trgName << " not found in triggerList." << std::endl;
        }
    }

    if (ptrigindex > 0) {
        passtrig_ = true;
    }
    return passtrig_;
}


bool Analysis::Trigger()
{
  bool trigpass = false;

  if (!TString(FileName_).Contains("Data")) {
      // MC: always apply trigName
      trigpass = SelTrigger(trigName);
      return trigpass;
  }

  /// Data Samples
  std::vector<std::string> seltrigName;
  std::vector<std::string> vetotrigName;

  if (RunPeriod.Contains("2018")) {
     if (TString(Decaymode).Contains("dimuon")) {
        if (TString(FileName_).Contains("Single")) {
           bool pass_single = SelTrigger(SLtrigName);
           bool pass_double = SelTrigger(DLtrigName);
           trigpass = pass_single && !pass_double;
           return trigpass;
        } else if (TString(FileName_).Contains("Double")) {
           bool pass_double = SelTrigger(DLtrigName);
           //bool pass_single = SelTrigger(SLtrigName);
           trigpass = pass_double;
           //trigpass = pass_single || pass_double;
           return trigpass;
        } else {
           logger_.Warning() << "[Trigger] Check FileName_ for dimuon in 2018" << std::endl;
           return false;
        }
     }

     else if (TString(Decaymode).Contains("muel")) {
        if (TString(FileName_).Contains("MuonEG")) {
           // MuonEG: Primary dataset for muel, use Double lepton triggers only
           bool pass_double = SelTrigger(DLtrigName);
           //bool pass_single = SelTrigger(SLtrigName);
           trigpass = pass_double;
           //trigpass = pass_single || pass_double;
           return trigpass;
        }
        else if (TString(FileName_).Contains("SingleMuon") || TString(FileName_).Contains("EGamma")) {
           // SingleMuon/EGamma: Use Single triggers only, veto Double to avoid overlap with MuonEG
           bool pass_single = SelTrigger(SLtrigName);
           bool pass_double = SelTrigger(DLtrigName);
           trigpass = pass_single && !pass_double;
           return trigpass;
        }
        else {
           logger_.Warning() << "[Trigger] Check FileName_ for muel in 2018" << std::endl;
           return false;
        }
     }

     else if (TString(Decaymode).Contains("dielec")) {
        if (TString(FileName_).Contains("EGamma")) {
           bool pass_single = SelTrigger(SLtrigName);
           bool pass_double = SelTrigger(DLtrigName);
           //trigpass = (pass_single && !pass_double) || (!pass_single && pass_double);
           trigpass = pass_single || pass_double;
           return trigpass;
        } else {
           logger_.Warning() << "[Trigger] Check FileName_ for dielec in 2018 : FileName_ : " << FileName_ << std::endl;
           return false;
        }
     }

     else {
        logger_.Warning() << "[Trigger] Check Decaymode for 2018" << std::endl;
        return false;
     }
  }

  /// 2016 / 2017
  else {
     if (TString(FileName_).Contains("Single")) {
        bool pass_single = SelTrigger(SLtrigName);
        bool pass_double = SelTrigger(DLtrigName);
        trigpass = pass_single && !pass_double;
        return trigpass;
     } else if (TString(FileName_).Contains("Double") || TString(FileName_).Contains("MuonEG")) {
        bool pass_double = SelTrigger(DLtrigName);
        trigpass = pass_double;
        return trigpass;
     } else {
        logger_.Warning() << "[Trigger] Check FileName_ for 2016/2017" << std::endl;
        return false;
     }
  }
}

TString Analysis::SetInputFileName(std::string inname)
{  
   TString inputName = inname;
   
   // Remove file extension
   if (inputName.Contains(".")) {
      Size_t dotIndex = inputName.Last('.');
      inputName.Remove(dotIndex, inputName.Length());
   }
   
   // Remove directory path if present
   if (inputName.Contains("/")) {
      Size_t slashIndex = inputName.Last('/');
      inputName.Remove(0, slashIndex+1);
   }
   
   // Remove numeric pattern after the last underscore
   Size_t last_underscore = inputName.Last('_');
   if (last_underscore != kNPOS) {
      TString suffix = inputName(last_underscore+1, inputName.Length()-last_underscore-1);
      // Check if the suffix consists only of digits
      bool is_numeric = true;
      for (Int_t i = 0; i < suffix.Length(); i++) {
         if (!isdigit(suffix[i])) {
            is_numeric = false;
            break;
         }
      }
      // If the suffix is numeric, remove it
      if (is_numeric) {
         inputName.Remove(last_underscore, inputName.Length()-last_underscore);
      }
   }
   
   logger_.Info() << "Original input: " << inname << ", processed name: " << inputName << std::endl;
   
   return inputName; 
}

void Analysis::MCSF()
{
    if (FileName_.Contains("Data")||FileName_.Contains("Single")||FileName_.Contains("EG")){ mc_sf_ = 1.; logger_.Info() << "mc_sf_ : " << mc_sf_ << std::endl; return; }
    FILE *xsecs_;
    char sampleName[1000];
    double xsec_ = -1.;
    double br_ = -1.; 
    int totalevt_  = -1.; 
    int positive_  = -1.; 
    int negative_  = -1.; 
    int posi_nega_ = -1.; 
    std::string xsec_dir= "./xsecAndsample/";
    std::string xsec_filePath = xsec_dir+ XsecTable_.Data();
    /// SampleName | TotalEvt | Positive+Negative | Xsection | Branching Fraction |
    xsecs_ = fopen(xsec_filePath.c_str(),"r");
    std::map<std::string, int> m_sam_totalevt;
    std::map<std::string, double> m_sam_xsec;
    std::map<std::string, double> m_sam_br;
    std::map<std::string, int> m_sam_positive;
    std::map<std::string, int> m_sam_negative;
    std::map<std::string, int> m_sam_posi_nega;
    if (xsecs_!=NULL) 
    { 
       while (fscanf(xsecs_, "%s %d %d %d %d %lf %lf\n", sampleName, &totalevt_, &positive_, &negative_, &posi_nega_, &xsec_, &br_ ) != EOF)
       {
          logger_.Info()
          << "sampleName : " << sampleName << " totalevt_ : " << totalevt_
          << " positive_ " << positive_ << " negative_ : " << negative_ 
          << " posi_nega_ " << posi_nega_ << " xsec_ : " << xsec_ 
          << " br_ " << br_ 
          << std::endl;
          m_sam_totalevt[sampleName] = totalevt_;
          m_sam_positive[sampleName] = positive_;
          m_sam_negative[sampleName] = negative_;
          m_sam_posi_nega[sampleName] = posi_nega_;
          m_sam_xsec[sampleName] = xsec_;
          m_sam_br[sampleName] = br_;
       }
       fclose(xsecs_);
    }
    else {logger_.Error() << "No xsec_filePath !!!" << xsec_filePath << std::endl;return;} 
    logger_.Info() << "Lumi : " << Lumi << std::endl;
    double lumi = Lumi/1000000;
    auto it = m_sam_xsec.find(FileName_.Data());
    if (it !=  m_sam_xsec.end()){
        logger_.Info() << "SK Key " << FileName_.Data() << " found in the std::map."<< std::endl;
        mc_sf_ = (m_sam_xsec[FileName_.Data()]*m_sam_br[FileName_.Data()]*lumi)/m_sam_posi_nega[FileName_.Data()];
        logger_.Info() << "mc_sf_ " << mc_sf_ << std::endl;
    }
    else {
        mc_sf_ =1.;
        logger_.Warning() << "Key " << FileName_.Data() << " not found in the std::map. mc sf is 1" << mc_sf_ << std::endl;
    }
    return;
}


void Analysis::MCSFApply()
{
    evt_weight_beforemcsf_ = 1;
    evt_weight_beforemcsf_ = evt_weight_;

    if ( !TString(FileName_).Contains( "Data") ){ evt_weight_ = evt_weight_*mc_sf_; }
    else {evt_weight_ = 1;}
}

void Analysis::L1PreFireApply()
{   
    evt_weight_beforeL1PreFire_ = evt_weight_; // keep event weight

    // Early return for Data or 2018+ (no L1 prefiring correction needed)
    if (isData || TString(RunPeriod).Contains("2018")) {return;}

    // Apply L1 prefiring correction for specific years
    bool needsL1PreFire = TString(RunPeriod).Contains("2016") || TString(RunPeriod).Contains("2017");
                         // Add Run3 years here if needed: || TString(RunPeriod).Contains("2022") || ...

    if (!needsL1PreFire) {
        return; // Skip L1 prefiring for years that don't need it
    }

    // Apply L1 prefiring weight for 2016/2017 MC only
    double l1prefire_ = 1.0;
    
    // Use const char* comparison for better performance
    const char* sys = L1PreFireSys.Data();

    if (strstr(sys, "central")) {
        l1prefire_ = branchReader_.RequiredFloatValue("L1PreFiringWeight_Nom");
    }
    else if (strstr(sys, "up")) {
        l1prefire_ = branchReader_.RequiredFloatValue("L1PreFiringWeight_Up");
    }
    else if (strstr(sys, "down")) {
        l1prefire_ = branchReader_.RequiredFloatValue("L1PreFiringWeight_Dn");
    }
    else if (!strstr(sys, "none")) {
        // Only print error if not "none"
        logger_.Warning() << "L1Prefiring sys Error ... Default is Weight_L1Prefiring ... : " << L1PreFireSys << std::endl;
        l1prefire_ = branchReader_.RequiredFloatValue("L1PreFiringWeight_Nom");
    }
    evt_weight_ *= l1prefire_;
}

// Safely create a TLorentzVector
TLorentzVector Analysis::createLorentzVector(float pt, float eta, float phi, float mass) {
    // Check validity of inputs
    if (pt < 0 || std::isnan(pt) || std::isnan(eta) || std::isnan(phi) || std::isnan(mass)) {
        throw std::runtime_error("Invalid inputs for TLorentzVector");
    }

    // Normalize phi to the range [-π, π]
    phi = TVector2::Phi_mpi_pi(phi);

    // Create the TLorentzVector
    TLorentzVector lv;// = new TLorentzVector();
    lv.SetPtEtaPhiM(static_cast<double>(pt),
                    static_cast<double>(eta),
                    static_cast<double>(phi),
                    static_cast<double>(mass));

    return lv;
}

void Analysis::LeptonSelector() {
    // Clear output vectors
    v_lepton_idx.clear();
    v_muon_idx.clear();
    v_electron_idx.clear();

    // Initialize vectors to store TLorentzVectors
    muons.clear();
    elecs.clear();

    // Prepare corrected lepton collections
    MakeMuonCollection();
    MakeElecCollection();

    // Lambda functions for common checks
    auto passKinematicCuts = [](float _pt, float _eta, float ptCut, float etaCut) {
      return _pt > ptCut && fabs(_eta) < etaCut;
    };

    auto passIsolation = [](float iso, float isoCut) {
      return iso <= isoCut;
    };

    auto passId = [](bool id) {
      return id;
    };

    auto elecSCBId = [](int id, int idcut) {
      return id >= idcut; 
    };
    
    auto elecCharge = [](int id) {
      return id > 0; // Electron_tightCharge check
    };

    // Dimuon channel
    if (TString(Decaymode).Contains("dimuon")) {
        Int_t nmu = muons_iso->GetSize();
        for (int i = 0; i < nmu; ++i) {
            // Skip muons that don't pass selection criteria
            if (!passIsolation(muons_iso->At(i), muon_isocut) ||
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i), muon_pt, muon_eta) ||
                !passId(muons_Id->At(i))) {
                continue;
            }
            
            if (v_lepton_idx.empty()) {
                // Add first selected muon
                v_lepton_idx.push_back(i);
                // Store its TLorentzVector in muons vector
                muons.push_back(pre_muons.at(i));
            } 
            else if (v_lepton_idx.size() == 1 && 
                     (*branchReader_.RequiredIntArray("Muon_charge"))[v_lepton_idx[0]] != (*branchReader_.RequiredIntArray("Muon_charge"))[i]) {
                // Add second selected muon with opposite charge
                v_lepton_idx.push_back(i);
                // Store its TLorentzVector in muons vector
                muons.push_back(pre_muons.at(i));
            }
        }
        // Assign final indices to v_muon_idx
        v_muon_idx = v_lepton_idx;
    }
    // Dielectron channel
    else if (TString(Decaymode).Contains("dielec")) {
        Int_t nel = elecs_pt->GetSize();
        for (int i = 0; i < nel; ++i) {
            // Skip electrons that don't pass selection criteria
            if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), elec_pt, elec_eta) ||
                !elecSCBId(branchReader_.GetIntArrayValue("Electron_cutBased", i), eleid_scbcut) ||
                (fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) < 1.566) ||
                !(elecCharge(branchReader_.GetIntArrayValue("Electron_tightCharge", i)))) {
                continue;
            }
            // Use same logic pattern as dimuon channel
            if (v_electron_idx.empty()) {
                // Add first selected electron
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i);
                // Store its TLorentzVector in electrons vector
                elecs.push_back(pre_elecs.at(i));
            } 
            else if (v_electron_idx.size() == 1 && 
                     (*branchReader_.RequiredIntArray("Electron_charge"))[v_electron_idx[0]] != (*branchReader_.RequiredIntArray("Electron_charge"))[i]) {
                // Add second selected electron with opposite charge
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i);
                // Store its TLorentzVector in electrons vector
                elecs.push_back(pre_elecs.at(i));
            }
        }
        v_electron_idx = v_lepton_idx;
    }
    // Muon-electron channel
    else if (TString(Decaymode).Contains("muel")) {
        // Process muons
        Int_t nmu = muons_iso->GetSize();
        for (int i = 0; i < nmu; ++i) {
            // Skip muons that don't pass selection criteria
            if (!passIsolation(muons_iso->At(i), muon_isocut) || 
                !passKinematicCuts(muons_pt->At(i), muons_eta->At(i), muon_pt, muon_eta) || 
                !passId(muons_Id->At(i))) {
                continue;
            }
    
            // Add selected muon
            v_muon_idx.push_back(i);
            // Store its TLorentzVector in muons vector
            muons.push_back(pre_muons.at(i));
        }

        // Process electrons
        Int_t nel = elecs_pt->GetSize();
        v_electron_idx.clear();

        // Only proceed if we have at least one muon
        if (!v_muon_idx.empty()) {
            for (int i = 0; i < nel; ++i) {
                // Skip electrons that don't pass selection criteria
                if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), elec_pt, elec_eta) ||
                    !elecSCBId(branchReader_.GetIntArrayValue("Electron_cutBased", i), eleid_scbcut) ||
                   (fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) > 1.4442 &&
                   fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) < 1.566) ||
                  !(elecCharge(branchReader_.GetIntArrayValue("Electron_tightCharge", i))) ||
                  !(*branchReader_.RequiredBoolArray("Electron_convVeto"))[i]) {
                    continue;
                }
        
                // Only select electrons with charge opposite to the first selected muon
                if ((*branchReader_.RequiredIntArray("Muon_charge"))[v_muon_idx[0]] != (*branchReader_.RequiredIntArray("Electron_charge"))[i]) {
                    // Add selected electron
                    v_electron_idx.push_back(i);
                    // Store its TLorentzVector in electrons vector
                    elecs.push_back(pre_elecs.at(i));
                }
            }
        }
    }
    else {
        logger_.Error() << "Lepton Selection error" << std::endl;
    }

    // Select veto leptons for jet cleaning & third lepton veto
    SelectVetoMuons();
    SelectVetoElectrons();
}


void Analysis::LeptonOrder() {
    Lep1.SetPxPyPzE(-999, -999, -999, -999);
    Lep2.SetPxPyPzE(-999, -999, -999, -999);
    Lep.SetPxPyPzE(-999, -999, -999, -999);
    AnLep.SetPxPyPzE(-999, -999, -999, -999);

    // Handle dimuon decay mode
    if (TString(Decaymode).Contains("dimuon")) {
        if (muons.size() >= 2) {
            Lep1 = muons[0];
            Lep2 = muons[1];

            // Set Lep & AnLep based on charge
            if ((*branchReader_.RequiredIntArray("Muon_charge"))[v_muon_idx[0]] < 0) {
                Lep = muons[0];
                AnLep = muons[1];
            } else {
                Lep = muons[1];
                AnLep = muons[0];
            }
        }
    }
    // Handle dielectron decay mode
    else if (TString(Decaymode).Contains("dielec")) {
        if (elecs.size() >= 2) {
            Lep1 = elecs[0];
            Lep2 = elecs[1];

            // Set Lep & AnLep based on charge
            if ((*branchReader_.RequiredIntArray("Electron_charge"))[v_electron_idx[0]] < 0) {
                Lep = elecs[0];
                AnLep = elecs[1];
            } else {
                Lep = elecs[1];
                AnLep = elecs[0];
            }
        }
    }
    // Handle muon-electron decay mode
    else if (TString(Decaymode).Contains("muel")) {
        if (muons.size() >= 1 && elecs.size() >= 1) {
            // Assign Lep1, Lep2 based on pT
            if (muons[0].Pt() > elecs[0].Pt()) {
                Lep1 = muons[0];
                Lep2 = elecs[0];
            } else {
                Lep1 = elecs[0];
                Lep2 = muons[0];
            }

            // Set Lep & AnLep based on muon charge
            if ((*branchReader_.RequiredIntArray("Muon_charge"))[v_muon_idx[0]] < 0) {
                Lep = muons[0];
                AnLep = elecs[0];
            } else {
                Lep = elecs[0];
                AnLep = muons[0];
            }
        }
    }
    // Handle invalid Decaymode
    else {
        logger_.Error() << "Lepton TLorentzVector Error: Decaymode = " << Decaymode << std::endl;
    }
}

void Analysis::MakeMuonCollection() {
    pre_muons.clear();

    for (int imu = 0; imu < muons_pt->GetSize(); ++imu) {
        TLorentzVector muon = createLorentzVector(muons_pt->At(imu), muons_eta->At(imu),
                                                 muons_phi->At(imu), muons_M->At(imu));

        // Apply Rochester correction if enabled
        if (applyRochester == "True") {
            double RoccoR = 1.0;

            // Check if Muon_charge exists
            auto muon_charge_it = branchReader_.intVectors.find("Muon_charge");
            if (muon_charge_it == branchReader_.intVectors.end() || !muon_charge_it->second) {
                logger_.Error() << "Muon_charge branch not available!" << std::endl;
                pre_muons.push_back(muon);
                continue;
            }

            int muon_charge = (*muon_charge_it->second)[imu];

            if (isData) {
                RoccoR = SSBCorr->RochesterCorrectionData(RunPeriod, muon_charge,
                                                        muon.Pt(), muon.Eta(), muon.Phi(), 0, 0);
            }
            else {
                // Muon_genPartIdx/Muon_nTrackerLayers are narrowed types in v15;
                // use branchReader_.GetIntArrayValue() rather than intVectors.
                auto GenPts = branchReader_.RequiredFloatArray("GenPart_pt");

                if (!branchReader_.BranchIsAvailable("Muon_genPartIdx") || !GenPts || !branchReader_.BranchIsAvailable("Muon_nTrackerLayers")) {
                    logger_.Error() << "Required MC branches for Rochester correction not available!" << std::endl;
                    logger_.Error() << "Muon_genPartIdx: " << (branchReader_.BranchIsAvailable("Muon_genPartIdx") ? "OK" : "NULL") << std::endl;
                    logger_.Error() << "GenPart_pt: " << (GenPts ? "OK" : "NULL") << std::endl;
                    logger_.Error() << "Muon_nTrackerLayers: " << (branchReader_.BranchIsAvailable("Muon_nTrackerLayers") ? "OK" : "NULL") << std::endl;
                    pre_muons.push_back(muon);
                    continue;
                }

                int GenID = static_cast<int>(branchReader_.GetIntArrayValue("Muon_genPartIdx", imu));
                double GenPt = (GenID >= 0) ? GenPts->At(GenID) : 0.0;
                int nLayers = static_cast<int>(branchReader_.GetIntArrayValue("Muon_nTrackerLayers", imu));

                RoccoR = SSBCorr->RochesterCorrectionMC(RunPeriod, muon_charge,
                                                      muon.Pt(), muon.Eta(), muon.Phi(),
                                                      GenID, GenPt, nLayers, 0, 0);
            }

            muon.SetPtEtaPhiM(muon.Pt() * RoccoR, muon.Eta(), muon.Phi(), muon.M());
        }

        pre_muons.push_back(muon);
    }
}


void Analysis::MakeElecCollection() {
    // Electron collection logic
    pre_elecs.clear();
    for (int iele = 0; iele < elecs_pt->GetSize(); ++iele){
        pre_elecs.push_back(createLorentzVector(elecs_pt->At(iele), elecs_eta->At(iele), elecs_phi->At(iele), elecs_M->At(iele) )); 
    }
 
    return;
}

void Analysis::MakeJetCollection() {
    pre_jets.clear();

    if (jets_pt == nullptr || jets_eta == nullptr || jets_phi == nullptr || jets_M == nullptr) {
        logger_.Error() << "Error: Some jet branch pointers are null in MakeJetCollection()" << std::endl;
        return;
    }

    Int_t nJets = jets_pt->GetSize();
    std::vector<TLorentzVector> rawJets;
    std::vector<float> rawFactors;
    std::vector<float> jetAreas;
    // Jet_muonSubtrFactor: needed (along with CorrT1METJet_*, read further
    // below) for a complete Type-1 MET recomputation - see
    // SSBCorrections::ApplyType1METWithCorrT1.
    std::vector<float> jetMuonSubtrFactors;
    std::vector<TLorentzVector> genJets;
    std::vector<int> genJetIndices;

    rawJets.reserve(nJets);
    rawFactors.reserve(nJets);
    jetAreas.reserve(nJets);
    jetMuonSubtrFactors.reserve(nJets);
    genJets.reserve(nJets);
    genJetIndices.reserve(nJets);

    // Every raw index gets an entry (real jet or, on a read failure below, a
    // dummy) so rawJets/rawFactors/jetAreas/... stay index-aligned with nJets.
    for (int ijet = 0; ijet < nJets; ++ijet) {
        try {
            TLorentzVector rawJet;
            rawJet.SetPtEtaPhiM(jets_pt->At(ijet), jets_eta->At(ijet), jets_phi->At(ijet), jets_M->At(ijet));
            rawJets.push_back(rawJet);

            float rawFactor = (branchReader_.OptionalFloatArray("Jet_rawFactor") && branchReader_.OptionalFloatArray("Jet_rawFactor")->GetSize() > ijet)
                                ? branchReader_.RequiredFloatArray("Jet_rawFactor")->At(ijet) : 0.0;
            rawFactors.push_back(rawFactor);

            float area = (branchReader_.OptionalFloatArray("Jet_area") && branchReader_.OptionalFloatArray("Jet_area")->GetSize() > ijet)
                           ? branchReader_.RequiredFloatArray("Jet_area")->At(ijet) : 0.5;
            jetAreas.push_back(area);

            float muonSubtr = (branchReader_.OptionalFloatArray("Jet_muonSubtrFactor") && branchReader_.OptionalFloatArray("Jet_muonSubtrFactor")->GetSize() > ijet)
                                ? branchReader_.RequiredFloatArray("Jet_muonSubtrFactor")->At(ijet) : 0.0;
            jetMuonSubtrFactors.push_back(muonSubtr);

            // Jet_genJetIdx is Short_t in NanoAODv15 (was Int_t in v9) - use the
            // version-agnostic accessor rather than assuming it's in branchReader_.intVectors.
            int genIdx = static_cast<int>(branchReader_.GetIntArrayValue("Jet_genJetIdx", ijet));
            genJetIndices.push_back(genIdx);
            
        } catch (const std::exception& e) {
            fallbackCounter_.RecordFallback("jet_raw_field_read");
            logger_.Error() << "Error reading jet at index " << ijet << ": " << e.what() << std::endl;
            logger_.Error() << "erro ! "  << std::endl;
            logger_.Error() << "Error reading jet at index " << ijet << ": " << e.what() << std::endl;
            logger_.Warning() << "Creating dummy jet at index " << ijet << std::endl;

            TLorentzVector dummyJet;
            dummyJet.SetPtEtaPhiM(-999.0, -999.0, -999.0, -999.0);  // Clearly invalid jet
            rawJets.push_back(dummyJet);
            rawFactors.push_back(0.0);
            jetAreas.push_back(0.5);
            jetMuonSubtrFactors.push_back(0.0);
            genJetIndices.push_back(-1);
        }
    }

    // Verify size consistency before proceeding
    if (rawJets.size() != static_cast<size_t>(nJets)) {
        logger_.Error() << "Size mismatch after jet reading: expected " << nJets
                  << ", got " << rawJets.size() << std::endl;
        pre_jets.clear();
        return;
    }

    if(!isData){
        if (gen_jets_pt && gen_jets_eta && gen_jets_phi && gen_jets_M) {
            Int_t nGenJets = gen_jets_pt->GetSize();
            for (int i = 0; i < nGenJets; ++i) {
                TLorentzVector gj;
                gj.SetPtEtaPhiM(gen_jets_pt->At(i), gen_jets_eta->At(i), gen_jets_phi->At(i), gen_jets_M->At(i));
                genJets.push_back(gj);
            }
        }
    }

    // fixedGridRhoFastjetAll -> Rho_fixedGridRhoFastjetAll and MET_pt/phi ->
    // PFMET_pt/phi were both renamed in NanoAODv15 - try both names.
    double rho = branchReader_.GetFloatSingleValueByAlias({"fixedGridRhoFastjetAll", "Rho_fixedGridRhoFastjetAll"});
    // Type-1 MET must start from raw (uncorrected) MET, not MET_pt/PFMET_pt/
    // PuppiMET_pt - those already have an old Type-1 correction baked in.
    double raw_met_pt; double raw_met_phi;
    if (METtype =="PF"){
        raw_met_pt  = branchReader_.GetFloatSingleValueByAlias({"RawMET_pt"});
        raw_met_phi = branchReader_.GetFloatSingleValueByAlias({"RawMET_phi"});
    }
    else if (METtype =="Puppi"){
        raw_met_pt  = branchReader_.GetFloatSingleValueByAlias({"RawPuppiMET_pt"});
        raw_met_phi = branchReader_.GetFloatSingleValueByAlias({"RawPuppiMET_phi"});
    }
    else {
	    logger_.Error() << "There is no MET type. Check out your Configuration! Default is PuppiMET " << std::endl;
        raw_met_pt  = branchReader_.GetFloatSingleValueByAlias({"RawPuppiMET_pt"});
        raw_met_phi = branchReader_.GetFloatSingleValueByAlias({"RawPuppiMET_phi"});
    }

    // GetFloatSingleValueByAlias silently returns 0.0 if not found, which
    // would look plausible but wrong here - so require Raw(Puppi)MET explicitly.
    bool haveRawMet = (METtype == "PF")
        ? (branchReader_.BranchIsAvailable("RawMET_pt") && branchReader_.BranchIsAvailable("RawMET_phi"))
        : (branchReader_.BranchIsAvailable("RawPuppiMET_pt") && branchReader_.BranchIsAvailable("RawPuppiMET_phi"));
    if (!haveRawMet) {
        logger_.Error() << "RawMET_pt/RawMET_phi (or RawPuppiMET_pt/RawPuppiMET_phi for Puppi) "
                  << "not available in this file - cannot compute the official Type-1 MET starting "
                  << "point (see NOTES.md, item 17)." << std::endl;
        throw std::runtime_error(
            "MakeJetCollection(): missing RawMET_pt/RawPuppiMET_pt - required as the Type-1 MET "
            "starting point (NanoAODv15).");
    }

    // Run number: needed as an extra evaluate() input for v15 Puppi DATA JEC.
    unsigned int run_number_for_jec = branchReader_.RequiredUIntValue("run");
    // Event number: seeds the deterministic per-jet JER smearing in SmearJER.
    unsigned long long event_number_for_jer = branchReader_.RequiredULongValue("event");
    // Builds jet pt+mass only; Type-1 MET is computed separately below.
    pre_jets = SSBCorr->ApplyJetCorrections(
        rawJets,
        rawFactors,
        jetAreas,
        rho,
        isData,
        dojes,
        dojer,
        genJets,
        genJetIndices,
        run_number_for_jec,
        event_number_for_jer,
        JERSys.Data()
    );

    // ApplyType1METWithCorrT1 is the sole MET calculation, matching the CMS
    // JERC tutorial's Applier::correctedMet() term-for-term.

    // CorrT1METJet_*/Jet_muonSubtrFactor are required for the official
    // Type-1 MET recipe; their absence is a hard error, not a silent fallback.
    if (branchReader_.BranchIsAvailable("CorrT1METJet_rawPt") &&
        branchReader_.BranchIsAvailable("Jet_muonSubtrFactor")) {
        auto *corrT1AreaArr      = branchReader_.OptionalFloatArray("CorrT1METJet_area");
        auto *corrT1EtaArr       = branchReader_.OptionalFloatArray("CorrT1METJet_eta");
        auto *corrT1PhiArr       = branchReader_.OptionalFloatArray("CorrT1METJet_phi");
        auto *corrT1RawPtArr     = branchReader_.RequiredFloatArray("CorrT1METJet_rawPt");
        auto *corrT1MuonSubtrArr = branchReader_.OptionalFloatArray("CorrT1METJet_muonSubtrFactor");

        Int_t nCorrT1 = corrT1RawPtArr ? corrT1RawPtArr->GetSize() : 0;
        std::vector<float> corrT1AreaVec, corrT1EtaVec, corrT1PhiVec, corrT1RawPtVec, corrT1MuonSubtrVec;
        corrT1AreaVec.reserve(nCorrT1); corrT1EtaVec.reserve(nCorrT1); corrT1PhiVec.reserve(nCorrT1);
        corrT1RawPtVec.reserve(nCorrT1); corrT1MuonSubtrVec.reserve(nCorrT1);
        for (int i = 0; i < nCorrT1; ++i) {
            corrT1AreaVec.push_back((corrT1AreaArr && corrT1AreaArr->GetSize() > i) ? corrT1AreaArr->At(i) : 0.5f);
            corrT1EtaVec.push_back((corrT1EtaArr && corrT1EtaArr->GetSize() > i) ? corrT1EtaArr->At(i) : 0.0f);
            corrT1PhiVec.push_back((corrT1PhiArr && corrT1PhiArr->GetSize() > i) ? corrT1PhiArr->At(i) : 0.0f);
            corrT1RawPtVec.push_back(corrT1RawPtArr->At(i));
            corrT1MuonSubtrVec.push_back((corrT1MuonSubtrArr && corrT1MuonSubtrArr->GetSize() > i) ? corrT1MuonSubtrArr->At(i) : 0.0f);
        }

        // Reused for the Type-1 MET selection cut (pt>15, |eta|<5.2, chEmEF+neEmEF<0.90).
        auto *jetChEmEFArr = branchReader_.OptionalFloatArray("Jet_chEmEF");
        auto *jetNeEmEFArr = branchReader_.OptionalFloatArray("Jet_neEmEF");
        std::vector<float> jetChEmEFVec, jetNeEmEFVec;
        jetChEmEFVec.reserve(nJets);
        jetNeEmEFVec.reserve(nJets);
        for (int i = 0; i < nJets; ++i) {
            jetChEmEFVec.push_back((jetChEmEFArr && jetChEmEFArr->GetSize() > i) ? jetChEmEFArr->At(i) : 0.0f);
            jetNeEmEFVec.push_back((jetNeEmEFArr && jetNeEmEFArr->GetSize() > i) ? jetNeEmEFArr->At(i) : 0.0f);
        }

        Met = SSBCorr->ApplyType1METWithCorrT1(
            raw_met_pt, raw_met_phi,
            rawJets, rawFactors, jetAreas, jetMuonSubtrFactors, jetChEmEFVec, jetNeEmEFVec,
            genJetIndices,
            corrT1RawPtVec, corrT1EtaVec, corrT1PhiVec, corrT1AreaVec, corrT1MuonSubtrVec,
            rho, isData, /*applyJES=*/dojes, /*applyJER=*/dojer,
            genJets, run_number_for_jec, event_number_for_jer,
            JERSys.Data()
        );
    } else {
        // Missing branches likely mean a wrong/older NanoAOD file was fed in.
        logger_.Error() << "CorrT1METJet_rawPt or Jet_muonSubtrFactor not available in this "
                  << "file - cannot compute the official Type-1 MET (see NOTES.md, item 16). "
                  << "This framework no longer falls back to the old simplified MET recipe."
                  << std::endl;
        throw std::runtime_error(
            "MakeJetCollection(): missing CorrT1METJet_rawPt or Jet_muonSubtrFactor - "
            "required for the official Type-1 MET recipe (NanoAODv15).");
    }

    TString yearForm;

    if (RunPeriod.Contains("2016Pre")) yearForm = "2016APV";
    else if (RunPeriod.Contains("2016Post")) yearForm = "2016nonAPV";
    else if (RunPeriod.Contains("2017")) yearForm = "2017";
    else if (RunPeriod.Contains("2018")) yearForm = "2018";
    else logger_.Warning() << "[Warning check the RunPeriod] : " << RunPeriod  << std::endl;
    bool isMC = !isData; bool isUL = true; bool isPuppi = METtype == "Puppi";


    if (applyMETXY == "True") Met = SSBCorr->METXYCorrection(Met, branchReader_.RequiredUIntValue("run"), yearForm, isMC, static_cast<int>(branchReader_.GetIntSingleValue("PV_npvsGood")), isUL, isPuppi);

    if (pre_jets.size() != static_cast<size_t>(nJets)) {
        logger_.Info() << "MakeJetCollection: Input jets=" << nJets
                  << ", Output corrected jets=" << pre_jets.size()
                  << " (some jets may have been filtered by JEC/JER)" << std::endl;
    }
}


bool Analysis::NumIsoLeptons(int nNLepsCut) // YOU SHOULD CALL THIS FUNCTION AFTER LEPTONSELETOR //
{
    bool numLeptons = true;
    
    if (TString(Decaymode).Contains("dimuon")) {
        // For dimuon channel, check using v_muon_idx size
        if (v_muon_idx.size() <= 1 || v_muon_idx.size() < nNLepsCut) {
            numLeptons = false;
        }
    }
    else if (TString(Decaymode).Contains("dielec")) {
        // For dielectron channel, check using v_electron_idx size
        if (v_electron_idx.size() <= 1 || v_electron_idx.size() < nNLepsCut) {
            numLeptons = false;
        }
    }
    else if (TString(Decaymode).Contains("muel")) {
        // For muon-electron channel, require at least one of each
        if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1 || 
            (v_muon_idx.size() + v_electron_idx.size()) < nNLepsCut) {
            numLeptons = false;
        }
    }
    else if (TString(Decaymode).Contains("muonJet")) {
        // For muon+jets channel, require exactly one lepton
        if (v_lepton_idx.size() != 1 || v_lepton_idx.size() < nNLepsCut) {
            numLeptons = false;
        }
    }
    else {
        logger_.Error() << "Error: Unrecognized decay mode in NumIsoLeptons()" << std::endl;
    }
    return numLeptons;
}

void Analysis::JetSelector() {
    // Pre-conditions check
    if (!object_variables_set_) {
        logger_.Error() << "JetSelector() called before SetObjectVariable()!" << std::endl;
        throw std::runtime_error("SetObjectVariable() must be called before JetSelector()");
    }

    // Reset state
    jets_selected_ = false;
    
    // Check if necessary pointers are initialized
    if (jets_pt == nullptr || jets_eta == nullptr || jets_phi == nullptr || jets_M == nullptr) {
        logger_.Error() << "Error: Basic jet variables (pt, eta, phi, M) are not initialized. Make sure SetObjectVariable() was called." << std::endl;
        return;
    }

    MakeJetCollection();
    CollectPUIDCandidates();  // for the weight calculation in PUIDSFApply()
    
    v_jet_idx.clear();
    jets.clear();
    
    // Lambda functions for common jet selection checks
    auto passKinematicCuts = [this](float pt, float eta) -> bool {
        return pt > jet_pt && fabs(eta) < jet_eta;
    };

    Int_t nJets = jets_pt->GetSize();

    for (int i = 0; i < nJets; i++) {
        TLorentzVector jetVec = pre_jets[i];
        float jetPt = jetVec.Pt();
        float jetEta = jetVec.Eta();
        
        // Apply kinematic cuts first
        if (!passKinematicCuts(jetPt, jetEta)) {
            continue;
        }

        // Check jet ID (NanoAODv15 Puppi jet ID, computed from energy fractions
        // since Jet_jetId no longer exists - see PassConfiguredJetId)
        if (!PassConfiguredJetId(i)) {
            continue;
        }

        // Jet cleaning
        if (!JetCleaning(&jetVec)) {
            continue;
        }

        // PUID selection only - the weight is computed separately in PUIDSFApply()
        if (apply_puid_ && jetPt <= puid_pt_threshold_) {
            int puId = (jets_puId != nullptr) ? jets_puId->At(i) : 0;
            bool passPUID = PassPileupID(jetPt, puId, puid_wp_);
            
            if (!passPUID) {
                continue;  // Skip this jet if PUID fails
            }
        }

	// Initialize jet veto event flag
        isjetveto_event_ = false;

        // Apply HEM15/16 veto based on configuration type
        bool should_apply_hem = false;
        if (RunPeriod.Contains("2018")) {
            if (isData) {
                if (branchReader_.BranchIsAvailable("run")) {
                    unsigned int run_num = branchReader_.RequiredUIntValue("run");
                    should_apply_hem = (run_num >= 319077);
                }
            } else {
                should_apply_hem = ((current_entry_ % 10000) < 6478);  // 64.78%
            }
        }

        // VetoChecker pre-selection requires (chEmEF+neEmEF)<0.90 before the veto map.
        double jetChEmEFForVeto = (branchReader_.OptionalFloatArray("Jet_chEmEF") && branchReader_.OptionalFloatArray("Jet_chEmEF")->GetSize() > i)
                                     ? branchReader_.RequiredFloatArray("Jet_chEmEF")->At(i) : 0.0;
        double jetNeEmEFForVeto = (branchReader_.OptionalFloatArray("Jet_neEmEF") && branchReader_.OptionalFloatArray("Jet_neEmEF")->GetSize() > i)
                                     ? branchReader_.RequiredFloatArray("Jet_neEmEF")->At(i) : 0.0;
        if (should_apply_hem && SSBCorr->ShouldVetoJet(jetVec, jetChEmEFForVeto, jetNeEmEFForVeto)) {
            if (SSBCorr->GetJetVetoType() == "jet") {
                continue; // Skip this jet
            } else { // "event" type
                isjetveto_event_ = true;
                return; // Exit JetSelector early
            }
        }
        // Add to selected jets
        v_jet_idx.push_back(i);
        jets.push_back(pre_jets[i]);
    }

    jets_selected_ = true;

    // NOTE: PUID weight calculation is now handled by PUIDSFApply() function
}


bool Analysis::JetCleaning(TLorentzVector* jet_)
{

    for (const auto& ele : elecsveto)
    {
        if (ele.DeltaR(*jet_) < 0.4)
            return false;
    }
    for (const auto& mu : muonsveto)
    {
        if (mu.DeltaR(*jet_) < 0.4)
            return false;
    }

    for (const auto& ele : elecs)
    {
        if (ele.DeltaR(*jet_) < 0.4)
            return false;
    }
    for (const auto& mu : muons)
    {
        if (mu.DeltaR(*jet_) < 0.4)
            return false;
    }
    return true;
}

TLorentzVector Analysis::JERSmearing(TLorentzVector* jet, int idx_, TString op_)
{
    double jerfrac_ = 1.0;

    if (op_.Contains("Norm")){
        //jerfrac_ = Jet_EnergyResolution_SF->at(idx_);
    }
    else if (op_ == "Up"){logger_.Debug() << "no op " << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFUp->at(idx_);
    }
    else if (op_ == "Down"){logger_.Debug() << "no op" << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFDown->at(idx_);
    }
    else{
        logger_.Warning() << "Check out your JERSmearing option !!" << "op_ : "<< op_<< std::endl;
    }
    if (dojer && !TString(FileName_).Contains("Data"))
    {
        return TLorentzVector(
            jet->Px() * jerfrac_,
            jet->Py() * jerfrac_,
            jet->Pz() * jerfrac_,
            jet->Energy() * jerfrac_
        );
    }

    return *jet;
}

bool Analysis::LeptonsPtAddtional()//YOU SHOULD REQUIRE THIS FUNCTION AFTER NumIsoLeptons //                                        
{     
    bool lepptadd = false;                                                                                                               
    if ( TString(Decaymode).Contains( "dimuon" ) ||
         TString(Decaymode).Contains( "dielec" ) ||                                                                                      
         TString(Decaymode).Contains( "muel" )      ){  if ( Lep1.Pt() > 25 && Lep2.Pt() > 20 ) {lepptadd=true;} } //Pt of Leading Lepton should be over than 25 GeV and Second Leading Lepton Pt should be over thand 20 GeV.
       
    else if ( TString(Decaymode).Contains( "muonJet" ) ){if( v_lepton_idx.size() == 1){lepptadd=true;} }
    else { logger_.Error() << "?? something wrong " << std::endl; }                                                                                      
    return lepptadd;
}
void Analysis::Start()
{
   //fout = new TFile(Form("output/%s",outfile),"RECREATE");
   if (strcmp(outdir.c_str(), "None") != 0 ) {
      fout = std::make_unique<TFile>(Form("gsidcap://cluster142.knu.ac.kr/%s/%s", outdir.c_str(), outfile.c_str()), "RECREATE");

   }
   else {
      //fout = new TFile(Form("output/%s",outfile),"RECREATE");
      fout = std::make_unique<TFile>(Form("output/%s", outfile.c_str()), "RECREATE");
   }
   logger_.Info() << "fout - getname : " << fout->GetName() << std::endl;
   fout->cd("");

   TDirectory *dir = gDirectory;
   dir->cd();

   DeclareHistos();
}

void Analysis::DeclareHistos()
{
 
   for (int i =0 ; i < 10 ; i++)
   {
      h_Lep1pt[i]  = new TH1D(Form("h_Lep1pt_%d" ,i), Form("Leading Lepton pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep1pt[i]->Sumw2(); 
      h_Lep2pt[i]  = new TH1D(Form("h_Lep2pt_%d" ,i), Form("Second Leading Lepton pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Lep2pt[i]->Sumw2();
      h_Lep1eta[i] = new TH1D(Form("h_Lep1eta_%d",i), Form("Leading Lepton Eta    %s"    ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep1eta[i]->Sumw2();
      h_Lep2eta[i] = new TH1D(Form("h_Lep2eta_%d",i), Form("Second Leading Lepton Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Lep2eta[i]->Sumw2();
      h_Lep1phi[i] = new TH1D(Form("h_Lep1phi_%d",i), Form("Leading Lepton Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Lep1phi[i]->Sumw2();
      h_Lep2phi[i] = new TH1D(Form("h_Lep2phi_%d",i), Form("Second Leading Lepton Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Lep2phi[i]->Sumw2();
   
      h_Jet1pt[i]  = new TH1D(Form("h_Jet1pt_%d", i), Form("Leading Jet pT %s"        ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet1pt[i]->Sumw2();
      h_Jet2pt[i]  = new TH1D(Form("h_Jet2pt_%d", i), Form("Second Leading Jet pT %s" ,cutflowName[i].Data()), 250, 0.0, 250); h_Jet2pt[i]->Sumw2();
      h_Jet1eta[i] = new TH1D(Form("h_Jet1eta_%d",i), Form("Leading Jet Eta %s"       ,cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet1eta[i]->Sumw2();
      h_Jet2eta[i] = new TH1D(Form("h_Jet2eta_%d",i), Form("Second Leading Jet Eta %s",cutflowName[i].Data()), 50, -2.5, 2.5); h_Jet2eta[i]->Sumw2();
      h_Jet1phi[i] = new TH1D(Form("h_Jet1phi_%d",i), Form("Leading Jet Phi %s"       ,cutflowName[i].Data()), 24, -1*pi, pi); h_Jet1phi[i]->Sumw2();
      h_Jet2phi[i] = new TH1D(Form("h_Jet2phi_%d",i), Form("Second Leading Jet Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_Jet2phi[i]->Sumw2();
    
     
      h_METpt[i]  = new TH1D(Form("h_METpt_%d",i), Form("MET pT %s" ,cutflowName[i].Data()), 200, 0.0, 200); h_METpt[i]->Sumw2();
      h_METphi[i] = new TH1D(Form("h_METphi_%d",i), Form("MET Phi %s",cutflowName[i].Data()), 24, -1*pi, pi); h_METphi[i]->Sumw2();

      h_DiLepMass[i] = new TH1D(Form("h_DiLepMass_%d",i),Form("Di-Lepton Invariant Mass %s",cutflowName[i].Data()), 300, 0.0, 300); h_DiLepMass[i]->Sumw2();
      h_Num_PV[i]    = new TH1D(Form("h_Num_PV_%d",i),     Form("Num of Primary Vertex after %s",cutflowName[i].Data()), 100, 0.0, 100); h_Num_PV[i]->Sumw2();
      h_Num_Jets[i]  = new TH1D(Form("h_Num_Jets_%d",i), Form("Num. of Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_Jets[i]->Sumw2();
      h_Num_bJets[i] = new TH1D(Form("h_Num_bJets_%d",i),Form("Num. of b Jets after %s",cutflowName[i].Data()), 20, 0.0, 20); h_Num_bJets[i]->Sumw2();
    }

    h_JetPUIDEvtWeight  = new TH1D(Form("h_JetPUIDEvtWeight"), Form("PUJetID Events Weight "), 200, 0.0, 2.0); h_JetPUIDEvtWeight->Sumw2(); 
    h_bTagEvtWeight  = new TH1D(Form("h_bTagEvtWeight"), Form("b-Tagging Events Weight "), 200, 0.0, 2.0); h_bTagEvtWeight->Sumw2(); 
    h_Top1Mass     = new TH1D(Form("h_Top1Mass"   ), Form("Top1 Mass"   ), 1000, 0.0, 1000); h_Top1Mass->Sumw2();
    h_Top1pt       = new TH1D(Form("h_Top1pt"   ), Form("Top1 pt"   ), 1000, 0.0, 1000); h_Top1pt->Sumw2();
    h_Top1Rapidity = new TH1D(Form("h_Top1Rapidity"   ), Form("Top1 Rapidity"   ), 100, -5, 5); h_Top1Rapidity->Sumw2();
    h_Top1phi      = new TH1D(Form("h_Top1phi"   ), Form("Top1 phi"   ), 24, -1*pi, pi); h_Top1phi->Sumw2();
    h_Top1Energy   = new TH1D(Form("h_Top1Energy"   ), Form("Top1 Energy"   ), 1000, 0.0, 1000); h_Top1Energy->Sumw2();
    h_Top2Mass     = new TH1D(Form("h_Top2Mass" ), Form("Top2 Mass" ), 1000, 0.0, 1000); h_Top2Mass->Sumw2();
    h_Top2pt       = new TH1D(Form("h_Top2pt"   ), Form("Top2 pt"   ), 1000, 0.0, 1000); h_Top2pt->Sumw2();
    h_Top2Rapidity = new TH1D(Form("h_Top2Rapidity"   ), Form("Top2 Rapidity"   ), 100, -5, 5); h_Top2Rapidity->Sumw2();
    h_Top2phi      = new TH1D(Form("h_Top2phi"   ), Form("Top2 phi"   ), 24, -1*pi, pi); h_Top2phi->Sumw2();
    h_Top2Energy = new TH1D(Form("h_Top2Energy"   ), Form("Top2 Energy"   ), 1000, 0.0, 1000); h_Top2Energy->Sumw2();
    
    h_TopMass       = new TH1D(Form("h_TopMass"   ), Form("Top Mass"   ), 1000, 0.0, 1000); h_TopMass->Sumw2();
    h_Toppt         = new TH1D(Form("h_Toppt"   ), Form("Top pt"   ), 1000, 0.0, 1000); h_Toppt->Sumw2();
    h_TopRapidity   = new TH1D(Form("h_TopRapidity"   ), Form("Top Rapidity"   ), 100, -5, 5); h_TopRapidity->Sumw2();
    h_Topphi        = new TH1D(Form("h_Topphi"   ), Form("Top phi"   ), 24, -1*pi, pi); h_Topphi->Sumw2();
    h_TopEnergy     = new TH1D(Form("h_TopEnergy"   ), Form("Top Energy"   ), 1000, 0.0, 1000); h_TopEnergy->Sumw2();
    h_AnTopMass     = new TH1D(Form("h_AnTopMass" ), Form("AnTop Mass" ), 1000, 0.0, 1000); h_AnTopMass->Sumw2();
    h_AnToppt       = new TH1D(Form("h_AnToppt"   ), Form("AnTop pt"   ), 1000, 0.0, 1000); h_AnToppt->Sumw2();
    h_AnTopRapidity = new TH1D(Form("h_AnTopRapidity"   ), Form("AnTop Rapidity"   ), 100, -5, 5); h_AnTopRapidity->Sumw2();
    h_AnTopphi      = new TH1D(Form("h_AnTopphi"   ), Form("AnTop phi"   ), 24, -1*pi, pi); h_AnTopphi->Sumw2();
    h_AnTopEnergy   = new TH1D(Form("h_AnTopEnergy"   ), Form("AnTop Energy"   ), 1000, 0.0, 1000); h_AnTopEnergy->Sumw2();

    h_W1Mass     = new TH1D(Form("h_W1Mass"  ), Form("W1 Mass" ), 300, 0.0, 300); h_W1Mass->Sumw2();
    h_W2Mass     = new TH1D(Form("h_W2Mass"  ), Form("W2 Mass" ), 300, 0.0, 300); h_W2Mass->Sumw2();
    
    h_W1Mt     = new TH1D(Form("h_W1Mt"  ), Form("W1 Transverse Mass" ), 300, 0.0, 300); h_W1Mt->Sumw2();
    h_W2Mt     = new TH1D(Form("h_W2Mt"  ), Form("W2 Transverse Mass" ), 300, 0.0, 300); h_W2Mt->Sumw2();
                      
    h_bJet1Energy = new TH1D(Form("h_bJet1Energy" ), Form("Leading bJet Energy"   ), 500, 0.0, 500); h_bJet1Energy->Sumw2();
    h_bJet2Energy = new TH1D(Form("h_bJet2Energy" ), Form("Second Leading bJet Energy" ), 500, 0.0, 500); h_bJet2Energy->Sumw2();
                      
    h_bJetEnergy   = new TH1D(Form("h_bJetEnergy" ), Form("bJet Energy"   ), 1000, 0.0, 1000); h_bJetEnergy->Sumw2();
    h_AnbJetEnergy = new TH1D(Form("h_AnbJetEnergy" ), Form("b-barJet Energy" ), 1000, 0.0, 1000); h_AnbJetEnergy->Sumw2();
    h_bJetPt   = new TH1D(Form("h_bJetPt" ), Form("bJet Pt"   ), 1000, 0.0, 1000); h_bJetPt->Sumw2();
    h_AnbJetPt = new TH1D(Form("h_AnbJetPt" ), Form("b-barJet Pt" ), 1000, 0.0, 1000); h_AnbJetPt->Sumw2();
                      
    h_Lep1Energy = new TH1D(Form("h_Lep1Energy" ), Form("Leading Lepton Energy"   ), 400, 0.0, 400); h_Lep1Energy->Sumw2();
    h_Lep2Energy = new TH1D(Form("h_Lep2Energy" ), Form("Second Leading Lepton Energy" ), 400, 0.0, 400); h_Lep2Energy->Sumw2();
                      
    h_LepEnergy   = new TH1D(Form("h_LepEnergy" ), Form("Lepton Energy"   ), 400, 0.0, 400); h_LepEnergy->Sumw2();
    h_AnLepEnergy = new TH1D(Form("h_AnLepEnergy" ), Form("Anti-Lepton Energy" ), 400, 0.0, 400); h_AnLepEnergy->Sumw2();
    
    h_Nu1Energy = new TH1D(Form("h_Nu1Energy" ), Form("Leading Nuetrino Energy"   ), 400, 0.0, 400); h_Nu1Energy->Sumw2();
    h_Nu2Energy = new TH1D(Form("h_Nu2Energy" ), Form("Second Leading Nuetrino Energy" ), 400, 0.0, 400); h_Nu2Energy->Sumw2();
    
    h_NuEnergy   = new TH1D(Form("h_NuEnergy" ), Form("Nuetrino Energy"   ), 400, 0.0, 400); h_NuEnergy->Sumw2();
    h_AnNuEnergy = new TH1D(Form("h_AnNuEnergy" ), Form("anti-Nuetrino Energy" ), 400, 0.0, 400); h_AnNuEnergy->Sumw2();

    for (int i =0; i < 13; ++i)
    {                 
        h_Reco_CPO_[i] = new TH1D(Form("h_Reco_CPO%d",i+1 ), Form("CPO%d",i+1   ), 200, -10, 10); h_Reco_CPO_[i]->Sumw2();
        h_Reco_CPO_ReRange_[i] = new TH1D(Form("h_Reco_CPO%d_ReRange",i+1 ), Form("CPO%d",i+1   ), 40, -2, 2); h_Reco_CPO_ReRange_[i]->Sumw2();
    }   
}


void Analysis::METDefiner()
{
    Met.SetPtEtaPhiM(static_cast<double>(**met_pt), 0, static_cast<double>(**met_phi), 0); // MET_phi MET_pt
}


void Analysis::JetOrder()
{
   Jet1.SetPtEtaPhiM(-999,-999,-999,-999);
   Jet2.SetPtEtaPhiM(-999,-999,-999,-999);
            
   if (jets.size() >=1){
      Jet1 = jets[0];
      if (jets.size() > 1){ Jet2 = jets[1]; } 
   }        
   return; 

}


bool Analysis::DiLeptonMassCut()
{        
   bool dimu_masscut = false;

   if ( ((Lep1)+(Lep2)).M() > 20 ){ dimu_masscut = true; }
         
   return dimu_masscut;
}        


void Analysis::SelectVetoMuons() {
    // Clear existing veto muon collection
    muonsveto.clear();
    v_vetomuon_idx.clear();
    
    // Lambda functions for common checks
    auto passKinematicCuts = [](float pt, float eta) {
        return pt > 20 && fabs(eta) < 2.4;
    };

    auto passIsolation = [](float iso, float isoCut) {
        return iso < isoCut;
    };

    auto passId = [](bool id) {
        return id;
    };
    
    // Channel-specific logic - same as before
    if (TString(Decaymode).Contains("dimuon")) {
        // For dimuon channel, check muons not in the selected pair
        Int_t nmu = muons_pt->GetSize();
        muonsveto.clear(); 
        for (int imu = 0; imu < nmu; ++imu) {
            // Skip if this muon is one of the selected muons
            if (std::find(v_muon_idx.begin(), v_muon_idx.end(), imu) != v_muon_idx.end()) {
                continue;
            }
            
            // Check criteria
            //if (passIsolation(muonsveto_iso->At(imu), veto_muoniso_cut) && 
            if (passIsolation(muons_iso->At(imu), muon_isocut) && 
                passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu)) && 
                passId(muonsveto_Id->At(imu))) {
                // Add to veto collection
                muonsveto.push_back(pre_muons.at(imu));
                v_vetomuon_idx.push_back(imu);
            }
        }
    }
    else if (TString(Decaymode).Contains("dielec")) {
        // For dielectron channel, check all muons
        Int_t nmu = muons_pt->GetSize();
        muonsveto.clear(); 
        for (int imu = 0; imu < nmu; ++imu) {
            // Check criteria
            if (passIsolation(muonsveto_iso->At(imu), veto_muoniso_cut) && 
                passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu)) && 
                passId(muonsveto_Id->At(imu))) {
                
                // Add to veto collection
                muonsveto.push_back(pre_muons.at(imu));
                v_vetomuon_idx.push_back(imu);
            }
        }
    }
    else if (TString(Decaymode).Contains("muel")) {
        muonsveto.clear(); 
        // For muon-electron channel, check muons not in the selected pair
        Int_t nmu = muons_pt->GetSize();
        muonsveto.clear();
        for (int imu = 0; imu < nmu; ++imu) {
            // Skip if this muon is the selected muon
            if (std::find(v_muon_idx.begin(), v_muon_idx.end(), imu) != v_muon_idx.end()) {
                continue;
            }
            
            // Check criteria
            if (passIsolation(muonsveto_iso->At(imu), veto_muoniso_cut) && 
                passKinematicCuts(muons_pt->At(imu), muons_eta->At(imu)) && 
                passId(muonsveto_Id->At(imu))) {
                
                // Add to veto collection
                muonsveto.push_back(pre_muons.at(imu));
                v_vetomuon_idx.push_back(imu);
            }
        }
    }
    else {
        logger_.Error() << "Error: Unrecognized decay mode in SelectVetoMuons!" << std::endl;
    }
    return;    
}



void Analysis::SelectVetoElectrons() {
    // Clear existing veto electron collection
    elecsveto.clear();
    
    // Lambda functions for common checks
    auto passKinematicCuts = [](float pt, float eta, float ptCut, float etaCut) {
        return pt >= ptCut && fabs(eta) <= etaCut;
    };

    auto passIsolation = [](float iso, float isoCut) {
        return iso <= isoCut;
    };

    auto elecSCBId = [](int id, int idcut) {
        return id >= idcut; 
    };
    
    auto elecCharge = [](int id) {
        return id > 0;
    };
    // Channel-specific logic
    if (TString(Decaymode).Contains("dimuon")) {
        // For dimuon channel, check all electrons
        //elecsveto.clear();    
        Int_t nel = elecs_pt->GetSize();
        for (int iel = 0; iel < nel; ++iel) {
            // Apply cuts
            if (!passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), elec_pt, elec_eta) ||
                !elecSCBId(branchReader_.GetIntArrayValue("Electron_cutBased", iel), elevetoid_scbcut) ||
                (fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[iel] + elecs_eta->At(iel)) > 1.4442 &&
                 fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[iel] + elecs_eta->At(iel)) < 1.566) ||
                !elecCharge(branchReader_.GetIntArrayValue("Electron_tightCharge", iel)) ||
                !(*branchReader_.RequiredBoolArray("Electron_convVeto"))[iel]
                ) {
                continue;
            }
            
            // Add to veto collection
            elecsveto.push_back(pre_elecs.at(iel));
        }
    }
    else if (TString(Decaymode).Contains("dielec")) {
        // For dielectron channel, check electrons not in the selected pair
        Int_t nel = elecs_pt->GetSize();
        for (int i = 0; i < nel; ++i) {
            // Skip if this electron is one of the selected electrons
            if (std::find(v_electron_idx.begin(), v_electron_idx.end(), i) != v_electron_idx.end()) {
                continue;
            }
            
            // Apply cuts (same logic as in LeptonSelector)
            if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), elec_pt, elec_eta) ||
                !elecSCBId(branchReader_.GetIntArrayValue("Electron_cutBased", i), elevetoid_scbcut) ||
                (fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) < 1.566) ||
                !elecCharge(branchReader_.GetIntArrayValue("Electron_tightCharge", i)) ||
                !(*branchReader_.RequiredBoolArray("Electron_convVeto"))[i]
                ) {
                continue;
            }

            // Add to veto collection
            elecsveto.push_back(pre_elecs.at(i));
        }
    }
    else if (TString(Decaymode).Contains("muel")) {
        // For muon-electron channel, check electrons not in the selected pair
        Int_t nel = elecs_pt->GetSize();
        for (int i = 0; i < nel; ++i) {
            // Skip if this electron is the selected electron
            if (std::find(v_electron_idx.begin(), v_electron_idx.end(), i) != v_electron_idx.end()) {
                continue;
            }
            
            // Apply cuts (same logic as in LeptonSelector)
            if (!passKinematicCuts(elecs_pt->At(i), elecs_eta->At(i), elec_pt, elec_eta) ||
                !elecSCBId(branchReader_.GetIntArrayValue("Electron_cutBased", i), elevetoid_scbcut) ||
                (fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[i] + elecs_eta->At(i)) < 1.566) ||
                !elecCharge(branchReader_.GetIntArrayValue("Electron_tightCharge", i)) ||
                !(*branchReader_.RequiredBoolArray("Electron_convVeto"))[i]) {
                continue;
            }
            
            // Add to veto collection
            elecsveto.push_back(pre_elecs.at(i));
        }
    }
    else {
        logger_.Error() << "Error: Unrecognized decay mode in SelectVetoElectrons!" << std::endl;
    }
}

bool Analysis::ThirdLeptonVeto()
{
    bool third_veto = true;

    if (TString(Decaymode).Contains("dimuon") || TString(Decaymode).Contains("dielec"))
    {
        // Check if we have at least 2 leptons
        if (TString(Decaymode).Contains("dimuon")) {
            // For dimuon channel
            if (v_muon_idx.size() <= 1) {
                third_veto = false;
            }

        }
        if (TString(Decaymode).Contains("dielec")) {
            // For dielectron channel
            if (v_electron_idx.size() <= 1) {
                third_veto = false;
            }

        }
        
        // Additional lepton veto - check the size of veto collections
        // If collections are not empty, fail the veto
        if (!muonsveto.empty() || !elecsveto.empty()) {
            third_veto = false;
        }
    }
    else if (TString(Decaymode).Contains("muel"))
    {
        // For muon-electron channel
        if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1) {
            third_veto = false;
        }
    //    else if ((*branchReader_.RequiredIntArray("Muon_charge"))[v_muon_idx.at(0)] == (*branchReader_.RequiredIntArray("Electron_charge"))[v_electron_idx.at(0)]) {
            // Require opposite sign
    //        third_veto = false;
    //    }
        
        // Additional lepton veto - check the size of veto collections
        // If collections are not empty, fail the veto
        if (!muonsveto.empty() || !elecsveto.empty()) {
            third_veto = false;
        }
    }
    else {
        logger_.Error() << "Error: Unrecognized decay mode in ThirdLeptonVeto function!" << std::endl;
    }
    
    return third_veto;
}
void Analysis::GenWeightApply()
{    
    double genweight = 1.0;
    if (!isData){
        if (branchReader_.RequiredFloatValue("Generator_weight") > 0.0){genweight =1;}
        else {genweight =-1;}
        evt_weight_ = evt_weight_*genweight;
    }  
    else {evt_weight_ = 1;}
}

void Analysis::bJetSelector() {
    // Clear any previous b-jet selections
    v_bjet_idx.clear();
    
    // Check if there are selected jets to work with
    if (v_jet_idx.empty() || jets_btag == nullptr) {
        return;
    }

    // Counter for number of b-tagged jets
    int nbtagged = 0;
            
    // Simple lambda function to determine if a jet is b-tagged
    auto isBTagged = [this](int jetIdx) -> bool {
        if (jets_btag != nullptr) {
            return (*jets_btag)[jetIdx] > bdisccut;
        }
        return false;
    };      
                
    // Loop over selected jets to find b-tagged ones
    for (const auto& jetIdx : v_jet_idx) {
        if (isBTagged(jetIdx)) {
            v_bjet_idx.push_back(jetIdx);
            nbtagged++;
        }
    }

    // Debug output
}

// ZVeto Cut : step 2
bool Analysis::ZVetoCut()
{        
   bool zvetocut = false;
         
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) )
   {        
      if ( ((Lep1)+(Lep2)).M() <= 76 || ((Lep1)+(Lep2)).M() >= 106 ){ zvetocut = true; }
   }        
   else if ( TString(Decaymode).Contains( "muel" ) ){ zvetocut = true; }
   else {logger_.Error() << "ZVeto Error !!" << std::endl;}
               
   return zvetocut;
}

// Num.Jet Cut : Step 3
bool Analysis::NumJetCut(std::vector<int> v_jets)
{           
   bool numjetcut = false;
   //if ( v_jets.size() == 2 ){ numjetcut = true; }
   if ( v_jets.size() >= 2 ){ numjetcut = true; }
   return numjetcut;
}

//MET Cut : step 4
bool Analysis::METCut(TLorentzVector met)        
{           
   bool metcut = false;
   if ( TString(Decaymode).Contains( "dimuon" ) || TString(Decaymode).Contains( "dielec" ) )
   {        
      if (met.Pt() > 40) { metcut =true; }
   }        
   else if ( TString(Decaymode).Contains( "muel" ) ){ metcut =true; }
   else {logger_.Error() << "METCut Error !!" << std::endl;}
   return metcut;
}

// Num.Jet Cut : Step 3
bool Analysis::NumbJetCut(std::vector<int> v_jets)
{           
   bool numbjetcut = false;
   //if ( v_jets.size() == 2 ){ numjetcut = true; }
   if ( v_jets.size() >= 1 ){ numbjetcut = true; }
   return numbjetcut;
}
void Analysis::SetUpKINObs()
{
   isKinSol=false;
   std::vector<double> jets_btag_vec;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   v_leptons_VLV.push_back(common::TLVtoLV(Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(AnLep));
   v_anlepidx_KIN.push_back(1);

   auto kinematicReconstruction = std::make_unique<KinematicReconstruction>(1, true);

   const LV met_LV = common::TLVtoLV(Met);
   
   // Create a map to translate from original jet indices to new indices in v_jets_VLV
   std::map<int, int> jet_idx_map;
   
   for (int i = 0; i < v_jet_idx.size(); ++i)
   {
      int idx_jet = v_jet_idx[i];
      v_jets_VLV.push_back(common::TLVtoLV(jets[i]));
      v_jetidx_KIN.push_back(i);  // Use sequential indices
      jets_btag_vec.push_back(static_cast<double>((*jets_btag)[idx_jet]));
      jet_idx_map[idx_jet] = i;  // Map original index to new index
   }

   // Translate b-jet indices to the new index system
   for (auto orig_idx : v_bjet_idx) {
      // Check if the original index exists in our jets vector
      if (jet_idx_map.find(orig_idx) != jet_idx_map.end()) {
         v_bjetidx_KIN.push_back(jet_idx_map[orig_idx]);
      }
   }
   

   // Only proceed if we have valid b-jet indices
   if (!v_bjetidx_KIN.empty() && v_jets_VLV.size() > 0) {
      KinematicReconstructionSolutions kinematicReconstructionSolutions = 
         kinematicReconstruction->solutions(v_lepidx_KIN, v_anlepidx_KIN, v_jetidx_KIN, 
                                         v_bjetidx_KIN, v_leptons_VLV, v_jets_VLV, 
                                         jets_btag_vec, met_LV);
                                         
      if (kinematicReconstructionSolutions.numberOfSolutions()) {
         isKinSol= true;
         LV top1 = kinematicReconstructionSolutions.solution().top();
         LV top2 = kinematicReconstructionSolutions.solution().antiTop();
         LV bjet1 = kinematicReconstructionSolutions.solution().bjet();
         LV bjet2 = kinematicReconstructionSolutions.solution().antiBjet();
         LV neutrino1 = kinematicReconstructionSolutions.solution().neutrino();
         LV neutrino2 = kinematicReconstructionSolutions.solution().antiNeutrino();
         
         Top       = common::LVtoTLV(top1);
         AnTop     = common::LVtoTLV(top2);
         bJet      = common::LVtoTLV(bjet1);
         AnbJet    = common::LVtoTLV(bjet2);
         Nu        = common::LVtoTLV(neutrino1);
         AnNu      = common::LVtoTLV(neutrino2);

         W1        = Lep  + AnNu;
         W2        = AnLep  + Nu;
      }
   } else {
      logger_.Warning() << "Not enough b-jets or jets for kinematic reconstruction" << std::endl;
   }
}

void Analysis::LeptonSFApply()
{
    lep_sf = 1.0;
    if (isData){
        evt_weight_ = 1.0;
        return;
    }
    else {
        if (TString(Decaymode).Contains("dimuon")) {
            lep_sf = SSBCorr->DoubleMuon_IDIsoEff(Lep1, Lep2, LepIdSFSys, LepIsoSFSys, LepTrackSFSys);
        }
        else if (TString(Decaymode).Contains("dielec")) {
            lep_sf = SSBCorr->DoubleElec_Eff(Lep1, Lep2,
                                            (*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[v_electron_idx[0]] + elecs_eta->At(v_electron_idx[0]),
                                            (*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[v_electron_idx[1]] + elecs_eta->At(v_electron_idx[1]),
                                            ElecId.Data(),LepIdSFSys.Data(), LepRecoSFSys.Data()); // LepRecoSFSys is for Electron
        }
        else if (TString(Decaymode).Contains("muel")) {
            // Calculate supercluster eta for electron
            double electron_sueta = (*branchReader_.RequiredFloatArray("Electron_deltaEtaSC"))[v_electron_idx[0]] + elecs_eta->At(v_electron_idx[0]);

            // Apply MuonElec_Eff function
            lep_sf = SSBCorr->MuonElec_Eff(Lep1, Lep2,  // Assuming Lep1=muon, Lep2=electron based on LeptonOrder()
                                          Lep1.Eta(),   // muon eta
                                          electron_sueta, // electron supercluster eta
                                          LepIdSFSys.Data(),   // muon ID systematic
                                          LepIsoSFSys.Data(),  // muon Iso systematic
                                          ElecId.Data(),       // electron ID working point
                                          LepIdSFSys.Data(),   // electron ID systematic (reuse muon ID sys)
                                          LepRecoSFSys.Data());// electron reco systematic
        }
        else {
           lep_sf = 1.0;
           logger_.Warning() << "Unknown decay mode: " << Decaymode
           << ", setting lep_sf = 1.0" << std::endl;

        }
        evt_weight_beforeLepsf_ = evt_weight_;
        evt_weight_ = lep_sf*evt_weight_;
    }
}

void Analysis::PUWeightApply()
{
     
    evt_weight_beforePileup_ = 1;
    evt_weight_beforePileup_ = evt_weight_; // keep event weight // 
    double puweight_ = 1.;
       
    if ( !TString(FileName_).Contains( "Data") )
    {  
       double pu_weight_central = SSBCorr->GetPUWeight( branchReader_.RequiredFloatValue("Pileup_nTrueInt") , PileUpSys.Data() );
       if (TString(PileUpSys).Contains("central") || TString(PileUpSys).Contains("nominal")  ) { puweight_   = pu_weight_central;}
       else {  
          logger_.Warning() << "PUWeightApply Error... Defalut is Weight_PileUp ... : " << PileUpSys << std::endl;
       }
       evt_weight_ = evt_weight_*puweight_;  // apply PileUpReweight //
    }  
    else {evt_weight_ = 1;}
}

void Analysis::TriggerSFApply()
{        
    double triggersf_ = 1.0;

    evt_weight_beforeTrigger_ = evt_weight_;  // Store weight before applying Trigger SF

    if (!TString(FileName_).Contains("Data")) {
        // Apply MC trigger scale factors based on decay mode
        if (TString(Decaymode).Contains("dielec")) {
            triggersf_ = SSBCorr->TrigDiElec_Eff(Lep1, Lep2, TrigSFSys);
        }
        else if (TString(Decaymode).Contains("muel")) {
            //triggersf_ = SSBCorr->TrigMuElec_Eff(Lep1, Lep2, TrigSFSys);
            triggersf_ = SSBCorr->TrigMuElec_Eff(elecs[0], muons[0], TrigSFSys);
        }
        else if (TString(Decaymode).Contains("dimu")) {
            triggersf_ = SSBCorr->TrigDiMuon_Eff(Lep1, Lep2, TrigSFSys);
        }
        else {
            logger_.Warning() << "[TriggerSFApply] Unknown Decaymode = " << Decaymode << std::endl;
        }

        evt_weight_ *= triggersf_;  // Apply the trigger scale factor
    }
    // For Data, no correction is applied
}

// PUID related functions implementation
bool Analysis::PassPileupID(float pt, int puId, const std::string& wp) const {
    // High-pT jets don't need PUID
    if (pt > 50.0) return true;

    // Puppi jets suppress pileup per-particle, so CMS ships no discrete PU-jet-ID
    // WP for them. Detect via Jet_puId availability rather than RunPeriod, since
    // v15 ships Puppi jets for older run periods too.
    if (!branchReader_.JetPuIdAvailable()) {
        return true;  // Puppi jets (or any file without Jet_puId): no PU-jet-ID to apply
    }

    // For 2016 - special handling due to potential bugs in NanoAODv9
    if (RunPeriod.Contains("2016")) {
        // Conservative approach: use bitmask check for 2016
        if (wp == "L") return (puId & 4) != 0;  // bit 2
        if (wp == "M") return (puId & 2) != 0;  // bit 1
        if (wp == "T") return (puId & 1) != 0;  // bit 0
        return false;
    }

    // For 2017, 2018 - use documented formula
    // puId = passlooseID*4 + passmediumID*2 + passtightID*1
    if (wp == "L") return puId >= 4;  // Loose: need at least 4 (100, 110, 111)
    if (wp == "M") return puId >= 6;  // Medium: need at least 6 (110, 111)
    if (wp == "T") return puId == 7;  // Tight: need exactly 7 (111)

    return false;
}

void Analysis::ApplyJetPUIDEventWeights() {
    // This function is now called automatically at the end of JetSelector()
    // Check if already applied
    if (jet_puid_weight_applied_) {
        return;  // Already applied in JetSelector()
    }

    logger_.Warning() << "ApplyJetPUIDEventWeights() called but PUID weights should be applied in JetSelector()!" << std::endl;
}

void Analysis::BTaggingSFApply() {
    evt_weight_beforeBtag_ = evt_weight_; // Store weight before btag SF
    btag_sf_weight_ = 1.0;
    
    // No b-tagging SF applied for Data
    if (isData) {
        return;
    }
    
    // Check if jets and btag info are available
    if (v_jet_idx.empty() || jets_btag == nullptr) {
        logger_.Info() << "No jets selected or btag info unavailable - B-tagging SF = 1.0" << std::endl;
        return;
    }
    
    // Prepare vectors for SF calculation
    std::vector<float> jet_pts, jet_etas;
    std::vector<int> jet_flavors;
    std::vector<bool> jet_isTagged;
    
    jet_pts.reserve(v_jet_idx.size());
    jet_etas.reserve(v_jet_idx.size());
    jet_flavors.reserve(v_jet_idx.size());
    jet_isTagged.reserve(v_jet_idx.size());
    
    // Jet_hadronFlavour is UChar_t in v15 (was Int_t in v9); use the version-agnostic accessor.
    bool hasHadronFlavour = branchReader_.BranchIsAvailable("Jet_hadronFlavour");

    if (!hasHadronFlavour) {
        logger_.Warning() << "Jet_hadronFlavour branch not available - using flavor=0 for all jets" << std::endl;
    }

    // Collect jet information
    for (size_t i = 0; i < v_jet_idx.size(); ++i) {
        int jet_idx = v_jet_idx[i];

        float jet_pt = jets[i].Pt();
        float jet_eta = jets[i].Eta();

        // Get jet flavor (MC only)
        int flavor = 0; // Default to light flavor
        if (hasHadronFlavour) {
            flavor = static_cast<int>(branchReader_.GetIntArrayValue("Jet_hadronFlavour", jet_idx));
            if (flavor < 0) flavor = 0; // branchReader_.GetIntArrayValue returns -999 if out of range
        }
        
        bool isTagged = (jets_btag->At(jet_idx) > bdisccut);
        
        // Fill vectors
        jet_pts.push_back(jet_pt);
        jet_etas.push_back(jet_eta);
        jet_flavors.push_back(flavor);
        jet_isTagged.push_back(isTagged);
        
    }
    

    // Set systematic variation - use "central" instead of "nominal"
    std::string syst_variation = "central";
    if (TString(BTagSFSys).Contains("up", TString::kIgnoreCase)) {
        syst_variation = "up";
    } else if (TString(BTagSFSys).Contains("down", TString::kIgnoreCase)) {
        syst_variation = "down";
    }
    
    // Calculate B-tagging event weight
    try {
        // BTagAlgoToString() is the canonical enum->string conversion, shared
        // with the constructor's parsing so the two can't drift apart.
        btag_sf_weight_ = SSBCorr->ComputeBTagEventWeight(
            jet_pts, jet_etas, jet_flavors, jet_isTagged,
            BTagAlgoToString(btag_algo_), btag_wp_, syst_variation
        );
        
        // Apply weight
        evt_weight_ *= btag_sf_weight_;
        
    } catch (const std::exception& e) {
        fallbackCounter_.RecordFallback("btag_sf");
        logger_.Error() << "Error in BTaggingSFApply: " << e.what() << std::endl;
        btag_sf_weight_ = 1.0; // Default on error
    }
    
    FillHisto(h_bTagEvtWeight, btag_sf_weight_);
}

bool Analysis::IsHardScatterJet(int jet_idx) const {
    // Jet_genJetIdx is Short_t in v15 (was Int_t in v9); version-agnostic accessor.
    if (!branchReader_.BranchIsAvailable("Jet_genJetIdx")) {
        return false;  // No gen matching info = consider as PileUp
    }
    int gen_jet_idx = static_cast<int>(branchReader_.GetIntArrayValue("Jet_genJetIdx", jet_idx));
    return gen_jet_idx >= 0;  // Gen matched = HardScatter, otherwise PileUp
}


void Analysis::CollectPUIDCandidates() {
    puid_hardscatter_jets_.clear();
    
    // Skip collection for Data or when PUID is disabled
    if (isData || !apply_puid_) {
        return;
    }
    
    // Check if jet collections are ready
    if (pre_jets.empty()) {
        // Normal case: event has no jets
        return;
    }
    
    if (jets_pt == nullptr) {
        logger_.Warning() << "CollectPUIDCandidates: jets_pt is nullptr" << std::endl;
        return;
    }
    
    Int_t nJets_original = jets_pt->GetSize();
    Int_t nJets_corrected = static_cast<Int_t>(pre_jets.size());
    
    // Use the smaller of the two sizes so indexing below never runs past
    // either array's end.
    Int_t nJets_safe = std::min(nJets_original, nJets_corrected);
    
    if (nJets_original != nJets_corrected) {
        logger_.Info() << "CollectPUIDCandidates: Size mismatch detected - "
                  << "original=" << nJets_original << ", corrected=" << nJets_corrected 
                  << ", using=" << nJets_safe << std::endl;
    }
    
    // Lambda for basic jet cuts (same as JetSelector)
    auto passBasicJetCuts = [this](int i, float jetPt, float jetEta) -> bool {
        // Kinematic cuts
        if (jetPt <= jet_pt || fabs(jetEta) >= jet_eta) return false;
        
        // Jet ID check (NanoAODv15 Puppi jet ID, see PassConfiguredJetId)
        if (!PassConfiguredJetId(i)) return false;
        
        // Jet cleaning (need TLorentzVector)
        if (i >= static_cast<int>(pre_jets.size())) return false;
        TLorentzVector jetVec = pre_jets[i];
        if (!JetCleaning(&jetVec)) return false;
        
        return true;
    };
    
    for (int i = 0; i < nJets_safe; i++) {
        // Skip invalid dummy jets (pT = -999 from error handling)
        if (pre_jets[i].Pt() < 0.0) continue;
        
        float jetPt = pre_jets[i].Pt();
        float jetEta = pre_jets[i].Eta();
        
        // Only consider jets that need PUID (CMS criteria)
        if (jetPt > puid_pt_threshold_) continue;
        
        // Skip Run3 PUPPI jets (no PUID needed)
        if (RunPeriod.Contains("2022") || RunPeriod.Contains("2023")) continue;
        
        // Must pass basic jet cuts to be PUID candidate
        if (!passBasicJetCuts(i, jetPt, jetEta)) continue;
        
        // Check if this is HardScatter jet (CMS criteria)
        bool is_hardscatter = IsHardScatterJet(i);
        if (!is_hardscatter) continue;  // Only collect HardScatter jets for weight calculation
        
        // Get PUID result with bounds checking
        int puId = 0;
        if (jets_puId != nullptr && i < jets_puId->GetSize()) {
            puId = jets_puId->At(i);
        }
        bool passes_puid = PassPileupID(jetPt, puId, puid_wp_);
        
        // Store information
        puid_hardscatter_jets_.emplace_back(i, jetPt, jetEta, passes_puid, is_hardscatter);
    }
}


void Analysis::PUIDSFApply() {
    evt_weight_beforePUID_ = evt_weight_;  // Store weight before PUID SF
    puid_sf_weight_ = 1.0;
    
    // Skip PUID SF for Data or when PUID is disabled
    if (isData || !apply_puid_) {
        return;
    }
    
    // Skip if no HardScatter candidate jets collected
    if (puid_hardscatter_jets_.empty()) {
        return;
    }
    
    // CMS Event Reweighting formula: P(DATA) / P(MC)
    float p_mc = 1.0;
    float p_data = 1.0;
    //          << " HardScatter jets" << std::endl;
    
    for (const auto& jet_info : puid_hardscatter_jets_) {
        try {
            // Get scale factor and efficiency using SSBCorrections
            float eff = SSBCorr->GetPUJetIDSFAndEff(
                jet_info.pt, jet_info.eta, jet_info.passes_puid, 
                true, puid_wp_, "MCEff", true);
                
            // Use systematic variation for scale factor
            std::string syst_variation = PUIDSFSys;
            if (syst_variation == "Central") syst_variation = "nominal";  // Handle legacy naming
            
            float sf = SSBCorr->GetPUJetIDSFAndEff(
                jet_info.pt, jet_info.eta, jet_info.passes_puid, 
                true, puid_wp_, syst_variation, false);
            
            // Validate efficiency values
            if (eff < 1e-5) eff = 1e-5;
            if (eff > 1.0 - 1e-5) eff = 1.0 - 1e-5;
            
            // Apply CMS Event Reweighting formula
            if (jet_info.passes_puid) {
                // Tagged: εᵢ vs SFᵢεᵢ
                p_mc *= eff;
                p_data *= (sf * eff);
            } else {
                // Not tagged: (1-εⱼ) vs (1-SFⱼεⱼ)  
                p_mc *= (1.0 - eff);
                p_data *= (1.0 - sf * eff);
            }
            
            
        } catch (const std::exception& e) {
            fallbackCounter_.RecordFallback("puid_sf");
            logger_.Error() << "PUID SF calculation failed for jet "
                      << jet_info.original_index << ": " << e.what() << std::endl;
            continue;  // Skip this jet on error
        }
    }
    
    // Calculate final event weight
    puid_sf_weight_ = (p_mc > 1e-10) ? p_data / p_mc : 1.0;
    
    // Apply weight to event
    evt_weight_ *= puid_sf_weight_;
    FillHisto(h_JetPUIDEvtWeight, puid_sf_weight_); 
}
