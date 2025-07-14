#include "../interface/Analysis.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// Constructor: initialize TTreeReader with TChain and branch list file
Analysis::Analysis(TChain *inputChain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt= -1)
    : chain(inputChain), fReader(inputChain), NumEvt(NumEvt), outdir(seDirName), outfile(outputName){
    if (!chain) {
        throw std::runtime_error("Error: Invalid TChain pointer!");
    }
    FileName_ = SetInputFileName(inputName);
    isData = TString(FileName_).Contains("Data");
    std::cout << "FileName_ : " << FileName_ << std::endl;
    // Load Configuration files //
    std::cout << "configFile : " << configFile << std::endl;  
    std::string confDir = "./configs/";
    std::string confpath = "";
    confpath = confDir+configFile;
    SSBConfReader = new TextReader();
    SSBConfReader->ReadFile(confpath);
    SSBConfReader->ReadVariables();
    SSBConfReader->PrintoutVariables();
    SSBCorr = new SSBCorrections(SSBConfReader, FileName_.Data());
    SSBCPVCal = new SSBCPVCalc();
    // Initialize branches based on branch list file
    InitBranches(branchListFile);
    cutflowName[0] = "Step_0";
    cutflowName[1] = "Step_1" ;
    cutflowName[2] = "Step_2";
    cutflowName[3] = "Step_3";
    cutflowName[4] = "Step_4";
    cutflowName[5] = "bTagged Jet >= 1";
    cutflowName[6] = "bTagged Jet >= 2";
    cutflowName[7] = "bTagged Jet == 2";
    cutflowName[8] = "Top-Recon.";
    cutflowName[9] = "Top-Pt-Rewight";
 
    pi = TMath::Pi();
    Start();
}

Analysis::~Analysis() {
    // Safely delete the TFile object

    if (fout) {
        fout->Write();
        std::cout << "[Info. output file Name] : " <<  fout->GetName() << std::endl;
        fout->Close(); // Close the file before deleting
        delete fout;
        fout = nullptr;
        std::cout << "fout successfully deleted." << std::endl;
    }

    // Safely delete the TextReader object
    if (SSBConfReader) {
        delete SSBConfReader;
        SSBConfReader = nullptr;
        std::cout << "SSBConfReader successfully deleted." << std::endl;
    }

    // Safely delete the TextReader object
    if (SSBCorr) {
        delete SSBCorr;
        SSBCorr = nullptr;
        std::cout << "SSBCorr successfully deleted." << std::endl;
    }
/*    else {
    std::cout << "Error !! SSBConfReader " << std::endl;
    }*/
    if (SSBCPVCal) {
        delete SSBCPVCal;
        SSBCPVCal = nullptr;
        std::cout << "SSBCPVCal successfully deleted." << std::endl;
    }
    // Add other cleanup as needed for dynamically allocated objects
    std::cout << "Analysis destructor completed." << std::endl;
}

void Analysis::InitBranches(const std::string &branchListFile) {
    std::cout << "branchListFile : " << branchListFile << std::endl;
    std::ifstream infile(branchListFile);
    if (!infile) {
        throw std::runtime_error("Error: Could not open branch list file " + branchListFile);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string branchName, objectType, dataType, varType;

        // Parse the line for branch name, object type, data type, and variable type (vector/single)
        if (std::getline(iss, branchName, ',') &&
            std::getline(iss, objectType, ',') &&
            std::getline(iss, dataType, ',') &&
            std::getline(iss, varType, ',')) {

            // Trim whitespace from each part (both ends)
            auto trim = [](std::string &s) {
                s.erase(s.find_last_not_of(" \t\n\r") + 1);
                s.erase(0, s.find_first_not_of(" \t\n\r"));
            };

            trim(branchName);
            trim(objectType);
            trim(dataType);
            trim(varType);
      	    if (isData && branchName.find("Gen") != std::string::npos) {continue;}
	          if (isData && branchName.find("gen") != std::string::npos) {continue;}
	          if (isData && branchName.find("Pileup_nTrueInt") != std::string::npos) {continue;}
	          if (isData && branchName.find("Jet_hadronFlavour") != std::string::npos) {continue;}
	          if (isData && branchName.find("L1PreFiringWeight") != std::string::npos) {continue;}
            // Initialize based on data type and variable type
            if (dataType == "Bool_t" && varType == "single") {
                boolSingles[branchName] = std::make_unique<TTreeReaderValue<Bool_t>>(fReader, branchName.c_str());
            } else if (dataType == "Int_t" && varType == "single") {
                intSingles[branchName] = std::make_unique<TTreeReaderValue<Int_t>>(fReader, branchName.c_str());
            } else if (dataType == "UInt_t" && varType == "single") {
                uintSingles[branchName] = std::make_unique<TTreeReaderValue<UInt_t>>(fReader, branchName.c_str());
            } else if (dataType == "Float_t" && varType == "single") {
                floatSingles[branchName] = std::make_unique<TTreeReaderValue<Float_t>>(fReader, branchName.c_str());
            } else if (dataType == "UChar_t" && varType == "single") { // New case for UChar_t single value
                ucharSingles[branchName] = std::make_unique<TTreeReaderValue<UChar_t>>(fReader, branchName.c_str());
            } else if (dataType == "Bool_t" && varType == "vector") {
                boolVectors[branchName] = std::make_unique<TTreeReaderArray<Bool_t>>(fReader, branchName.c_str());
            } else if (dataType == "Int_t" && varType == "vector") {
                intVectors[branchName] = std::make_unique<TTreeReaderArray<Int_t>>(fReader, branchName.c_str());
            } else if (dataType == "UInt_t" && varType == "vector") {
                uintVectors[branchName] = std::make_unique<TTreeReaderArray<UInt_t>>(fReader, branchName.c_str());
            } else if (dataType == "Float_t" && varType == "vector") {
                floatVectors[branchName] = std::make_unique<TTreeReaderArray<Float_t>>(fReader, branchName.c_str());
            } else if (dataType == "UChar_t" && varType == "vector") { // New case for UChar_t vector
                ucharVectors[branchName] = std::make_unique<TTreeReaderArray<UChar_t>>(fReader, branchName.c_str());
            } else {
                std::cerr << "Warning: Unsupported data type or format in branch list: "
                          << dataType << ", " << varType << " for branch " << branchName << std::endl;
            }

        } else {
            std::cerr << "Warning: Invalid branch format in list: " << line << std::endl;
        }
    }
}



void Analysis::SetVariables() {
    std::cout << "Set varibles " << std::endl;

    Lumi = SSBConfReader->GetNumber( "Luminosity" );
    RunPeriod = SSBConfReader->GetText( "RunRange" );
    Decaymode = SSBConfReader->GetText( "Channel" ); // Channel
    XsecTable_ = SSBConfReader->GetText( "XSecTablesName" );
    /// Set Trigger List ///
    num_dleptrig = SSBConfReader->Size( "dileptrigger" );
    num_sleptrig = SSBConfReader->Size( "singleleptrigger" );

    ///

    for(int i =0; i < num_dleptrig; ++i)
    {
        std::cout << SSBConfReader->GetText("dileptrigger",i+1) << std::endl;
        std::string tmptrg = SSBConfReader->GetText("dileptrigger",i+1);
        DLtrigName.push_back( removeSubstring( tmptrg, "_v") );
        trigName.push_back( removeSubstring( tmptrg, "_v")  );
    }

    for(int i =0; i < num_sleptrig; ++i)
    {
        std::cout << SSBConfReader->GetText("singleleptrigger",i+1) << std::endl;
        std::string tmptrg = SSBConfReader->GetText("singleleptrigger",i+1);
        SLtrigName.push_back( removeSubstring(tmptrg,"_v") );
        trigName.push_back( removeSubstring(tmptrg,"_v") );
    }

    for (int i = 0; i < trigName.size(); ++i){
        std::cout << "test trigName[i] " << trigName[i] << std::endl;
        triggerList[trigName[i]] =DeepCopy<bool>( boolSingles[trigName[i]]);
    }

    /// Set Noise filter(MET, events filter) ///
    for (int i = 0; i < SSBConfReader->Size("METFilters"); ++i){
       //std::cout << "METFilters: " << SSBConfReader->GetText("METFilters",i+1) << std::endl;
       std::string tmpnoisefl = SSBConfReader->GetText("METFilters",i+1);
       //std::cout << "METFilters: " <<  tmpnoisefl << std::endl;
       //std::cout << boolSingles[tmpnoisefl] << std::endl;
       if(boolSingles[tmpnoisefl] ==NULL) {std::cout << "Error!!!" << tmpnoisefl << std::endl;}
       noiseFilters[tmpnoisefl] = DeepCopy<bool>( boolSingles[tmpnoisefl]);
       //std::cout << "test "<< std::endl;
    }
    // Kinematic cut variables for Object 
    muon_pt     = SSBConfReader->GetNumber( "MuonPt_cut"     );// Muon pT Cut 
    muon_eta    = SSBConfReader->GetNumber( "MuonEta_cut"    );// Muon Eta Cut 
    muon_isocut  = SSBConfReader->GetNumber( "MuonIso_cut" );// 

    elec_pt     = SSBConfReader->GetNumber( "ElecPt_cut"     );
    elec_eta    = SSBConfReader->GetNumber( "ElecEta_cut"    );
    elec_isocut  = SSBConfReader->GetNumber( "ElecIso_cut" );


    // Muon Infor. ID  ISO // 
    MuonIsoType = SSBConfReader->GetText( "MuonIso_type" );
    MuonId      = SSBConfReader->GetText( "Muon_ID"      );
    
    // Electron ID ISO // 
    ElecIsoType = SSBConfReader->GetText( "ElecIso_type" );
    ElecId      = SSBConfReader->GetText( "Elec_ID"      );  
    // Jet Infor
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

    dojer = SSBConfReader->GetBool("DoJER");
    
    PileUpSys    = SSBConfReader->GetText("PileupSys");
    L1PreFireSys = SSBConfReader->GetText("L1PreFireSys");
    TrigSFSys    = SSBConfReader->GetText("TrigSFSys");


    // ============================================================================
    // Step 1: Add PUID configuration (newly added section)
    // ============================================================================
    
    // PUID configuration
    puid_wp_ = SSBConfReader->GetText("PUIDWorkingPoint");
    PUIDSFSys = SSBConfReader->GetText("PUIDSFSys"); 
    apply_puid_ = SSBConfReader->GetBool("ApplyPUID");
    puid_pt_threshold_ = SSBConfReader->GetNumber("PUIDPtThreshold");

    // Validate PUID settings
    if (puid_wp_ != "L" && puid_wp_ != "M" && puid_wp_ != "T") {
        std::cerr << "[WARNING] Invalid PUID working point: " << puid_wp_ 
                  << ". Using default 'L'." << std::endl;
        puid_wp_ = "L";
    }
    
    if (puid_pt_threshold_ <= 0 || puid_pt_threshold_ > 200) {
        std::cerr << "[WARNING] Invalid PUID pT threshold: " << puid_pt_threshold_ 
                  << ". Using default 50.0 GeV." << std::endl;
        puid_pt_threshold_ = 50.0;
    }
    
    // Print PUID configuration
    std::cout << "PUID Configuration:" << std::endl;
    std::cout << "  Working Point: " << puid_wp_ << std::endl;
    std::cout << "  Systematic: " << PUIDSFSys << std::endl;
    std::cout << "  Apply PUID: " << (apply_puid_ ? "true" : "false") << std::endl;
    std::cout << "  pT Threshold: " << puid_pt_threshold_ << " GeV" << std::endl;


    ///
    //std::cout << "triggerList : " << triggerList.size()<< std::endl;
    //SetObjectVariable();
}

void Analysis::SetObjectVariable() {
    //std::cout << "!!! SetObjectVariable start!!!" << std::endl;
    // Leptons //
    // muon //
    muons_pt  = floatVectors["Muon_pt"].get(); 
    muons_eta = floatVectors["Muon_eta"].get(); 
    muons_phi = floatVectors["Muon_phi"].get();
    muons_M   = floatVectors["Muon_mass"].get();
    
    muons_Id  = boolVectors["Muon_looseId"].get();
    muons_iso = floatVectors["Muon_pfRelIso03_all"].get();

    /// Muon ID
    if      (TString(MuonId).Contains( "Loose"  ) )  { muons_Id = boolVectors["Muon_looseId"].get(); }
    else if (TString(MuonId).Contains( "Medium"  ) ) { muons_Id = boolVectors["Muon_mediumId"].get(); }
    else if (TString(MuonId).Contains( "Tight"  ) )  { muons_Id = boolVectors["Muon_tightId"].get(); }
    else { std::cout << "Muon ID Error" << std::endl; }

    if (TString(MuonIsoType).Contains("PFIsodbeta03")) {
        if (floatVectors["Muon_pfRelIso03_all"] == nullptr) {
            std::cerr << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
            return;
        }
        muons_iso = floatVectors["Muon_pfRelIso03_all"].get();
    
    } else if (TString(MuonIsoType).Contains("PFIsodbeta04")) {
        if (floatVectors["Muon_pfRelIso04_all"] == nullptr) {
            std::cerr << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
            return;
        }
        muons_iso = floatVectors["Muon_pfRelIso04_all"].get();
    } else {
            std::cerr << "Muon Iso type Error" << std::endl;
            return;
    }

    // electron //
    elecs_pt  = floatVectors["Electron_pt"].get();
    elecs_eta = floatVectors["Electron_eta"].get(); 
    elecs_phi = floatVectors["Electron_phi"].get();
    elecs_M   = floatVectors["Electron_mass"].get();

    //elecs_iso = 
    if (TString(ElecIsoType).Contains("PFIsoRho03")) {
        elecs_iso = floatVectors["Electron_pfRelIso03_all"].get();
    } else if (TString(ElecIsoType).Contains("PFIsoRho04")) {
        elecs_iso = floatVectors["Electron_pfRelIso03_all"].get();
        std::cerr << "No PFIsoRho04 in NanoAOD..." << std::endl;
    } else {
        std::cerr << "Electron Iso type Error" << std::endl;
    }

    /// Electron iso type

    /// Electron ID
    elecs_scbId = intVectors["Electron_cutBased"].get();
    if (TString(ElecId).Contains("SCBLoose")) {
        //elecIdVariant = intVectors["Electron_cutBased"].get();
        eleid_scbcut = 2;
    } else if (TString(ElecId).Contains("SCBMedium")) {
        eleid_scbcut = 3;
    } else if (TString(ElecId).Contains("SCBTight")) {
        eleid_scbcut = 4;
    } else if (TString(ElecId).Contains("SCBVeto")) {
        eleid_scbcut = 1;
    } else if (TString(ElecId).Contains("MVALoose")) {
        elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        eleid_scbcut = 2;
    } else if (TString(ElecId).Contains("MVAMedium")) {
        elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
        eleid_scbcut = 3;
    } else if (TString(ElecId).Contains("MVATight")) {
        elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
        eleid_scbcut = 4;
    } else if (TString(ElecId).Contains("MVAVeto")) {
        elecs_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        eleid_scbcut = 1;
    } else {
        std::cerr << "Electron ID Error" << std::endl;
    }


    //////////////////////////////////////////////////////////////////////
    /// Set leptons for veto lepton (jet cleaning & third lepton veto) ///
    //////////////////////////////////////////////////////////////////////
    ////////////////////////
    /// Muon information ///
    ////////////////////////
    if (TString(veto_muoniso_type).Contains("PFIsodbeta03")) {
        if (floatVectors["Muon_pfRelIso03_all"] == nullptr) {
            std::cerr << "Error: Muon_pfRelIso03_all branch not initialized!" << std::endl;
            return;
        }
        muonsveto_iso = floatVectors["Muon_pfRelIso03_all"].get();
    } else if (TString(veto_muoniso_type).Contains("PFIsodbeta04")) {
        if (floatVectors["Muon_pfRelIso04_all"] == nullptr) {
            std::cerr << "Error: Muon_pfRelIso04_all branch not initialized!" << std::endl;
            return;
        }
        muonsveto_iso = floatVectors["Muon_pfRelIso04_all"].get();
    } else {
        std::cerr << "Muon Iso type Error" << std::endl;
        return;
    }

    
    if (muonsveto_iso != nullptr) {
        //std::cout << "muonsveto_iso size: " << muonsveto_iso->GetSize() << std::endl;
    } else {
        std::cerr << "Error: muonsveto_iso is null after assignment!" << std::endl;
    }


    /// Muon ID
    if      (TString(veto_muonid).Contains( "Loose"  ) )  { muonsveto_Id = boolVectors["Muon_looseId"].get(); }
    else if (TString(veto_muonid).Contains( "Medium"  ) ) { muonsveto_Id = boolVectors["Muon_mediumId"].get();}
    else if (TString(veto_muonid).Contains( "Tight"  ) )  { muonsveto_Id = boolVectors["Muon_tightId"].get(); }
    else { std::cout << "Muon ID Error" << std::endl; }
    
    if (muonsveto_Id == nullptr) {
        std::cerr << "Error: muonsveto_Id is null after assignment!" << std::endl;
    }


    ////////////////////////////
    /// Electron information ///
    ////////////////////////////
    elecsveto_scbId = intVectors["Electron_cutBased"].get();
    if (!intVectors["Electron_cutBased"].get()) {std::cerr <<"Error:  Electron_cutBased !! is nullptr! ub Electron information" <<std::endl;}
    if (!elecsveto_scbId) std::cerr << "Error: elecsveto_scbId is nullptr! ub Electron information " << std::endl;

    if (TString(veto_elecid).Contains("SCBLoose")) {
        //elecIdVariant = intVectors["Electron_cutBased"].get();
        elevetoid_scbcut = 2;
    } else if (TString(veto_elecid).Contains("SCBMedium")) {
        elevetoid_scbcut = 3;
    } else if (TString(veto_elecid).Contains("SCVTight")) {
        elevetoid_scbcut = 4;
    } else if (TString(veto_elecid).Contains("SCBVeto")) {
        elevetoid_scbcut = 1;
    } else if (TString(veto_elecid).Contains("MVALoose")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        elevetoid_scbcut = 2;
    } else if (TString(veto_elecid).Contains("MVAMedium")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP90"].get();
        elevetoid_scbcut = 3;
    } else if (TString(veto_elecid).Contains("MVATight")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WP80"].get();
        elevetoid_scbcut = 4;
    } else if (TString(veto_elecid).Contains("MVAVeto")) {
        elecsveto_mvaId = boolVectors["Electron_mvaFall17V2Iso_WPL"].get();
        elevetoid_scbcut = 1;
    } else {
        std::cerr << "Electron ID in veto selection Error " << veto_elecid << std::endl;
    }

    if (TString(veto_eleciso_type).Contains("PFIsoRho03")) {
        if (floatVectors.find("Electron_pfRelIso03_all") != floatVectors.end()) {
            auto* ptr = floatVectors["Electron_pfRelIso03_all"].get();
            if (ptr) {
                elecsveto_iso = ptr;
                //std::cout << "in sele.... elecsveto_iso .. " << std::endl;
            } else {
                std::cerr << "Error: 'Electron_pfRelIso03_all' is a null unique_ptr." << std::endl;
            }
        } else {
            std::cerr << "Error: 'Electron_pfRelIso03_all' not found in floatVectors." << std::endl;
        }
    } else if (TString(veto_eleciso_type).Contains("PFIsoRho04")) {
        std::cerr << "No PFIsoRho04 in NanoAOD..." << std::endl;
    } else {
        std::cerr << "Electron Iso type Error" << std::endl;
    }

    if (elecsveto_iso == nullptr) {
        std::cerr << "Error: elecsveto_iso is null after assignment!" << std::endl;
    }


    //////////////////////
    /// Set Jet object ///
    //////////////////////
     
    jets_pt  = floatVectors["Jet_pt"].get(); 
    jets_eta = floatVectors["Jet_eta"].get(); 
    jets_phi = floatVectors["Jet_phi"].get();
    jets_M   = floatVectors["Jet_mass"].get();
    jets_Id = intVectors["Jet_jetId"].get();
    jets_puId = intVectors["Jet_puId"].get();// for Run 2 50 GeV Jets have to pass PUID //

    if (!isData){ 
        gen_jets_pt  = floatVectors["Gen_Jet_pt"].get();
        gen_jets_eta = floatVectors["Gen_Jet_eta"].get();
        gen_jets_phi = floatVectors["Gen_Jet_phi"].get();
        gen_jets_M   = floatVectors["Gen_Jet_mass"].get();
    }

    if (RunPeriod.Contains("2016")) {
        if      (JetId == "PFLoose") { jet_id = 1; } // Loose ID
        else if (JetId == "PFTight") { jet_id = 3; } // Loose + Tight ID
        else if (JetId == "PFLooseLepVeto") { jet_id = 5; } // Loose + Tight + TightLeptonVeto
        else if (JetId == "PFTightLepVeto") { jet_id = 7; } // Loose + Tight + TightLeptonVeto
        else { std::cout << "Jet condition error for 2016!" << std::endl; }
    } 
    else if (RunPeriod.Contains("2017") || RunPeriod.Contains("2018")) {
        if      (JetId == "PFTight") { jet_id = 2; } // Tight ID
        else if (JetId == "PFTightLepVeto") { jet_id = 6; } // Tight + TightLeptonVeto
        else { std::cout << "Jet condition error for 2017/2018!" << std::endl; }
    } 
    else if (RunPeriod.Contains("2022") || RunPeriod.Contains("2023")) {
        if      (JetId == "PFTight") { jet_id = 2; } // Tight ID
        else if (JetId == "PFTightLepVeto") { jet_id = 6; } // Tight + TightLeptonVeto
        else { std::cout << "Jet condition error for 2022/2023!" << std::endl; }
    } 
    else {
        std::cout << "RunPeriod not recognized!" << std::endl;
    }


    // Set b-tag discriminator and threshold based on algorithm, WP and RunPeriod
    
    // First set the appropriate b-tag discriminator variable
    if (TString(JetbTag).Contains("deepCSV")) {
        jets_btag = floatVectors["Jet_btagDeepB"].get();
    }
    else if (TString(JetbTag).Contains("deepJet")) {
        jets_btag = floatVectors["Jet_btagDeepFlavB"].get();
    }
    else if (TString(JetbTag).Contains("pfCSVV2")) {
        jets_btag = floatVectors["Jet_btagCSVV2"].get();
    }
    else if (TString(JetbTag).Contains("CSV") && !TString(JetbTag).Contains("deep")) {
        jets_btag = floatVectors["Jet_btagCSV"].get(); // Modify with appropriate variable name if needed
    }
    else if (TString(JetbTag).Contains("CISV")) {
        jets_btag = floatVectors["Jet_btagCISV"].get(); // Modify with appropriate variable name if needed
    }
    else {
        std::cout << "Error: Unknown b-tagging algorithm in " << JetbTag << std::endl;
        jets_btag = nullptr;
        bdisccut = -1.0;
        return;
    }
   

    // Parse B-tagging algorithm from JetbTag
    btag_algo_ = "DeepCSV";  // Default
    if (TString(JetbTag).Contains("deepCSV")) {
        btag_algo_ = "DeepCSV";
    } else if (TString(JetbTag).Contains("deepJet")) {
        btag_algo_ = "DeepJet";
    } else if (TString(JetbTag).Contains("pfCSVV2")) {
        btag_algo_ = "CSVv2";
    } else {
        std::cerr << "Unknown b-tagging algorithm in JetbTag: " << JetbTag << std::endl;
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
        std::cerr << "Unknown working point in JetbTag: " << JetbTag << std::endl;
    }
    
    // Debug output (optional)
    //std::cout << "Parsed B-tag config: " << btag_algo_ 
    //          << " " << btag_wp_ << std::endl;
 
    // Set appropriate bdisccut value based on RunPeriod, algorithm, and WP
    if (TString(RunPeriod).Contains("2016PreVFP")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.2027;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.6001;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.8819;
            else std::cout << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0508;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.2598;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.6502;
            else std::cout << "Unknown deepJet working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("pfCSVV2")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.5426;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.8484;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.9535;
            else std::cout << "Unknown pfCSVV2 working point!" << std::endl;
        }
    }
    else if (TString(RunPeriod).Contains("2016PostVFP")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.1918;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.5847;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.8767;
            else std::cout << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0480;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.2489;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.6377;
            else std::cout << "Unknown deepJet working point!" << std::endl;
        }
    }
    else if (TString(RunPeriod).Contains("2017")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.1355;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.4506;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7738;
            else std::cout << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0532;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.3040;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7476;
            else std::cout << "Unknown deepJet working point!" << std::endl;
        }
    }
    else if (TString(RunPeriod).Contains("2018")) {
        if (TString(JetbTag).Contains("deepCSV")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.1208;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.4168;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7665;
            else std::cout << "Unknown deepCSV working point!" << std::endl;
        }
        else if (TString(JetbTag).Contains("deepJet")) {
            if (TString(JetbTag).Contains("L")) bdisccut = 0.0490;
            else if (TString(JetbTag).Contains("M")) bdisccut = 0.2783;
            else if (TString(JetbTag).Contains("T")) bdisccut = 0.7100;
            else std::cout << "Unknown deepJet working point!" << std::endl;
        }
    }
    else {
        std::cout << "Error: Unsupported run period: " << RunPeriod << std::endl;
        jets_btag = nullptr;
        bdisccut = -1.0;
    }
    
    /*std::cout << "Setup B-tagging: " << JetbTag << " for period " << RunPeriod 
              << " with cut value " << bdisccut << std::endl;*/
    /// MET ///
    met_pt  = floatSingles["PuppiMET_pt"].get(); 
    met_phi  = floatSingles["PuppiMET_phi"].get();
    object_variables_set_ = true; 
//    std::cout << " met_pt : " << met_pt << std::endl; 
//    std::cout << " met_phi : " << met_phi << std::endl; 

    //std::cout << "End of SetObjectVariable !" << std::endl;
}



// Event loop function
void Analysis::Loop() {
    if (!chain) {
        std::cerr << "Error: TChain is null!" << std::endl;
        return;
    }

    Long64_t nEntries = chain->GetEntries();
    fReader.Restart(); // Reset the reader to the beginning

    MCSF();

    for (Long64_t ientry = 0; ientry < ((NumEvt == -1 || NumEvt > nEntries) ? nEntries : NumEvt); ++ientry) {


        if (!fReader.Next()) { // Move to the next entry
            std::cerr << "Error: Failed to read entry " << ientry << std::endl;
            break;
        }
        evt_weight_ = 1.;
        MCSFApply();

        GenWeightApply();
        PUWeightApply();
        L1PreFireApply();

        SetObjectVariable(); //
        LeptonSelector(); 

        LeptonOrder();
        JetSelector();
        PUIDSFApply();  // Apply PUID event weight using collected information
 
        JetOrder(); 
        bJetSelector(); 
        //METDefiner();

        //if (i > 100) break; //%lld supports Long64_t
        if (ientry % 10000 == 0) {
            printf("Event %lld\n", ientry); //%lld supports Long64_t
        }


        // Noise (MET) Filter //
        if ( METFilterAPP() == false ) {continue;}
        // Trigger Requirement //
        if ( Trigger() == false ) {continue;}
  
        // Good Primary Vertex Selection //
        if (intSingles["PV_npvsGood"] && **intSingles["PV_npvsGood"] < 1) {
            continue;
        }

        //std::cout << "**intSingles[PV_npvsGood] : " << **intSingles["PV_npvsGood"] <<std::endl;
        //std::cout << "evt_weight_ : " << evt_weight_ << std::endl;
        FillHisto( h_Num_PV[0]    , **intSingles["PV_npvsGood"] , evt_weight_ );
        if ( NumIsoLeptons(2) == false ) {continue;}

        if (ThirdLeptonVeto() == false ) {continue;}
        //std::cout << "after ThirdLeptons : " << std::endl;
        if (LeptonsPtAddtional() == false ) {continue;}
        if (DiLeptonMassCut() == false) {continue;}

        LeptonSFApply();
        TriggerSFApply();

        /// Step 1 ///
        //std::cout << "? evet weight " << evt_weight_ << std::endl;
        num_pv =  **intSingles["PV_npvsGood"]; 
        FillHisto( h_DiLepMass[1], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Num_PV[1],     num_pv, evt_weight_ );
        FillHisto( h_Lep1pt[1] ,    (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[1],    (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[1],    (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[1] ,    (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[1],    (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[1],    (Lep2).Phi() , evt_weight_ );
        FillHisto( h_METpt[1]   ,   Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[1]  ,   Met.Phi()  , evt_weight_ );
        FillHisto( h_Num_Jets[1]  , v_jet_idx.size(), evt_weight_ );
        //FillHisto( h_Num_bJets[1], nbtagged, evt_weight_ );

        if (ZVetoCut() == false) {continue;}

        FillHisto( h_DiLepMass[2], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Num_PV[2],     num_pv, evt_weight_ );
        FillHisto( h_Lep1pt[2] ,    (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[2],    (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[2],    (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[2] ,    (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[2],    (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[2],    (Lep2).Phi() , evt_weight_ );
        FillHisto( h_METpt[2]   ,   Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[2]  ,   Met.Phi()  , evt_weight_ );
        FillHisto( h_Num_Jets[2]  , v_jet_idx.size(), evt_weight_ );      

        if (NumJetCut(v_jet_idx) == false) {continue;}

        FillHisto( h_DiLepMass[3], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Num_PV[3],     num_pv, evt_weight_ );

        FillHisto( h_Lep1pt[3] ,    (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[3],    (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[3],    (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[3] ,    (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[3],    (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[3],    (Lep2).Phi() , evt_weight_ );

        //if ((Jet2).Pt() < 30.) printf("jet2pt %lf eta %lf \n",Jet2.Pt(), Jet2.Eta());//std::cout << "Wrong! " << Jet2.Pt() << std::endl;
        FillHisto( h_Jet1pt[3] ,    (Jet1).Pt()  , evt_weight_ );
        FillHisto( h_Jet1eta[3],    (Jet1).Eta() , evt_weight_ );
        FillHisto( h_Jet1phi[3],    (Jet1).Phi() , evt_weight_ );
        FillHisto( h_Jet2pt[3] ,    (Jet2).Pt()  , evt_weight_ );
        FillHisto( h_Jet2eta[3],    (Jet2).Eta() , evt_weight_ );
        FillHisto( h_Jet2phi[3],    (Jet2).Phi() , evt_weight_ );
        FillHisto( h_Num_Jets[3]  , v_jet_idx.size(), evt_weight_ );      

        FillHisto( h_METpt[3]   ,   Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[3]  ,   Met.Phi()  , evt_weight_ );

        if (METCut(Met) == false) {continue;}

        FillHisto( h_DiLepMass[4], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Num_PV[4],     num_pv, evt_weight_ );
        FillHisto( h_Lep1pt[4] ,    (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[4],    (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[4],    (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[4] ,    (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[4],    (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[4],    (Lep2).Phi() , evt_weight_ );

        FillHisto( h_Jet1pt[4] ,    (Jet1).Pt()  , evt_weight_ );
        FillHisto( h_Jet1eta[4],    (Jet1).Eta() , evt_weight_ );
        FillHisto( h_Jet1phi[4],    (Jet1).Phi() , evt_weight_ );
        FillHisto( h_Jet2pt[4] ,    (Jet2).Pt()  , evt_weight_ );
        FillHisto( h_Jet2eta[4],    (Jet2).Eta() , evt_weight_ );
        FillHisto( h_Jet2phi[4],    (Jet2).Phi() , evt_weight_ );
        FillHisto( h_Num_Jets[4]  , v_jet_idx.size(), evt_weight_ );      

        FillHisto( h_METpt[4]   ,   Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[4]  ,   Met.Phi()  , evt_weight_ );

        BTaggingSFApply();

        if (NumbJetCut(v_bjet_idx) == false) {continue;}

        FillHisto( h_DiLepMass[5], ( (Lep1)+(Lep2) ).M(), evt_weight_ );
        FillHisto( h_Num_PV[5],     num_pv, evt_weight_ );
        FillHisto( h_Lep1pt[5] ,    (Lep1).Pt()  , evt_weight_ );
        FillHisto( h_Lep1eta[5],    (Lep1).Eta() , evt_weight_ );
        FillHisto( h_Lep1phi[5],    (Lep1).Phi() , evt_weight_ );
        FillHisto( h_Lep2pt[5] ,    (Lep2).Pt()  , evt_weight_ );
        FillHisto( h_Lep2eta[5],    (Lep2).Eta() , evt_weight_ );
        FillHisto( h_Lep2phi[5],    (Lep2).Phi() , evt_weight_ );

        FillHisto( h_Jet1pt[5] ,    (Jet1).Pt()  , evt_weight_ );
        FillHisto( h_Jet1eta[5],    (Jet1).Eta() , evt_weight_ );
        FillHisto( h_Jet1phi[5],    (Jet1).Phi() , evt_weight_ );
        FillHisto( h_Jet2pt[5] ,    (Jet2).Pt()  , evt_weight_ );
        FillHisto( h_Jet2eta[5],    (Jet2).Eta() , evt_weight_ );
        FillHisto( h_Jet2phi[5],    (Jet2).Phi() , evt_weight_ );
        FillHisto( h_Num_Jets[5]  , v_jet_idx.size(), evt_weight_ );      

        FillHisto( h_METpt[5]   ,   Met.Pt()  , evt_weight_ );
        FillHisto( h_METphi[5]  ,   Met.Phi()  , evt_weight_ );
        SetUpKINObs();
        if (isKinSol)
        {
            FillHisto( h_Lep1pt[8] , Lep1.Pt() , evt_weight_ );
            FillHisto( h_Lep2pt[8] , Lep2.Pt() , evt_weight_ );
            FillHisto( h_Lep1eta[8], Lep1.Eta(), evt_weight_ );
            FillHisto( h_Lep2eta[8], Lep2.Eta(), evt_weight_ );
            FillHisto( h_Lep1phi[8], Lep1.Phi(), evt_weight_ );
            FillHisto( h_Lep2phi[8], Lep2.Phi(), evt_weight_ );

            FillHisto( h_Jet1pt[8] , Jet1.Pt() , evt_weight_ );
            FillHisto( h_Jet2pt[8] , Jet2.Pt() , evt_weight_ );
            FillHisto( h_Jet1eta[8], Jet1.Eta(), evt_weight_ );
            FillHisto( h_Jet2eta[8], Jet2.Eta(), evt_weight_ );
            FillHisto( h_Jet1phi[8], Jet1.Phi(), evt_weight_ );
            FillHisto( h_Jet2phi[8], Jet2.Phi(), evt_weight_ );
            FillHisto( h_METpt[8]  , Met.Pt()  , evt_weight_ );
            FillHisto( h_METphi[8] , Met.Phi() , evt_weight_ );
            //FillHisto( h_HT[8]     , AllJetpt   , evt_weight_);
            
            FillHisto( h_DiLepMass[8], ( Lep1+Lep2 ).M(), evt_weight_ );
            
            FillHisto( h_Num_PV[8]   , num_pv          ,  evt_weight_ );
            FillHisto( h_Num_Jets[8] , v_jet_idx.size(),  evt_weight_ );
            FillHisto( h_Num_bJets[8], v_bjet_idx.size(), evt_weight_ );
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

            for (int i = 0; i < v_recocp_O.size(); ++ i)
            {
               FillHisto( h_Reco_CPO_[i], v_recocp_O[i] , evt_weight_ );
               FillHisto( h_Reco_CPO_ReRange_[i], v_recocp_O[i] , evt_weight_ );
            } 

                  
 
        }


    }// end of event iteration //
    std::cout << "End Loop !!" << std::endl;
}


template <typename T>
std::unique_ptr<TTreeReaderValue<T>> Analysis::DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src) {
    if (src) {
        // If source pointer is valid, make a deep copy
        return std::make_unique<TTreeReaderValue<T>>(*src);
    } else {
	std::cerr << "no source !!! "  << std::endl;
        // If source is nullptr, return nullptr
        return nullptr;
    }
}



std::string Analysis::removeSubstring(std::string &str, const std::string &keyword) {
    size_t pos = str.find(keyword);  // find specific keyword ! 
    if (pos != std::string::npos) {
        str.erase(pos);  // 
    }
    std::cout << "str : " << str << std::endl;
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
            std::cerr << "Error: Trigger " << trgName << " not found in triggerList." << std::endl;
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
            trigpass = pass_double;
            return trigpass;
         } else {
            std::cout << "[Trigger] Check FileName_ for dimuon in 2018" << std::endl;
            return false;
         }
      }

      else if (TString(Decaymode).Contains("muel")) {
         if (TString(FileName_).Contains("MuonEG")) {
            bool pass_double = SelTrigger(DLtrigName);
            trigpass = pass_double;
            return trigpass;
         } else if (TString(FileName_).Contains("EGamma")) {
            bool pass_single = SelTrigger(SLtrigName);
            trigpass = pass_single;
            return trigpass;
         } else {
            std::cout << "[Trigger] Check FileName_ for muel in 2018" << std::endl;
            return false;
         }
      }

      else if (TString(Decaymode).Contains("dielec")) {
         if (TString(FileName_).Contains("EGamma")) {
            bool pass_single = SelTrigger(SLtrigName);
            bool pass_double = SelTrigger(DLtrigName);
            trigpass = (pass_single && !pass_double) || (!pass_single && pass_double);
            return trigpass;
         } else {
            std::cout << "[Trigger] Check FileName_ for dielec in 2018" << std::endl;
            return false;
         }
      }

      else {
         std::cout << "[Trigger] Check Decaymode for 2018" << std::endl;
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
         std::cout << "[Trigger] Check FileName_ for 2016/2017" << std::endl;
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
   
   std::cout << "Original input: " << inname << ", processed name: " << inputName << std::endl;
   
   return inputName; 
}

void Analysis::MCSF()
{
    if (FileName_.Contains("Data")||FileName_.Contains("Single")||FileName_.Contains("EG")){ mc_sf_ = 1.; return; }
    /// Open Xsec Tables ///
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
    //cout << "xsec_filePath : " << xsec_filePath << std::endl;
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
       //cout << "Load Xsection Table!" << std::endl;
       while (fscanf(xsecs_, "%s %d %d %d %d %lf %lf\n", sampleName, &totalevt_, &positive_, &negative_, &posi_nega_, &xsec_, &br_ ) != EOF)
       {
          std::cout 
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
    else {std::cout << "No xsec_filePath !!!" << xsec_filePath << std::endl;return;} 
    std::cout << "Lumi : " << Lumi << std::endl;
    double lumi = Lumi/1000000;
    auto it = m_sam_xsec.find(FileName_.Data());
    if (it !=  m_sam_xsec.end()){
        std::cout << "SK Key " << FileName_.Data() << " found in the std::map."<< std::endl;
        mc_sf_ = (m_sam_xsec[FileName_.Data()]*m_sam_br[FileName_.Data()]*lumi)/m_sam_posi_nega[FileName_.Data()];
        std::cout << "mc_sf_ " << mc_sf_ << std::endl;
    }
    else {
        mc_sf_ =1.;
        std::cout << "Key " << FileName_.Data() << " not found in the std::map. mc sf is 1" << mc_sf_ << std::endl;
    }
    return;
}


// Apply MC SF To Event //
void Analysis::MCSFApply()
{
    evt_weight_beforemcsf_ =1; // Initailize evt_weight_beforemcsf_ //
    evt_weight_beforemcsf_ = evt_weight_; // keep event weight //
     
    if ( !TString(FileName_).Contains( "Data") ){ evt_weight_ = evt_weight_*mc_sf_; } // apply MC scale factor // 
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
        l1prefire_ = **floatSingles["L1PreFiringWeight_Nom"];
    }
    else if (strstr(sys, "up")) {
        l1prefire_ = **floatSingles["L1PreFiringWeight_Up"];
    }
    else if (strstr(sys, "down")) {
        l1prefire_ = **floatSingles["L1PreFiringWeight_Dn"];
    }
    else if (!strstr(sys, "none")) {
        // Only print error if not "none"
        std::cout << "L1Prefiring sys Error ... Default is Weight_L1Prefiring ... : " << L1PreFireSys << std::endl;
        l1prefire_ = **floatSingles["L1PreFiringWeight_Nom"];
    }
    //std::cout << "l1prefire_ : " << l1prefire_ << std::endl;
    evt_weight_ *= l1prefire_;
}

// Safely create a TLorentzVector
TLorentzVector Analysis::createLorentzVector(float pt, float eta, float phi, float mass) {
    //std::cout << "createLorentzVector " << std::endl;
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

    // Ensure E^2 >= P^2 (physical validity check)
    /*if (lv.E() < lv.P()) {
        std::cout << "pt : " << pt << " eta: " << eta << " phi :" << phi << " mass " << mass << std::endl;
        throw std::runtime_error("Invalid LorentzVector: E < P");
    }*/

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
                     (*intVectors["Muon_charge"])[v_lepton_idx[0]] != (*intVectors["Muon_charge"])[i]) {
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
                !elecSCBId(elecs_scbId->At(i), eleid_scbcut) ||
                (fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                !(elecCharge((*intVectors["Electron_tightCharge"])[i]))) {
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
                     (*intVectors["Electron_charge"])[v_electron_idx[0]] != (*intVectors["Electron_charge"])[i]) {
                // Add second selected electron with opposite charge
                v_lepton_idx.push_back(i);
                v_electron_idx.push_back(i);
                // Store its TLorentzVector in electrons vector
                elecs.push_back(pre_elecs.at(i));
            }
        }
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
                    !elecSCBId(elecs_scbId->At(i), eleid_scbcut) ||
                   (fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                   fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                  !(elecCharge((*intVectors["Electron_tightCharge"])[i])) ||
                  !(*boolVectors["Electron_convVeto"])[i]) {
                    continue;
                }
        
                // Only select electrons with charge opposite to the first selected muon
                if ((*intVectors["Muon_charge"])[v_muon_idx[0]] != (*intVectors["Electron_charge"])[i]) {
                    // Add selected electron
                    v_electron_idx.push_back(i);
                    // Store its TLorentzVector in electrons vector
                    elecs.push_back(pre_elecs.at(i));
                }
            }
        }
    }
    else {
        std::cerr << "Lepton Selection error" << std::endl;
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
    // Helper function for assigning leptons and debugging
    auto assignLeptons = [&](const std::vector<TLorentzVector>& leptons,
                             const std::string& chargeKey,
                             const std::vector<int>& indices,
                             int idx1, int idx2) {
        if (indices.size() <= std::max(idx1, idx2)) {
            std::cerr << "Error: indices size (" << indices.size() 
                      << ") is smaller than required index." << std::endl;
            return;
        }

        Lep1 = leptons.at(indices[idx1]);
        Lep2 = leptons.at(indices[idx2]);

        // Set Lep & AnLep //
        if ((*intVectors[chargeKey])[indices[idx1]] < 0) {
            Lep = leptons.at(indices[idx1]);
            AnLep = leptons.at(indices[idx2]);
        } else {
            Lep = leptons.at(indices[idx2]);
            AnLep = leptons.at(indices[idx1]);
        }

    }; // end of assignLeptons //

    //std::cout << "sk1 " << std::endl;

    // Handle dimuon decay mode
    if (TString(Decaymode).Contains("dimuon")) {
    //std::cout << "sk2 " << std::endl;

        if (v_muon_idx.size() > 1) {
    //std::cout << "sk3 " << std::endl;
            assignLeptons(pre_muons, "Muon_charge", v_muon_idx, 0, 1);
        } /*else {
            std::cerr << "Lepton TLorentzVector Error: v_muon_idx is empty or too small for Decaymode = dimuon" 
                      << " (size: " << v_muon_idx.size() << ")" << std::endl;
            return;
        }*/
    //std::cout << "sk4 " << std::endl;
    }
    else if (TString(Decaymode).Contains("dielec")) {
    // Handle dielectron decay mode
        if (v_electron_idx.size() > 1) {
            assignLeptons(pre_elecs, "Electron_charge", v_electron_idx, 0, 1);
        } else {
            std::cerr << "Lepton TLorentzVector Error: v_electron_idx is empty or too small for Decaymode = dielec" 
                      << " (size: " << v_electron_idx.size() << ")" << std::endl;
            return;
        }
    }
    else if (TString(Decaymode).Contains("muel")) {
    // Handle muon-electron decay mode
        if (v_muon_idx.size() > 0 && v_electron_idx.size() > 0) {
            if (pre_muons.at(v_muon_idx.at(0)).Pt() > pre_elecs.at(v_electron_idx.at(0)).Pt()) {
                Lep1 = pre_muons.at(v_muon_idx.at(0));
                Lep2 = pre_elecs.at(v_electron_idx.at(0));
            } else {
                Lep1 = pre_elecs.at(v_electron_idx.at(0));
                Lep2 = pre_muons.at(v_muon_idx.at(0));
            }

            if ((*intVectors["Muon_charge"])[v_muon_idx[0]] < 0) {
                Lep = pre_muons.at(v_muon_idx.at(0));
                AnLep = pre_elecs.at(v_electron_idx.at(0));
            } else {
                Lep = pre_elecs.at(v_electron_idx.at(0));
                AnLep = pre_muons.at(v_muon_idx.at(0));
            }
        } else {
         //   std::cerr << "Lepton TLorentzVector Error: v_muon_idx size = " 
         //             << v_muon_idx.size() << ", v_electron_idx size = " 
         //             << v_electron_idx.size() << " for Decaymode = muel" << std::endl;
            return;
        }
    }
    // Handle invalid Decaymode
    else {
        std::cerr << "Lepton TLorentzVector Error: Decaymode = " << Decaymode << std::endl;
    }

}


void Analysis::MakeMuonCollection() {
    // make muon collection tlorentzvectors//
    //std::cout << "Muon Correction !!" << std::endl;
    pre_muons.clear();
    for (int imu = 0; imu < muons_pt->GetSize(); ++imu){
        //std::cout << "imu " << imu << " muons_pt->At(imu) : " << muons_pt->At(imu) << " muons_eta->At(imu) : " << muons_eta->At(imu) << " muons_phi->At(imu) : " << muons_phi->At(imu) << " muons_M->At(imu)  : "<< muons_M->At(imu)  << std::endl; 
        pre_muons.push_back(createLorentzVector(muons_pt->At(imu), muons_eta->At(imu), muons_phi->At(imu), muons_M->At(imu) )); 
    }
//    std::cout << "end of MakeMuonCollection !" << std::endl;
//    std::cout << "size of pre_muons.size : " << muons.size() << std::endl; 
    return;
}


void Analysis::MakeElecCollection() {
    // Electron collection logic
    pre_elecs.clear();
    for (int iele = 0; iele < elecs_pt->GetSize(); ++iele){
        pre_elecs.push_back(createLorentzVector(elecs_pt->At(iele), elecs_eta->At(iele), elecs_phi->At(iele), elecs_M->At(iele) )); 
    }
 
//    std::cout << "end of MakeElecCollection !" << std::endl;
//    std::cout << "size of pre_elecs.size : " << pre_elecs.size() << std::endl; 
    return;
}
/*
void Analysis::MakeJetCollection() {
    pre_jets.clear();

    if (jets_pt == nullptr || jets_eta == nullptr || jets_phi == nullptr || jets_M == nullptr) {
        std::cerr << "Error: Some jet branch pointers are null in MakeJetCollection()" << std::endl;
        return;
    }

    Int_t nJets = jets_pt->GetSize();
    std::vector<TLorentzVector> rawJets;
    std::vector<float> rawFactors;
    std::vector<float> jetAreas;
    std::vector<TLorentzVector> genJets;
    std::vector<int> genJetIndices;

    rawJets.reserve(nJets);
    rawFactors.reserve(nJets);
    jetAreas.reserve(nJets);
    genJets.reserve(nJets);
    genJetIndices.reserve(nJets);

    for (int ijet = 0; ijet < nJets; ++ijet) {
        try {
            TLorentzVector rawJet;
            rawJet.SetPtEtaPhiM(jets_pt->At(ijet), jets_eta->At(ijet), jets_phi->At(ijet), jets_M->At(ijet));
            rawJets.push_back(rawJet);

            float rawFactor = (floatVectors.count("Jet_rawFactor") && floatVectors["Jet_rawFactor"]->GetSize() > ijet)
                                ? floatVectors["Jet_rawFactor"]->At(ijet) : 0.0;
            rawFactors.push_back(rawFactor);

            float area = (floatVectors.count("Jet_area") && floatVectors["Jet_area"]->GetSize() > ijet)
                           ? floatVectors["Jet_area"]->At(ijet) : 0.5;
            jetAreas.push_back(area);

            int genIdx = (intVectors.count("Jet_genJetIdx") && intVectors["Jet_genJetIdx"]->GetSize() > ijet)
                           ? intVectors["Jet_genJetIdx"]->At(ijet) : -1;
            genJetIndices.push_back(genIdx);
        } catch (const std::exception& e) {
            std::cerr << "Error reading jet at index " << ijet << ": " << e.what() << std::endl;
        }
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

    double rho = (floatSingles.count("fixedGridRhoFastjetAll") > 0) ? **floatSingles["fixedGridRhoFastjetAll"] : 0.0;
    double raw_met_pt = (floatSingles.count("MET_pt") > 0)  ? **floatSingles["MET_pt"]  : 0.0;
    double raw_met_phi = (floatSingles.count("MET_phi") > 0) ? **floatSingles["MET_phi"] : 0.0;
    //bool isData = TString(FileName_).Contains("Data");

    JetCorrectionOutput corr_output = SSBCorr->ApplyJetCorrectionsWithMET(
        rawJets,
        rawFactors,
        jetAreas,
        rho,
        isData,
        true,
        true,
        raw_met_pt,
        raw_met_phi,
        genJets,
        genJetIndices
    );

    pre_jets = corr_output.corrected_jets;
    Met = corr_output.corrected_met;
}*/

void Analysis::MakeJetCollection() {
    pre_jets.clear();

    if (jets_pt == nullptr || jets_eta == nullptr || jets_phi == nullptr || jets_M == nullptr) {
        std::cerr << "Error: Some jet branch pointers are null in MakeJetCollection()" << std::endl;
        return;
    }

    Int_t nJets = jets_pt->GetSize();
    std::vector<TLorentzVector> rawJets;
    std::vector<float> rawFactors;
    std::vector<float> jetAreas;
    std::vector<TLorentzVector> genJets;
    std::vector<int> genJetIndices;

    rawJets.reserve(nJets);
    rawFactors.reserve(nJets);
    jetAreas.reserve(nJets);
    genJets.reserve(nJets);
    genJetIndices.reserve(nJets);

    // ============================================================================
    // SOLUTION 1: Ensure all jets are included, even with invalid data
    // ============================================================================
    for (int ijet = 0; ijet < nJets; ++ijet) {
        try {
            TLorentzVector rawJet;
            rawJet.SetPtEtaPhiM(jets_pt->At(ijet), jets_eta->At(ijet), jets_phi->At(ijet), jets_M->At(ijet));
            rawJets.push_back(rawJet);

            float rawFactor = (floatVectors.count("Jet_rawFactor") && floatVectors["Jet_rawFactor"]->GetSize() > ijet)
                                ? floatVectors["Jet_rawFactor"]->At(ijet) : 0.0;
            rawFactors.push_back(rawFactor);

            float area = (floatVectors.count("Jet_area") && floatVectors["Jet_area"]->GetSize() > ijet)
                           ? floatVectors["Jet_area"]->At(ijet) : 0.5;
            jetAreas.push_back(area);

            int genIdx = (intVectors.count("Jet_genJetIdx") && intVectors["Jet_genJetIdx"]->GetSize() > ijet)
                           ? intVectors["Jet_genJetIdx"]->At(ijet) : -1;
            genJetIndices.push_back(genIdx);
            
        } catch (const std::exception& e) {
            std::cerr << "Error reading jet at index " << ijet << ": " << e.what() << std::endl;
            std::cout << "erro ! "  << std::endl;
            std::cerr << "Error reading jet at index " << ijet << ": " << e.what() << std::endl;
            std::cout << ">>> CREATING DUMMY JET AT INDEX " << ijet << " <<<" << std::endl;  // 추가
 
            // ============================================================================
            // FIX: Add invalid dummy jet with -999 values to maintain size consistency
            // ============================================================================
            TLorentzVector dummyJet;
            dummyJet.SetPtEtaPhiM(-999.0, -999.0, -999.0, -999.0);  // Clearly invalid jet
            rawJets.push_back(dummyJet);
            rawFactors.push_back(0.0);
            jetAreas.push_back(0.5);
            genJetIndices.push_back(-1);
        }
    }

    // Verify size consistency before proceeding
    if (rawJets.size() != static_cast<size_t>(nJets)) {
        std::cerr << "[ERROR] Size mismatch after jet reading: expected " << nJets 
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

    double rho = (floatSingles.count("fixedGridRhoFastjetAll") > 0) ? **floatSingles["fixedGridRhoFastjetAll"] : 0.0;
    double raw_met_pt = (floatSingles.count("PuppiMET_pt") > 0)  ? **floatSingles["PuppiMET_pt"]  : 0.0;
    double raw_met_phi = (floatSingles.count("PuppiMET_phi") > 0) ? **floatSingles["PuppiMET_phi"] : 0.0;

    JetCorrectionOutput corr_output = SSBCorr->ApplyJetCorrectionsWithMET(
        rawJets,
        rawFactors,
        jetAreas,
        rho,
        isData,
        true,
        true,
        raw_met_pt,
        raw_met_phi,
        genJets,
        genJetIndices
    );

    pre_jets = corr_output.corrected_jets;
    Met = corr_output.corrected_met;
    
    // ============================================================================
    // VERIFICATION: Check final size consistency
    // ============================================================================
    if (pre_jets.size() != static_cast<size_t>(nJets)) {
        std::cout << "[INFO] MakeJetCollection: Input jets=" << nJets 
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
        std::cerr << "Error: Unrecognized decay mode in NumIsoLeptons()" << std::endl;
    }
    //std::cout << "v_muon_idx.size() "<< v_muon_idx.size() << "numLeptons " << numLeptons << std::endl; 
    return numLeptons;
}

void Analysis::JetSelector() {
    // Pre-conditions check
    if (!object_variables_set_) {
        std::cerr << "ERROR: JetSelector() called before SetObjectVariable()!" << std::endl;
        throw std::runtime_error("SetObjectVariable() must be called before JetSelector()");
    }

    // Reset state
    jets_selected_ = false;
    
    // Check if necessary pointers are initialized
    if (jets_pt == nullptr || jets_eta == nullptr || jets_phi == nullptr || jets_M == nullptr) {
        std::cerr << "Error: Basic jet variables (pt, eta, phi, M) are not initialized. Make sure SetObjectVariable() was called." << std::endl;
        return;
    }

    MakeJetCollection();
    
    // ============================================================================
    // Step 4: Collect PUID candidate information (for weight calculation in PUIDSFApply)
    // ============================================================================
    CollectPUIDCandidates();
    
    v_jet_idx.clear();
    jets.clear();
    
    // Lambda functions for common jet selection checks
    auto passKinematicCuts = [this](float pt, float eta) -> bool {
        return pt > jet_pt && fabs(eta) < jet_eta;
    };

    Int_t nJets = jets_pt->GetSize();
    
    // ============================================================================
    // Step 4: CLEAN JET SELECTION LOOP - No PUID weight calculation
    // ============================================================================
    for (int i = 0; i < nJets; i++) {
        TLorentzVector jetVec = pre_jets[i];
        float jetPt = jetVec.Pt();
        float jetEta = jetVec.Eta();
        
        // Apply kinematic cuts first
        if (!passKinematicCuts(jetPt, jetEta)) {
            continue;
        }

        // Check jet ID
        if (jets_Id != nullptr && jets_Id->At(i) < jet_id) {
            continue;
        }

        // Jet cleaning
        if (!JetCleaning(&jetVec)) {
            continue;
        }

        // ============================================================================
        // Step 4: PUID selection ONLY (no weight calculation)
        // ============================================================================
        if (apply_puid_ && jetPt <= puid_pt_threshold_) {
            int puId = (jets_puId != nullptr) ? jets_puId->At(i) : 0;
            bool passPUID = PassPileupID(jetPt, puId, puid_wp_);
            
            if (!passPUID) {
                continue;  // Skip this jet if PUID fails
            }
        }

        // Add to selected jets
        v_jet_idx.push_back(i);
        jets.push_back(pre_jets[i]);
    }

    // ============================================================================
    // Step 4: Clean completion - No PUID weight application here
    // ============================================================================
    jets_selected_ = true;
    
    //std::cout << "JetSelector completed: " << v_jet_idx.size() << " jets selected" << std::endl;
    
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
    else if (op_ == "Up"){std::cout << "no op " << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFUp->at(idx_);
    }
    else if (op_ == "Down"){std::cout << "no op" << std::endl;
        //jerfrac_ = Jet_EnergyResolution_SFDown->at(idx_);
    }
    else{
        std::cout << "Check out your JERSmearing option !!" << "op_ : "<< op_<< std::endl;
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
    else { std::cerr << "?? something wrong " << std::endl; }                                                                                      
    return lepptadd;
}
void Analysis::Start()
{
   //fout = new TFile(Form("output/%s",outfile),"RECREATE");
   if (strcmp(outdir.c_str(), "None") != 0 ) {
      fout = new TFile(Form("gsidcap://cluster142.knu.ac.kr/%s/%s", outdir.c_str(), outfile.c_str()), "RECREATE");

   }
   else {
      //fout = new TFile(Form("output/%s",outfile),"RECREATE");
      fout = new TFile(Form("output/%s", outfile.c_str()), "RECREATE");
   }
   std::cout << "fout - getname : " << fout->GetName() << std::endl;
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
//   Met.SetPtEtaPhiM(,0,,0);//MET_phi MET_pt // 
    //Met.SetPtEtaPhiM(met_pt, 0, met_phi, 0); // MET_phi MET_pt
    //Met.SetPtEtaPhiM(*met_pt, 0, *met_phi, 0); // MET_phi MET_pt
    //Met.SetPtEtaPhiM(static_cast<double>(*met_pt), 0, static_cast<double>(*met_phi), 0); // MET_phi MET_pt
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
                /*std::cout << "In SelectVetoMuons : muonsveto_iso : " << muonsveto_iso->At(imu) 
                          << " muons_pt->At(imu) : " << muons_pt->At(imu) 
                          << " muons_eta->At(imu) : " << muons_eta->At(imu) 
                          << " muonsveto_Id->At(imu) : " << muonsveto_Id->At(imu) 
                          << std::endl; */
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
        std::cerr << "Error: Unrecognized decay mode in SelectVetoMuons!" << std::endl;
    }
    //std::cout << "IN SelectVetoMuons : muonsveto : " << muonsveto.size() << std::endl;
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
        //std::cout << "nel :" << nel << std::endl;
        for (int iel = 0; iel < nel; ++iel) {
            // Apply cuts
            if (!passKinematicCuts(elecs_pt->At(iel), elecs_eta->At(iel), elec_pt, elec_eta) ||
                !elecSCBId(elecsveto_scbId->At(iel), elevetoid_scbcut) ||
                (fabs((*floatVectors["Electron_deltaEtaSC"])[iel] + elecs_eta->At(iel)) > 1.4442 &&
                 fabs((*floatVectors["Electron_deltaEtaSC"])[iel] + elecs_eta->At(iel)) < 1.566) ||
                !elecCharge((*intVectors["Electron_tightCharge"])[iel]) ||
                !(*boolVectors["Electron_convVeto"])[iel]
                ) {
                continue;
            }
                 /*std::cout << "In Veto Elec. elecsveto_scbId : " << elecsveto_scbId->At(iel) 
                          << " elevetoid_scbcut : " << elevetoid_scbcut 
                          << " elecs_pt->At(i) : " << elecs_pt->At(iel) 
                          << " elecCharge : " << (*intVectors["Electron_tightCharge"])[iel] 
                          << " Electron_convVeto : " << (*boolVectors["Electron_convVeto"])[iel]
                          << std::endl;*/
            
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
                !elecSCBId(elecsveto_scbId->At(i), elevetoid_scbcut) ||
                (fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                !elecCharge((*intVectors["Electron_tightCharge"])[i]) ||
                !(*boolVectors["Electron_convVeto"])[i]
                ) {
                continue;
            }
            
            // Add to veto collection
            elecsveto.push_back(pre_elecs.at(i));
        }
        //std::cout << "elecsveto size in SelectVetoElectrons : "<< elecsveto.size() << std::endl;
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
                !elecSCBId(elecsveto_scbId->At(i), elevetoid_scbcut) ||
                (fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) > 1.4442 &&
                 fabs((*floatVectors["Electron_deltaEtaSC"])[i] + elecs_eta->At(i)) < 1.566) ||
                !elecCharge((*intVectors["Electron_tightCharge"])[i]) ||
                !(*boolVectors["Electron_convVeto"])[i]) {
                continue;
            }
            
            // Add to veto collection
            elecsveto.push_back(pre_elecs.at(i));
        }
    }
    else {
        std::cerr << "Error: Unrecognized decay mode in SelectVetoElectrons!" << std::endl;
    }
    //std::cout << "elecsveto : " << elecsveto.size()<< std::endl;
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
//        std::cout << " 11 --22 third_veto : " << third_veto << std::endl;
    }
    else if (TString(Decaymode).Contains("muel"))
    {
        // For muon-electron channel
        if (v_muon_idx.size() < 1 || v_electron_idx.size() < 1) {
            third_veto = false;
        }
    //    else if ((*intVectors["Muon_charge"])[v_muon_idx.at(0)] == (*intVectors["Electron_charge"])[v_electron_idx.at(0)]) {
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
        std::cerr << "Error: Unrecognized decay mode in ThirdLeptonVeto function!" << std::endl;
    }
    
    return third_veto;
}
void Analysis::GenWeightApply()
{    
    //std::cout << "start GenWieghtApply !!" << std::endl; 
    double genweight = 1.0;
    if (!isData){
        if (**floatSingles["Generator_weight"] > 0.0){genweight =1;}
        else {genweight =-1;}
        evt_weight_ = evt_weight_*genweight;
    }  
    else {evt_weight_ = 1;}
    //std::cout << "End GenWieghtApply !!" << std::endl; 
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
    //std::cout << "Selected " << nbtagged << " b-tagged jets out of " << v_jet_idx.size() << " jets" << std::endl;
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
   else {std::cerr << "ZVeto Error !!" << std::endl;}
               
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
   else {std::cerr << "METCut Error !!" << std::endl;}
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
   //std::cout << "SetUpKINObs Start ! " << std::endl;
   isKinSol=false;
   std::vector<double> jets_btag_vec;
   v_leptons_VLV.clear();
   v_jets_VLV.clear();
   v_lepidx_KIN.clear();
   v_anlepidx_KIN.clear();
   v_jetidx_KIN.clear();
   v_bjetidx_KIN.clear();
   v_btagging_KIN.clear();
   /// lepton ///
   v_leptons_VLV.push_back(common::TLVtoLV(Lep));
   v_lepidx_KIN.push_back(0);
   v_leptons_VLV.push_back(common::TLVtoLV(AnLep));
   v_anlepidx_KIN.push_back(1);

   const KinematicReconstruction* kinematicReconstruction(0); 
   kinematicReconstruction = new KinematicReconstruction(1, true);

   const LV met_LV = common::TLVtoLV(Met);
   //std::cout << "v_jet_idx.size() : " << v_jet_idx.size() << std::endl;
   
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
   
   // Debug output
   //std::cout << "v_bjetidx_KIN size: " << v_bjetidx_KIN.size() << std::endl;
   /*for (auto idx : v_bjetidx_KIN) {
      std::cout << "b-jet index: " << idx << " (valid range: 0-" << (v_jets_VLV.size()-1) << ")" << std::endl;
   }*/

   // Only proceed if we have valid b-jet indices
   if (!v_bjetidx_KIN.empty() && v_jets_VLV.size() > 0) {
      KinematicReconstructionSolutions kinematicReconstructionSolutions = 
         kinematicReconstruction->solutions(v_lepidx_KIN, v_anlepidx_KIN, v_jetidx_KIN, 
                                         v_bjetidx_KIN, v_leptons_VLV, v_jets_VLV, 
                                         jets_btag_vec, met_LV);
                                         
      if (kinematicReconstructionSolutions.numberOfSolutions()) {
         //std::cout << "Num Sol : " << kinematicReconstructionSolutions.numberOfSolutions() << std::endl;
         //std::cout << "MET ? " << met_LV.pt() << std::endl;
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
      std::cout << "Not enough b-jets or jets for kinematic reconstruction" << std::endl;
   }
   
   // Clean up
   delete kinematicReconstruction;
}
void Analysis::LeptonSFApply()
{
    //std::cout << "start ! LeptonSFApply  " << std::endl;
    lep_sf = 1.0;
    
    if (TString(Decaymode).Contains("dimuon")) {
        //std::cout << "dimuon case !" << std::endl;
        lep_sf = SSBCorr->DoubleMuon_IDIsoEff(Lep1, Lep2, LepIdSFSys, LepIsoSFSys, LepTrackSFSys);
    //    std::cout << "lep_sf " << lep_sf << std::endl;
    }
    else if (TString(Decaymode).Contains("dielec")) {
    //   lep_sf = SSBEffcal->DoubleElec_EffROOT(Lep1, Lep2,
    //                   Elec_Supercluster_Eta->at(v_lepton_idx[0]),
    //                   Elec_Supercluster_Eta->at(v_lepton_idx[1]),
    //                   LepIdSFSys, LepRecoSFSys); // LepRecoSFSys is for Electron
    }
    else if (TString(Decaymode).Contains("muel")) {
    //   lep_sf = SSBEffcal->MuonElec_EffROOT(TMuon, TElectron,
    //                   Elec_Supercluster_Eta->at(v_electron_idx[0]),
    //                   LepIdSFSys, LepIsoSFSys, LepTrackSFSys, LepRecoSFSys); // RecoSF is for Electron
    }
    else {
       lep_sf = 1.0;
       std::cout << "LeptonGetSF error !!!!" << std::endl;
    }
}
void Analysis::PUWeightApply()
{
    //std::cout << "start ! PUSFApply  " << std::endl;
     
    evt_weight_beforePileup_ = 1;
    evt_weight_beforePileup_ = evt_weight_; // keep event weight // 
    double puweight_ = 1.;
       
    if ( !TString(FileName_).Contains( "Data") )
    {  
       double pu_weight_central = SSBCorr->GetPUWeight( **floatSingles["Pileup_nTrueInt"] , PileUpSys.Data() );
       if (TString(PileUpSys).Contains("central") || TString(PileUpSys).Contains("nominal")  ) { puweight_   = pu_weight_central;}
       else {  
          std::cerr << "PUWeightApply Error... Defalut is Weight_PileUp ... : " << PileUpSys << std::endl;
       }
       evt_weight_ = evt_weight_*puweight_;  // apply PileUpReweight //
    }  
    else {evt_weight_ = 1;}
    //std::cout << "PileUp evt : " << evt_weight_ << " pu weight " << puweight_ << std::endl;
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
            triggersf_ = SSBCorr->TrigMuElec_Eff(Lep1, Lep2, TrigSFSys);
        }
        else if (TString(Decaymode).Contains("dimu")) {
            triggersf_ = SSBCorr->TrigDiMuon_Eff(Lep1, Lep2, TrigSFSys);
        }
        else {
            std::cerr << "[TriggerSFApply] WARNING: Unknown Decaymode = " << Decaymode << std::endl;
        }

        evt_weight_ *= triggersf_;  // Apply the trigger scale factor
    }
    // For Data, no correction is applied
}

// PUID related functions implementation
bool Analysis::PassPileupID(float pt, int puId, const std::string& wp) const {
    // High-pT jets don't need PUID
    if (pt > 50.0) return true;

    // For Run3 (2022, 2023) - no PUID branch exists for PUPPI jets
    if (RunPeriod.Contains("2022") || RunPeriod.Contains("2023")) {
        return true;  // No PUID needed for PUPPI jets
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

    std::cerr << "WARNING: ApplyJetPUIDEventWeights() called but PUID weights should be applied in JetSelector()!" << std::endl;
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
        std::cout << "No jets selected or btag info unavailable - B-tagging SF = 1.0" << std::endl;
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
    
    // Check if hadron flavor branch exists (MC only)
    bool hasHadronFlavour = (intVectors.count("Jet_hadronFlavour") && 
                             intVectors["Jet_hadronFlavour"] && 
                             intVectors["Jet_hadronFlavour"]->GetSize() > 0);
    
    if (!hasHadronFlavour) {
        std::cout << "WARNING: Jet_hadronFlavour branch not available - using flavor=0 for all jets" << std::endl;
    }
    
    // Collect jet information
    for (size_t i = 0; i < v_jet_idx.size(); ++i) {
        int jet_idx = v_jet_idx[i];
        
        float jet_pt = jets[i].Pt();
        float jet_eta = jets[i].Eta();
        
        // Get jet flavor (MC only)
        int flavor = 0; // Default to light flavor
        if (hasHadronFlavour && jet_idx < intVectors["Jet_hadronFlavour"]->GetSize()) {
            flavor = intVectors["Jet_hadronFlavour"]->At(jet_idx);
        }
        
        bool isTagged = (jets_btag->At(jet_idx) > bdisccut);
        
        // Fill vectors
        jet_pts.push_back(jet_pt);
        jet_etas.push_back(jet_eta);
        jet_flavors.push_back(flavor);
        jet_isTagged.push_back(isTagged);
        
        /*std::cout << "  Jet " << i << ": pT=" << jet_pt 
                  << ", eta=" << jet_eta 
                  << ", flavor=" << flavor
                  << ", isTagged=" << isTagged 
                  << " (btag_score=" << jets_btag->At(jet_idx) 
                  << ", cut=" << bdisccut << ")" << std::endl;*/
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
        btag_sf_weight_ = SSBCorr->ComputeBTagEventWeight(
            jet_pts, jet_etas, jet_flavors, jet_isTagged,
            //btag_algo, btag_wp, syst_variation
            btag_algo_, btag_wp_, syst_variation
        );
        
        // Apply weight
        evt_weight_ *= btag_sf_weight_;
        
        /*std::cout << "B-tagging SF computed: " << btag_sf_weight_ 
                  << " (total evt_weight: " << evt_weight_ << ")" << std::endl;*/
                  
    } catch (const std::exception& e) {
        std::cerr << "Error in BTaggingSFApply: " << e.what() << std::endl;
        btag_sf_weight_ = 1.0; // Default on error
    }
    
    FillHisto(h_bTagEvtWeight, btag_sf_weight_);
}

bool Analysis::IsHardScatterJet(int jet_idx) const {
    // Check if jet is matched to gen jet (CMS criteria: ΔR < 0.4)
    auto it = intVectors.find("Jet_genJetIdx");
    if (it == intVectors.end() || !it->second || 
        jet_idx >= it->second->GetSize()) {
        return false;  // No gen matching info = consider as PileUp
    }
    
    int gen_jet_idx = it->second->At(jet_idx);
    return gen_jet_idx >= 0;  // Gen matched = HardScatter, otherwise PileUp
}

/*void Analysis::CollectPUIDCandidates() {
    puid_hardscatter_jets_.clear();
    
    // Skip collection for Data or when PUID is disabled
    if (isData || !apply_puid_) {
        return;
    }
    
    // Check if jet collections are ready
    if (pre_jets.empty() || jets_pt == nullptr) {
        std::cerr << "[WARNING] CollectPUIDCandidates: Jet collections not ready" << std::endl;
        return;
    }
    
    Int_t nJets = jets_pt->GetSize();
    
    // Lambda for basic jet cuts (same as JetSelector)
    auto passBasicJetCuts = [this](int i, float jetPt, float jetEta) -> bool {
        // Kinematic cuts
        if (jetPt <= jet_pt || fabs(jetEta) >= jet_eta) return false;
        
        // Jet ID
        if (jets_Id != nullptr && jets_Id->At(i) < jet_id) return false;
        
        // Jet cleaning (need TLorentzVector)
        TLorentzVector jetVec = pre_jets[i];
        if (!JetCleaning(&jetVec)) return false;
        
        return true;
    };
    
    for (int i = 0; i < nJets; i++) {
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
        
        // Get PUID result
        int puId = (jets_puId != nullptr && i < jets_puId->GetSize()) ? jets_puId->At(i) : 0;
        bool passes_puid = PassPileupID(jetPt, puId, puid_wp_);
        
        // Store information
        puid_hardscatter_jets_.emplace_back(i, jetPt, jetEta, passes_puid, is_hardscatter);
    }
    
    //std::cout << "[PUID] Collected " << puid_hardscatter_jets_.size() 
    //          << " HardScatter candidate jets for weight calculation" << std::endl;
}*/

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
        std::cerr << "[WARNING] CollectPUIDCandidates: jets_pt is nullptr" << std::endl;
        return;
    }
    
    Int_t nJets_original = jets_pt->GetSize();
    Int_t nJets_corrected = static_cast<Int_t>(pre_jets.size());
    
    // ============================================================================
    // HANDLE SIZE MISMATCH: Use the smaller size for safety
    // ============================================================================
    Int_t nJets_safe = std::min(nJets_original, nJets_corrected);
    
    if (nJets_original != nJets_corrected) {
        std::cout << "[INFO] CollectPUIDCandidates: Size mismatch detected - "
                  << "original=" << nJets_original << ", corrected=" << nJets_corrected 
                  << ", using=" << nJets_safe << std::endl;
    }
    
    // Lambda for basic jet cuts (same as JetSelector)
    auto passBasicJetCuts = [this](int i, float jetPt, float jetEta) -> bool {
        // Kinematic cuts
        if (jetPt <= jet_pt || fabs(jetEta) >= jet_eta) return false;
        
        // Jet ID check with bounds checking
        if (jets_Id != nullptr && i < jets_Id->GetSize() && jets_Id->At(i) < jet_id) return false;
        
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
        //std::cout << "[PUID] No HardScatter candidate jets - PUID SF = 1.0" << std::endl;
        return;
    }
    
    // CMS Event Reweighting formula: P(DATA) / P(MC)
    float p_mc = 1.0;
    float p_data = 1.0;
    
    //std::cout << "[PUID] Computing event weight for " << puid_hardscatter_jets_.size() 
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
            
            // Debug output (can be removed for production)
            /*std::cout << "  Jet " << jet_info.original_index 
                      << ": pT=" << jet_info.pt << ", eta=" << jet_info.eta
                      << ", passes=" << jet_info.passes_puid 
                      << ", SF=" << sf << ", Eff=" << eff << std::endl;*/
            
        } catch (const std::exception& e) {
            std::cerr << "[ERROR] PUID SF calculation failed for jet " 
                      << jet_info.original_index << ": " << e.what() << std::endl;
            continue;  // Skip this jet on error
        }
    }
    
    // Calculate final event weight
    puid_sf_weight_ = (p_mc > 1e-10) ? p_data / p_mc : 1.0;
    
    // Apply weight to event
    evt_weight_ *= puid_sf_weight_;
    FillHisto(h_JetPUIDEvtWeight, puid_sf_weight_); 
    //std::cout << "[PUID] Event weight: " << puid_sf_weight_ 
    //          << " (p_mc=" << p_mc << ", p_data=" << p_data << ")" << std::endl;
}

