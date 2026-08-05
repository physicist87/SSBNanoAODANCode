#include "../interface/SSBCorrections.h"
#include "../TextReader/TextReader.hpp"
#include <fstream>
#include <cstdlib>
#include "correction.h"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <vector>
#include "TLorentzVector.h"



SSBCorrections::SSBCorrections(TextReader* reader, const std::string inputfileName) {
    std::cout << "TextReader in SSBCorrections ! " << std::endl;
    std::cout << "Current directory: " << std::filesystem::current_path() << std::endl;
    reader->PrintoutVariables();


    // JSON correction file root, in priority order:
    //   1. New CMS "CAT" cvmfs distribution (cms-griddata.cern.ch), organized as
    //      <POG>/<campaign>/latest/<file>.json.gz - e.g. for NanoAODv15 Puppi
    //      jets: JME/Run2-2018-UL-NanoAODv15/latest/jet_jerc.json.gz. This is the
    //      one that actually has Puppi-era (AK4PFPuppi) JEC/JER corrections.
    //   2. The older jsonpog-integration cvmfs distribution (kept for v9/PFchs
    //      running on this same code, or if the new mount isn't available).
    //   3. A local jsonpog-integration/ copy next to the executable, for laptops
    //      without cvmfs.
    // Config *Path values (JECPath, PUWeightPath, etc.) are relative to whichever
    // root is picked here, so they need to match that root's directory layout -
    // see the configs for the "<POG>/<campaign>/latest/<file>" paths used with
    // the new CAT layout. NOTE: the <campaign> folder name is NOT uniform
    // across POGs - confirmed via `ls LUM/Run*`: LUM keeps
    // "Run2-2018-UL-NanoAODv9" even when running a v15 analysis (pileup
    // profiles don't depend on the Jet collection), while JME needed a new
    // "Run2-2018-UL-NanoAODv15" folder specifically for the Puppi-era
    // corrections. Don't assume one suffix applies to every POG's path.
    std::string jsonDir;
    const std::string catCvmfsPath  = "/cvmfs/cms-griddata.cern.ch/cat/metadata/";
    const std::string oldCvmfsPath  = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";
    if (std::filesystem::exists(catCvmfsPath)) {
        std::cout << "[INFO] Using new CAT CVMFS path for JSONs: " << catCvmfsPath << std::endl;
        jsonDir = catCvmfsPath;
    } else if (std::filesystem::exists(oldCvmfsPath)) {
        std::cout << "[INFO] New CAT CVMFS path not found - using legacy jsonpog-integration CVMFS path: "
                  << oldCvmfsPath << std::endl;
        jsonDir = oldCvmfsPath;
    } else {
        std::string localPath = std::filesystem::current_path().string() + "/jsonpog-integration/POG/";
        std::cout << "[WARNING] No CVMFS JSON path found. Falling back to local path: " << localPath << std::endl;
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
        std::cerr << "[SSBCorrections] Unknown RunPeriod: " << RunPeriod << std::endl;
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

    std::cout << "era : " << era <<  std::endl;
    std::cout << "jec_name: " << jec_name << std::endl;
    // JetAlgoTag/JECLevel are optional config overrides; if absent, ExpandJECName
    // falls back to its NanoAODv15 Puppi-jet defaults ("AK4PFPuppi"/"L1L2L3Res").
    // Set JetAlgoTag="AK4PFchs" in the config to run this same code over v9/PFchs
    // input instead.
    std::string jetAlgoTag = reader->Check("JetAlgoTag") ? reader->GetText("JetAlgoTag") : "AK4PFPuppi";
    std::string jecLevelTag = reader->Check("JECLevel") ? reader->GetText("JECLevel") : "L1L2L3Res";
    // Keep the un-expanded base tag (e.g. "Summer20UL18NanoV15_V1") around so
    // we can also derive the L1FastJet-only correction name below (needed as
    // the Type-1 MET baseline - see ApplyType1METWithCorrT1) with a second
    // ExpandJECName call using jecLevel="L1FastJet" instead - jec_name itself
    // gets overwritten with the full L1L2L3Res-level name on the next line.
    std::string jecBaseNameForL1 = jec_name;
    jec_name = ExpandJECName(jec_name, RunPeriod, era, is_data, jetAlgoTag, jecLevelTag);
    std::cout << "after ExpandJECName jec_name: " << jec_name << std::endl;
    std::string jecL1Name = ExpandJECName(jecBaseNameForL1, RunPeriod, era, is_data, jetAlgoTag, "L1FastJet");
    


    std::cout << "jsonDir + puw_path : " << jsonDir + puw_path << std::endl;
    auto puset = correction::CorrectionSet::from_file(jsonDir + puw_path);
    
    pu_weight_ = puset->at(puJson);


    auto jec_set = correction::CorrectionSet::from_file(jsonDir + jec_path);
    jec_ = jec_set->compound().at(jec_name);

    // Detect whether this compound correction needs a 5th "run" input, by
    // probing with dummy values. This was discovered from the CAT
    // jet_jerc.json.gz for 2018 Puppi: the DATA compound correction
    // ("..._V1_DATA_L1L2L3Res_AK4PFPuppi") declares 5 inputs (area, eta, pt,
    // rho, run) and appears to resolve era-dependent behavior internally from
    // the run number, whereas the old v9/PFchs scheme (and the v15 MC
    // compound, "..._V1_MC_L1L2L3Res_AK4PFPuppi") only take 4 (area, eta, pt,
    // rho) and instead bake the era into the correction *name* (see
    // ExpandJECName's per-era "_RunA/_RunB/..." logic). Probing at startup
    // (rather than hardcoding is_data==true implies 5 inputs) means this
    // keeps working even if a future json goes back to the old per-era-name
    // scheme, or some other campaign/POG does something different again.
    jec_needs_run_ = false;
    try {
        jec_->evaluate({1.0, 0.0, 30.0, 10.0});
    } catch (const std::exception&) {
        try {
            jec_->evaluate({1.0, 0.0, 30.0, 10.0, 1});
            jec_needs_run_ = true;
            std::cout << "[INFO] JEC compound correction '" << jec_name
                      << "' expects a 5th 'run' input (detected automatically) - "
                      << "will pass the event run number to it." << std::endl;
        } catch (const std::exception& e2) {
            std::cerr << "[WARNING] JEC compound correction '" << jec_name
                      << "' probe failed for both the 4-input (area,eta,pt,rho) and "
                      << "5-input (+run) signatures: " << e2.what()
                      << ". JEC corrections may not evaluate correctly - check jec_name "
                      << "against your jet_jerc.json.gz." << std::endl;
        }
    }

    // L1FastJet-only correction (single, not compound) - needed as the Type-1
    // MET baseline (see ApplyType1METWithCorrT1/GetL1CorrectedJetPt). Same
    // file as jec_, just a plain (non-compound) correction lookup.
    try {
        jec_l1_ = jec_set->at(jecL1Name);
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] Could not load L1FastJet-only correction '" << jecL1Name
                  << "' from " << jec_path << ": " << e.what()
                  << ". Type-1 MET recomputation will fall back to a (fullyCorrected-raw) "
                  << "delta instead of the CMS-recommended (fullyCorrected-L1only) one - "
                  << "this is a known-less-accurate fallback, not the standard recipe."
                  << std::endl;
        jec_l1_ = nullptr;
    }

    auto jer_set = correction::CorrectionSet::from_file(jsonDir + jer_path);
    jer_ = jer_set->at(jer_res_name);

    auto jer_sf_set = correction::CorrectionSet::from_file(jsonDir + jer_sf_path);
    jer_sf_ = jer_sf_set->at(jer_name);

    // JER SF uncertainty - confirmed (by directly inspecting jet_jerc.json.gz)
    // to be a SEPARATE correction from jer_sf_ itself, not a "systematic"
    // input to it (jer_sf_'s own schema is evaluate({eta, pt}), no
    // string/tag input at all) - the CMS JERC tutorial combines them as
    // sf*(1 +/- sfUnc) for up/down. Optional: some campaigns may not publish
    // this correction, in which case up/down variations just aren't
    // available (SmearJER only uses "nominal" today regardless).
    try {
        std::string jerSfUncName = jer_name;
        size_t pos = jerSfUncName.find("ScaleFactor");
        if (pos != std::string::npos) {
            jerSfUncName.replace(pos, std::string("ScaleFactor").size(), "SFUncertainty");
        }
        jer_sfunc_ = jer_sf_set->at(jerSfUncName);
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] Could not load JER SF uncertainty correction (derived from JERName "
                  << "by replacing 'ScaleFactor' with 'SFUncertainty'): " << e.what()
                  << ". JER SF up/down variations will not be available." << std::endl;
        jer_sfunc_ = nullptr;
    }

    // CMS JME's official "JERSmear" correctionlib tool (jer_smear.json.gz) -
    // does the hybrid-method/stochastic JER smearing decision and the actual
    // random smearing internally (a `hashprng` node keyed off EventID etc.),
    // per the CMS JERC ApplicationTutorial. Path/name are config-driven since
    // this is commonly a separate file from jet_jerc.json.gz. If either key
    // is missing from the config or the file/correction can't be loaded,
    // SmearJER() falls back to an in-house TRandom3-based implementation
    // (see there) - functional, but not the officially recommended approach.
    std::string jerSmearPath = reader->Check("JERSmearPath") ? reader->GetText("JERSmearPath") : "";
    std::string jerSmearName = reader->Check("JERSmearName") ? reader->GetText("JERSmearName") : "JERSmear";
    if (!jerSmearPath.empty()) {
        try {
            auto jerSmearSet = correction::CorrectionSet::from_file(jsonDir + jerSmearPath);
            jer_smear_ = jerSmearSet->at(jerSmearName);
            std::cout << "[INFO] Loaded JER smearing tool '" << jerSmearName << "' from "
                      << jerSmearPath << " - using it for JER smearing (CMS-recommended)." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "[WARNING] Could not load JERSmearPath=" << jerSmearPath << ": " << e.what()
                      << ". Falling back to an in-house JER smearing implementation." << std::endl;
            jer_smear_ = nullptr;
        }
    } else {
        std::cerr << "[WARNING] JERSmearPath not set in config - falling back to an in-house "
                  << "JER smearing implementation instead of CMS's official JERSmear tool. "
                  << "Set JERSmearPath (e.g. 'JME/Run2-2018-UL-NanoAODv15/latest/jer_smear.json.gz') "
                  << "to use the recommended approach." << std::endl;
        jer_smear_ = nullptr;
    }

    // JMAR (PU jet ID SF) - was NOT wrapped in try/catch, unlike every other
    // optional correction load in this constructor (jec_l1_, jer_sfunc_,
    // jer_smear_, jetvetomap_ below). correction::CorrectionSet::from_file()
    // throws std::runtime_error if the file doesn't exist, and with no
    // try/catch that propagates all the way out of the constructor uncaught -
    // std::terminate()/abort(), crashing the whole job, instead of a graceful
    // [WARNING] + fallback like everything else here. GetPUJetIDSFAndEff()
    // already null-checks pujetid_sf_ and returns 1.0 (no PU-ID reweighting)
    // if it's not loaded, so there was already a safe fallback available -
    // this constructor just never used it. Found from a real run where
    // jmar.json.gz didn't exist at the expected CAT cvmfs path (2018 UL JMAR
    // apparently doesn't have a NanoAODv15 folder there, or is named/located
    // differently - worth an `ls` on your cvmfs to confirm the real path).
    try {
        auto jmar_sf_set = correction::CorrectionSet::from_file(jsonDir + jmar_path);
        pujetid_sf_ = jmar_sf_set->at("PUJetID_eff");
        std::cout << "[INFO] Loaded PU jet ID SF from " << jmar_path << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] Could not load JMARPath=" << jmar_path << ": " << e.what()
                  << ". PU jet ID SF will be unavailable (GetPUJetIDSFAndEff() returns 1.0 - "
                  << "no PU-ID reweighting applied)." << std::endl;
        pujetid_sf_ = nullptr;
    }

    jveto_name_ = jveto_name;
    jveto_key_ = jveto_map_key;
    jveto_type_ = jveto_type;	

    if (!jveto_path.empty() && !jveto_name.empty()) {
       try {
           auto jetveto_set = correction::CorrectionSet::from_file(jsonDir + jveto_path);
           jetvetomap_ = jetveto_set->at(jveto_name);  // jveto_name 사용
           std::cout << "[INFO] Loaded jet veto map: " << jveto_name
                     << " with key: " << jveto_map_key
                     << " type: " << jveto_type << std::endl;
       } catch (const std::exception& e) {
           std::cerr << "[WARNING] Failed to load jet veto map: " << e.what() << std::endl;
           jetvetomap_ = nullptr;
       }
    } else {
       std::cout << "[INFO] Jet veto map not configured, skipping..." << std::endl;
       jetvetomap_ = nullptr;
    }

    if (btag_sf_type_ != "comb" && btag_sf_type_ != "mujets") {
        std::cerr << "[WARNING] Invalid BTagSFType: " << btag_sf_type_ 
                  << ". Using default 'comb'." << std::endl;
        btag_sf_type_ = "comb";
    }

    // Analysis type specific recommendation
    if (btag_sf_type_ == "comb") {
        std::cout << "[INFO] Using 'comb' SF (QCD + ttbar enriched regions)" << std::endl;
    } else {
        std::cout << "[INFO] Using 'mujets' SF (QCD enriched regions, bias avoidance)" << std::endl;
    }


    std::string btag_algo = "";
    std::string btag_wp = "";

    if (jet_btag_conf.find("deepCSV") != std::string::npos) {
        btag_algo = "DeepCSV";
    } else if (jet_btag_conf.find("deepJet") != std::string::npos) {
        btag_algo = "DeepJet";
    } else if (jet_btag_conf.find("UParT") != std::string::npos) {
        // UParT (Unified Particle Transformer) AK4 b-tagger, new in NanoAODv15.
        btag_algo = "UParT";
    } else if (jet_btag_conf.find("pfCSVV2") != std::string::npos) {
        btag_algo = "CSVv2";
    } else {
        std::cerr << "[WARNING] Unknown b-tag algorithm in Jet_btag: " << jet_btag_conf << std::endl;
    }

    char last = jet_btag_conf.back();
    if (last == 'L' || last == 'l') btag_wp = "Loose";
    else if (last == 'M' || last == 'm') btag_wp = "Medium";
    else if (last == 'T' || last == 't') btag_wp = "Tight";
    else {
        std::cerr << "[WARNING] Unknown WP in Jet_btag: " << jet_btag_conf << std::endl;
    }

    // Path scheme: config BTagSFJsonPath overrides everything (safest, since the
    // new CAT cvmfs layout for BTV hasn't been directly confirmed - only JME's
    // "Run2-<year>-UL-NanoAODv15/latest/" convention was confirmed from the
    // user's cvmfs listing). Otherwise, guess the same campaign-folder scheme
    // JME uses when running against the new CAT jsonDir, and fall back to the
    // legacy "<year>_UL/" scheme for the old jsonpog-integration path.
    std::string btag_sf_json;
    if (reader->Check("BTagSFJsonPath")) {
        btag_sf_json = reader->GetText("BTagSFJsonPath");
    } else if (jsonDir == catCvmfsPath) {
        btag_sf_json = "BTV/Run2-" + year_ + "-UL-NanoAODv15/latest/btagging.json.gz";
    } else {
        btag_sf_json = "BTV/" + year_ + "_UL/btagging.json.gz";
    }

    // correctionlib correction-name string for the b-tag SF, keyed off btag_sf_type_
    // ("comb"/"mujets" - only meaningful for DeepJet/DeepCSV in the CMS BTV scheme).
    // UParT's correction name in btagging.json.gz is NOT verified here - it may not
    // even follow the "<tagger>_<sftype>" pattern used by DeepJet/DeepCSV. Override
    // via config key BTagTaggerName if "particleTransformerAK4_<btag_sf_type_>" turns
    // out to be wrong; check with:
    //   python3 -c "import correctionlib; print(list(correctionlib.CorrectionSet.from_file('btagging.json.gz')))"
    std::string btag_tagger;
    if (reader->Check("BTagTaggerName")) {
        btag_tagger = reader->GetText("BTagTaggerName");
    } else if (btag_algo == "DeepJet") {
        btag_tagger = "deepJet_" + btag_sf_type_;
    } else if (btag_algo == "UParT") {
        // Confirmed from the user's actual btagging.json.gz: the corrections
        // are named "UParTAK4_comb" (b/c-jet SF) and "UParTAK4_light"
        // (light-jet SF) - NOT "particleTransformerAK4_<type>" as originally
        // guessed. Also unlike DeepJet/DeepCSV, UParT only ships a single
        // "comb"-type heavy-flavor SF (no separate "..._mujets" correction),
        // per its description: "Should not be used by ttbar measurements
        // with large overlap of SF derivation phase spaces (emu, >=2j and
        // e/mu==4j)" - i.e. CMS BTV's own docs flag exactly this dilepton
        // ttbar analysis as a case to be careful with when using 'comb'.
        // BTagSFType is therefore not applied to the tagger name here (see
        // InitBtagSFCorrection, which special-cases the UParTAK4_* naming).
        if (btag_sf_type_ != "comb") {
            std::cout << "[INFO] btag_algo=UParT: BTagSFType='" << btag_sf_type_
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
        // File-naming ONLY fix: your actual efficiency ROOT files are named
        // "btagEff_UParTAK4.root", not "btagEff_UParT.root" (confirmed via
        // `ls CorrectionFiles/BTag/UL2018/*/`). Using a separate
        // btag_eff_filename_algo here (rather than changing btag_algo
        // itself) intentionally, because btag_algo is ALSO used as the
        // in-ROOT-file histogram name prefix ("eff_" + btag_algo + "_b_...")
        // via LoadMCBtagEfficiencies/GetMCBtagEfficiency/
        // ComputeBTagEventWeight (and Analysis.cpp's separate btag_algo_
        // member, which matches this same "UParT" convention) - that's an
        // internal, self-consistent naming scheme independent of the
        // external .root filename, and I have no way to confirm here
        // whether the histograms INSIDE btagEff_UParTAK4.root are named
        // "eff_UParT_b_Medium" or "eff_UParTAK4_b_Medium". If you see
        // "[WARNING] Histogram not found: eff_UParT_..." after this fix,
        // the histograms themselves use the AK4 suffix too - tell me and
        // I'll change btag_algo (both here and in Analysis.cpp) to match.
        std::string btag_eff_filename_algo = (btag_algo == "UParT") ? "UParTAK4" : btag_algo;
        btag_eff_path = "CorrectionFiles/BTag/UL" + RunPeriod
                      + "/" + process_subdir
                      + "/btagEff_" + btag_eff_filename_algo + ".root";
        std::cout << "[INFO] btag_eff_path: " << btag_eff_path << std::endl;
    }

    InitBtagSFCorrection(jsonDir + btag_sf_json, btag_tagger);
    if (!is_data) {
        //LoadMCBtagEfficiencies(btag_eff_path, btag_algo);
        LoadMCBtagEfficiencies(btag_eff_path, btag_algo, btag_wp);

    }


    // Load muon SF
    //auto muon_set = CorrectionSet::from_file(jsonDir+muon_path);
    std::cout << "jsonDir+muon_path " << jsonDir+muon_path << std::endl; 
    auto muon_set = correction::CorrectionSet::from_file(jsonDir + muon_path);

    std::cout << "muon_id_corName : " << muon_id_corName << std::endl;
    std::cout << "muon_iso_corName : " << muon_iso_corName << std::endl;

    muon_id_   = muon_set->at(muon_id_corName);
    muon_iso_  = muon_set->at(muon_iso_corName);

    auto cset = correction::CorrectionSet::from_file(jsonDir + elec_path);
    //std::cout << "sk ele 1 ele_sf_name_ : " << ele_sf_name_ << std::endl;
    ele_sf_ = cset->at(ele_sf_name_);
    //std::cout << "sk ele 2 ele_reco_sf_name_: " << ele_reco_sf_name_ << std::endl;
    std::string TrigSFPath = std::filesystem::current_path().string() + "/CorrectionFiles/Trig/";
    TFile *f_trg     = new TFile((TrigSFPath+Trig_sf_name_).c_str());
    H_trig = (TH2D*) f_trg->Get(Trig_sf_histname_.c_str());

    /// Muon Rochester correction ///
    std::string rochesterCorrFile;
    std::cout << "RunPeriod : " << RunPeriod << std::endl;
    if (RunPeriod.find("2016Pre")  != std::string::npos) {    rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2016aUL.txt";}
    else if (RunPeriod.find("2016Post") != std::string::npos) { rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2016bUL.txt";}
    else if (RunPeriod.find("2017")     != std::string::npos) { rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2017UL.txt";}
    else if (RunPeriod.find("2018")     != std::string::npos) { rochesterCorrFile =  "./CorrectionFiles/Rochester/RoccoR2018UL.txt";}
    else {
        std::cerr << "[SSBCorrections] Unknown RunPeriod in rochesterCorrFile : " << RunPeriod << std::endl;
    }

    std::cout << "[Debug!!] in SSBCorrections !! rochesterCorrFile  : " << rochesterCorrFile<< std::endl;
    std::ifstream file(rochesterCorrFile);
    if(!file.good()){
            std::cerr << "ERROR:  Rochester file is not found:" << rochesterCorrFile << std::endl;
            std::exit(1);
    }

    rc.init(rochesterCorrFile);
}

// Destructor ///
//
SSBCorrections::~SSBCorrections() {
    std::cout << "[SSBCorrections] Destructor called - cleaning up memory..." << std::endl;
    
    // 1. Safely delete H_trig histogram
    if (H_trig) {
        std::cout << "[SSBCorrections] Deleting H_trig histogram..." << std::endl;
        delete H_trig;
        H_trig = nullptr;
    }
    
    // 2. Safely delete all TH2D pointers in eff_histograms_ map
    std::cout << "[SSBCorrections] Cleaning up " << eff_histograms_.size() 
              << " efficiency histograms..." << std::endl;
              
    for (auto& pair : eff_histograms_) {
        if (pair.second) {
            std::cout << "[SSBCorrections] Deleting histogram: " << pair.first << std::endl;
            delete pair.second;
            pair.second = nullptr;
        }
    }
    
    // Clear the map itself
    eff_histograms_.clear();
    
    std::cout << "[SSBCorrections] Destructor completed successfully." << std::endl;
}
/*
SSBCorrections::~SSBCorrections() {
    std::cout << "[SSBCorrections] Destructor called - cleaning up memory..." << std::endl;

    if (H_trig) {
        // Check if ROOT object is still valid
        if (gROOT && gROOT->FindObject(H_trig)) {
            std::cout << "[SSBCorrections] Deleting H_trig histogram..." << std::endl;
            delete H_trig;
        }
        H_trig = nullptr;
    }

    std::cout << "[SSBCorrections] Cleaning up " << eff_histograms_.size()
              << " efficiency histograms..." << std::endl;

    for (auto& pair : eff_histograms_) {
        if (pair.second) {
            // Check if ROOT object is still valid
            if (gROOT && gROOT->FindObject(pair.second)) {
                std::cout << "[SSBCorrections] Deleting histogram: " << pair.first << std::endl;
                delete pair.second;
            }
            pair.second = nullptr;
        }
    }

    // Clear the map itself
    eff_histograms_.clear();

    std::cout << "[SSBCorrections] Destructor completed successfully." << std::endl;
}
*/

double SSBCorrections::GetCorrectedJetPt(double raw_pt, double eta, double area, double rho, unsigned int run_number) const {
    //std::cout << "in GetCorrectedJetPt raw_pt : " << raw_pt
    //          << " eta " << eta << " area : " << area << " rho : " << rho << std::endl;

    double sf = jec_needs_run_
        ? jec_->evaluate({area, eta, raw_pt, rho, static_cast<int>(run_number)})
        : jec_->evaluate({area, eta, raw_pt, rho});
    //std::cout << "sf in GetCorrectedJetPt : " << sf << std::endl;

    return raw_pt * sf;
}

double SSBCorrections::GetL1CorrectedJetPt(double raw_pt, double eta, double area, double rho) const {
    if (!jec_l1_) {
        // No L1-only correction loaded (see constructor warning) - returning
        // raw_pt here makes ApplyType1METWithCorrT1's (corrected-L1only)
        // delta fall back to (corrected-raw), i.e. the old/less-accurate
        // behavior, rather than silently producing a wrong number some other way.
        return raw_pt;
    }
    double c1 = jec_l1_->evaluate({area, eta, raw_pt, rho});
    return raw_pt * c1;
}

double SSBCorrections::GetCorrectedJetMass(double raw_mass, double raw_pt, double eta, double area, double rho, unsigned int run_number) const {
    //double sf = jec_->evaluate({eta, raw_pt, area});
    double sf = jec_needs_run_
        ? jec_->evaluate({area, eta, raw_pt, rho, static_cast<int>(run_number)})
        : jec_->evaluate({area, eta, raw_pt, rho});
    return raw_mass * sf;
}

double SSBCorrections::GetJER(double eta, double pt, double rho) const {
    return jer_->evaluate({eta, pt, rho});
}

namespace {
// Deterministic per-(event, jet) seed for JER stochastic smearing, so that:
//  (a) results are reproducible across job re-runs/resubmissions, and
//  (b) if JER-up/down systematic passes are added later, they can reuse the
//      exact same random draw per jet (only the resolution/SF value differs
//      between nominal/up/down), which is the standard way to isolate the
//      systematic from jet-to-jet statistical fluctuation.
// This mixes the event number with the jet's eta/phi (scaled to integers) so
// different jets in the same event get different but reproducible seeds.
// NOTE: this is a reasonable in-house deterministic scheme, not necessarily
// bit-identical to CMSSW/nanoAOD-tools' internal JER smearer seed formula -
// if bit-exact reproducibility with an existing official recipe is required,
// substitute that exact formula here instead.
UInt_t ComputeJerSeed(ULong64_t event, double eta, double phi) {
    ULong64_t etaBits = static_cast<ULong64_t>(std::llround(eta * 1.0e4));
    ULong64_t phiBits = static_cast<ULong64_t>(std::llround(phi * 1.0e4));
    ULong64_t mixed = (event * 2654435761ULL)
                     ^ (etaBits * 40503ULL)
                     ^ (phiBits * 2246822519ULL);
    // Avoid seed 0 (TRandom3 treats 0 as "seed from clock/entropy", which
    // would defeat determinism).
    UInt_t seed = static_cast<UInt_t>(mixed & 0xFFFFFFFFULL);
    return seed == 0 ? 1u : seed;
}
}

double SSBCorrections::SmearJER(double reco_pt, double gen_pt, double eta, double phi, double rho,
                                 ULong64_t event, const std::string& jer_tag) const {
    // jer_sf_/jer_sfunc_ evaluate({eta, pt}) - confirmed by directly
    // inspecting jet_jerc.json.gz: this is a plain 2-real-input schema, NOT
    // {eta, systematic_string} as originally assumed here. Up/down variations
    // come from a SEPARATE SFUncertainty correction, combined arithmetically
    // (sf*(1+-unc)) - matches the CMS JERC tutorial's Applier::jerFactor.
    double sf = jer_sf_->evaluate({eta, reco_pt});
    if (jer_sfunc_ && jer_tag != "nominal") {
        double sfUnc = jer_sfunc_->evaluate({eta, reco_pt});
        if (jer_tag == "up")        sf = sf * (1.0 + sfUnc);
        else if (jer_tag == "down") sf = sf * (1.0 - sfUnc);
    }
    double resolution = jer_->evaluate({eta, reco_pt, rho});

    if (jer_smear_) {
        // CMS's official "JERSmear" correctionlib tool - handles the hybrid
        // method decision (scaling vs. stochastic, including its own
        // 3-sigma gen-pt consistency check) and the actual random smearing
        // internally via a deterministic `hashprng` node. Input order per
        // the CMS JERC ApplicationTutorial's Applier::jerFactor: JetPt,
        // JetEta, GenPt (-1 if no match), Rho, EventID, JER (resolution), JERSF.
        double smear = jer_smear_->evaluate({reco_pt, eta, gen_pt, rho,
                                              static_cast<int>(event), resolution, sf});
        double corr = (std::isfinite(smear) && smear > 0.0) ? smear : 1.0;
        return std::max(0.0, reco_pt * corr);
    }

    // Fallback (jer_smear_ not loaded - see constructor warning): an in-house
    // reimplementation of the same hybrid method. Not the officially
    // recommended approach (no access to correctionlib's own hashprng
    // internals), but functional and deterministic via a local per-(event,
    // jet) seeded TRandom3 rather than the shared global gRandom.
    bool useScaling = (gen_pt >= 0.0) &&
                      (std::abs(reco_pt - gen_pt) < 3.0 * resolution * reco_pt);

    if (useScaling) {
        double delta_pt = reco_pt - gen_pt;
        double smeared_pt = gen_pt + sf * delta_pt;
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
        std::cerr << "[SSBCorrections::GetPUJetIDSFAndEff] PUJetID correction not loaded." << std::endl;
        return 1.0;
    }
    
    if (pt >= 50.0) return 1.0;
    
    try {
        std::string eval_type = getEff ? "MCEff" : syst;  // Efficiency or Scale Factor
        std::variant<double, std::vector<double>> val = pujetid_sf_->evaluate({eta, pt, eval_type, wp});
        return std::get<double>(val);
    } catch (const std::exception& e) {
        std::cerr << "[SSBCorrections::GetPUJetIDSFAndEff] Evaluation failed: " << e.what() << std::endl;
        return 1.0;
    }
}


double SSBCorrections::GetMuonRecoSF(double pt, double eta) const {
    //auto val = muon_reco_->evaluate({pt, eta});
    std::variant<double, std::vector<double>> val = muon_reco_->evaluate({pt, eta});
    //return std::get<double>(val);
    return std::get<double>(val);
}

double SSBCorrections::GetMuonIDSF(double pt, double eta, const std::string& syst_tag) const {
    //auto val = muon_id_->evaluate({eta, pt, syst_tag});
    std::variant<double, std::vector<double>> val = muon_id_->evaluate({eta, pt, syst_tag});
    return std::get<double>(val);
}

double SSBCorrections::GetMuonIsoSF(double pt, double eta, const std::string& syst_tag) const {
    //auto val = muon_iso_->evaluate({eta, pt, syst_tag});
    std::variant<double, std::vector<double>> val = muon_iso_->evaluate({eta, pt, syst_tag});
    return std::get<double>(val);
}


double SSBCorrections::DoubleMuon_IDIsoEff(TLorentzVector lep1, TLorentzVector lep2,
                                           TString muidsys, TString muisosys, TString tracksys) const {
    float pt1 = std::min(lep1.Pt(), 119.999);
    float pt2 = std::min(lep2.Pt(), 119.999);
    float abseta1 = std::abs(lep1.Eta());
    float abseta2 = std::abs(lep2.Eta());

    //std::cout << "muidsys : " << muidsys << std::endl;
    //std::cout << "muisosys : " << muisosys << std::endl;

    std::string IDSyst = "nominal", IsoSyst = "nominal";

    if (muidsys.Contains("up", TString::kIgnoreCase)) IDSyst = "up";
    else if (muidsys.Contains("down", TString::kIgnoreCase)) IDSyst = "down";

    if (muisosys.Contains("up", TString::kIgnoreCase)) IsoSyst = "up";
    else if (muisosys.Contains("down", TString::kIgnoreCase)) IsoSyst = "down";
    //std::cout << "IDSyst : " << IDSyst << std::endl;
    double mu1id  = GetMuonIDSF(pt1, abseta1, IDSyst);
    //std::cout << "After mu1id  " << mu1id << std::endl;
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
    
    // Step 1: Apply pt/eta limits based on JSON ranges
    // pT range in JSON: 10.0 to Infinity, but clamp very high values
    float lep1pt = std::clamp(static_cast<float>(lep1.Pt()), 10.0f, 999.0f);
    float lep2pt = std::clamp(static_cast<float>(lep2.Pt()), 10.0f, 999.0f);
    
    // Eta range in JSON: -Infinity to +Infinity, but avoid extreme values  
    // Keep within reasonable detector range
    float lep1sueta_clamped = std::clamp(static_cast<float>(ele1sueta), -3.0f, 3.0f);
    float lep2sueta_clamped = std::clamp(static_cast<float>(ele2sueta), -3.0f, 3.0f);
    
    // Step 2: Calculate ID SF for each electron using GetElectronSF
    // Check if id_wp is empty and set default
    std::string actual_id_wp = id_wp;
    if (actual_id_wp.empty()) {
        actual_id_wp = "Tight";  // Set default working point
        std::cout << "[WARNING] Empty ID working point, using default: " << actual_id_wp << std::endl;
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
    
    // Debug output (can be removed in production)
    /*
    std::cout << "[DoubleElec_Eff] Debug info:" << std::endl;
    std::cout << "  Ele1: pt=" << lep1pt << ", eta=" << lep1sueta_clamped 
              << ", ID_SF=" << ele1id << ", Reco_SF=" << ele1reco << std::endl;
    std::cout << "  Ele2: pt=" << lep2pt << ", eta=" << lep2sueta_clamped 
              << ", ID_SF=" << ele2id << ", Reco_SF=" << ele2reco << std::endl;
    std::cout << "  Total efficiency: " << doubleEleff << std::endl;
    */ 
    
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
        std::cout << "[WARNING] Empty electron ID working point, using default: " << actual_ele_id_wp << std::endl;
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

    // Debug output (can be removed in production)
    /*
    std::cout << "[MuonElec_Eff] Debug info:" << std::endl;
    std::cout << "  Muon: pt=" << muon_pt << ", abseta=" << muon_abseta 
              << ", ID_SF=" << mu_id << ", Iso_SF=" << mu_iso << ", Trk_SF=" << mu_trk << std::endl;
    std::cout << "  Electron: pt=" << electron_pt << ", sueta=" << electron_sueta_clamped 
              << ", ID_SF=" << ele_id << ", Reco_SF=" << ele_reco << ", Iso_SF=" << ele_iso << std::endl;
    std::cout << "  Total efficiency: " << muonelec_eff << std::endl;
    */

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
        std::cerr << "[WARNING] Unknown systematic: '" << syst << "', using nominal (sf)" << std::endl;
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
            std::cerr << "[WARNING] Unknown SCB working point: " << sf_type << ", using Tight" << std::endl;
            working_point = "Tight";
        }
    } else if (sf_type.find("MVA") == 0) {
        // Convert MVA format to JSON format if needed
        if (sf_type == "MVALoose") working_point = "wp90noiso";
        else if (sf_type == "MVAMedium") working_point = "wp80noiso";
        else if (sf_type == "MVATight") working_point = "wp90iso";
        else {
            std::cerr << "[WARNING] Unknown MVA working point: " << sf_type << ", using wp90iso" << std::endl;
            working_point = "wp90iso";
        }
    }
    // For other cases (direct JSON names like "Tight", "Medium"), use as-is

    try {
        if (!ele_sf_) {
            std::cerr << "[ERROR] ele_sf_ pointer is null!" << std::endl;
            return 1.0;
        }

        // Validate parameters before evaluation
        if (year_.empty()) {
            std::cerr << "[ERROR] Year is empty!" << std::endl;
            return 1.0;
        }

        if (working_point.empty()) {
            std::cerr << "[ERROR] Working point is empty!" << std::endl;
            return 1.0;
        }

        // Use the correct order from JSON: {year, ValType, WorkingPoint, eta, pt}
        float result = ele_sf_->evaluate({year_, valtype, working_point, eta, pt});

        // Sanity check on result
        if (result <= 0.0 || result > 10.0) {
            std::cerr << "[WARNING] Unusual SF value: " << result
                      << " for parameters: " << working_point << ", eta=" << eta << ", pt=" << pt << std::endl;
        }

        return result;

    } catch (const std::exception& e) {
        std::cerr << "[GetElectronSF] Evaluation failed: " << e.what() << std::endl;
        std::cerr << "  Parameters: year='" << year_ << "', valtype='" << valtype
                  << "', working_point='" << working_point << "', eta=" << eta << ", pt=" << pt << std::endl;
        return 1.0;
    }
}

float SSBCorrections::GetPUWeight(float nTrueInt, const std::string& systTag) const {
    std::string variation = systTag;
    
    if (!pu_weight_) {
        std::cerr << "[SSBCorrections::GetPUWeight] PU correction not loaded.\n";
        return 1.0;
    }
    
    // If systTag is "Central" or empty, treat as "nominal"
    if (systTag == "Central" || systTag == "") {
        variation = "nominal";
    }
    
    //std::cout << "variation in PU " << variation << std::endl;

    // If the variation is not valid, fallback to "nominal"
    if (variation != "nominal" && variation != "up" && variation != "down") {
        std::cerr << "PileUp sys Error ... Defaulting to Weight_PileUp ... Original value: " << variation << std::endl;
        variation = "nominal";
    }
    
    //std::cout << "nTrueInt : " << nTrueInt << std::endl;
    
    try {
        std::variant<double, std::vector<double>> val = pu_weight_->evaluate({ nTrueInt, variation });
        double weight = std::get<double>(val);
        //std::cout << "PileUp weight: " << weight << std::endl;  
        return weight;
    } catch (const std::exception& e) {
        std::cerr << "[SSBCorrections::GetPUWeight] Evaluation failed: " << e.what() << std::endl;
        return 1.0;
    }
}


float SSBCorrections::MatchGenPt(const TLorentzVector& reco_jet,
                                  const std::vector<TLorentzVector>& gen_jets,
                                  float maxDR) const {
    float minDR = maxDR;
    float matched_genpt = -1.0;

    for (const auto& gen_jet : gen_jets) {
        float dR = reco_jet.DeltaR(gen_jet);
        if (dR < minDR) {
            minDR = dR;
            matched_genpt = gen_jet.Pt();
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
    ULong64_t event_number
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

        // Snapshot the pt BEFORE JER smearing so the mass-rescaling below only
        // folds in the JER-specific ratio, not the JEC factor a second time.
        // (Bug found 2026-08: this used to divide by raw_pt instead of
        // pt_before_jer, so corrected_mass ended up carrying jec_sf^2 *
        // jer_ratio instead of jec_sf * jer_ratio - i.e. the JEC scale factor
        // was silently squared into the jet mass whenever applyJES was on,
        // pre-dating this session's other JEC/JER fixes.)
        double pt_before_jer = corrected_pt;

        if (!isData && applyJER) {
            // Always smear (previously this block only ran when a valid gen
            // match existed, silently skipping JER entirely otherwise). Pass
            // -1.0 as the "no gen match" sentinel when Jet_genJetIdx has no
            // valid entry for this jet - SmearJER() falls back to stochastic
            // smearing in that case rather than leaving the jet unsmeared.
            float matched_genpt = -1.0;
            if (genJetIndices.size() > i && genJetIndices[i] >= 0 &&
                static_cast<size_t>(genJetIndices[i]) < genJets.size()) {
                matched_genpt = genJets[genJetIndices[i]].Pt();
            }
            corrected_pt = SmearJER(corrected_pt, matched_genpt, eta, phi, rho, event_number, "nominal");
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
    ULong64_t event_number
) const {
    double met_px = raw_met_pt * std::cos(raw_met_phi);
    double met_py = raw_met_pt * std::sin(raw_met_phi);

    // Shared per-jet logic, confirmed against the CMS JERC tutorial's
    // Applier::correctedMet(): muon-subtract the raw pt, apply JES/JER
    // (reusing GetCorrectedJetPt/GetL1CorrectedJetPt/SmearJER so this can't
    // silently drift out of sync with the physics-jet path), require the
    // Type-1 selection (corrected pt>15, |eta|<5.2, chEmEF+neEmEF<0.90), and
    // accumulate the (L1only_no_mu - corrected_no_mu) differential into MET -
    // NOT (raw_no_mu - corrected_no_mu), which was the bug in an earlier
    // version of this function (and was also present in the old fused
    // ApplyJetCorrectionsWithMET/RecomputeMET path, now removed entirely -
    // this function is the sole MET computation).
    auto accumulate = [&](double eta, double phi, double rawPtRaw, double area,
                          double muonSubtrFactor, double chEmEF, double neEmEF) {
        double rawPtNoMu = rawPtRaw * (1.0 - muonSubtrFactor);
        if (rawPtNoMu <= 0.0) return;

        double l1PtNoMu = rawPtNoMu;
        double corrPtNoMu = rawPtNoMu;
        if (applyJES) {
            l1PtNoMu   = GetL1CorrectedJetPt(rawPtNoMu, eta, area, rho);
            corrPtNoMu = GetCorrectedJetPt(rawPtNoMu, eta, area, rho, run_number);
        }
        if (!isData && applyJER) {
            TLorentzVector recoJet4v;
            recoJet4v.SetPtEtaPhiM(rawPtNoMu, eta, phi, 0.0);
            // -1.0 sentinel handled inside SmearJER as "no gen match" (falls
            // back to stochastic smearing) if MatchGenPt finds nothing within maxDR.
            // NOTE (from the tutorial's own README to-do list): CorrT1METJet_
            // has no genJetIdx branch, so this dR-based rematch is the same
            // acknowledged-imperfect approach the official tutorial itself
            // uses for those jets - not something specific to our code.
            float matchedGenPt = MatchGenPt(recoJet4v, genJets, 0.2f);
            corrPtNoMu = SmearJER(corrPtNoMu, matchedGenPt, eta, phi, rho, event_number, "nominal");
        }

        // Type-1 selection cut (CMS JERC tutorial's JvmApplication-adjacent
        // correctedMet() logic) - only jets passing this contribute to MET.
        bool passSel = (corrPtNoMu > 15.0) && (std::abs(eta) < 5.2) && ((chEmEF + neEmEF) < 0.90);
        if (!passSel) return;

        met_px += (l1PtNoMu - corrPtNoMu) * std::cos(phi);
        met_py += (l1PtNoMu - corrPtNoMu) * std::sin(phi);
    };

    // Regular Jet_ branch jets: jetsAsStored holds the NanoAOD-stored
    // (already centrally-corrected) pt, so undo rawFactor first to get the
    // true raw pt - same convention as ApplyJetCorrections.
    for (size_t i = 0; i < jetsAsStored.size(); ++i) {
        double rawFactor = (i < jetRawFactors.size()) ? jetRawFactors[i] : 0.0f;
        double area       = (i < jetAreas.size())       ? jetAreas[i]       : 0.5f;
        double muonSubtr  = (i < jetMuonSubtrFactors.size()) ? jetMuonSubtrFactors[i] : 0.0f;
        double chEmEF     = (i < jetChEmEF.size())      ? jetChEmEF[i]      : 0.0f;
        double neEmEF     = (i < jetNeEmEF.size())      ? jetNeEmEF[i]      : 0.0f;
        double rawPtRaw   = jetsAsStored[i].Pt() * (1.0 - rawFactor);
        accumulate(jetsAsStored[i].Eta(), jetsAsStored[i].Phi(), rawPtRaw, area, muonSubtr, chEmEF, neEmEF);
    }

    // Low-pT CorrT1METJet_ branch jets: rawPt is already raw (no rawFactor to
    // undo). No EM-fraction branches exist for these (too low-pT to have
    // them produced) - pass 0/0, same as the CMS JERC tutorial's
    // CollectJetMet.hpp ("EM fractions are not provided for CorrT1METJet;
    // set to zero explicitly").
    for (size_t i = 0; i < corrT1RawPt.size(); ++i) {
        double area      = (i < corrT1Area.size())               ? corrT1Area[i]               : 0.5f;
        double muonSubtr = (i < corrT1MuonSubtrFactor.size())    ? corrT1MuonSubtrFactor[i]    : 0.0f;
        accumulate(corrT1Eta[i], corrT1Phi[i], corrT1RawPt[i], area, muonSubtr, 0.0, 0.0);
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
        std::cerr << "[WARNING] Jet veto map not loaded, skipping HEM veto" << std::endl;
        return false;
    }

    // Pre-selection cut confirmed from the CMS JERC tutorial's
    // JvmApplication::VetoChecker (kMaxEmFrac = 0.90) - pt/jetId pre-selection
    // are already applied at the only call site (Analysis::JetSelector)
    // before ShouldVetoJet() is invoked, so only the EM-fraction cut is
    // checked here.
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
        
        // Debug output to see what values we're getting
//        std::cout << "[DEBUG] Jet eta=" << eta << ", phi=" << phi 
//                  << " -> veto_flag=" << veto_flag << std::endl;
        
        // Return true if jet should be vetoed (non-zero value)
        bool should_veto = (veto_flag > 0.0);
        
//        if (should_veto) {
//            std::cout << "[DEBUG] -> VETOED (flag=" << veto_flag << ")" << std::endl;
//        }
        
        return should_veto;
        
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] Jet veto map evaluation failed for jet (eta=" 
                  << eta << ", phi=" << phi << ") with key=" << jveto_key_ 
                  << ": " << e.what() << std::endl;
        return false;
    }
}


void SSBCorrections::InitBtagSFCorrection(const std::string& json_path, 
                                          const std::string& tagger_name) {
    std::cout << "[SSBCorrections] Loading b-tagging SF from JSON: " << json_path << std::endl;
    
    auto cset = correction::CorrectionSet::from_file(json_path);

    // UParT's btagging.json.gz (confirmed by directly inspecting the file)
    // only ships ONE heavy-flavor SF, named "UParTAK4_comb" (no per-BTagSFType
    // "..._mujets" variant the way DeepJet/DeepCSV do), and names its
    // light-flavor SF "UParTAK4_light" instead of "..._incl". Detected here by
    // the tagger_name prefix, since that already fully identifies which
    // scheme is in play (constructor already built it as exactly
    // "UParTAK4_comb" for this case - see there).
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
            std::cerr << "[ERROR] Expected 'comb' in tagger_name: " << tagger_name << std::endl;
            return;
        }
    } else {
        // Heavy flavor correction name (config-based selection)
        size_t pos = heavy_flavor_name.find("comb");
        if (pos != std::string::npos) {
            heavy_flavor_name.replace(pos, 4, btag_sf_type_);
        } else {
            std::cerr << "[ERROR] Expected 'comb' in tagger_name: " << tagger_name << std::endl;
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
    
    std::cout << "[SSBCorrections] Loaded corrections: " 
              << heavy_flavor_name << " and " << light_flavor_name << std::endl;
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
        std::cerr << "[ERROR] B-tag correction not found: " << corr_name << std::endl;
        return 1.0;
    }

    try {
        return it->second->evaluate({syst, wp, flav, eta_abs, pt_clamped});
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] GetBtagSF failed: " << e.what() << std::endl;
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
        std::cout << "Input vectors have different sizes" << std::endl;
        return 1.0;
    }
    //std::cout << "start ! " << std::endl; 
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
        
        /*std::cout << "  Jet " << i << ": pt=" << pt << ", eta=" << eta 
                  << ", flav=" << flav << ", tagged=" << tagged 
                  << ", SF=" << sf << ", Eff=" << eff 
                  << ", p_mc=" << p_mc << ", p_data=" << p_data << std::endl;*/
    }
    
    // Calculate final event weight
    float btag_evt_weight = (p_mc > 1e-10) ? p_data / p_mc : 1.0;
    
    return btag_evt_weight;
}


void SSBCorrections::LoadMCBtagEfficiencies(const std::string& filepath,
                                             const std::string& algo,
                                             const std::string& wp) {
    TFile* f = TFile::Open(filepath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Failed to open efficiency file: " << filepath << std::endl;
        return;
    }

    for (const std::string& flav : {"b", "c", "l"}) {
        std::string name = "eff_" + algo + "_" + flav + "_" + wp;
        TH2D* hist = (TH2D*)f->Get(name.c_str());
        if (hist) {
            TH2D* hist_copy = (TH2D*)hist->Clone((name + "_copy").c_str());
            hist_copy->SetDirectory(0);
            auto it = eff_histograms_.find(algo + "_" + flav + "_" + wp);
            if (it != eff_histograms_.end() && it->second) {
                delete it->second;
            }
            eff_histograms_[algo + "_" + flav + "_" + wp] = hist_copy;
            std::cout << "[INFO] Loaded and copied hist: " << name << std::endl;
        } else {
            std::cerr << "[WARNING] Histogram not found: " << name << std::endl;
        }
    }
    f->Close();
    delete f;
}

/*
void SSBCorrections::LoadMCBtagEfficiencies(const std::string& filepath, const std::string& algo) {
    TFile* f = TFile::Open(filepath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Failed to open efficiency file: " << filepath << std::endl;
        return;
    }

    for (const std::string& flav : {"b", "c", "l"}) {
        for (const std::string& wp : {"Loose", "Medium", "Tight"}) {
            std::string name = "eff_" + algo + "_" + flav + "_" + wp;
	    //std::cout << "name in LoadMCBtagEfficiencies : " << name << std::endl;
            TH2D* hist = (TH2D*)f->Get(name.c_str());
            if (hist) {
                hist->SetDirectory(0); // Detach histogram from file
                eff_histograms_[algo + "_" + flav + "_" + wp] = hist;
                //std::cout << "[INFO] Loaded hist: " << name << std::endl;
            } else {
                std::cerr << "[WARNING] Histogram not found: " << name << std::endl;
            }
        }
    }
    f->Close();
}
*/
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
    //std::cout << "[Info] key in GetMCBtagEfficiency " << key << std::endl;
    auto it = eff_histograms_.find(key);
    if (it == eff_histograms_.end()) {
        std::cerr << "[WARNING] Efficiency hist not found: " << key << std::endl;
        return 1.0;
    }

    TH2D* hist = it->second;
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
        std::cerr << "[SSBCorrections::METXYCorrection] MET correction object not loaded.\n";
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
            std::cerr << "[METXYCorrection] Warning: unexpected type for x correction\n";
        }

        if (std::holds_alternative<double>(val_y)) {
            corr_y = std::get<double>(val_y);
        } else {
            std::cerr << "[METXYCorrection] Warning: unexpected type for y correction\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "[SSBCorrections::METXYCorrection] Evaluation failed: " << e.what() << std::endl;
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
        std::cerr << "[ERROR] Invalid jec_name format: " << base_jec_name << std::endl;
        return base_jec_name; // fallback
    }

    std::string prefix = base_jec_name.substr(0, pos);             // "Summer19UL16"
    std::string version = base_jec_name.substr(pos + 2);           // "7"
    // Was hardcoded to "_L1L2L3Res_AK4PFchs" - NanoAODv15's Jet_ collection is
    // AK4 Puppi, not PFchs, so this now comes from the caller (config-driven:
    // JetAlgoTag/JECLevel), defaulting to the Puppi values. See the CMS CAT
    // jet_jerc.json.gz correction names (e.g. "..._L1FastJet_AK4PFPuppi",
    // "..._L2Relative_AK4PFPuppi") - if the compound correction there isn't
    // literally named "..._L1L2L3Res_AK4PFPuppi", set JECLevel in the config
    // to whatever the actual compound entry is called instead.
    std::string suffix = "_" + jecLevel + "_" + jetAlgo;

    std::string expanded = prefix;  // prefix

    // Confirmed from the user's jet_jerc.json.gz screenshot: the v15 Puppi
    // compound corrections are named e.g.
    // "Summer20UL18NanoV15_V1_DATA_L1L2L3Res_AK4PFPuppi" - a SINGLE name for
    // all of 2018 data, with no per-era ("_RunA"/"_RunB"/...) component. The
    // DATA compound's declared inputs also included "run" (a 5th evaluate()
    // input beyond area/eta/pt/rho - see jec_needs_run_), meaning era
    // resolution now happens *inside* the correction from the run number,
    // not via the correction *name* as in the old v9/PFchs scheme below. So:
    // skip the whole per-era name-building chain for Puppi jets; only do it
    // for PFchs (kept for v9 compatibility / in case some other campaign's
    // Puppi json still uses the old per-era-name convention).
    bool useEraInName = (jetAlgo.find("Puppi") == std::string::npos);

    // Year-specific logic (only for the old v9/PFchs per-era-name convention -
    // see useEraInName above)
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

    std::cerr << "[WARNING] GetProcessSubDir: no match for '" << inputfileName
              << "', falling back to TTbar_Signal" << std::endl;
    return "TTbar_Signal";
}
