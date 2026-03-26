#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>

#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector2.h>
#include <TH1D.h>
#include <TFile.h>

// Top reconstruction //TOP-17-014
#include "./../KinSolv/analysisUtils.h"                                                                                       
#include "./../KinSolv/KinematicReconstruction.h"                                                                             
#include "./../KinSolv/KinematicReconstructionSolution.h"

// Common tools 
#include "./../CommonTools.hpp"

// Correction tools
#include "./SSBCorrections.h"

// Textreader
#include "./../TextReader/TextReader.hpp"

// CPVObservables Calculator 
#include "./../interface/SSBCPVCalc.h"


class Analysis {
public:
    // Constructor and destructor
Analysis(TChain *chain, std::string inputName, std::string seDirName, std::string outputName, const std::string &branchListFile, const std::string &configFile, int NumEvt);

    //~Analysis() = default;
    ~Analysis();

    // Set variable function
    void SetVariables();
    // Event loop function
    void Loop();
    // create TLorentzVector //
    TLorentzVector createLorentzVector(float pt, float eta, float phi, float mass);
    void Start();
    void DeclareHistos();
private:
    // TTreeReader and TChain
    TChain *chain;
    std::string outfile;
    std::string outdir;
    TTreeReader fReader;
    TFile *fout;

    Long64_t current_entry_;
    bool isjetveto_event_;

    bool isData;
    //TextReader from Jaehoon.
    TextReader *SSBConfReader;
    SSBCorrections *SSBCorr;
    SSBCPVCalc *SSBCPVCal;
    int NumEvt; //

    double Lumi;
    TString FileName_;
    TString RunPeriod;
    TString Decaymode;
    TString XsecTable_;
    TString METtype;
    TString applyMETXY;
    TString applyRochester;


    TString cutflowName[10];

    int num_pv;


    // PUID related state variables
    bool jets_selected_ = false;
    bool jet_puid_weight_applied_ = false;
    bool object_variables_set_ = false;

    // PUID related configuration
    bool apply_puid_ = true;     // whether to apply PUID


    // PUID related variables (from Step 1)
    std::string puid_wp_;
    std::string PUIDSFSys;
    float puid_pt_threshold_;
    double evt_weight_beforePUID_;
    double puid_sf_weight_;
    
    // ============================================================================
    // Step 2: NEW - PUID candidate jets information structure
    // ============================================================================
    struct PUIDJetInfo {
        int original_index;     // Original jet index in jets collection
        float pt, eta;         // Jet kinematics
        bool passes_puid;      // Whether jet passes PUID
        bool is_hardscatter;   // Whether jet is matched to gen jet (HardScatter)
        
        // Constructor for easy initialization
        PUIDJetInfo(int idx, float p, float e, bool pass, bool hard) 
            : original_index(idx), pt(p), eta(e), passes_puid(pass), is_hardscatter(hard) {}
    };
    
    std::vector<PUIDJetInfo> puid_hardscatter_jets_;  // HardScatter jets for weight calculation
    
    // ============================================================================
    // Step 2: NEW - Function declarations
    // ============================================================================
    void CollectPUIDCandidates();
    bool IsHardScatterJet(int jet_idx) const;



    // Maps for dynamic branch storage
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> boolSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Int_t>>> intSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t>>> uintSingles; // Add uintSingles for UInt_t types
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UChar_t>>> ucharSingles; // Add uintSingles for UInt_t types
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Float_t>>> floatSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>> floatVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Bool_t>>> boolVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Int_t>>> intVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UInt_t>>> uintVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UChar_t>>> ucharVectors;

    // Initialize TTreeReader with branches from the list
    void InitBranches(const std::string &branchListFile);
    std::string removeSubstring(std::string &str, const std::string &keyword);
    //std::unique_ptr<TTreeReaderValue<bool>> DeepCopy(const std::unique_ptr<TTreeReaderValue<bool>>& src);
    template <typename T>std::unique_ptr<TTreeReaderValue<T>> DeepCopy(const std::unique_ptr<TTreeReaderValue<T>>& src);
    bool METFilterAPP();
    bool Trigger();
    bool SelTrigger(std::vector<std::string> v_sel);
    //bool ;
    //TString SetInputFileName( char *inname );
    TString SetInputFileName( std::string inname );
    void SetObjectVariable();
    void MCSF();
    void MCSFApply();
    void GenWeightApply();
    void PUWeightApply();
    void L1PreFireApply();
    void NumPVCount();
    void LeptonSelector(); //(muon & electron)
    void LeptonOrder(); // Lepton (muon & electron)
    
    void SelectVetoMuons();
    void SelectVetoElectrons();

    //void CorrectedMuonCollection(); // Muon  
    void MakeMuonCollection(); // Muon  
    void MakeElecCollection(); // Electron 
    void JetSelector(); // Jet

    // PUID related functions
    bool PassPileupID(float pt, int puId, const std::string& wp = "L") const;
    void ApplyJetPUIDEventWeights();

    void MakeJetCollection(); // Jet 
    //void MakeGenJetCollection(); // Jet 
    bool JetCleaning(TLorentzVector* jet_);
    TLorentzVector JERSmearing(TLorentzVector* jet, int idx_, TString op_);
    void bJetSelector();
    void METDefiner();
    void JetOrder();

    //
    bool NumIsoLeptons(int nNLepsCut);
    bool ThirdLeptonVeto(); // Lepton (muon & electron)
    bool MuonVeto(); // Lepton (muon & electron)
    bool ElectronVeto(); // Lepton (muon & electron)
    bool LeptonsPtAddtional();
    bool DiLeptonMassCut();
    bool ZVetoCut();
    bool NumJetCut(std::vector<int> v_jets);
    bool METCut(TLorentzVector met);
    bool NumbJetCut(std::vector<int> v_jets);

    // Event Weight //
    void LeptonSFApply();
    void TriggerSFApply();

    ////////////////////////////
    /// New Kinematic Solver ///
    ////////////////////////////
    void SetUpKINObs();
    bool isKinSol;
    VLV v_leptons_VLV; 
    VLV v_jets_VLV; 
    VLV v_bjets_VLV; 
    std::vector<int> v_lepidx_KIN; 
    std::vector<int> v_anlepidx_KIN; 
    std::vector<int> v_jetidx_KIN; 
    std::vector<int> v_bjetidx_KIN; 
    std::vector<double> v_btagging_KIN; 


    ////////////////////////////
    /// Trigger & MET Filter ///
    ////////////////////////////
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> triggerList;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> noiseFilters;

    TLorentzVector Lep1; // 
    TLorentzVector Lep2; //
    TLorentzVector Lep; // 
    TLorentzVector AnLep; //
    TLorentzVector Muon; // only for muon-electron channel
    TLorentzVector Elec; // only for muon-electron channel
    TLorentzVector Met; // 
    TLorentzVector Jet1; // 
    TLorentzVector Jet2; // 

    // Generator //
    TTreeReaderArray<Float_t>* Gen_Weight;
    // muons //
    TTreeReaderArray<Float_t>* muons_pt;
    TTreeReaderArray<Float_t>* muons_eta;
    TTreeReaderArray<Float_t>* muons_phi;
    TTreeReaderArray<Float_t>* muons_M;
    TTreeReaderArray<float>* muons_iso;
    TTreeReaderArray<float>* muonsveto_iso;
    TTreeReaderArray<bool>   *muons_Id;
    TTreeReaderArray<bool>   *muonsveto_Id;
    std::vector<TLorentzVector> pre_muons;
    std::vector<TLorentzVector> muons;
    std::vector<TLorentzVector> muonsveto;
    TLorentzVector TMuon;
    TLorentzVector TElectron;

    // electrons //
    TTreeReaderArray<Float_t>* elecs_pt;
    TTreeReaderArray<Float_t>* elecs_eta;
    TTreeReaderArray<Float_t>* elecs_phi;
    TTreeReaderArray<Float_t>* elecs_M;
    TTreeReaderArray<float>*   elecs_iso;
    TTreeReaderArray<float>*   elecsveto_iso;
    TTreeReaderArray<Int_t>*   elecs_scbId;
    TTreeReaderArray<Bool_t>*  elecs_mvaId;
    TTreeReaderArray<Int_t>*   elecsveto_scbId;
    TTreeReaderArray<Bool_t>*  elecsveto_mvaId;
    std::vector<TLorentzVector> pre_elecs;
    std::vector<TLorentzVector> elecs;
    std::vector<TLorentzVector> elecsveto;

    // jets //
    TTreeReaderArray<Float_t>* jets_pt;
    TTreeReaderArray<Float_t>* jets_eta;
    TTreeReaderArray<Float_t>* jets_phi;
    TTreeReaderArray<Float_t>* jets_M;
    //TTreeReaderArray<Bool_t>*  jets_Id;
    TTreeReaderArray<Int_t>*  jets_Id;
    TTreeReaderArray<Int_t>*  jets_puId;
    TTreeReaderArray<Float_t>* jets_btag;
    std::vector<TLorentzVector> pre_jets;
    std::vector<TLorentzVector> jets;
    // gen jets // 
    TTreeReaderArray<Float_t>* gen_jets_pt;
    TTreeReaderArray<Float_t>* gen_jets_eta;
    TTreeReaderArray<Float_t>* gen_jets_phi;
    TTreeReaderArray<Float_t>* gen_jets_M;
    //std::vector<TLorentzVector> genJets;
    bool dojer;

    // MET //
    TTreeReaderValue<Float_t>* met_pt;
    TTreeReaderValue<Float_t>* met_phi;

    // Top Recontruction //
    TLorentzVector Top;
    TLorentzVector AnTop;
    TLorentzVector Top1;
    TLorentzVector Top2;
    TLorentzVector bJet;
    TLorentzVector AnbJet;
    TLorentzVector Nu;
    TLorentzVector AnNu;
    TLorentzVector W1;
    TLorentzVector W2;


    //LepIdType   v_muon_Id;
    //LepIdType   v_elec_Id;

    // index vector for object (muon, electron, jet)
    std::vector<int> v_lepton_idx; // Indecies of lepton
    std::vector<int> v_muon_idx; // Indecies of muon
    std::vector<int> v_vetomuon_idx; // Indecies of veto muon
    std::vector<int> v_electron_idx; // Indecies of electron
    std::vector<int> v_vetoelectron_idx; // Indecies of veto muon
    std::vector<int> v_jet_idx;  // Indecies of jet
    std::vector<int> v_bjet_idx;  // Indecies of jet

    //JetCleanning information from config.
    TString veto_muoniso_type;
    TString veto_muonid;
    TString veto_eleciso_type;
    TString veto_elecid;

    double pi;

    //to get trigger information from config.
    int num_dleptrig;
    int num_sleptrig;
    std::vector<std::string> DLtrigName; // trigger 
    std::vector<std::string> SLtrigName; // trigger 
    std::vector<std::string> trigName; // trigger 

    // lepton variables

    //Muon Type for MuEl - channel.
    TString MuonIsoType;
    TString MuonId;
    TString ElecId;
    TString ElecIsoType;
    TString JetId;
    TString JetbTag;
    TString JetEnSys;
    TString JetResSys;
    TString BTagSFSys;
    TString BTagEffSys;
    // MET Systematic ...
    TString MetSys;
    
    // Lepton Systematic
    TString LepIdSFSys;
    TString LepIsoSFSys;
    TString LepRecoSFSys;
    TString LepTrackSFSys;
    
    // PileUp Systematic ...
    std::string PileUpMCFile;
    std::string PileUpDATAFile;
    TString PileUpSys;
    TString L1PreFireSys;
    TString TrigSFSys;

    // PDF Systematic ... Getting the Num of PDF Set
    int PDFSys;
    // FactReno Systematic ... Getting the Num of FactReno Set
    int FactRenoSys;
    // Fragment Systematic ... Getting the Num of Fragment Set
    TString FragmentSys;
    // DecayTable Systematic ... Getting the Num of DecayTable Set
    TString DecayTableSys;
    // TopPtReweight Sys
    TString TopPtSys;


    int eleid_scbcut;
    int elevetoid_scbcut;

    /// Event weights ...
    double mc_sf_;
    double evt_weight_; 
    double evt_weight_beforemcsf_;
    double evt_weight_beforegenweight_;
    double evt_weight_beforefactrenoweight_;
    double evt_weight_beforedectabweight_;
    double evt_weight_beforefragmentweight_;
    double evt_weight_beforepdfweight_;
    double evt_weight_beforePileup_;
    double evt_weight_beforeL1PreFire_;
    double evt_weight_beforeTrigger_;
    double evt_weight_beforeLepsf_;
    // B-tagging SF related variables
    double evt_weight_beforeBtag_;
    double btag_sf_weight_;
   
    void PUIDSFApply(); 
    // B-tagging SF application function
    void BTaggingSFApply();

    double lep_sf;

    //Kinetic Variables                                                                                                               
    int    n_elep_Id;                                                                                                                 
    double el_id_1;                                                                                                                   
    double el_id_2;                                                                                                                   
    
    double muon_pt;
    double muon_eta;
    double muon_isocut;                                                                                                                
    double veto_muoniso_cut;                                                                                                                
    double elec_pt;
    double elec_eta;
    double elec_isocut;
    int    n_elec_Id;


    double jet_pt;
    double jet_eta;
    int    jet_id;
    float  bdisccut;
    double met_cut;

    //Jetcleannig Variables
    double mu_pt_jetcl;
    double mu_eta_jetcl;
    double mu_iso_jetcl;

    double el_pt_jetcl;
    double el_eta_jetcl;
    double el_iso_jetcl;
    int    n_el_id_jetcl;
    double el_id_jetcl_1;
    double el_id_jetcl_2;

    /// BTag -- Variables ///
    std::string btag_algo_;     // "DeepCSV", "DeepJet", "CSVv2"
    std::string btag_wp_; // "L", "M", "T"


    TH1D *h_JetPUIDEvtWeight;
    TH1D *h_bTagEvtWeight;
    TH1D *h_Lep1pt[10];
    TH1D *h_Lep2pt[10];
    TH1D *h_Lep1eta[10];
    TH1D *h_Lep2eta[10];
    TH1D *h_Lep1phi[10];
    TH1D *h_Lep2phi[10];
    TH1D *h_Jet1pt[10];
    TH1D *h_Jet2pt[10];
    TH1D *h_Jet1eta[10];
    TH1D *h_Jet2eta[10];
    TH1D *h_Jet1phi[10];
    TH1D *h_Jet2phi[10];
    TH1D *h_HT[10]; 
    TH1D *h_METpt[10]; 
    TH1D *h_METphi[10]; 
    TH1D *h_DiLepMass[10];
    TH1D *h_cf_DiLepMass[10];

    TH1D *h_NoMETCut_Lep1pt[10];
    TH1D *h_NoMETCut_Lep2pt[10];
    TH1D *h_NoMETCut_Lep1eta[10];
    TH1D *h_NoMETCut_Lep2eta[10];
    TH1D *h_NoMETCut_Lep1phi[10];
    TH1D *h_NoMETCut_Lep2phi[10];
    TH1D *h_NoMETCut_Jet1pt[10];
    TH1D *h_NoMETCut_Jet2pt[10];
    TH1D *h_NoMETCut_Jet1eta[10];
    TH1D *h_NoMETCut_Jet2eta[10];
    TH1D *h_NoMETCut_Jet1phi[10];
    TH1D *h_NoMETCut_Jet2phi[10];
    TH1D *h_NoMETCut_HT[10];
    TH1D *h_NoMETCut_METpt[10];
    TH1D *h_NoMETCut_METphi[10];
    TH1D *h_NoMETCut_DiLepMass[10];
    TH1D *h_cf_NoMETCut_DiLepMass[10];
    TH1D *h_NoMETCut_Num_PV[10];
    TH1D *h_NoMETCut_Num_Jets[10];
    TH1D *h_NoMETCut_Num_bJets[10];



    TH1D *h_Num_PV[10];
    TH1D *h_Num_Jets[10];
    TH1D *h_Num_bJets[10];

    TH1D *h_Reco_CPO_[13];
    TH1D *h_Reco_CPO_ReRange_[13];

    TH1D* h_Top1Mass;
    TH1D* h_Top1pt;
    TH1D* h_Top1Rapidity;
    TH1D* h_Top1phi;
    TH1D* h_Top1Energy;
    TH1D* h_Top2Mass;
    TH1D* h_Top2pt;
    TH1D* h_Top2Rapidity;
    TH1D* h_Top2phi;
    TH1D* h_Top2Energy;
    
    TH1D* h_TopMass;
    TH1D* h_Toppt;
    TH1D* h_TopRapidity;
    TH1D* h_Topphi;
    TH1D* h_TopEnergy;
    TH1D* h_AnTopMass;
    TH1D* h_AnToppt;
    TH1D* h_AnTopRapidity;
    TH1D* h_AnTopphi;
    TH1D* h_AnTopEnergy;

    TH1D* h_W1Mass;
    TH1D* h_W2Mass;
    TH1D* h_W1Mt;
    TH1D* h_W2Mt;

    TH1D* h_bJet1Energy;
    TH1D* h_bJet2Energy;

    TH1D* h_bJetEnergy;
    TH1D* h_AnbJetEnergy;
    TH1D* h_bJetPt;
    TH1D* h_AnbJetPt;

    TH1D* h_Lep1Energy;
    TH1D* h_Lep2Energy;
    
    TH1D* h_LepEnergy;
    TH1D* h_AnLepEnergy;

    TH1D* h_Nu1Energy;
    TH1D* h_Nu2Energy;

    TH1D* h_NuEnergy;
    TH1D* h_AnNuEnergy;

};

#endif // ANALYSIS_H
