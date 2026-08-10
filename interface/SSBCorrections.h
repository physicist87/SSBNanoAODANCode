#ifndef SSBCORRECTIONS_H
#define SSBCORRECTIONS_H

#include <memory>
#include <string>
#include <correction.h>
#include "TLorentzVector.h"
#include <variant> 
#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include "../CorrectionFiles/Rochester/RoccoR.h"
#include "../CorrectionFiles/METXY/XYMETCorrection_withUL17andUL18andUL16.h"
#include "DebugTools.h"

using correction::CorrectionSet;
using correction::CompoundCorrection;

// Forward declarations
class TextReader;

// Canonical b-tag algorithm enum, shared by SSBCorrections.cpp and
// Analysis.cpp so their Jet_btag string parsing can't drift apart again (it
// did once - see NOTES.md).
enum class BTagAlgo { DeepCSV, DeepJet, UParTAK4, CSVv2, Unknown };

// Parses the algorithm portion of a Jet_btag config string (e.g. "UParTM"),
// ignoring the trailing working-point letter.
BTagAlgo ParseBTagAlgo(const std::string& jetBtagConfig);

// String form of a BTagAlgo, for filenames/histogram names only - compare
// the enum directly with == elsewhere.
std::string BTagAlgoToString(BTagAlgo algo);

// Loads and applies correctionlib-based JEC, JER, and lepton scale factors.
class SSBCorrections {
public:
    explicit SSBCorrections(TextReader* reader, const std::string inputfileName);
    ~SSBCorrections();
    
    SSBCorrections(const SSBCorrections&) = delete;
    SSBCorrections& operator=(const SSBCorrections&) = delete;
    
    // jetAlgo/jecLevel default to NanoAODv15 Puppi but are config-driven so
    // v9/PFchs input or a different CAT naming scheme still work.
    std::string ExpandJECName(const std::string& base_jec_name, const std::string runPeriod, const std::string& era, bool is_data,
                               const std::string& jetAlgo = "AK4PFPuppi", const std::string& jecLevel = "L1L2L3Res");
    float GetPUWeight(float nTrueInt, const std::string& variation = "nominal") const;

    double GetJEC(double eta, double pt, double rho) const;

    // Unused elsewhere - kept 3-input (eta, pt, rho) to match jet_jerc.json.gz's
    // PtResolution schema.
    double GetJER(double eta, double pt, double rho) const;

    // Delegates to CMS JME's "JERSmear" tool when loaded, else falls back in-house.
    // Pass gen_pt=-1.0 for "no match" (eta/phi unused then).
    double SmearJER(double reco_pt, double gen_pt, double gen_eta, double gen_phi,
                     double eta, double phi, double rho,
                     ULong64_t event, const std::string& jer_tag = "nominal") const;

    // Reco/gen matching for JER hybrid smearing; out_eta/out_phi receive the
    // matched gen jet's eta/phi if found (check the return value, not these).
    float MatchGenPt(const TLorentzVector& reco_jet,
                     const std::vector<TLorentzVector>& gen_jets,
                     float maxDR = 0.2,
                     float* out_eta = nullptr,
                     float* out_phi = nullptr) const;

    // Builds the physics jet collection (pt+mass) via JES/JER only - MET
    // correction is a separate, independent computation (ApplyType1METWithCorrT1),
    // per the CMS JERC tutorial's own split. run_number is a no-op unless the
    // loaded JEC needs it (jec_needs_run_); jerSysTag is "nominal"/"up"/"down".
    std::vector<TLorentzVector> ApplyJetCorrections(
        const std::vector<TLorentzVector>& rawJets,
        const std::vector<float>& rawFactors,
        const std::vector<float>& areas,
        float rho,
        bool isData,
        bool applyJES,
        bool applyJER,
        const std::vector<TLorentzVector>& genJets,
        const std::vector<int>& genJetIndices,
        unsigned int run_number = 0,
        ULong64_t event_number = 0,
        const std::string& jerSysTag = "nominal"
    ) const;

    // Type-1 MET including CorrT1METJet_ and muon-subtraction. Returns MET
    // only, not corrected jets (use ApplyJetCorrections for those).
    TLorentzVector ApplyType1METWithCorrT1(
        double raw_met_pt,
        double raw_met_phi,
        const std::vector<TLorentzVector>& jetsAsStored,     // Jet_pt/eta/phi/mass, as stored (pre-undo of rawFactor)
        const std::vector<float>& jetRawFactors,             // Jet_rawFactor
        const std::vector<float>& jetAreas,                  // Jet_area
        const std::vector<float>& jetMuonSubtrFactors,       // Jet_muonSubtrFactor
        const std::vector<float>& jetChEmEF,                 // Jet_chEmEF (Type-1 selection cut)
        const std::vector<float>& jetNeEmEF,                 // Jet_neEmEF (Type-1 selection cut)
        const std::vector<int>& jetGenJetIndices,             // Jet_genJetIdx - same convention as ApplyJetCorrections
        const std::vector<float>& corrT1RawPt,               // CorrT1METJet_rawPt (already raw)
        const std::vector<float>& corrT1Eta,                 // CorrT1METJet_eta
        const std::vector<float>& corrT1Phi,                 // CorrT1METJet_phi
        const std::vector<float>& corrT1Area,                // CorrT1METJet_area
        const std::vector<float>& corrT1MuonSubtrFactor,      // CorrT1METJet_muonSubtrFactor
        float rho,
        bool isData,
        bool applyJES,
        bool applyJER,
        const std::vector<TLorentzVector>& genJets,
        unsigned int run_number = 0,
        ULong64_t event_number = 0,
        // Same tag as ApplyJetCorrections' jerSysTag - applies to both the
        // regular-Jet_ and CorrT1METJet_ loops, no per-collection mixing.
        const std::string& jerSysTag = "nominal"
    ) const;

    double GetMuonRecoSF(double pt, double eta) const;
    double GetMuonIDSF(double pt, double eta, const std::string& tag = "nominal") const;
    double GetMuonIsoSF(double pt, double eta, const std::string& tag = "nominal") const;
    double DoubleMuon_IDIsoEff(TLorentzVector lep1, TLorentzVector lep2, TString muidsys, TString muisosys, TString tracksys) const;
    float GetElectronSF(const std::string& sf_type, float eta, float pt, const std::string& syst = "sf") const;
    double DoubleElec_Eff(
        const TLorentzVector& lep1, const TLorentzVector& lep2,
        double ele1sueta, double ele2sueta,
        const std::string& id_wp = "Tight",       // Default: Tight working point
        const std::string& id_syst = "nominal",   // Default: nominal (no systematic)
        const std::string& reco_syst = "nominal"  // Default: nominal (no systematic)
    ) const;

    double MuonElec_Eff(const TLorentzVector& muon, const TLorentzVector& electron,
                        double muon_eta, double electron_sueta,
                        const std::string& mu_id_syst = "nominal",
                        const std::string& mu_iso_syst = "nominal",
                        const std::string& ele_id_wp = "Tight",
                        const std::string& ele_id_syst = "nominal",
                        const std::string& ele_reco_syst = "nominal") const;



    double TrigDiMuon_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_ = "nominal");
    double TrigDiElec_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_ = "nominal");
    double TrigMuElec_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_ = "nominal");

    // L1FastJet-only corrected pt (needed as the Type-1 MET baseline, used by
    // ApplyType1METWithCorrT1 - the sole MET computation now; RecomputeMET,
    // the old fused-into-ApplyJetCorrectionsWithMET MET calculation, was
    // removed since it was never the official recipe and its output had
    // stopped being used). Loaded from the same jet_jerc.json.gz as jec_, via
    // the same base-name/era/data expansion as the full compound correction
    // but with jecLevel="L1FastJet" - see the constructor.
    // Confirmed (both from the actual file and from the CMS JERC tutorial's
    // own code) that L1FastJet alone never needs a "run" input, unlike the
    // full L1L2L3Res compound for data - so no run_number parameter here.
    double GetL1CorrectedJetPt(double raw_pt, double eta, double area, double rho) const;

    // run_number is a no-op unless jec_needs_run_ is true (see cpp for detection).
    double GetCorrectedJetPt(double raw_pt, double eta, double area, double rho, unsigned int run_number = 0) const;
    double GetCorrectedJetMass(double raw_mass, double raw_pt, double eta, double area, double rho, unsigned int run_number = 0) const;
    TLorentzVector GetCorrectedJet(const TLorentzVector& raw_jet,
                                float rawFactor, float eta, float area, float rho,
                                int jet_index, int seed, bool isData,
                                bool applyJES, bool applyJER) const;

    TLorentzVector METXYCorrection(const TLorentzVector& type1_met,
                                               int runnb, TString year, bool isMC, int npv, bool isUL, bool ispuppi
                                               ) const;

    TLorentzVector METXYCorrection_corrlib(const TLorentzVector& type1_met,
                                const std::string& era,
                                bool isData,
                                int npv) const;

    float GetPUJetIDSFAndEff(float pt, float eta, bool passPU, bool genMatched, const std::string& wp, const std::string& syst, bool getEff = false) const;
    double RochesterCorrectionData(TString year, int Q, double pt, double eta, double phi, int s,int m) const;
    double RochesterCorrectionMC(TString year, int Q, double pt, double eta,double phi,int genID,double genPt,int nl, int s,int m) const;
    // EM-fraction veto cut only (kMaxEmFrac=0.90) - pt/jetId pre-selection
    // is already applied by the only caller, Analysis::JetSelector.
    bool ShouldVetoJet(const TLorentzVector& jet, double chEmEF = 0.0, double neEmEF = 0.0) const;
    std::string GetJetVetoType() const;
    void InitBtagSFCorrection(const std::string& json_path, const std::string& tagger_name);
    float GetBtagSF(float pt, float eta, int flav, const std::string& wp, const std::string& syst = "nominal") const;
    // Numeric discriminant cut for a WP from btagging.json.gz's wp_values
    // correction. Returns -1.0 (fails every cut) if unavailable - caller
    // should fall back to an explicit BTagDiscCut, not treat -1.0 as valid.
    double GetBtagWPCut(const std::string& wp) const;
    //void LoadMCBtagEfficiencies(const std::string& filepath, const std::string& algo);
    void LoadMCBtagEfficiencies(const std::string& filepath, const std::string& algo, const std::string& wp);
    float GetMCBtagEfficiency(float pt, float eta, int flav, const std::string& algo, const std::string& wp) const;
    float ComputeBTagEventWeight(const std::vector<float>& pts,
                             const std::vector<float>& etas,
                             const std::vector<int>& flavs,
                             const std::vector<bool>& isTagged,
                             const std::string& algo,
                             const std::string& wp,
                             const std::string& syst = "nominal") const;

private:
    // mutable: most logging call sites are inside const member functions,
    // and Logger's stream methods aren't const.
    mutable Logger logger_;

    std::string year_;
    std::string jveto_name_; // correction name (e.g., "Summer19UL18_V1") 
    std::string jveto_key_;  // veto map key (e.g., "jetvetomap", "hem1516")
    std::string jveto_type_;   // "jet" or "event"  
    std::string btag_sf_type_; // Added: "comb" or "mujets" 

    // Trigger variables
    double GetTrgEff(double pt1, double pt2, TString Sys_);
    // Owned; cloned out of the trigger SF file and detached via
    // SetDirectory(nullptr) at load time (see constructor), so it outlives
    // the TFile it was read from and is cleaned up automatically. TH2 (not
    // TH2D) because the actual file stores TH2F - GetBinContent/GetXaxis/
    // FindBin/GetBinError are all virtual on the common TH1/TH2 base, so
    // this works regardless of the concrete histogram type.
    std::unique_ptr<TH2> H_trig;
    
    // Rochester correction
    RoccoR rc;

    // Correction objects
    std::shared_ptr<const correction::Correction> pu_weight_;
    std::shared_ptr<const correction::Correction> muon_id_sf_;
    std::shared_ptr<const correction::Correction> muon_iso_sf_;
    std::shared_ptr<const correction::Correction> muon_reco_;
    std::shared_ptr<const correction::Correction> muon_id_;
    std::shared_ptr<const correction::Correction> muon_iso_;
    std::shared_ptr<const correction::Correction> ele_sf_;
    std::shared_ptr<const correction::Correction> ele_reco_sf_;
    std::shared_ptr<const correction::Correction> jetvetomap_;
    std::shared_ptr<const correction::CompoundCorrection> jec_;
    // Whether the loaded JEC compound takes a 5th "run" input (v15 Puppi
    // DATA does; detected once at construction - see cpp).
    bool jec_needs_run_ = false;
    // L1FastJet-only correction (Type-1 MET baseline); same base as jec_,
    // jecLevel="L1FastJet".
    std::shared_ptr<const correction::Correction> jec_l1_;
    std::shared_ptr<const correction::Correction> jer_;
    std::shared_ptr<const correction::Correction> jer_sf_; // JER Scale factor - evaluate({eta, pt}), confirmed 2-input (not {eta, syst})
    std::shared_ptr<const correction::Correction> jer_sfunc_; // JER SF uncertainty - evaluate({eta, pt}); combined arithmetically as sf*(1+-unc) for up/down, not a separate correction per variation
    // CMS's official JER smearing tool; optional - SmearJER() falls back
    // in-house (with a warning) if unavailable.
    std::shared_ptr<const correction::Correction> jer_smear_;
    // JES full-uncertainty-set systematic, resolved from JESSys/JESSysDir at
    // construction. nullptr = no JES systematic; applied inside
    // GetCorrectedJetPt() so both the jet and MET paths pick it up.
    std::shared_ptr<const correction::Correction> jes_unc_source_;
    std::string jes_sys_dir_; // "up" or "down"; only meaningful if jes_unc_source_ is loaded
    std::shared_ptr<const correction::Correction> pujetid_sf_; // PU JetID SF

    // B-tagging corrections map
    std::map<std::string, std::shared_ptr<const correction::Correction>> btag_corrections_;
    // "<tagger>_wp_values" correction (currently only loaded for UParTAK4 -
    // see InitBtagSFCorrection) - optional, nullptr if not present/loadable.
    std::shared_ptr<const correction::Correction> btag_wp_values_;

    // Helper function to get appropriate correction name
    std::string getBtagCorrectionName(int flavor) const;
    std::string GetProcessSubDir(const std::string& inputfileName) const;
    
    // Owned MC b-tag efficiency histograms; unique_ptr means no manual
    // delete loop in the destructor and no dangling/double-free risk.
    // TH2 (not TH2D) for the same reason as H_trig above.
    std::map<std::string, std::unique_ptr<TH2>> eff_histograms_;
    
    std::shared_ptr<const correction::Correction> metphi_corr_;
};

#endif  // SSBCORRECTIONS_H
