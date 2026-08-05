#ifndef SSBCORRECTIONS_H
#define SSBCORRECTIONS_H

#include <memory>
#include <string>
#include <correction.h>
#include "TLorentzVector.h"
#include <variant> 
#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include "../CorrectionFiles/Rochester/RoccoR.h"
#include "../CorrectionFiles/METXY/XYMETCorrection_withUL17andUL18andUL16.h"

using correction::CorrectionSet;
using correction::CompoundCorrection;

// Forward declarations
class TextReader;

// A utility class for loading and applying correctionlib-based
// JEC, JER, and muon scale factors using configuration
class SSBCorrections {
public:
    // Constructor using a configuration reader (TextReader)
    explicit SSBCorrections(TextReader* reader, const std::string inputfileName);
    
    // Explicit destructor declaration to fix memory management issues
    ~SSBCorrections();
    
    SSBCorrections(const SSBCorrections&) = delete;
    SSBCorrections& operator=(const SSBCorrections&) = delete;
    
    // jetAlgo/jecLevel default to the NanoAODv15 Puppi-jet values but are
    // config-driven (JetAlgoTag/JECLevel) so this keeps working if CMS JME's
    // exact correction-name scheme for the new CAT json turns out to differ,
    // or if this code is ever pointed at PFchs (v9) input again.
    std::string ExpandJECName(const std::string& base_jec_name, const std::string runPeriod, const std::string& era, bool is_data,
                               const std::string& jetAlgo = "AK4PFPuppi", const std::string& jecLevel = "L1L2L3Res");
    // Get PU weight for a given nTrueInt and variation
    float GetPUWeight(float nTrueInt, const std::string& variation = "nominal") const;

    // Jet Energy Correction (JEC) factor
    double GetJEC(double eta, double pt, double rho) const;

    // Jet Energy Resolution (JER) sigma. Currently unused elsewhere in this
    // codebase - fixed to the correct 3-input (eta, pt, rho) evaluate() call
    // to match jet_jerc.json.gz's PtResolution schema (it previously omitted
    // rho, unlike SmearJER's internal call, which was already correct).
    double GetJER(double eta, double pt, double rho) const;

    // Smear JER for a MC jet. Delegates the actual hybrid-method/stochastic
    // decision and random smearing to CMS JME's own "JERSmear" correctionlib
    // tool (jer_smear.json.gz, a correctionlib `hashprng` node) instead of a
    // hand-rolled implementation - this is the officially recommended
    // approach (see the CMS JERC ApplicationTutorial's JecApplication.cpp,
    // Applier::jerFactor): it takes (JetPt, JetEta, GenPt-or-(-1), Rho,
    // EventID, JER resolution, JER SF) and returns a multiplicative smearing
    // factor directly - reproducibility/seeding is handled internally by
    // correctionlib, not by us. gen_pt < 0 is the "no gen match" sentinel
    // (previously, no gen match meant JER was skipped entirely for that jet -
    // now it's passed through and the JERSmear tool decides scaling vs
    // stochastic itself, including its own 3-sigma consistency check).
    double SmearJER(double reco_pt, double gen_pt, double eta, double phi, double rho,
                     ULong64_t event, const std::string& jer_tag = "nominal") const;

    // Reco/gen matching for JER hybrid smearing
    float MatchGenPt(const TLorentzVector& reco_jet,
                     const std::vector<TLorentzVector>& gen_jets,
                     float maxDR = 0.2) const;

    // Apply JES/JER corrections to build the physics jet collection (pt +
    // mass). Confirmed against the CMS JERC ApplicationTutorial's
    // JecApplication.h: the tutorial's Applier class keeps jet-pt/mass
    // correction (jesFactorNominal()/jerFactor(), which the caller multiplies
    // into Jet_pt/mass itself) and Type-1 MET correction (correctedMet(), a
    // fully independent computation from raw quantities) as two SEPARATE
    // computations - they are not fused into one function/return value.
    // This function used to also return a "corrected_met" (via an internal
    // RecomputeMET() call using a non-official raw-vs-corrected delta, plus
    // -MET propagation missing muon subtraction and the Type-1 jet selection
    // cut) - that MET value is no longer computed here at all now that MET
    // always comes from ApplyType1METWithCorrT1 (the function that actually
    // matches the tutorial's correctedMet() term-for-term). Renamed from
    // ApplyJetCorrectionsWithMET to reflect that this only returns jets now.
    // run_number is only actually used when the loaded JEC compound correction
    // turns out to need it (see jec_needs_run_) - pass the event's "run" branch
    // value here regardless; it's a no-op for corrections that don't need it.
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
        ULong64_t event_number = 0
    ) const;

    // Type-1 MET recomputation that additionally includes NanoAOD's
    // CorrT1METJet_ branch (low-pT jets below the Jet_ collection's storage
    // threshold that NanoAOD stores separately, but which still need to
    // enter Type-1 MET propagation for it to be complete/correct - see the
    // CMS JERC tutorial's CollectJetMet::collectInputsForType1Met) and
    // Jet_/CorrT1METJet_muonSubtrFactor (avoids double-counting a muon's
    // momentum in MET, since muons are already handled directly in MET
    // rather than via jet energy).
    //
    // This does NOT return corrected physics jets (use
    // ApplyJetCorrections for that, unaffected by muon subtraction -
    // muon-subtracting the physics jet collection used for jet
    // selection/counting would be wrong) - it only returns the recomputed
    // MET, built as raw_met + sum_over_selected_T1_jets(L1only_no_mu - corrected_no_mu),
    // confirmed against the CMS JERC ApplicationTutorial's
    // JecApplication::Applier::correctedMet(): the delta is (L1-only minus
    // fully-corrected), not (raw minus fully-corrected) - the L1 (pileup-only)
    // correction is treated as already implicitly reflected in the raw MET.
    // Only jets passing the Type-1 selection (corrected pt > 15, |eta| < 5.2,
    // chEmEF+neEmEF < 0.90) contribute - CorrT1METJet_ jets have no EM
    // fraction branches (definitionally 0, always pass that part).
    TLorentzVector ApplyType1METWithCorrT1(
        double raw_met_pt,
        double raw_met_phi,
        const std::vector<TLorentzVector>& jetsAsStored,     // Jet_pt/eta/phi/mass, as stored (pre-undo of rawFactor)
        const std::vector<float>& jetRawFactors,             // Jet_rawFactor
        const std::vector<float>& jetAreas,                  // Jet_area
        const std::vector<float>& jetMuonSubtrFactors,       // Jet_muonSubtrFactor
        const std::vector<float>& jetChEmEF,                 // Jet_chEmEF (Type-1 selection cut)
        const std::vector<float>& jetNeEmEF,                 // Jet_neEmEF (Type-1 selection cut)
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
        ULong64_t event_number = 0
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
    // chEmEF/neEmEF added per the CMS JERC tutorial's JvmApplication::VetoChecker
    // (kMaxEmFrac = 0.90) - jets with (chEmEF+neEmEF) >= 0.90 should not be
    // evaluated against the veto map at all. pt/jetId pre-selection (the
    // tutorial's kMinPt=15, kMinJetId=6/TightLepVeto) are NOT re-checked here
    // since the only call site (Analysis::JetSelector) already applies a
    // tighter pt cut and PassConfiguredJetId() before calling this.
    bool ShouldVetoJet(const TLorentzVector& jet, double chEmEF = 0.0, double neEmEF = 0.0) const;
    std::string GetJetVetoType() const;
    void InitBtagSFCorrection(const std::string& json_path, const std::string& tagger_name);
    float GetBtagSF(float pt, float eta, int flav, const std::string& wp, const std::string& syst = "nominal") const;
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
    std::string year_;
    std::string jveto_name_; // correction name (e.g., "Summer19UL18_V1") 
    std::string jveto_key_;  // veto map key (e.g., "jetvetomap", "hem1516")
    std::string jveto_type_;   // "jet" or "event"  
    std::string btag_sf_type_; // Added: "comb" or "mujets" 

    // Trigger variables
    double GetTrgEff(double pt1, double pt2, TString Sys_);
    TH2D* H_trig;    
    
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
    // Detected once at construction (see cpp): the v15 Puppi-jet DATA compound
    // JEC correction (per the CAT jet_jerc.json.gz) appears to take a 5th
    // "run" input in addition to (area, eta, pt, rho) - unlike the old v9/
    // PFchs scheme, which encoded the era directly in the correction name
    // instead. jec_needs_run_ lets GetCorrectedJetPt/Mass build the right
    // number of evaluate() arguments without hardcoding which schema is in use.
    bool jec_needs_run_ = false;
    // L1FastJet-only correction (single, not compound) - needed as the Type-1
    // MET baseline (see RecomputeMET/GetL1CorrectedJetPt). Loaded from the
    // same jet_jerc.json.gz/jec_name base as jec_, just with jecLevel
    // ="L1FastJet" instead of "L1L2L3Res".
    std::shared_ptr<const correction::Correction> jec_l1_;
    std::shared_ptr<const correction::Correction> jer_;
    std::shared_ptr<const correction::Correction> jer_sf_; // JER Scale factor - evaluate({eta, pt}), confirmed 2-input (not {eta, syst})
    std::shared_ptr<const correction::Correction> jer_sfunc_; // JER SF uncertainty - evaluate({eta, pt}); combined arithmetically as sf*(1+-unc) for up/down, not a separate correction per variation
    // CMS JME's official "JERSmear" correctionlib tool (jer_smear.json.gz) -
    // does the actual hybrid-method/stochastic JER smearing internally via a
    // `hashprng` node, given (pt, eta, genPt-or-(-1), rho, eventID,
    // resolution, sf). Optional/nullable: if this file isn't available for a
    // given campaign, SmearJER() falls back to an in-house implementation
    // (see cpp) - loudly warned once, since that fallback is not the
    // officially recommended approach.
    std::shared_ptr<const correction::Correction> jer_smear_;
    std::shared_ptr<const correction::Correction> pujetid_sf_; // PU JetID SF

    // B-tagging corrections map
    std::map<std::string, std::shared_ptr<const correction::Correction>> btag_corrections_;

    // Helper function to get appropriate correction name
    std::string getBtagCorrectionName(int flavor) const;
    std::string GetProcessSubDir(const std::string& inputfileName) const;
    
    // Primary cause of segmentation fault - TH2D pointers need manual cleanup
    std::map<std::string, TH2D*> eff_histograms_;
    
    std::shared_ptr<const correction::Correction> metphi_corr_;
};

#endif  // SSBCORRECTIONS_H
