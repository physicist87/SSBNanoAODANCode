#ifndef KINRECO_INPUT_MAKER_H
#define KINRECO_INPUT_MAKER_H

// Produces the input histograms read by KinSolv/KinematicReconstruction
// (loadData()) from ttbar dilepton MC:
//
//   KinReco_W_mass_step0       gen W mass (both W's, last copy, from top)
//   KinReco_mbl_true_step0     gen m(l+ b), m(l- bbar): final-state (status 1)
//                              prompt lepton from W, b quark = direct top
//                              daughter (parton, no GenJet)
//   KinReco_fE_lep_step7       E_gen / E_reco for reco leptons matched to the
//                              gen final-state leptons (same charge, dR<0.1)
//   KinReco_d_angle_lep_step7  3D angle(reco, gen) for the same leptons
//   KinReco_fE_jet_step7       E_bquark / E_jet for selected jets matched to
//                              the top b quarks (dR<0.4, unambiguous)
//   KinReco_d_angle_jet_step7  3D angle(jet, b quark) for the same jets
//
// step0 = every gen dilepton event that reaches the job (after the gen
// channel split), step7 = after the full dilepton selection incl. >=1 b-tag
// (same objects that are handed to the solver). Binning is identical to the
// old KinSolv/Kin_input.root so the solver code needs no changes.
//
// Filled with the sign of Generator_weight only (unit weights), so the bin
// contents stay raw-count-like - the solver uses h_mbl bin contents directly
// as a weight (/1e8), and resolution shapes do not need SFs.

#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>

class NanoAODBranchReader;

class KinRecoInputMaker {
public:
    struct GenTtbar {
        bool hasLep = false, hasAnLep = false, hasB = false, hasBbar = false;
        bool hasWplus = false, hasWminus = false;
        TLorentzVector lep, anLep;   // l- (pdgId>0), l+ (pdgId<0); final state
        TLorentzVector b, bbar;      // direct top / antitop daughters
        TLorentzVector wPlus, wMinus;
        bool Complete() const { return hasLep && hasAnLep && hasB && hasBbar; }
    };

    // Books the histograms in the current gDirectory (i.e. the job's fout).
    KinRecoInputMaker();

    // Reads the NanoAOD GenPart collection. Needs GenPart_{pt,eta,phi,mass,
    // pdgId,status,statusFlags,genPartIdxMother} in the branch list.
    static GenTtbar FindGenTtbar(const NanoAODBranchReader& reader);

    void FillGen(const GenTtbar& gen, double w);
    void FillReco(const GenTtbar& gen,
                  const TLorentzVector& recoLep, const TLorentzVector& recoAnLep,
                  const std::vector<TLorentzVector>& jets, double w);

    void PrintSummary() const;

    static constexpr double kLepMatchDR = 0.1;
    static constexpr double kJetMatchDR = 0.4;

private:
    TH1D* h_wmass_;
    TH1D* h_mbl_;
    TH1D* h_fE_lep_;
    TH1D* h_angle_lep_;
    TH1D* h_fE_jet_;
    TH1D* h_angle_jet_;
    // diagnostics (not read by the solver)
    TH2D* h_fE_jet_vs_pt_;
    TH2D* h_fE_lep_vs_pt_;
    TH1D* h_dR_lep_;
    TH1D* h_dR_jet_;
    TH1D* h_count_;   // bookkeeping, see bin labels

    long nGenAll_ = 0, nGenComplete_ = 0, nReco_ = 0, nRecoLepMatch_ = 0, nRecoJetMatch_ = 0;
};

#endif
