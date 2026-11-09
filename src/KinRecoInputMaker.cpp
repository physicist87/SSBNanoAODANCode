#include "../interface/KinRecoInputMaker.h"
#include "../interface/NanoAODBranchReader.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <initializer_list>

namespace {
// NanoAOD GenPart_statusFlags bits
constexpr int kIsPrompt          = 0;
constexpr int kIsTauDecayProduct = 2;
constexpr int kFromHardProcess   = 8;
constexpr int kIsLastCopy        = 13;

inline bool flag(Long64_t flags, int bit) { return (flags >> bit) & 1; }
}

KinRecoInputMaker::KinRecoInputMaker()
{
    // --- read by KinematicReconstruction::loadData() (same binning as the old Kin_input.root)
    h_wmass_     = new TH1D("KinReco_W_mass_step0",      "gen W mass;m_{W} [GeV]",                250, 30.0, 130.0);
    h_mbl_       = new TH1D("KinReco_mbl_true_step0",    "mbl_true;m_{lb}^{gen} [GeV]",           100,  0.0, 180.0);
    h_angle_jet_ = new TH1D("KinReco_d_angle_jet_step7", "angle(jet, b quark);#alpha [rad]",      200,  0.0,   0.5);
    h_angle_lep_ = new TH1D("KinReco_d_angle_lep_step7", "angle(reco lep, gen lep);#alpha [rad]", 200,  0.0,   0.2);
    h_fE_jet_    = new TH1D("KinReco_fE_jet_step7",      "Jet Energy Correction Factor;E_{b}/E_{jet}",        100, 0.0, 4.0);
    h_fE_lep_    = new TH1D("KinReco_fE_lep_step7",      "Lepton Energy Correction Factor;E_{gen}/E_{reco}", 200, 0.5, 2.5);

    // --- diagnostics
    h_fE_jet_vs_pt_ = new TH2D("KinRecoDiag_fE_jet_vs_pt_step7", ";jet p_{T} [GeV];E_{b}/E_{jet}",      40, 0, 400, 100, 0.0, 4.0);
    h_fE_lep_vs_pt_ = new TH2D("KinRecoDiag_fE_lep_vs_pt_step7", ";lep p_{T} [GeV];E_{gen}/E_{reco}",   40, 0, 400, 200, 0.5, 2.5);
    h_dR_lep_       = new TH1D("KinRecoDiag_dR_lep_step7", "min #DeltaR(reco lep, gen lep);#DeltaR", 100, 0, 0.5);
    h_dR_jet_       = new TH1D("KinRecoDiag_dR_jet_step7", "min #DeltaR(jet, b quark);#DeltaR",       100, 0, 1.0);
    h_count_        = new TH1D("KinRecoDiag_count", "bookkeeping (unweighted)", 8, 0, 8);
    const char* labels[8] = {"gen evts", "gen complete", "gen W+W- found", "reco evts",
                             "lep matched", "anlep matched", "b matched", "bbar matched"};
    for (int i = 0; i < 8; ++i) h_count_->GetXaxis()->SetBinLabel(i + 1, labels[i]);

    for (TH1* h : std::initializer_list<TH1*>{h_wmass_, h_mbl_, h_angle_jet_, h_angle_lep_, h_fE_jet_, h_fE_lep_,
                                               h_fE_jet_vs_pt_, h_fE_lep_vs_pt_, h_dR_lep_, h_dR_jet_})
        h->Sumw2();
}

KinRecoInputMaker::GenTtbar KinRecoInputMaker::FindGenTtbar(const NanoAODBranchReader& reader)
{
    GenTtbar g;
    auto* pt   = reader.RequiredFloatArray("GenPart_pt");
    auto* eta  = reader.RequiredFloatArray("GenPart_eta");
    auto* phi  = reader.RequiredFloatArray("GenPart_phi");
    auto* mass = reader.RequiredFloatArray("GenPart_mass");
    auto* pdg  = reader.RequiredIntArray("GenPart_pdgId");
    const int n = pt->GetSize();

    auto mother = [&](int i) { return static_cast<int>(reader.GetIntArrayValue("GenPart_genPartIdxMother", i)); };
    auto p4 = [&](int i) {
        TLorentzVector v;
        v.SetPtEtaPhiM((*pt)[i], (*eta)[i], (*phi)[i], (*mass)[i]);
        return v;
    };
    // First ancestor whose |pdgId| differs from the particle's own |pdgId|
    // (skips the copy chain t->t->t, W->W, l->l (FSR), ...). Returns index or -1.
    auto firstDifferentAncestor = [&](int i) {
        const int self = std::abs((*pdg)[i]);
        int m = mother(i);
        int guard = 0;
        while (m >= 0 && m < n && std::abs((*pdg)[m]) == self && guard++ < 100) m = mother(m);
        return (m >= 0 && m < n) ? m : -1;
    };

    double lepPt = -1, anLepPt = -1;
    for (int i = 0; i < n; ++i) {
        const int id = (*pdg)[i];
        const int aid = std::abs(id);

        // b quark: direct daughter of the top (first b in the chain).
        if (aid == 5) {
            const int m = mother(i);
            if (m >= 0 && m < n && (*pdg)[m] == 6 * (id > 0 ? 1 : -1)) {
                if (id > 0 && !g.hasB)    { g.b    = p4(i); g.hasB = true; }
                if (id < 0 && !g.hasBbar) { g.bbar = p4(i); g.hasBbar = true; }
            }
            continue;
        }

        if (aid == 24) {
            const Long64_t flags = reader.GetIntArrayValue("GenPart_statusFlags", i);
            if (!flag(flags, kIsLastCopy)) continue;
            const int a = firstDifferentAncestor(i);
            if (a < 0 || (*pdg)[a] != 6 * (id > 0 ? 1 : -1)) continue;   // W+ <- t, W- <- tbar
            if (id > 0) { g.wPlus  = p4(i); g.hasWplus  = true; }
            else        { g.wMinus = p4(i); g.hasWminus = true; }
            continue;
        }

        if (aid == 11 || aid == 13) {
            if (reader.GetIntArrayValue("GenPart_status", i) != 1) continue;    // final state
            const Long64_t flags = reader.GetIntArrayValue("GenPart_statusFlags", i);
            if (!flag(flags, kIsPrompt) || !flag(flags, kFromHardProcess) || flag(flags, kIsTauDecayProduct)) continue;
            const int w = firstDifferentAncestor(i);
            // l- (pdgId>0) <- W-,  l+ (pdgId<0) <- W+
            if (w < 0 || (*pdg)[w] != (id > 0 ? -24 : 24)) continue;
            const double ptv = (*pt)[i];
            if (id > 0 && ptv > lepPt)   { g.lep   = p4(i); g.hasLep   = true; lepPt   = ptv; }
            if (id < 0 && ptv > anLepPt) { g.anLep = p4(i); g.hasAnLep = true; anLepPt = ptv; }
        }
    }
    return g;
}

void KinRecoInputMaker::FillGen(const GenTtbar& gen, double w)
{
    ++nGenAll_;
    h_count_->Fill(0);
    if (gen.hasWplus && gen.hasWminus) {
        h_count_->Fill(2);
        h_wmass_->Fill(gen.wPlus.M(),  w);
        h_wmass_->Fill(gen.wMinus.M(), w);
    }
    if (!gen.Complete()) return;
    ++nGenComplete_;
    h_count_->Fill(1);
    h_mbl_->Fill((gen.anLep + gen.b).M(),  w);
    h_mbl_->Fill((gen.lep + gen.bbar).M(), w);
}

void KinRecoInputMaker::FillReco(const GenTtbar& gen,
                                 const TLorentzVector& recoLep, const TLorentzVector& recoAnLep,
                                 const std::vector<TLorentzVector>& jets, double w)
{
    ++nReco_;
    h_count_->Fill(3);

    // ---- leptons: charge-ordered reco lepton vs same-charge gen lepton
    auto fillLep = [&](const TLorentzVector& reco, const TLorentzVector& genl, int bin) {
        if (reco.E() <= 0) return false;
        const double dR = reco.DeltaR(genl);
        h_dR_lep_->Fill(dR, w);
        if (dR >= kLepMatchDR) return false;
        h_fE_lep_->Fill(genl.E() / reco.E(), w);
        h_angle_lep_->Fill(reco.Angle(genl.Vect()), w);
        h_fE_lep_vs_pt_->Fill(reco.Pt(), genl.E() / reco.E(), w);
        h_count_->Fill(bin);
        return true;
    };
    bool okL = gen.hasLep   && fillLep(recoLep,   gen.lep,   4);
    bool okA = gen.hasAnLep && fillLep(recoAnLep, gen.anLep, 5);
    if (okL && okA) ++nRecoLepMatch_;

    // ---- jets: closest selected jet to each top b quark, dR<0.4, and the two
    // b quarks must not pick the same jet.
    auto closest = [&](const TLorentzVector& q, double& dRmin) {
        int best = -1; dRmin = 999.;
        for (size_t j = 0; j < jets.size(); ++j) {
            const double dR = jets[j].DeltaR(q);
            if (dR < dRmin) { dRmin = dR; best = static_cast<int>(j); }
        }
        return best;
    };
    double dRb = 999., dRbb = 999.;
    int jb  = gen.hasB    ? closest(gen.b,    dRb)  : -1;
    int jbb = gen.hasBbar ? closest(gen.bbar, dRbb) : -1;
    if (jb  >= 0) h_dR_jet_->Fill(dRb,  w);
    if (jbb >= 0) h_dR_jet_->Fill(dRbb, w);
    if (jb >= 0 && jb == jbb) return;     // ambiguous, skip both
    int nJetMatched = 0;
    auto fillJet = [&](int j, double dR, const TLorentzVector& q, int bin) {
        if (j < 0 || dR >= kJetMatchDR || jets[j].E() <= 0) return;
        h_fE_jet_->Fill(q.E() / jets[j].E(), w);
        h_angle_jet_->Fill(jets[j].Angle(q.Vect()), w);
        h_fE_jet_vs_pt_->Fill(jets[j].Pt(), q.E() / jets[j].E(), w);
        h_count_->Fill(bin);
        ++nJetMatched;
    };
    fillJet(jb,  dRb,  gen.b,    6);
    fillJet(jbb, dRbb, gen.bbar, 7);
    if (nJetMatched == 2) ++nRecoJetMatch_;
}

void KinRecoInputMaker::PrintSummary() const
{
    std::cout << "[KinRecoInputMaker] gen events: " << nGenAll_
              << ", gen complete (l+ l- b bbar): " << nGenComplete_
              << " | reco events: " << nReco_
              << ", both leptons matched: " << nRecoLepMatch_
              << ", both b jets matched: " << nRecoJetMatch_ << std::endl;
}
