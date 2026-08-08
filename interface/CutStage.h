#ifndef CUTSTAGE_H
#define CUTSTAGE_H

// TODO(agreed 2026-08-06): rename CutStage -> SelectionStage (file, enum,
// cutIndex()) in a later, separate step. Left as CutStage for now so this
// doesn't get bundled into the same diff as unrelated work - purely a
// rename, no behavior change, but still deserves its own small step and
// regression check like everything else here.

#include <cstddef>

// Named selection-stage indices for the per-stage control-plot histogram
// arrays (h_Lep1pt[10], h_Num_PV[10], etc. - see Analysis::DeclareHistos()/
// cutflowName[10] in interface/Analysis.h and src/Analysis.cpp).
//
// This enum is a DIRECT TRANSCRIPTION of what Analysis::Loop() actually
// fills today - it does not invent or "clean up" any stage. Traced from
// src/Analysis.cpp (both the cutflowName[i] assignments in SetVariables()
// and the FillHisto(h_*[i], ...) call sites in Loop()):
//
//   cutflowName[0] = "Step_0"              -> filled right after the
//       Trigger()+PV_npvsGood>=1 requirement, BEFORE the 2-lepton cut.
//   cutflowName[1] = "Step_1"              -> filled after NumIsoLeptons(2)
//       + ThirdLeptonVeto() + LeptonsPtAddtional() + DiLeptonMassCut(),
//       i.e. right after LeptonSFApply()/TriggerSFApply().
//   cutflowName[2] = "Step_2"              -> filled after ZVetoCut().
//   cutflowName[3] = "Step_3"              -> filled after NumJetCut().
//   cutflowName[4] = "Step_4"              -> filled after METCut().
//   cutflowName[5] = "bTagged Jet >= 1"    -> filled after BTaggingSFApply()
//       + NumbJetCut(). (Name documents the *typical* config, not a literal
//       check done at this exact point - NumbJetCut()'s actual threshold is
//       config-driven.)
//   cutflowName[6] = "bTagged Jet >= 2"    -> NOT currently filled anywhere
//       in Loop() (checked - no h_*[6] FillHisto call exists as of this
//       writing). Declared/named but effectively dead. Kept here as-is
//       rather than silently dropped, so this doesn't quietly change output
//       structure - flag for cleanup once confirmed genuinely unused.
//   cutflowName[7] = "bTagged Jet == 2"    -> same as [6]: declared, not
//       filled.
//   cutflowName[8] = "Top-Recon."          -> filled only inside
//       `if (isKinSol)` - i.e. only for events where KinSolv produced a
//       solution.
//   cutflowName[9] = "Top-Pt-Rewight"      -> same as [6]/[7]: declared,
//       not filled.
//
// If a later step DOES start filling stage 6/7/9 (or removes them), update
// this file and the comment above in the same change - that's the whole
// point of centralizing this mapping instead of leaving raw integers
// scattered through Loop().
enum class CutStage : std::size_t {
    AfterTriggerAndPV   = 0,  // "Step_0"
    AfterDileptonCuts    = 1,  // "Step_1" (2 leptons, 3rd-lepton veto, lepton pT, OS, dilepton mass)
    AfterZVeto           = 2,  // "Step_2"
    AfterJetMultiplicity = 3,  // "Step_3"
    AfterMETCut          = 4,  // "Step_4"
    AfterBTagMultiplicity = 5, // "bTagged Jet >= 1"
    Reserved6            = 6,  // "bTagged Jet >= 2" - declared, not filled (see above)
    Reserved7            = 7,  // "bTagged Jet == 2" - declared, not filled (see above)
    AfterTopReconstruction = 8, // "Top-Recon." (isKinSol == true)
    Reserved9            = 9,  // "Top-Pt-Rewight" - declared, not filled (see above)
    Count                = 10
};

constexpr std::size_t kNumCutStages = static_cast<std::size_t>(CutStage::Count);

// Use this instead of a raw static_cast at every call site - e.g.
//   h_Lep1pt[cutIndex(CutStage::AfterZVeto)]->Fill(...)
// instead of
//   h_Lep1pt[2]->Fill(...)
constexpr std::size_t cutIndex(CutStage stage) {
    return static_cast<std::size_t>(stage);
}

#endif  // CUTSTAGE_H
