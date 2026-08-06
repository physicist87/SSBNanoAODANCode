#ifndef JET_H
#define JET_H

#include <TLorentzVector.h>

// One physics jet's raw (pre-JEC/JER), per-event NanoAOD inputs, read
// together in a single per-jet pass (see Analysis::BuildRawJets(), called
// from MakeJetCollection()) instead of the same bounds-checked "branch
// present? array long enough? else default" idiom being repeated separately
// for each field across the function. Consolidating into one struct/one
// loop removes a real (if currently harmless) latent risk: before this
// refactor, MakeJetCollection() read Jet_rawFactor/area/muonSubtrFactor/
// genJetIdx in one loop over nJets and Jet_chEmEF/neEmEF in a SEPARATE
// second loop, also over nJets - nothing structurally guaranteed those two
// loops stayed in sync if either field list changed later.
//
// idx indexes the raw NanoAOD Jet_ collection - same convention as
// jets_pt/jets_eta/etc and JetID's idx. This struct intentionally does NOT
// replace SSBCorrections::ApplyJetCorrections()/ApplyType1METWithCorrT1()'s
// existing parallel-vector signatures (that verified, tutorial-cross-checked
// code is out of scope for this step - see NOTES.md item 28) - Analysis
// projects vector<RawJet> back into those parameter vectors right before
// calling them.
struct RawJet {
    TLorentzVector p4;             // as stored (pre-undo of rawFactor) - Jet_pt/eta/phi/mass
    float rawFactor = 0.0f;        // Jet_rawFactor
    float area = 0.5f;             // Jet_area
    float muonSubtrFactor = 0.0f;  // Jet_muonSubtrFactor
    float chEmEF = 0.0f;           // Jet_chEmEF (Type-1 MET selection cut)
    float neEmEF = 0.0f;           // Jet_neEmEF (Type-1 MET selection cut)
    int genJetIdx = -1;            // Jet_genJetIdx (Short_t in v15, widened) - -1 = no gen match
};

// One CorrT1METJet_ collection jet's raw inputs - NanoAOD's separate,
// lower-pT-threshold jet collection needed to complete Type-1 MET
// propagation (see SSBCorrections::ApplyType1METWithCorrT1). rawPt is
// already raw (no rawFactor to undo, unlike RawJet::p4).
struct CorrT1Jet {
    float rawPt = 0.0f;
    float eta = 0.0f;
    float phi = 0.0f;
    float area = 0.5f;
    float muonSubtrFactor = 0.0f;
};

#endif  // JET_H
