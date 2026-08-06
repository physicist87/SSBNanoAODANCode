#ifndef JETID_H
#define JETID_H

#include <TString.h>
#include "NanoAODBranchReader.h"

// NanoAODv15 PUPPI jet ID (Jet_jetId branch removed upstream; POG pseudocode
// reimplemented directly from the jet energy fractions/multiplicities - see
// CMS NanoAODv15 documentation for the source pseudocode this mirrors).
//
// Extracted out of Analysis.cpp (same precedent as the earlier
// NanoAODBranchReader extraction - see NOTES.md "OO refactor" entries):
// pulling a self-contained piece of physics logic out of the monolithic
// event-loop file so it can be read, reasoned about, and (eventually) unit
// tested independently of the rest of Analysis.
//
// idx indexes the raw NanoAOD Jet collection (same indexing convention
// Analysis uses for jets_pt/jets_eta/etc, not the post-selection v_jet_idx).
class JetID {
public:
    // branchReader must outlive this object. runPeriod is captured by value
    // (a TString copy, cheap) rather than by reference: callers construct a
    // JetID fresh per lookup (see Analysis::PassConfiguredJetId) instead of
    // keeping one around for the whole job, which sidesteps any question of
    // whether Analysis's own RunPeriod member is populated yet at the point
    // a longer-lived JetID would have been constructed (it's read from
    // config in SetVariables(), not at Analysis's own construction time).
    JetID(const NanoAODBranchReader& branchReader, TString runPeriod);

    // POG tight working point, from jet energy fractions/multiplicities
    // (era-dependent eta binning/thresholds - see .cpp).
    bool PassTight(int idx) const;

    // Tight + a muon/charged-EM veto in the barrel-ish region (era-dependent
    // eta boundary) - see .cpp.
    bool PassTightLepVeto(int idx) const;

    // Dispatches on a Jet_ID config string ("PFTight"/"PFTightLepVeto").
    // Unknown values warn and fall back to PFTightLepVeto - matches the
    // previous Analysis::PassConfiguredJetId behavior exactly.
    bool PassConfigured(int idx, const TString& jetIdConfig) const;

private:
    const NanoAODBranchReader& branchReader_;
    TString runPeriod_;

    // Bounds-checked float-array lookup, matching the getF lambda this class
    // was extracted from - branchReader_.floatVectors is a public map by
    // design (see NanoAODBranchReader.h), this is just the small
    // out-of-range/missing-branch guard around it.
    float GetFloat(const char* name, int idx) const;
};

#endif  // JETID_H
