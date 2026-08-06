#include "../interface/JetID.h"
#include <cmath>
#include <iostream>

JetID::JetID(const NanoAODBranchReader& branchReader, TString runPeriod)
    : branchReader_(branchReader), runPeriod_(std::move(runPeriod)) {}

float JetID::GetFloat(const char* name, int idx) const {
    auto it = branchReader_.floatVectors.find(name);
    if (it == branchReader_.floatVectors.end() || !it->second || idx >= it->second->GetSize()) return 0.0f;
    return it->second->At(idx);
}

bool JetID::PassTight(int idx) const {
    // chMultiplicity / neMultiplicity are UChar_t in NanoAODv15
    int chMult = static_cast<int>(branchReader_.GetIntArrayValue("Jet_chMultiplicity", idx));
    int neMult = static_cast<int>(branchReader_.GetIntArrayValue("Jet_neMultiplicity", idx));

    float eta    = GetFloat("Jet_eta", idx);
    float neHEF  = GetFloat("Jet_neHEF", idx);
    float neEmEF = GetFloat("Jet_neEmEF", idx);
    float chHEF  = GetFloat("Jet_chHEF", idx);

    float absEta = fabs(eta);
    bool passTight = false;

    if (runPeriod_.Contains("2016")) {
        if (absEta <= 2.4) {
            passTight = (neHEF < 0.9) && (neEmEF < 0.9) && ((chMult + neMult) > 1) && (chHEF > 0.0) && (chMult > 0);
        } else if (absEta <= 2.7) {
            passTight = (neHEF < 0.98) && (neEmEF < 0.99);
        } else if (absEta <= 3.0) {
            passTight = (neMult >= 1);
        } else {
            passTight = (neMult > 2) && (neEmEF < 0.9);
        }
    } else if (runPeriod_.Contains("2017") || runPeriod_.Contains("2018")) {
        if (absEta <= 2.6) {
            passTight = (neHEF < 0.9) && (neEmEF < 0.9) && ((chMult + neMult) > 1) && (chHEF > 0.0) && (chMult > 0);
        } else if (absEta <= 2.7) {
            passTight = (neHEF < 0.90) && (neEmEF < 0.99);
        } else if (absEta <= 3.0) {
            passTight = (neHEF < 0.9999);
        } else {
            passTight = (neMult > 2) && (neEmEF < 0.9);
        }
    } else {
        std::cerr << "[WARNING] JetID::PassTight: no PUPPI jet ID formula defined for RunPeriod "
                  << runPeriod_ << " - treating jet as failing tight ID." << std::endl;
        return false;
    }

    return passTight;
}

bool JetID::PassTightLepVeto(int idx) const {
    if (!PassTight(idx)) return false;

    float eta    = GetFloat("Jet_eta", idx);
    float muEF   = GetFloat("Jet_muEF", idx);
    float chEmEF = GetFloat("Jet_chEmEF", idx);
    float absEta = fabs(eta);

    // 2016: TightLepVeto only differs from Tight below |eta| <= 2.4
    // 2017/2018: below |eta| <= 2.7 (per the pasted pseudocode)
    float etaBoundary = runPeriod_.Contains("2016") ? 2.4 : 2.7;

    if (absEta <= etaBoundary) {
        return (muEF < 0.8) && (chEmEF < 0.8);
    }
    return true;  // outside the boundary, TightLepVeto == Tight
}

bool JetID::PassConfigured(int idx, const TString& jetIdConfig) const {
    if (jetIdConfig == "PFTightLepVeto") return PassTightLepVeto(idx);
    if (jetIdConfig == "PFTight")        return PassTight(idx);
    // "PFLoose"/"PFLooseLepVeto" (2016 v9-era working points) no longer exist
    // for PUPPI jets - the pasted CMS pseudocode only defines Tight/TightLepVeto.
    std::cerr << "[WARNING] JetID::PassConfigured: Jet_ID='" << jetIdConfig
              << "' has no PUPPI jet ID implementation - update configs to "
              << "PFTight or PFTightLepVeto. Defaulting to PFTightLepVeto." << std::endl;
    return PassTightLepVeto(idx);
}
