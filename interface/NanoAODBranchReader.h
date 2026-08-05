#ifndef NANOAOD_BRANCH_READER_H
#define NANOAOD_BRANCH_READER_H

#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

// Owns the dynamic branch-name -> TTreeReaderValue/Array bindings for a
// NanoAOD TChain, and knows how to read them back out in a way that is
// resilient to NanoAOD storage-type changes between versions (Int_t <->
// UChar_t/Short_t/UInt_t) and to a handful of outright branch renames
// (e.g. MET_pt -> PFMET_pt, fixedGridRhoFastjetAll -> Rho_fixedGridRhoFastjetAll).
//
// This is the version-adapter boundary: Analysis (the physics selection code)
// should not need to know or care which underlying map a given branch's
// value lives in, or whether its name changed between NanoAOD versions -
// it just asks this class for a value by branch name.
class NanoAODBranchReader {
public:
    NanoAODBranchReader(TChain *chain, TTreeReader &reader);

    // Reads branchListFile (name, category, dataType-hint, varType-hint per
    // line) and, for each branch listed, asks the TChain directly what its
    // real type and single/vector structure are (TBranch/TLeaf), then binds
    // a TTreeReaderValue/Array of the right type into the matching map below.
    // dataType/varType from the file are only used to print a mismatch
    // warning, never trusted blindly - this is what lets the same branch
    // list & code work across NanoAOD versions without editing this
    // function every time a branch's storage width changes.
    // isData is passed in (rather than read from a member) so this class
    // doesn't need to know anything about the analysis beyond "which
    // branches to bind".
    void InitBranches(const std::string &branchListFile, bool isData);

    // Version-agnostic integer accessors: look a branch up across whichever
    // integer-like map it actually landed in (Int_t/UInt_t/UChar_t/Short_t)
    // and widen the result, instead of the caller needing to know/guess.
    Long64_t GetIntArrayValue(const std::string &branchName, int idx) const;
    Long64_t GetIntSingleValue(const std::string &branchName) const;
    bool BranchIsAvailable(const std::string &branchName) const;

    // For branches renamed outright between versions (type auto-detection
    // can't help with a name change) - tries each candidate name in order
    // and returns the first one found.
    double GetFloatSingleValueByAlias(const std::vector<std::string> &candidateNames) const;
    TTreeReaderValue<Float_t>* GetFloatSinglePtrByAlias(const std::vector<std::string> &candidateNames) const;

    // True once InitBranches() has determined Jet_puId is present in the
    // input file. False for NanoAODv15 Puppi jets, where only the continuous
    // Jet_puIdDisc remains (no correctionlib SF yet, not wired in).
    bool JetPuIdAvailable() const { return puidBranchAvailable_; }

    // Maps for dynamic branch storage. Public by design: this class's whole
    // job is to be a thin, type-aware wrapper around TTreeReader bindings -
    // hiding them behind another layer of getters/setters that just forward
    // to the same map lookup would add indirection without adding safety.
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Bool_t>>> boolSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Int_t>>> intSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UInt_t>>> uintSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<ULong64_t>>> ulongSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<UChar_t>>> ucharSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderValue<Float_t>>> floatSingles;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Float_t>>> floatVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Bool_t>>> boolVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Int_t>>> intVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UInt_t>>> uintVectors;
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<UChar_t>>> ucharVectors;
    // NanoAODv15 narrows several index branches (Jet_genJetIdx, Muon_jetIdx,
    // Muon_genPartIdx, Muon_fsrPhotonIdx, ...) from Int_t to Short_t.
    std::unordered_map<std::string, std::unique_ptr<TTreeReaderArray<Short_t>>> shortVectors;

private:
    TChain *chain_;
    TTreeReader &reader_;
    bool puidBranchAvailable_ = true;
};

#endif // NANOAOD_BRANCH_READER_H
