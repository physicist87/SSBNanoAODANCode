#include "../interface/NanoAODBranchReader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <TBranch.h>
#include <TLeaf.h>

NanoAODBranchReader::NanoAODBranchReader(TChain *chain, TTreeReader &reader)
    : chain_(chain), reader_(reader) {}

void NanoAODBranchReader::InitBranches(const std::string &branchListFile, bool isData) {
    std::cout << "branchListFile : " << branchListFile << std::endl;
    std::ifstream infile(branchListFile);
    if (!infile) {
        throw std::runtime_error("Error: Could not open branch list file " + branchListFile);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines and comment lines (leading '#', ignoring leading
        // whitespace) - the branch list's own header/inline comments were
        // previously falling through to the "Invalid branch format" warning
        // below instead of being recognized as comments, spamming the log
        // with one warning per comment line on every single run.
        if (line.empty()) continue;
        {
            size_t firstNonWs = line.find_first_not_of(" \t\r\n");
            if (firstNonWs == std::string::npos || line[firstNonWs] == '#') continue;
        }

        std::istringstream iss(line);
        std::string branchName, objectType, dataType, varType;

        // Parse the line for branch name, object type, data type, and variable type (vector/single)
        // NOTE: dataType/varType are now only used as a fallback hint (printed on mismatch)
        // -- the branch list file's job is to say WHICH branches to read, not their type.
        if (std::getline(iss, branchName, ',') &&
            std::getline(iss, objectType, ',') &&
            std::getline(iss, dataType, ',') &&
            std::getline(iss, varType, ',')) {

            // Trim whitespace from each part (both ends)
            auto trim = [](std::string &s) {
                s.erase(s.find_last_not_of(" \t\n\r") + 1);
                s.erase(0, s.find_first_not_of(" \t\n\r"));
            };

            trim(branchName);
            trim(objectType);
            trim(dataType);
            trim(varType);
      	    if (isData && branchName.find("Gen") != std::string::npos) {continue;}
	    if (isData && branchName.find("gen") != std::string::npos) {continue;}
	    if (isData && branchName.find("Pileup_nTrueInt") != std::string::npos) {continue;}
	    if (isData && branchName.find("Jet_hadronFlavour") != std::string::npos) {continue;}
	    if (isData && branchName.find("L1PreFiringWeight") != std::string::npos) {continue;}

            // ----------------------------------------------------------------
            // Auto-detect the real type & structure of this branch from the
            // TChain itself (TLeaf), instead of trusting the branch list's
            // dataType/varType columns. This is what lets the same code read
            // NanoAODv9 and NanoAODv15 files: branch storage widths change
            // between versions (e.g. Int_t -> UChar_t / Short_t, UInt_t ->
            // Int_t) but we no longer need to update this function or the
            // branch list every time that happens.
            // ----------------------------------------------------------------
            if (!chain_) {
                std::cerr << "Warning: no TChain available while initializing branch " << branchName << std::endl;
                continue;
            }
            TBranch *br = chain_->GetBranch(branchName.c_str());
            if (!br) {
                std::cerr << "Warning: branch '" << branchName
                          << "' not found in input tree (removed/renamed in this NanoAOD version?) - skipping."
                          << std::endl;
                continue;
            }
            TLeaf *leaf = br->GetLeaf(branchName.c_str());
            if (!leaf) {
                std::cerr << "Warning: could not resolve leaf for branch '" << branchName << "' - skipping." << std::endl;
                continue;
            }

            std::string actualType = leaf->GetTypeName();
            // Arrays (e.g. Jet_pt[nJet]) have a leaf-count leaf pointing at the
            // counter branch (nJet); plain scalars (e.g. MET_pt) do not.
            bool isVector = (leaf->GetLeafCount() != nullptr);

            if (!dataType.empty() && !varType.empty()) {
                std::string expectedStruct = isVector ? "vector" : "single";
                if (dataType != actualType || varType != expectedStruct) {
                    std::cout << "[InitBranches] " << branchName
                              << ": branch list said (" << dataType << ", " << varType
                              << ") but tree actually has (" << actualType << ", " << expectedStruct
                              << ") - using the tree's actual type." << std::endl;
                }
            }

            // Initialize based on the auto-detected type and structure
            if (actualType == "Bool_t" && !isVector) {
                boolSingles[branchName] = std::make_unique<TTreeReaderValue<Bool_t>>(reader_, branchName.c_str());
            } else if (actualType == "Int_t" && !isVector) {
                intSingles[branchName] = std::make_unique<TTreeReaderValue<Int_t>>(reader_, branchName.c_str());
            } else if (actualType == "UInt_t" && !isVector) {
                uintSingles[branchName] = std::make_unique<TTreeReaderValue<UInt_t>>(reader_, branchName.c_str());
            } else if (actualType == "ULong64_t" && !isVector) {
                ulongSingles[branchName] = std::make_unique<TTreeReaderValue<ULong64_t>>(reader_, branchName.c_str());
            } else if (actualType == "Float_t" && !isVector) {
                floatSingles[branchName] = std::make_unique<TTreeReaderValue<Float_t>>(reader_, branchName.c_str());
            } else if (actualType == "UChar_t" && !isVector) {
                ucharSingles[branchName] = std::make_unique<TTreeReaderValue<UChar_t>>(reader_, branchName.c_str());
            } else if (actualType == "Bool_t" && isVector) {
                boolVectors[branchName] = std::make_unique<TTreeReaderArray<Bool_t>>(reader_, branchName.c_str());
            } else if (actualType == "Int_t" && isVector) {
                intVectors[branchName] = std::make_unique<TTreeReaderArray<Int_t>>(reader_, branchName.c_str());
            } else if (actualType == "UInt_t" && isVector) {
                uintVectors[branchName] = std::make_unique<TTreeReaderArray<UInt_t>>(reader_, branchName.c_str());
            } else if (actualType == "Float_t" && isVector) {
                floatVectors[branchName] = std::make_unique<TTreeReaderArray<Float_t>>(reader_, branchName.c_str());
            } else if (actualType == "UChar_t" && isVector) {
                ucharVectors[branchName] = std::make_unique<TTreeReaderArray<UChar_t>>(reader_, branchName.c_str());
            } else if (actualType == "Short_t" && isVector) {
                shortVectors[branchName] = std::make_unique<TTreeReaderArray<Short_t>>(reader_, branchName.c_str());
            } else {
                std::cerr << "Warning: Unsupported/unhandled leaf type for branch '" << branchName
                          << "': type=" << actualType << ", vector=" << isVector
                          << " - add a case in InitBranches() if this branch is actually needed." << std::endl;
            }

        } else {
            std::cerr << "Warning: Invalid branch format in list: " << line << std::endl;
        }
    }

    // Jet_puId (the categorical PU-jet-ID working point) does not exist in
    // NanoAODv15 -- only the continuous Jet_puIdDisc remains, and there is no
    // correctionlib SF for it yet (CMS JME recommendation is still pending).
    // Detect that here so PUID isn't silently mis-applied (e.g. rejecting all
    // low-pT jets because puId defaults to 0).
    puidBranchAvailable_ = BranchIsAvailable("Jet_puId");
    if (!puidBranchAvailable_) {
        std::cerr << "[WARNING] Jet_puId branch not found in input tree - "
                  << "disabling PU-jet-ID selection/weighting for this job "
                  << "(no replacement wired in yet; see Jet_puIdDisc)." << std::endl;
    }
}

// ------------------------------------------------------------------
// Version-agnostic integer accessors (see Analysis.h for rationale).
// ------------------------------------------------------------------
Long64_t NanoAODBranchReader::GetIntArrayValue(const std::string &branchName, int idx) const {
    if (idx < 0) return -999;

    auto itI = intVectors.find(branchName);
    if (itI != intVectors.end() && itI->second && idx < itI->second->GetSize()) {
        return static_cast<Long64_t>(itI->second->At(idx));
    }
    auto itS = shortVectors.find(branchName);
    if (itS != shortVectors.end() && itS->second && idx < itS->second->GetSize()) {
        return static_cast<Long64_t>(itS->second->At(idx));
    }
    auto itU = ucharVectors.find(branchName);
    if (itU != ucharVectors.end() && itU->second && idx < itU->second->GetSize()) {
        return static_cast<Long64_t>(itU->second->At(idx));
    }
    auto itUI = uintVectors.find(branchName);
    if (itUI != uintVectors.end() && itUI->second && idx < itUI->second->GetSize()) {
        return static_cast<Long64_t>(itUI->second->At(idx));
    }
    return -999; // branch missing / index out of range
}

Long64_t NanoAODBranchReader::GetIntSingleValue(const std::string &branchName) const {
    auto itI = intSingles.find(branchName);
    if (itI != intSingles.end() && itI->second) return static_cast<Long64_t>(**itI->second);
    auto itU = ucharSingles.find(branchName);
    if (itU != ucharSingles.end() && itU->second) return static_cast<Long64_t>(**itU->second);
    auto itUI = uintSingles.find(branchName);
    if (itUI != uintSingles.end() && itUI->second) return static_cast<Long64_t>(**itUI->second);
    auto itUL = ulongSingles.find(branchName);
    if (itUL != ulongSingles.end() && itUL->second) return static_cast<Long64_t>(**itUL->second);
    return -999;
}

double NanoAODBranchReader::GetFloatSingleValueByAlias(const std::vector<std::string> &candidateNames) const {
    for (const auto &name : candidateNames) {
        auto it = floatSingles.find(name);
        if (it != floatSingles.end() && it->second) return static_cast<double>(**it->second);
    }
    return 0.0; // none of the candidate names were found
}

TTreeReaderValue<Float_t>* NanoAODBranchReader::GetFloatSinglePtrByAlias(const std::vector<std::string> &candidateNames) const {
    for (const auto &name : candidateNames) {
        auto it = floatSingles.find(name);
        if (it != floatSingles.end() && it->second) return it->second.get();
    }
    return nullptr;
}

bool NanoAODBranchReader::BranchIsAvailable(const std::string &branchName) const {
    return intVectors.count(branchName)   || shortVectors.count(branchName) ||
           ucharVectors.count(branchName) || uintVectors.count(branchName)  ||
           floatVectors.count(branchName) || boolVectors.count(branchName)  ||
           intSingles.count(branchName)   || ucharSingles.count(branchName) ||
           uintSingles.count(branchName)  || ulongSingles.count(branchName) ||
           floatSingles.count(branchName) || boolSingles.count(branchName);
}
