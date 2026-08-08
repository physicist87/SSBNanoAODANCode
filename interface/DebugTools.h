#ifndef DEBUGTOOLS_H
#define DEBUGTOOLS_H

#include <cstdint>
#include <iostream>
#include <map>
#include <string>

// ============================================================================
// Phase 1 debug/logging infrastructure - additive only.
//
// Nothing in Analysis/SSBCorrections currently includes or calls into this
// file. It exists so the *next* small step (migrating a handful of the
// existing 140+ std::cout call sites in Analysis.cpp, one function at a
// time) has something to migrate to, without this step itself touching the
// event loop. Keeping infrastructure-addition and infrastructure-adoption as
// two separate diffs means a regression failure after this step can't be
// caused by this file (nothing calls it yet), and a regression failure after
// the adoption step is easy to bisect to "which cout did we replace".
//
// Deliberately NOT a third-party logging library - this is a grid/cvmfs job
// environment, so a new dependency costs more than the convenience it buys.
// Everything here is std::cerr/std::cout underneath, same as the code today.
// ============================================================================

// ----------------------------------------------------------------------------
// LogLevel / Logger
// ----------------------------------------------------------------------------
enum class LogLevel {
    Error = 0,
    Warning = 1,
    Info = 2,
    Debug = 3,
    Trace = 4
};

inline const char* ToString(LogLevel lvl) {
    switch (lvl) {
        case LogLevel::Error:   return "ERROR";
        case LogLevel::Warning: return "WARN";
        case LogLevel::Info:    return "INFO";
        case LogLevel::Debug:   return "DEBUG";
        case LogLevel::Trace:   return "TRACE";
    }
    return "?";
}

// One Logger per Analysis job (own it as a member, don't make it a global -
// keeps this testable and avoids static-init-order issues with other
// globals in this codebase, e.g. the TROOT root(...) in main_ssb.cpp).
//
// Usage (once call sites are migrated - see note above):
//   logger_.Warning() << "[BTaggingSFApply] correctionlib evaluate() failed: " << e.what();
//   logger_.Debug()   << "nJet=" << v_jet_idx.size() << " nBJet=" << v_bjet_idx.size();
//
// Messages above the configured level are cheaply dropped (the `if` guards
// the stream formatting itself, not just the final write) - Trace-level
// per-event prints won't cost anything in a normal production job run at
// the default Info level.
class Logger {
public:
    explicit Logger(LogLevel level = LogLevel::Info) : level_(level) {}

    void SetLevel(LogLevel level) { level_ = level; }
    LogLevel Level() const { return level_; }

    // Returns std::cerr/std::cout prefixed with the level tag if the level
    // passes the threshold, otherwise a null stream that discards everything
    // written to it - callers don't need their own `if (level_ >= X)` guard.
    std::ostream& Stream(LogLevel lvl) {
        if (static_cast<int>(lvl) > static_cast<int>(level_)) {
            return NullStream();
        }
        std::ostream& out = (lvl == LogLevel::Error || lvl == LogLevel::Warning) ? std::cerr : std::cout;
        out << "[" << ToString(lvl) << "] ";
        return out;
    }

    std::ostream& Error()   { return Stream(LogLevel::Error); }
    std::ostream& Warning() { return Stream(LogLevel::Warning); }
    std::ostream& Info()    { return Stream(LogLevel::Info); }
    std::ostream& Debug()   { return Stream(LogLevel::Debug); }
    std::ostream& Trace()   { return Stream(LogLevel::Trace); }

private:
    LogLevel level_;

    static std::ostream& NullStream() {
        // A single shared discard sink - not thread-safe, but neither is
        // anything else in this single-threaded event loop.
        static struct NullBuf : std::streambuf {
            int overflow(int c) override { return c; }
        } buf;
        static std::ostream null_stream(&buf);
        return null_stream;
    }
};

// ----------------------------------------------------------------------------
// DebugEventFilter - "trace exactly this run:lumi:event through every stage"
// ----------------------------------------------------------------------------
// Config keys (read the usual way via TextReader, see SSBConfReader in
// Analysis - not wired up yet, this class only does the matching):
//   DebugMode  = true/false
//   DebugRun   = 315257
//   DebugLumi  = 112
//   DebugEvent = 18452731
//
// A field left at 0 (the default) matches any run/lumi/event, e.g. leaving
// DebugLumi=0 traces every lumi block of the given run+event... but since
// (run,lumi,event) should be unique per event in practice, ordinarily all
// three are set together to pin exactly one event.
class DebugEventFilter {
public:
    DebugEventFilter() = default;

    void Configure(bool enabled, unsigned int run, unsigned int lumi, unsigned long long event) {
        enabled_ = enabled;
        run_ = run;
        lumi_ = lumi;
        event_ = event;
    }

    bool Enabled() const { return enabled_; }

    bool Matches(unsigned int run, unsigned int lumi, unsigned long long event) const {
        if (!enabled_) return false;
        if (run_ != 0 && run_ != run) return false;
        if (lumi_ != 0 && lumi_ != lumi) return false;
        if (event_ != 0 && event_ != event) return false;
        return true;
    }

private:
    bool enabled_ = false;
    unsigned int run_ = 0;
    unsigned int lumi_ = 0;
    unsigned long long event_ = 0;
};

// ----------------------------------------------------------------------------
// CorrectionFallbackCounter - "did any correctionlib lookup silently fall
// back to a default value, and for how many events"
// ----------------------------------------------------------------------------
// Motivating example (src/Analysis.cpp, BTaggingSFApply()):
//   catch (const std::exception& e) {
//       std::cerr << "Error in BTaggingSFApply: " << e.what() << std::endl;
//       btag_sf_weight_ = 1.0; // Default on error
//   }
// This prints to cerr today, but it's one line among 140+ other std::cout/
// std::cerr calls in a grid job's log - easy to miss. Recording it here
// instead (or in addition) means a job-end summary can answer "did this run
// silently default any weights" in one line instead of grepping megabytes
// of log.
//
// Not wired into SSBCorrections/Analysis yet - same reasoning as Logger
// above (keep addition and adoption separate).
class CorrectionFallbackCounter {
public:
    void RecordFallback(const std::string& correctionName) {
        ++counts_[correctionName];
    }

    long Count(const std::string& correctionName) const {
        auto it = counts_.find(correctionName);
        return it == counts_.end() ? 0 : it->second;
    }

    bool HasAnyFallback() const { return !counts_.empty(); }

    void PrintSummary(std::ostream& out = std::cout) const {
        if (counts_.empty()) {
            out << "[CorrectionFallbackCounter] no fallbacks recorded.\n";
            return;
        }
        out << "[CorrectionFallbackCounter] fallback summary:\n";
        for (const auto& kv : counts_) {
            out << "    " << kv.first << " : " << kv.second << " event(s)\n";
        }
    }

private:
    std::map<std::string, long> counts_;
};

#endif  // DEBUGTOOLS_H
