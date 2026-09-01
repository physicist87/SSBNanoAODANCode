#ifndef DEBUGTOOLS_H
#define DEBUGTOOLS_H

#include <cstdint>
#include <iostream>
#include <map>
#include <string>

// Logging/debug infrastructure: Logger (leveled cout/cerr wrapper),
// DebugEventFilter (trace one run:lumi:event), CorrectionFallbackCounter
// (counts correctionlib fallbacks). No third-party logging lib, by design.

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

// One Logger per Analysis job, owned as a member (not a global) to avoid
// static-init-order issues. Messages above the configured level are dropped
// cheaply via a null stream, so per-event Trace prints cost nothing by default.
class Logger {
public:
    explicit Logger(LogLevel level = LogLevel::Info) : level_(level) {}

    void SetLevel(LogLevel level) { level_ = level; }
    LogLevel Level() const { return level_; }

    // Returns std::cerr/std::cout prefixed with the level tag if the level
    // passes the threshold, otherwise a null stream that discards everything
    // written to it - callers don't need their own `if (level_ >= X)` guard.
    // const: logging doesn't mutate the Logger's logical state (level_ is
    // only read here), so this can be called from const member functions
    // like Analysis::ChannelIndex() without needing a mutable logger_.
    std::ostream& Stream(LogLevel lvl) const {
        if (static_cast<int>(lvl) > static_cast<int>(level_)) {
            return NullStream();
        }
        std::ostream& out = (lvl == LogLevel::Error || lvl == LogLevel::Warning) ? std::cerr : std::cout;
        out << "[" << ToString(lvl) << "] ";
        return out;
    }

    std::ostream& Error()   const { return Stream(LogLevel::Error); }
    std::ostream& Warning() const { return Stream(LogLevel::Warning); }
    std::ostream& Info()    const { return Stream(LogLevel::Info); }
    std::ostream& Debug()   const { return Stream(LogLevel::Debug); }
    std::ostream& Trace()   const { return Stream(LogLevel::Trace); }

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

// DebugEventFilter: traces exactly one run:lumi:event through every stage,
// configured via DebugMode/DebugRun/DebugLumi/DebugEvent. A field left at 0
// (default) matches any value for that field.
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

// CorrectionFallbackCounter: counts correctionlib evaluate() failures that
// fell back to a default, printed once as a summary at end of job.
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
