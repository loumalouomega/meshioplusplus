//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// External includes
#include <gtest/gtest.h>

// System includes
#include <iostream>
#include <sstream>
#include <string>

// Project includes
#include "meshioplusplus/detail/source_location_compat.hpp"
#include "meshioplusplus/log.hpp"

using namespace meshioplusplus;

namespace {

// RAII redirect of std::cerr into a string buffer for the duration of a scope.
class CerrCapture {
public:
    CerrCapture() : mOld(std::cerr.rdbuf(mBuf.rdbuf())) {}
    ~CerrCapture() { std::cerr.rdbuf(mOld); }
    std::string str() const { return mBuf.str(); }

private:
    std::ostringstream mBuf;
    std::streambuf* mOld;
};

}  // namespace

// The threshold is memoized from MESHIOPLUSPLUS_LOG_LEVEL on first use, so these
// tests avoid asserting an absolute level; they check the internal consistency
// of enabled()/threshold() and the exact line format produced by write().

TEST(Log, LevelsAreOrdered) {
    EXPECT_LT(static_cast<int>(log::Level::Debug), static_cast<int>(log::Level::Info));
    EXPECT_LT(static_cast<int>(log::Level::Info), static_cast<int>(log::Level::Warn));
    EXPECT_LT(static_cast<int>(log::Level::Warn), static_cast<int>(log::Level::Error));
    EXPECT_LT(static_cast<int>(log::Level::Error), static_cast<int>(log::Level::Off));
}

TEST(Log, EnabledIsConsistentWithThreshold) {
    const int thr = static_cast<int>(log::threshold());
    for (auto lvl : {log::Level::Debug, log::Level::Info, log::Level::Warn, log::Level::Error}) {
        EXPECT_EQ(log::enabled(lvl), static_cast<int>(lvl) >= thr)
            << "level " << static_cast<int>(lvl) << " vs threshold " << thr;
    }
    // The Off level is never a message level, but by construction it is >= any
    // threshold below Off.
    EXPECT_TRUE(log::enabled(log::Level::Off));
}

TEST(Log, WriteEmitsFormattedLine) {
    // write() emits unconditionally (the level gate lives in debug/info/...), so
    // this is deterministic regardless of MESHIOPLUSPLUS_LOG_LEVEL.
    std::string out;
    {
        CerrCapture cap;
        log::write(log::Level::Warn, "payload 42", detail::source_location::current());
        out = cap.str();
    }
    EXPECT_NE(out.find("meshio warning"), std::string::npos) << out;
    EXPECT_NE(out.find("payload 42"), std::string::npos) << out;
    // The call site basename (this file) and a line number are included.
    EXPECT_NE(out.find("test_log.cpp:"), std::string::npos) << out;
    EXPECT_EQ(out.back(), '\n');
}

TEST(Log, WriteLabelsEachLevel) {
    struct Case {
        log::Level level;
        const char* label;
    };
    for (const Case& c :
         {Case{log::Level::Debug, "meshio debug"}, Case{log::Level::Info, "meshio info"},
          Case{log::Level::Warn, "meshio warning"}, Case{log::Level::Error, "meshio error"}}) {
        std::string out;
        {
            CerrCapture cap;
            log::write(c.level, "msg", detail::source_location::current());
            out = cap.str();
        }
        EXPECT_NE(out.find(c.label), std::string::npos) << out;
    }
}
