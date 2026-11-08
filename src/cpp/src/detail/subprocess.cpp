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
// Running an installed program as a separate process. See subprocess.hpp.

#include "subprocess.hpp"

// System includes
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <system_error>

// Project includes
#include "meshioplusplus/exceptions.hpp"

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>

#include <atomic>
#include <chrono>
#else
#include <stdlib.h>  // mkdtemp
#endif
#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace meshioplusplus::detail {

namespace {

bool sp_is_file(const std::filesystem::path& rPath) {
    std::error_code ec;
    return std::filesystem::is_regular_file(rPath, ec);
}

}  // namespace

#if defined(__EMSCRIPTEN__)

bool subprocess_available() noexcept {
    return false;
}

std::string find_executable(const std::string&, const char*) {
    return {};
}

int run_process(const std::vector<std::string>&, const std::string&, const std::string&) {
    throw Unsupported("meshio++: running an external program is not available in WebAssembly");
}

#else

bool subprocess_available() noexcept {
    return true;
}

std::string find_executable(const std::string& rName, const char* pOverride) {
    if (pOverride) {
        const char* p_value = std::getenv(pOverride);
        if (p_value && *p_value)
            return sp_is_file(p_value) ? std::string(p_value) : std::string();
    }
    const char* p_path = std::getenv("PATH");
    if (!p_path)
        return {};
#if defined(_WIN32)
    constexpr char kSeparator = ';';
    const std::string name = rName + ".exe";
#else
    constexpr char kSeparator = ':';
    const std::string& name = rName;
#endif
    const std::string path(p_path);
    std::size_t start = 0;
    while (start <= path.size()) {
        std::size_t end = path.find(kSeparator, start);
        if (end == std::string::npos)
            end = path.size();
        const std::string dir = path.substr(start, end - start);
        start = end + 1;
        if (dir.empty())
            continue;
        const std::filesystem::path candidate = std::filesystem::path(dir) / name;
        if (!sp_is_file(candidate))
            continue;
#if !defined(_WIN32)
        if (access(candidate.c_str(), X_OK) != 0)
            continue;
#endif
        return candidate.string();
    }
    return {};
}

#if defined(_WIN32)

namespace {

std::wstring sp_wide(const std::string& rText) {
    const int n = MultiByteToWideChar(CP_UTF8, 0, rText.c_str(), -1, nullptr, 0);
    std::wstring wide(static_cast<std::size_t>(n > 0 ? n : 1), L'\0');
    if (n > 0)
        MultiByteToWideChar(CP_UTF8, 0, rText.c_str(), -1, wide.data(), n);
    wide.resize(wide.size() - 1);
    return wide;
}

// One argument quoted as CommandLineToArgvW reads it back.
std::wstring sp_quote(const std::wstring& rArg) {
    if (!rArg.empty() && rArg.find_first_of(L" \t\n\v\"") == std::wstring::npos)
        return rArg;
    std::wstring out = L"\"";
    for (std::size_t i = 0;; ++i) {
        std::size_t slashes = 0;
        while (i < rArg.size() && rArg[i] == L'\\') {
            ++slashes;
            ++i;
        }
        if (i == rArg.size()) {
            out.append(slashes * 2, L'\\');
            break;
        }
        if (rArg[i] == L'"')
            out.append(slashes * 2 + 1, L'\\');
        else
            out.append(slashes, L'\\');
        out.push_back(rArg[i]);
    }
    out.push_back(L'"');
    return out;
}

}  // namespace

int run_process(const std::vector<std::string>& rArgv, const std::string& rWorkDir,
                const std::string& rLogPath) {
    if (rArgv.empty())
        throw std::runtime_error("meshio++: run_process: no program");
    std::wstring command;
    for (const std::string& r_arg : rArgv)
        command += (command.empty() ? L"" : L" ") + sp_quote(sp_wide(r_arg));
    SECURITY_ATTRIBUTES sa{};
    sa.nLength = sizeof(sa);
    sa.bInheritHandle = TRUE;
    HANDLE log = CreateFileW(sp_wide(rLogPath).c_str(), FILE_APPEND_DATA, FILE_SHARE_READ, &sa,
                             OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    STARTUPINFOW si{};
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = INVALID_HANDLE_VALUE;
    si.hStdOutput = log;
    si.hStdError = log;
    PROCESS_INFORMATION pi{};
    const std::wstring dir = sp_wide(rWorkDir);
    const BOOL ok = CreateProcessW(sp_wide(rArgv[0]).c_str(), command.data(), nullptr, nullptr,
                                   TRUE, CREATE_NO_WINDOW, nullptr, dir.c_str(), &si, &pi);
    if (log != INVALID_HANDLE_VALUE)
        CloseHandle(log);
    if (!ok)
        throw std::runtime_error("meshio++: cannot start '" + rArgv[0] + "'");
    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD code = 0;
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return static_cast<int>(code);
}

#else

int run_process(const std::vector<std::string>& rArgv, const std::string& rWorkDir,
                const std::string& rLogPath) {
    if (rArgv.empty())
        throw std::runtime_error("meshio++: run_process: no program");
    // fork and exec rather than posix_spawn, which has no portable
    // working-directory action: the child changes directory, then runs the
    // program directly (no shell), with its output in the log.
    std::vector<std::string> args = rArgv;
    std::vector<char*> argv;
    for (std::string& r_arg : args)
        argv.push_back(r_arg.data());
    argv.push_back(nullptr);
    const pid_t pid = fork();
    if (pid < 0)
        throw std::runtime_error("meshio++: cannot start '" + rArgv[0] + "'");
    if (pid == 0) {
        // The child: only async-signal-safe calls from here.
        const int null = open("/dev/null", O_RDONLY);
        if (null >= 0)
            dup2(null, STDIN_FILENO);
        // 0600: the log is the caller's; nobody else needs to read it.
        const int log = open(rLogPath.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0600);
        if (log >= 0) {
            dup2(log, STDOUT_FILENO);
            dup2(log, STDERR_FILENO);
        }
        if (chdir(rWorkDir.c_str()) != 0)
            _exit(126);
        execv(argv[0], argv.data());
        _exit(127);
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR)
            throw std::runtime_error("meshio++: lost track of '" + rArgv[0] + "'");
    }
    if (WIFEXITED(status))
        return WEXITSTATUS(status);
    return -1;
}

#endif

#endif  // __EMSCRIPTEN__

std::string make_temp_dir(const std::string& rTag) {
    std::error_code ec;
    // The directory is created exclusively and private below, so another user
    // can neither place it first nor read what a run leaves in it.
    const std::filesystem::path base =
        std::filesystem::temp_directory_path(ec);  // NOSONAR(cpp:S5443): see above
#if !defined(_WIN32)
    // mkdtemp: an exclusive create with mode 0700.
    std::string pattern = (base / ("meshioplusplus_" + rTag + "_XXXXXX")).string();
    if (mkdtemp(pattern.data()) != nullptr)
        return pattern;
#else
    // The temporary directory is the user's own on Windows; create_directory
    // fails on an existing name, so the loop takes a fresh one.
    static std::atomic<unsigned long> s_counter{0};
    const auto stamp = static_cast<unsigned long long>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    for (int attempt = 0; attempt < 100; ++attempt) {
        const std::filesystem::path dir =
            base / ("meshioplusplus_" + rTag + "_" + std::to_string(stamp) + "_" +
                    std::to_string(s_counter.fetch_add(1)));
        if (std::filesystem::create_directory(dir, ec))
            return dir.string();
    }
#endif
    throw std::runtime_error("meshio++: cannot create a temporary directory for " + rTag);
}

}  // namespace meshioplusplus::detail
