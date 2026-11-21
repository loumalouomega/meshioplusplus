# Security policy

## Supported versions

Security fixes go into the latest release (the newest minor version of the current major version). Older releases are not patched; upgrade to the latest release first. The release history is in [`CHANGELOG.md`](CHANGELOG.md).

## Reporting a vulnerability

Please report a vulnerability **privately**, through GitHub's security advisories: open the repository's **Security** tab and choose **Report a vulnerability** ([how it works](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing-information-about-vulnerabilities/privately-reporting-a-security-vulnerability)). Do not open a public issue or pull request for it.

A useful report names the version, the platform and the entry point (Python, C API, CLI, browser build or MCP server), and attaches the smallest input that reproduces it.

You can expect an acknowledgement within 7 days and an assessment within 30 days. A confirmed issue is fixed in the latest release, and the advisory credits the reporter unless they prefer not to be named.

## What is in scope

meshio++ parses files from sources its users do not control, so these are in scope:

- **The readers.** Any input file that makes a native or Python reader crash, hang, exhaust memory, read or write out of bounds, or do anything other than return a mesh or raise `ReadError`. The native readers are fuzzed ([doc/fuzzing.md](doc/fuzzing.md)); a reproducer found outside the fuzzer is as welcome.
- **The writers**, when they write outside the path they were given, or write data the input did not contain.
- **The MCP server's path sandbox and the HTTP companion** ([doc/mcp.md](doc/mcp.md)): reading or writing outside the configured root, or reaching the server without its token.
- **The C API, the browser (WebAssembly) build and the language bindings**, for memory-safety issues reachable from a documented call.

Out of scope: vulnerabilities in third-party libraries meshio++ links (report those to the library; tell us if a meshio++ release should pin a fixed version), and denial of service through a mesh that is legitimately very large.

## Fuzzing

The readers run under AddressSanitizer and UndefinedBehaviorSanitizer on every pull request, and a libFuzzer campaign runs weekly; see [doc/fuzzing.md](doc/fuzzing.md). Every crash it finds is fixed with the minimized input committed as a regression test.
