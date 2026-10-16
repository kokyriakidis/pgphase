#include "pgcat_command.hpp"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

// pgcat (tools/pgcat, `make pgcat`) is a Rust program on the GBZ library; the
// pgphase commands that need the graph run it with the same arguments.

namespace {

// pgcat next to this executable, as `make pgcat` and the release bundle put it.
std::string pgcat_beside_self(const char* argv0) {
    std::string self;
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0) {
        self.assign(buf, static_cast<size_t>(n));
    } else if (argv0 != nullptr && std::strchr(argv0, '/') != nullptr) {
        self = argv0;
    }
    const size_t slash = self.rfind('/');
    if (slash == std::string::npos) return {};
    const std::string candidate = self.substr(0, slash + 1) + "pgcat";
    return access(candidate.c_str(), X_OK) == 0 ? candidate : std::string();
}

} // namespace

int run_pgcat(const char* subcommand, int argc, char* argv[], const char* argv0) {
    std::string pgcat;
    if (const char* env = std::getenv("PGCAT"); env != nullptr && *env != '\0') pgcat = env;
    if (pgcat.empty()) pgcat = pgcat_beside_self(argv0);
    if (pgcat.empty()) pgcat = "pgcat";  // from PATH

    std::vector<char*> args;
    args.push_back(const_cast<char*>(pgcat.c_str()));
    args.push_back(const_cast<char*>(subcommand));
    for (int i = 1; i < argc; ++i) args.push_back(argv[i]);
    args.push_back(nullptr);
    execvp(pgcat.c_str(), args.data());
    std::cerr << "error: cannot run " << pgcat << " (" << std::strerror(errno) << ").\n"
              << "pgphase runs pgcat " << subcommand << ": build pgcat with `make pgcat` (needs Rust),\n"
              << "or point PGCAT at the binary.\n";
    return 127;
}
