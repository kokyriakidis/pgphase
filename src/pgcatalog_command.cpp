#include "pgcatalog_command.hpp"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

// pgcatalog (tools/pgcatalog, `make pgcatalog`) is a Rust program on the GBZ
// library; the pgphase commands that need the graph run it with the same
// arguments.

std::string find_tool(const char* name, const char* env, const char* argv0) {
    if (const char* set = std::getenv(env); set != nullptr && *set != '\0') return set;
    std::string self;
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0) {
        self.assign(buf, static_cast<size_t>(n));
    } else if (argv0 != nullptr && std::strchr(argv0, '/') != nullptr) {
        self = argv0;
    }
    if (const size_t slash = self.rfind('/'); slash != std::string::npos) {
        const std::string beside = self.substr(0, slash + 1) + name;
        if (access(beside.c_str(), X_OK) == 0) return beside;
    }
    return name;  // from PATH
}

int run_pgcatalog(const char* subcommand, int argc, char* argv[], const char* argv0) {
    const std::string pgcatalog = find_tool("pgcatalog", "PGCATALOG", argv0);
    std::vector<char*> args;
    args.push_back(const_cast<char*>(pgcatalog.c_str()));
    args.push_back(const_cast<char*>(subcommand));
    for (int i = 1; i < argc; ++i) args.push_back(argv[i]);
    args.push_back(nullptr);
    execvp(pgcatalog.c_str(), args.data());
    std::cerr << "error: cannot run " << pgcatalog << " (" << std::strerror(errno) << ").\n"
              << "pgphase runs pgcatalog " << subcommand << ": build pgcatalog with `make pgcatalog` (needs Rust),\n"
              << "or point PGCATALOG at the binary.\n";
    return 127;
}
