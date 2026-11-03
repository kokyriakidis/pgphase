#ifndef PGPHASE_PGCATALOG_COMMAND_HPP
#define PGPHASE_PGCATALOG_COMMAND_HPP

#include <string>

/// The program `name`: from the environment variable `env` if set, else beside
/// this executable (as `make` and the release bundle put it), else `name` itself
/// (found on PATH when run). argv0 is the pgphase executable.
std::string find_tool(const char* name, const char* env, const char* argv0);

/// Runs `pgcatalog SUBCOMMAND ARGS`, where argv[1..] are ARGS (argv[0] is the
/// pgphase command name), replacing this process. pgcatalog builds what pgphase
/// reads from the graph: `build` the site catalogs (pgphase
/// build-snarl-catalog). Returns only if pgcatalog cannot be run.
int run_pgcatalog(const char* subcommand, int argc, char* argv[], const char* argv0);

/// pgphase index: the reads' graph alignments (GAM or GAF) indexed by
/// reference position, by `pgcatalog ref-map` (once per graph) and `dromos
/// index`. argv[0] is the command name. Returns the exit status.
int run_index(int argc, char* argv[], const char* argv0);

#endif
