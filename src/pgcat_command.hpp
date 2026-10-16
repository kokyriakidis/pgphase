#ifndef PGPHASE_PGCAT_COMMAND_HPP
#define PGPHASE_PGCAT_COMMAND_HPP

/// Runs `pgcat SUBCOMMAND ARGS`, where argv[1..] are ARGS (argv[0] is the
/// pgphase command name), replacing this process. pgcat builds what pgphase
/// reads from the graph: `build` the site catalogs (pgphase
/// build-snarl-catalog), `index-gaf` the reads indexed by reference position
/// (pgphase index-gaf). argv0 is the pgphase executable, to find pgcat beside
/// it. Returns only if pgcat cannot be run.
int run_pgcat(const char* subcommand, int argc, char* argv[], const char* argv0);

#endif
