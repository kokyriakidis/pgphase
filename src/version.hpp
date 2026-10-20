#ifndef PGPHASE_VERSION_HPP
#define PGPHASE_VERSION_HPP

// Release version. A build may override it, e.g. with the git description:
//   make CXXFLAGS+='-DPGPHASE_VERSION=\"$(git describe --tags --dirty)\"'
#ifndef PGPHASE_VERSION
#define PGPHASE_VERSION "0.2.0"
#endif

#endif
