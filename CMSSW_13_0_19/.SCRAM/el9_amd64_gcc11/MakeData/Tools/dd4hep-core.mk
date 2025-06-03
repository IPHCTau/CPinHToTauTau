ALL_TOOLS      += dd4hep-core
dd4hep-core_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/dd4hep/v01-23x-e2cf61cff4a911ed9506a3153f3788a1/include
dd4hep-core_EX_LIB := DDCore DDParsers
dd4hep-core_EX_USE := root_cxxdefaults root boost xerces-c clhep
dd4hep-core_EX_FLAGS_CPPDEFINES  := -DDD4HEP_USE_GEANT4_UNITS=1
dd4hep-core_EX_FLAGS_LISTCOMPONENTS  := $(DD4HEP_CORE_BASE)/bin/listcomponents_dd4hep

