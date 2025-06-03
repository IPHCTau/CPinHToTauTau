ALL_TOOLS      += geant4_interface
geant4_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/geant4/10.7.2-0e62b7ff14b69147f66862fcd6cd565b/include/Geant4 /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/geant4/10.7.2-0e62b7ff14b69147f66862fcd6cd565b/include
geant4_interface_EX_USE := clhep vecgeom zlib expat xerces-c root_cxxdefaults
geant4_interface_EX_FLAGS_CPPDEFINES  := -DGNU_GCC -DG4V9
geant4_interface_EX_FLAGS_CXXFLAGS  := -ftls-model=global-dynamic -pthread -DG4GEOM_USE_USOLIDS

