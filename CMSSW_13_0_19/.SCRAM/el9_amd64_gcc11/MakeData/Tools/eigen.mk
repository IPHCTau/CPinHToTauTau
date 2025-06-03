ALL_TOOLS      += eigen
eigen_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-9ac4aed18ac60d0189693c592862694d/include/eigen3
eigen_EX_FLAGS_CPPDEFINES  := -DEIGEN_DONT_PARALLELIZE
eigen_EX_FLAGS_CUDA_FLAGS  := --diag-suppress 20014

