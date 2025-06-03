ALL_TOOLS      += rocm
rocm_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/rocm/5.4.3-d077dcce910eb4e8562ad8848b8fcd3a/include
rocm_EX_LIB := amdhip64
rocm_EX_FLAGS_CPPDEFINES  := -D__HIP_PLATFORM_HCC__ -D__HIP_PLATFORM_AMD__
rocm_EX_FLAGS_ROCM_FLAGS  := -fgpu-rdc --offload-arch=gfx900 --offload-arch=gfx906 --offload-arch=gfx908 --offload-arch=gfx90a --target=x86_64-redhat-linux-gnu --gcc-toolchain=$(COMPILER_PATH)
rocm_EX_FLAGS_SYSTEM_INCLUDE  := 1

