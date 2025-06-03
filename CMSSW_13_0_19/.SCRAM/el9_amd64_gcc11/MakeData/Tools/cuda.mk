ALL_TOOLS      += cuda
ALL_LIB_TYPES += CUDA_LIB
cuda_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cuda/11.5.2-91421d5831fe9d02c7a55c817b17889e/include
cuda_EX_LIB := cudart cudadevrt nvToolsExt
cuda_EX_CUDA_LIB := cudadevrt
cuda_EX_USE := cuda-stubs
cuda_EX_FLAGS_CUDA_FLAGS  := -std=c++17 -O3 --generate-line-info --source-in-ptx --display-error-number --expt-relaxed-constexpr --extended-lambda -gencode arch=compute_60,code=[sm_60,compute_60] -gencode arch=compute_70,code=[sm_70,compute_70] -gencode arch=compute_75,code=[sm_75,compute_75] -Wno-deprecated-gpu-targets -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored --cudart shared
cuda_EX_FLAGS_CUDA_HOST_CXXFLAGS  := -std=c++17
cuda_EX_FLAGS_REM_CUDA_HOST_CXXFLAGS  := -std=% %potentially-evaluated-expression

