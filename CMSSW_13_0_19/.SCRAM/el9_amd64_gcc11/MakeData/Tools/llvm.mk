ALL_TOOLS      += llvm
llvm_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/llvm/12.0.1-34ef258eea96cb5b111081adc51a827d/include
llvm_EX_LIB := clang
llvm_EX_FLAGS_CXXFLAGS  := -D_DEBUG -D_GNU_SOURCE -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -O3 -fomit-frame-pointer -fPIC -Wno-enum-compare -Wno-strict-aliasing -fno-rtti
llvm_EX_FLAGS_LDFLAGS  := -Wl,-undefined -Wl,suppress

