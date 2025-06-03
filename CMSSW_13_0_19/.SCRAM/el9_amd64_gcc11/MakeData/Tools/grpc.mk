ALL_TOOLS      += grpc
grpc_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/grpc/1.35.0-4ae66c6d769e6f806f3b64de826cbe10/include
grpc_EX_LIB := grpc grpc++ grpc++_reflection
grpc_EX_USE := protobuf openssl pcre abseil-cpp c-ares re2
grpc_EX_FLAGS_SYSTEM_INCLUDE  := 1

