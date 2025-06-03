ALL_TOOLS      += abseil-cpp
abseil-cpp_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/abseil-cpp/20210324.2-036b208600e4ad74d685693991a468c1/include
abseil-cpp_EX_LIB := absl_bad_any_cast_impl absl_bad_optional_access absl_bad_variant_access absl_base absl_city absl_civil_time absl_cord absl_debugging_internal absl_demangle_internal absl_examine_stack absl_exponential_biased absl_failure_signal_handler absl_flags_commandlineflag_internal absl_flags_commandlineflag absl_flags_config absl_flags_internal absl_flags_marshalling absl_flags_parse absl_flags_private_handle_accessor absl_flags_program_name absl_flags_reflection absl_flags absl_flags_usage_internal absl_flags_usage absl_graphcycles_internal absl_hash absl_hashtablez_sampler absl_int128 absl_leak_check_disable absl_leak_check absl_log_severity absl_malloc_internal absl_periodic_sampler absl_random_distributions absl_random_internal_distribution_test_util absl_random_internal_platform absl_random_internal_pool_urbg absl_random_internal_randen_hwaes_impl absl_random_internal_randen_hwaes absl_random_internal_randen_slow absl_random_internal_randen absl_random_internal_seed_material absl_random_seed_gen_exception absl_random_seed_sequences absl_raw_hash_set absl_raw_logging_internal absl_scoped_set_env absl_spinlock_wait absl_stacktrace absl_statusor absl_status absl_strerror absl_str_format_internal absl_strings_internal absl_strings absl_symbolize absl_synchronization absl_throw_delegate absl_time absl_time_zone absl_wyhash

ALL_TOOLS      += alpaka-cuda
alpaka-cuda_EX_USE := alpaka cuda
alpaka-cuda_EX_FLAGS_CUDA_FLAGS  := -DALPAKA_ACC_GPU_CUDA_ENABLED -UALPAKA_HOST_ONLY
alpaka-cuda_EX_FLAGS_CXXFLAGS  := -DALPAKA_ACC_GPU_CUDA_ENABLED -DALPAKA_HOST_ONLY
alpaka-cuda_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DALPAKA_ACC_GPU_CUDA_ENABLED -DALPAKA_HOST_ONLY

ALL_TOOLS      += alpaka
alpaka_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/alpaka/develop-20230621-c24a22ec5c5d0e67ae2ff41d221065df/include
alpaka_EX_USE := boost
alpaka_EX_FLAGS_CUDA_FLAGS  := -DALPAKA_DEFAULT_HOST_MEMORY_ALIGNMENT=128
alpaka_EX_FLAGS_CXXFLAGS  := -DALPAKA_DEFAULT_HOST_MEMORY_ALIGNMENT=128
alpaka_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DALPAKA_DEFAULT_HOST_MEMORY_ALIGNMENT=128

ALL_TOOLS      += alpaka-rocm
alpaka-rocm_EX_USE := alpaka rocm rocm-rocrand
alpaka-rocm_EX_FLAGS_CXXFLAGS  := -DALPAKA_ACC_GPU_HIP_ENABLED -DALPAKA_HOST_ONLY
alpaka-rocm_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DALPAKA_ACC_GPU_HIP_ENABLED -DALPAKA_HOST_ONLY
alpaka-rocm_EX_FLAGS_ROCM_FLAGS  := -DALPAKA_ACC_GPU_HIP_ENABLED -UALPAKA_HOST_ONLY

ALL_TOOLS      += alpaka-serial
alpaka-serial_EX_USE := alpaka
alpaka-serial_EX_FLAGS_CXXFLAGS  := -DALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED
alpaka-serial_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DALPAKA_ACC_CPU_B_SEQ_T_SEQ_ENABLED

ALL_TOOLS      += alpaka-tbb
alpaka-tbb_EX_USE := alpaka tbb
alpaka-tbb_EX_FLAGS_CXXFLAGS  := -DALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED
alpaka-tbb_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DALPAKA_ACC_CPU_B_TBB_T_SEQ_ENABLED

ALL_TOOLS      += alpgen

ALL_TOOLS      += blackhat
blackhat_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/blackhat/0.9.9-7aab172897dfb2b49eb67ba086524dd5/include
blackhat_EX_LIB := Ampl_eval BG BH BHcore CutPart Cut_wCI Cuteval Integrals Interface OLA RatPart Rateval Spinors assembly ratext
blackhat_EX_USE := qd

ALL_TOOLS      += boost_chrono
boost_chrono_EX_LIB := boost_chrono
boost_chrono_EX_USE := boost_system boost

ALL_TOOLS      += boost_filesystem
boost_filesystem_EX_LIB := boost_filesystem
boost_filesystem_EX_USE := boost_system boost

ALL_TOOLS      += boost_header
boost_header_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.80.0-1dda32caa80f65c9270532a6ef741559/include
boost_header_EX_USE := sockets root_cxxdefaults
boost_header_EX_FLAGS_CPPDEFINES  := -DBOOST_SPIRIT_THREADSAFE -DPHOENIX_THREADSAFE -DBOOST_MATH_DISABLE_STD_FPCLASSIFY -DBOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
boost_header_EX_FLAGS_CXXFLAGS  := -Wno-error=unused-variable
boost_header_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += boost_iostreams
boost_iostreams_EX_LIB := boost_iostreams
boost_iostreams_EX_USE := boost

ALL_TOOLS      += boost
boost_EX_LIB := boost_thread boost_date_time
boost_EX_USE := boost_header

ALL_TOOLS      += boost_mpi
boost_mpi_EX_LIB := boost_mpi
boost_mpi_EX_USE := boost boost_serialization

ALL_TOOLS      += boost_program_options
boost_program_options_EX_LIB := boost_program_options
boost_program_options_EX_USE := boost

ALL_TOOLS      += boost_python
boost_python_EX_LIB := boost_python39
boost_python_EX_USE := boost_header python3

ALL_TOOLS      += boost_regex
boost_regex_EX_LIB := boost_regex
boost_regex_EX_USE := boost

ALL_TOOLS      += boost_serialization
boost_serialization_EX_LIB := boost_serialization
boost_serialization_EX_USE := boost

ALL_TOOLS      += boost_system
boost_system_EX_LIB := boost_system
boost_system_EX_USE := boost

ALL_TOOLS      += boost_test
boost_test_EX_LIB := boost_unit_test_framework
boost_test_EX_USE := boost

ALL_TOOLS      += bz2lib
bz2lib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/bz2lib/1.0.6-2c1f18484cb66c30aba7929f2be5e7d4/include
bz2lib_EX_LIB := bz2
bz2lib_EX_USE := root_cxxdefaults

ALL_TOOLS      += c-ares
c-ares_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/c-ares/1.15.0-4766aa7a4eb7a186cf2e3b08a98f821c/include
c-ares_EX_LIB := cares

ALL_TOOLS      += catch2
catch2_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/catch2/2.13.6-84fde6927c396f2973a03635196e7a49/include

ALL_TOOLS      += cepgen
cepgen_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cepgen/1.2.3-ccbc6efda2caa8731758de3deb5b48db/include
cepgen_EX_LIB := CepGen CepGenHepMC2 CepGenHepMC3 CepGenLHAPDF CepGenProcesses CepGenPythia6
cepgen_EX_USE := gsl openblas hepmc hepmc3 lhapdf pythia6

ALL_TOOLS      += classlib
classlib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/classlib/3.1.3-52a2fe4d34a5f733d734f75ee05cb886/include
classlib_EX_LIB := classlib
classlib_EX_USE := pcre root_cxxdefaults
classlib_EX_FLAGS_CPPDEFINES  := -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS

ALL_TOOLS      += clhepheader
clhepheader_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/clhep/2.4.5.1-a18b92988540ef0a928de6c03833f73b/include
clhepheader_EX_USE := root_cxxdefaults
clhepheader_EX_FLAGS_ROOTCLING_ARGS  := -moduleMapFile=$(CLHEP_BASE)/include/module.modulemap
clhepheader_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += clhep
clhep_EX_LIB := CLHEP
clhep_EX_USE := clhepheader

ALL_TOOLS      += cmsswdata
cmsswdata_EX_FLAGS_CMSSW_DATA_PACKAGE  := Alignment/OfflineValidation=V00-03-00 CalibCalorimetry/CaloMiscalibTools=V01-00-00 CalibCalorimetry/EcalTrivialCondModules=V00-03-00 CalibPPS/ESProducers=V01-04-00 CalibTracker/SiPixelESProducers=V02-03-00 CalibTracker/SiStripDCS=V01-01-00 Calibration/Tools=V01-00-00 CondFormats/JetMETObjects=V01-00-03 CondTools/SiPhase2Tracker=V00-02-00 CondTools/SiStrip=V00-02-00 Configuration/Generator=V01-08-00 DQM/DTMonitorClient=V00-01-00 DQM/EcalMonitorClient=V00-03-00 DQM/Integration=V00-02-00 DQM/PhysicsHWW=V01-00-00 DQM/SiStripMonitorClient=V01-01-00 DataFormats/PatCandidates=V01-01-00 DetectorDescription/Schema=V02-03-00 EgammaAnalysis/ElectronTools=V00-03-01 EventFilter/L1TRawToDigi=V01-00-00 FWCore/Modules=V00-01-00 FastSimulation/MaterialEffects=V05-00-00 FastSimulation/TrackingRecHitProducer=V01-00-03 Fireworks/Geometry=V07-06-00 GeneratorInterface/EvtGenInterface=V02-09-00 GeneratorInterface/ReggeGribovPartonMCInterface=V00-00-02 Geometry/DTGeometryBuilder=V00-01-00 Geometry/TestReference=V00-10-00 HLTrigger/JetMET=V01-00-00 HeterogeneousCore/SonicTriton=V00-01-00 IOPool/Input=V00-01-00 L1Trigger/CSCTriggerPrimitives=V00-12-00 L1Trigger/DTTriggerPhase2=V00-02-00 L1Trigger/L1TCalorimeter=V01-01-00 L1Trigger/L1TGlobal=V00-04-00 L1Trigger/L1THGCal=V01-07-00 L1Trigger/L1TMuon=V01-08-00 L1Trigger/Phase2L1ParticleFlow=V00-04-00 L1Trigger/RPCTrigger=V00-15-00 L1Trigger/TrackFindingTMTT=V00-02-00 L1Trigger/TrackFindingTracklet=V00-03-00 L1Trigger/TrackTrigger=V00-02-00 L1TriggerConfig/L1TConfigProducers=V00-01-00 MagneticField/Engine=V00-01-00 MagneticField/Interpolation=V01-02-00 PhysicsTools/NanoAOD=V01-03-00 PhysicsTools/PatUtils=V00-05-00 RecoBTag/CTagging=V01-00-03 RecoBTag/Combined=V01-17-00 RecoBTag/SecondaryVertex=V02-00-04 RecoBTag/SoftLepton=V01-00-01 RecoCTPPS/TotemRPLocal=V00-02-00 RecoEcal/EgammaClusterProducers=V00-02-00 RecoEgamma/EgammaPhotonProducers=V00-01-00 RecoEgamma/ElectronIdentification=V01-12-00 RecoEgamma/PhotonIdentification=V01-06-00 RecoHGCal/TICL=V00-03-00 RecoHI/HiJetAlgos=V01-00-01 RecoJets/JetProducers=V05-14-00 RecoLocalCalo/EcalDeadChannelRecoveryAlgos=V01-01-00 RecoMET/METPUSubtraction=V01-02-00 RecoMTD/TimingIDTools=V00-01-00 RecoMuon/MuonIdentification=V01-14-00 RecoMuon/TrackerSeedGenerator=V00-05-00 RecoParticleFlow/PFBlockProducer=V02-04-02 RecoParticleFlow/PFProducer=V16-02-00 RecoParticleFlow/PFTracking=V13-01-00 RecoTauTag/TrainingFiles=V00-07-00 RecoTracker/FinalTrackSelectors=V01-04-00 RecoTracker/MkFit=V00-12-00 RecoTracker/TkSeedGenerator=V00-02-00 SLHCUpgradeSimulations/Geometry=V01-00-10 SimG4CMS/Calo=V03-04-00 SimG4CMS/Forward=V02-04-00 SimG4CMS/HGCalTestBeam=V01-00-00 SimPPS/PPSPixelDigiProducer=V00-00-02 SimTracker/SiStripDigitizer=V01-01-00 SimTransport/HectorProducer=V01-00-01 SimTransport/PPSProtonTransport=V00-02-00 SimTransport/TotemRPProtonTransportParametrization=V00-01-00 Validation/Geometry=V00-07-00 Validation/HGCalValidation=V00-05-00

ALL_TOOLS      += coral
ALL_SCRAM_PROJECTS += coral
coral_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/coral/CORAL_2_3_21-35c87e33fa3474474e315a875fff39ef/include/LCG
coral_EX_USE := root_cxxdefaults
coral_ORDER := 98000

ALL_TOOLS      += correctionlib
correctionlib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-correctionlib/2.1.0-48c0dd2db01e7bf30cbe3a7382a656ce/lib/python3.9/site-packages/correctionlib/include
correctionlib_EX_LIB := correctionlib

ALL_TOOLS      += cppunit
cppunit_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cppunit/1.15.x-af807824c6eacdb03fd692e08bfbc1c4/include
cppunit_EX_LIB := cppunit
cppunit_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += cpu_features
cpu_features_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cpu_features/0.7.0-9eb21eea97f43a6d1e119932a8c2df90/include
cpu_features_EX_LIB := cpu_features

ALL_TOOLS      += csctrackfinderemulation
csctrackfinderemulation_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/CSCTrackFinderEmulation/1.2-f743171ba99b52f2b8933c74be392a12/include
csctrackfinderemulation_EX_LIB := CSCTrackFinderEmulation

ALL_TOOLS      += cuda-compatible-runtime
cuda-compatible-runtime_EX_USE := cuda cuda-stubs

ALL_TOOLS      += cuda-cublas
cuda-cublas_EX_LIB := cublas cublasLt
cuda-cublas_EX_USE := cuda

ALL_TOOLS      += cuda-cufft
cuda-cufft_EX_LIB := cufft cufftw
cuda-cufft_EX_USE := cuda

ALL_TOOLS      += cuda-curand
cuda-curand_EX_LIB := curand
cuda-curand_EX_USE := cuda

ALL_TOOLS      += cuda-cusolver
cuda-cusolver_EX_LIB := cusolver cusolverMg
cuda-cusolver_EX_USE := cuda

ALL_TOOLS      += cuda-cusparse
cuda-cusparse_EX_LIB := cusparse
cuda-cusparse_EX_USE := cuda

ALL_TOOLS      += cuda-gcc-support

ALL_TOOLS      += cuda
ALL_LIB_TYPES += CUDA_LIB
cuda_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cuda/11.5.2-91421d5831fe9d02c7a55c817b17889e/include
cuda_EX_LIB := cudart cudadevrt nvToolsExt
cuda_EX_CUDA_LIB := cudadevrt
cuda_EX_USE := cuda-stubs
cuda_EX_FLAGS_CUDA_FLAGS  := -std=c++17 -O3 --generate-line-info --source-in-ptx --display-error-number --expt-relaxed-constexpr --extended-lambda -gencode arch=compute_60,code=[sm_60,compute_60] -gencode arch=compute_70,code=[sm_70,compute_70] -gencode arch=compute_75,code=[sm_75,compute_75] -Wno-deprecated-gpu-targets -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored --cudart shared
cuda_EX_FLAGS_CUDA_HOST_CXXFLAGS  := -std=c++17
cuda_EX_FLAGS_REM_CUDA_HOST_CXXFLAGS  := -std=% %potentially-evaluated-expression

ALL_TOOLS      += cuda-npp
cuda-npp_EX_LIB := nppial nppicc nppidei nppif nppig nppim nppist nppisu nppitc npps nppc
cuda-npp_EX_USE := cuda

ALL_TOOLS      += cuda-nvgraph
cuda-nvgraph_EX_LIB := nvgraph
cuda-nvgraph_EX_USE := cuda

ALL_TOOLS      += cuda-nvjpeg
cuda-nvjpeg_EX_LIB := nvjpeg
cuda-nvjpeg_EX_USE := cuda

ALL_TOOLS      += cuda-nvml
cuda-nvml_EX_LIB := nvidia-ml
cuda-nvml_EX_USE := cuda-stubs

ALL_TOOLS      += cuda-nvrtc
cuda-nvrtc_EX_LIB := nvrtc
cuda-nvrtc_EX_USE := cuda

ALL_TOOLS      += cuda-stubs
cuda-stubs_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cuda/11.5.2-91421d5831fe9d02c7a55c817b17889e/include
cuda-stubs_EX_LIB := cuda
cuda-stubs_EX_LIBDIR := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cuda/11.5.2-91421d5831fe9d02c7a55c817b17889e/lib64/stubs
cuda-stubs_EX_FLAGS_SKIP_TOOL_SYMLINKS  := 1

ALL_TOOLS      += cudnn
cudnn_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cudnn/8.3.3.40-1fb182f10266436d1bce396b5f64d822/include
cudnn_EX_LIB := cudnn
cudnn_EX_USE := cuda

ALL_TOOLS      += cupti
cupti_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cuda/11.5.2-91421d5831fe9d02c7a55c817b17889e/include
cupti_EX_LIB := cupti

ALL_TOOLS      += curl
curl_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/curl/7.79.0-5e48e0bf013ba13376a33ec8da72dabc/include
curl_EX_LIB := curl
curl_EX_USE := root_cxxdefaults

ALL_TOOLS      += dablooms
dablooms_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/dablooms/0.9.1-dfdd93a51a36797f9aeb5dfcea5af9ab/include
dablooms_EX_LIB := dablooms

ALL_TOOLS      += das_client

ALL_TOOLS      += davix
davix_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/davix/0.8.1-5a38d8a24145ab7f5b977e62088ce313/include/davix
davix_EX_LIB := davix
davix_EX_USE := boost_system libxml2

ALL_TOOLS      += db6
db6_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/db6/6.2.32-16f729842d558005f3c577e0de97893d/include
db6_EX_LIB := db

ALL_TOOLS      += dcap
dcap_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/dcap/2.47.12-84a42bcf47b5916931bf051c5d3eea7d/include
dcap_EX_LIB := dcap
dcap_EX_USE := root_cxxdefaults

ALL_TOOLS      += dd4hep-core
dd4hep-core_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/dd4hep/v01-23x-e2cf61cff4a911ed9506a3153f3788a1/include
dd4hep-core_EX_LIB := DDCore DDParsers
dd4hep-core_EX_USE := root_cxxdefaults root boost xerces-c clhep
dd4hep-core_EX_FLAGS_CPPDEFINES  := -DDD4HEP_USE_GEANT4_UNITS=1
dd4hep-core_EX_FLAGS_LISTCOMPONENTS  := $(DD4HEP_CORE_BASE)/bin/listcomponents_dd4hep

ALL_TOOLS      += dd4hep-geant4
dd4hep-geant4_EX_LIB := DDG4-static
dd4hep-geant4_EX_USE := geant4core dd4hep-core

ALL_TOOLS      += dd4hep
dd4hep_EX_LIB := DDAlign DDCond
dd4hep_EX_USE := dd4hep-core

ALL_TOOLS      += dip_interface
dip_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/dip/8693f00cc422b4a15858fcd84249acaeb07b6316-08f74240b433fac5e3b5bcaf3edb45a0/include

ALL_TOOLS      += dip
dip_EX_LIB := dip
dip_EX_USE := dip-platform-dependent log4cplus

ALL_TOOLS      += dip-platform-dependent
dip-platform-dependent_EX_LIB := platform-dependent
dip-platform-dependent_EX_USE := dip_interface

ALL_TOOLS      += eigen
eigen_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/eigen/82dd3710dac619448f50331c1d6a35da673f764a-9ac4aed18ac60d0189693c592862694d/include/eigen3
eigen_EX_FLAGS_CPPDEFINES  := -DEIGEN_DONT_PARALLELIZE
eigen_EX_FLAGS_CUDA_FLAGS  := --diag-suppress 20014

ALL_TOOLS      += evtgen
evtgen_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/evtgen/2.0.0-e462ba602a2435b89fa768ce18660219/include
evtgen_EX_LIB := EvtGen EvtGenExternal
evtgen_EX_USE := hepmc pythia8 tauolapp photospp

ALL_TOOLS      += expat
expat_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/expat/2.1.0-5f6457b4c04e97afec6079bd7d2db998/include
expat_EX_LIB := expat
expat_EX_USE := root_cxxdefaults

ALL_TOOLS      += fastjet-contrib-archive
fastjet-contrib-archive_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fastjet-contrib/1.044-9f8288cf11100a284a45f857da9ac46e/include
fastjet-contrib-archive_EX_LIB := EnergyCorrelator GenericSubtractor JetCleanser JetFFMoments Nsubjettiness ScJet SubjetCounting VariableR

ALL_TOOLS      += fastjet-contrib
fastjet-contrib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fastjet-contrib/1.044-9f8288cf11100a284a45f857da9ac46e/include
fastjet-contrib_EX_LIB := fastjetcontribfragile

ALL_TOOLS      += fastjet
fastjet_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fastjet/3.4.0-07dc8fd99576cb469cf9a5928b5c31be/include
fastjet_EX_LIB := fastjetplugins fastjettools siscone siscone_spherical fastjet
fastjet_EX_USE := root_cxxdefaults

ALL_TOOLS      += fftjet
fftjet_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fftjet/1.5.0-1cd7b7c71b42be65b840b68281156a1e/include
fftjet_EX_LIB := fftjet
fftjet_EX_USE := root_cxxdefaults

ALL_TOOLS      += fftw3
fftw3_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fftw3/3.3.8-5f403c3dc6c4147769a1f10565ec7e26/include
fftw3_EX_LIB := fftw3
fftw3_EX_USE := root_cxxdefaults

ALL_TOOLS      += flatbuffers
flatbuffers_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/flatbuffers/1.12.0-11947d9da619f98dfdaa5a32d32902d7/include
flatbuffers_EX_LIB := flatbuffers

ALL_TOOLS      += fmt
fmt_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/fmt/8.0.1-43b841663c2a0d6622910a1ad66d228d/include
fmt_EX_LIB := fmt

ALL_TOOLS      += freetype
freetype_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/freetype/2.10.0-1c5d3293d887fb2fec3e2e99c7ba52fe/include
freetype_EX_LIB := freetype-cms
freetype_EX_USE := root_cxxdefaults

ALL_TOOLS      += frontier_client
frontier_client_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/frontier_client/2.9.1-7eb5a7e2f18435d9d07e1bb911c3b2fd/include
frontier_client_EX_LIB := frontier_client
frontier_client_EX_USE := root_cxxdefaults zlib expat

ALL_TOOLS      += g4hepemcore
g4hepemcore_EX_LIB := g4HepEmData g4HepEmInit g4HepEmRun g4HepEm
g4hepemcore_EX_USE := g4hepem_interface

ALL_TOOLS      += g4hepem_interface
g4hepem_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/g4hepem/20221014-7267d0044c03f73d979697ff84064193/include/G4HepEm

ALL_TOOLS      += g4hepemstatic
g4hepemstatic_EX_LIB := g4hepem-static
g4hepemstatic_EX_USE := g4hepem_interface

ALL_TOOLS      += gbl
gbl_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gbl/V02-04-01-bf18851db250ca8ca2851f6bef12b0b3/include
gbl_EX_LIB := GBL
gbl_EX_USE := eigen

ALL_TOOLS      += gcc-atomic
gcc-atomic_EX_LIB := atomic

ALL_TOOLS      += gcc-ccompiler
gcc-ccompiler_EX_FLAGS_CFLAGS  := -O2 -pthread
gcc-ccompiler_EX_FLAGS_CSHAREDOBJECTFLAGS  := -fPIC

ALL_TOOLS      += gcc-cxxcompiler
gcc-cxxcompiler_EX_FLAGS_CPPDEFINES  := -DGNU_GCC -D_GNU_SOURCE
gcc-cxxcompiler_EX_FLAGS_CXXFLAGS  := -O2 -pthread -pipe -Werror=main -Werror=pointer-arith -Werror=overlength-strings -Wno-vla -Werror=overflow -std=c++1z -ftree-vectorize -Werror=array-bounds -Werror=format-contains-nul -Werror=type-limits -fvisibility-inlines-hidden -fno-math-errno --param vect-max-version-for-alias-checks=50 -Xassembler --compress-debug-sections -fuse-ld=bfd -msse3 -felide-constructors -fmessage-length=0 -Wall -Wno-non-template-friend -Wno-long-long -Wreturn-type -Wextra -Wpessimizing-move -Wclass-memaccess -Wno-cast-function-type -Wno-unused-but-set-parameter -Wno-ignored-qualifiers -Wno-deprecated-copy -Wno-unused-parameter -Wunused -Wparentheses -Wno-deprecated -Werror=return-type -Werror=missing-braces -Werror=unused-value -Werror=unused-label -Werror=address -Werror=format -Werror=sign-compare -Werror=write-strings -Werror=delete-non-virtual-dtor -Werror=strict-aliasing -Werror=narrowing -Werror=unused-but-set-variable -Werror=reorder -Werror=unused-variable -Werror=conversion-null -Werror=return-local-addr -Wnon-virtual-dtor -Werror=switch -fdiagnostics-show-option -Wno-unused-local-typedefs -Wno-attributes -Wno-psabi
gcc-cxxcompiler_EX_FLAGS_CXXSHAREDFLAGS  := -shared -Wl,-E
gcc-cxxcompiler_EX_FLAGS_CXXSHAREDOBJECTFLAGS  := -fPIC
gcc-cxxcompiler_EX_FLAGS_LDFLAGS  := -Wl,-E -Wl,--hash-style=gnu
gcc-cxxcompiler_EX_FLAGS_LD_UNIT  := -r -z muldefs
gcc-cxxcompiler_EX_FLAGS_LTO_FLAGS  := -flto -fipa-icf -flto-odr-type-merging -fno-fat-lto-objects -Wodr

ALL_TOOLS      += gcc-f77compiler
gcc-f77compiler_EX_LIB := gfortran m
gcc-f77compiler_EX_FLAGS_FFLAGS  := -fno-second-underscore -Wunused -Wuninitialized -O2 -cpp -std=legacy
gcc-f77compiler_EX_FLAGS_FOPTIMISEDFLAGS  := -O2
gcc-f77compiler_EX_FLAGS_FSHAREDOBJECTFLAGS  := -fPIC

ALL_TOOLS      += gcc-plugin
gcc-plugin_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gcc/11.2.1-f9b9dfdd886f71cd63f5538223d8f161/bin/../lib/gcc/x86_64-redhat-linux-gnu/11.2.1/plugin/include
gcc-plugin_EX_LIB := cc1plugin cp1plugin

ALL_TOOLS      += gdb

ALL_TOOLS      += gdbm
gdbm_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gdbm/1.10-94fd72446cd6c73834b291fb1d1c6f46/include
gdbm_EX_LIB := gdbm
gdbm_EX_USE := root_cxxdefaults

ALL_TOOLS      += gdrcopy
gdrcopy_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gdrcopy/2.3-dabcd07891cf3cc84481d5fa41e5dcda/include
gdrcopy_EX_LIB := gdrapi

ALL_TOOLS      += geant4core
geant4core_EX_LIB := G4digits_hits G4error_propagation G4event G4geometry G4global G4graphics_reps G4intercoms G4interfaces G4materials G4parmodels G4particles G4persistency G4physicslists G4processes G4readout G4run G4tracking G4track G4analysis
geant4core_EX_USE := geant4_interface
geant4core_EX_FLAGS_CPPDEFINES  := -DGNU_GCC -DG4V9

ALL_TOOLS      += geant4data

ALL_TOOLS      += geant4_interface
geant4_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/geant4/10.7.2-0e62b7ff14b69147f66862fcd6cd565b/include/Geant4 /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/geant4/10.7.2-0e62b7ff14b69147f66862fcd6cd565b/include
geant4_interface_EX_USE := clhep vecgeom zlib expat xerces-c root_cxxdefaults
geant4_interface_EX_FLAGS_CPPDEFINES  := -DGNU_GCC -DG4V9
geant4_interface_EX_FLAGS_CXXFLAGS  := -ftls-model=global-dynamic -pthread -DG4GEOM_USE_USOLIDS

ALL_TOOLS      += geant4
geant4_EX_USE := geant4core geant4vis

ALL_TOOLS      += geant4-parfullcms

ALL_TOOLS      += geant4static
geant4static_EX_LIB := geant4-static
geant4static_EX_USE := geant4_interface

ALL_TOOLS      += geant4vis
geant4vis_EX_LIB := G4FR G4modeling G4RayTracer G4Tree G4visHepRep G4vis_management G4visXXX G4VRML G4GMocren
geant4vis_EX_USE := geant4core

ALL_TOOLS      += giflib
giflib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/giflib/5.2.0-a87f472cd9febee3811c0f3e8b1546ec/include
giflib_EX_LIB := gif
giflib_EX_USE := root_cxxdefaults

ALL_TOOLS      += git

ALL_TOOLS      += glimpse

ALL_TOOLS      += gmake

ALL_TOOLS      += gmp
gmp_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gmp-static/6.2.1-71a50667e82274ff697be2437a3fd6fa/include
gmp_EX_LIB := gmp
gmp_EX_USE := mpfr

ALL_TOOLS      += gnuplot

ALL_TOOLS      += google-benchmark-main
google-benchmark-main_EX_LIB := benchmark_main
google-benchmark-main_EX_USE := google-benchmark

ALL_TOOLS      += google-benchmark
google-benchmark_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/google-benchmark/1.4.x-c77e5f08281e2b8b61ac86d962c74824/include
google-benchmark_EX_LIB := benchmark
google-benchmark_EX_USE := sockets

ALL_TOOLS      += gosamcontrib
gosamcontrib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gosamcontrib/2.0-20150803-001c255025c150fdbe081040c1f49536/include

ALL_TOOLS      += gosam

ALL_TOOLS      += gperf
gperf_EX_LIB := profiler
gperf_EX_USE := libunwind

ALL_TOOLS      += grpc
grpc_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/grpc/1.35.0-4ae66c6d769e6f806f3b64de826cbe10/include
grpc_EX_LIB := grpc grpc++ grpc++_reflection
grpc_EX_USE := protobuf openssl pcre abseil-cpp c-ares re2
grpc_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += gsl
gsl_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/gsl/2.6-fcf47bcbedd800ca8386c7e2920fa474/include
gsl_EX_LIB := gsl
gsl_EX_USE := openblas root_cxxdefaults

ALL_TOOLS      += hdf5
hdf5_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hdf5/1.10.6-9692e139de3c800f61c237bec8aca294/include
hdf5_EX_LIB := hdf5 hdf5_hl
hdf5_EX_USE := openmpi

ALL_TOOLS      += hector
hector_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hector/1.3.4_patch1-d4676c710c031ed4217494ab52c63da8/include
hector_EX_LIB := Hector
hector_EX_USE := root_cxxdefaults

ALL_TOOLS      += hepmc3
hepmc3_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hepmc3/3.2.5-a8c3d3d81a1b670617f2caa48c1801ef/include
hepmc3_EX_LIB := HepMC3 HepMC3search

ALL_TOOLS      += hepmc_headers
hepmc_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hepmc/2.06.10-8d052e96205063f65fddfc4cbf40f506/include
hepmc_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += hepmc
hepmc_EX_LIB := HepMCfio HepMC
hepmc_EX_USE := hepmc_headers

ALL_TOOLS      += heppdt
heppdt_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/heppdt/3.04.01-720968ba28a92b392511de93b515db3b/include
heppdt_EX_LIB := HepPDT HepPID
heppdt_EX_USE := root_cxxdefaults

ALL_TOOLS      += herwig7
herwig7_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/herwig7/7.2.2-c154d39c8c55adbe3ab379346e7cb558/include
herwig7_EX_LIB := HerwigAPI
herwig7_EX_USE := root_cxxdefaults lhapdf thepeg madgraph5amcatnlo openloops

ALL_TOOLS      += highfive
highfive_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/highfive/2.3.1-3e6d2d063cbb3255a45420653e2cc6bd/include
highfive_EX_USE := boost hdf5

ALL_TOOLS      += histfactory
histfactory_EX_LIB := HistFactory
histfactory_EX_USE := roofitcore roofit roostats rootcore roothistmatrix rootgpad rootxml rootfoam

ALL_TOOLS      += hls
hls_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hls/2019.08-fd724004387c2a6770dc3517446d30d9/include
hls_EX_USE := root_cxxdefaults
hls_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += hwloc
hwloc_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hwloc/2.9.0-8abc90483f9c655cad9e31cc4d119b92/include/hwloc
hwloc_EX_LIB := hwloc

ALL_TOOLS      += hydjet2
hydjet2_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/hydjet2/2.4.3-2ee8d225ebc9db4420b10c71404db389/include
hydjet2_EX_LIB := hydjet2
hydjet2_EX_USE := pyquen pythia6 lhapdf root
hydjet2_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += hydjet
hydjet_EX_LIB := hydjet
hydjet_EX_USE := pyquen pythia6 lhapdf

ALL_TOOLS      += igprof

ALL_TOOLS      += intel-license

ALL_TOOLS      += ittnotify
ittnotify_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/ittnotify/16.06.18-55dda586bd53ea172280a0e952fdff95/include
ittnotify_EX_LIB := ittnotify

ALL_TOOLS      += iwyu-cxxcompiler
iwyu-cxxcompiler_EX_USE := llvm-cxxcompiler

ALL_TOOLS      += jemalloc
jemalloc_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/jemalloc/5.3.0-0e4eec620afb3820f8b34db85bac3dc0/include
jemalloc_EX_LIB := jemalloc
jemalloc_EX_USE := root_cxxdefaults

ALL_TOOLS      += json
json_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/json/3.10.2-a6d86565b09ec3d0e02bf7b52c31bbfc/include

ALL_TOOLS      += ktjet
ktjet_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/ktjet/1.06-2553ebaf812ddb9f45fb09136cae448d/include
ktjet_EX_LIB := KtEvent
ktjet_EX_USE := root_cxxdefaults
ktjet_EX_FLAGS_CPPDEFINES  := -DKTDOUBLEPRECISION

ALL_TOOLS      += lcov

ALL_TOOLS      += lhapdf
lhapdf_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/lhapdf/6.4.0-576e285c9e352d7f1106d0eaafdd3acb/include
lhapdf_EX_LIB := LHAPDF
lhapdf_EX_USE := pythia6 root_cxxdefaults

ALL_TOOLS      += libffi
libffi_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libffi/3.4.2-8d81640cbf7bb880677ab7baa97b0341/include
libffi_EX_LIB := ffi

ALL_TOOLS      += libibumad
libibumad_EX_LIB := ibumad
libibumad_EX_USE := rdma-core

ALL_TOOLS      += libibverbs
libibverbs_EX_LIB := ibverbs
libibverbs_EX_USE := rdma-core

ALL_TOOLS      += libjpeg-turbo
libjpeg-turbo_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libjpeg-turbo/2.0.2-470e29b709f28ed0f6cede9a044313f1/include
libjpeg-turbo_EX_LIB := jpeg turbojpeg
libjpeg-turbo_EX_USE := root_cxxdefaults

ALL_TOOLS      += libpciaccess
libpciaccess_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libpciaccess/0.16-9394bc1511efecfa229a0ff28678a57e/include
libpciaccess_EX_LIB := libpciaccess

ALL_TOOLS      += libpng
libpng_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libpng/1.6.37-fa0ad33da9354d525ff104b6029da553/include
libpng_EX_LIB := png
libpng_EX_USE := root_cxxdefaults zlib

ALL_TOOLS      += librdmacm
librdmacm_EX_LIB := rdmacm
librdmacm_EX_USE := rdma-core libibverbs

ALL_TOOLS      += libtiff
libtiff_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libtiff/4.0.10-74c83ecf4126bf7e57595c36393df2f2/include
libtiff_EX_LIB := tiff
libtiff_EX_USE := root_cxxdefaults libjpeg-turbo zlib

ALL_TOOLS      += libungif
libungif_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libungif/4.1.4-6b4f3e28660c3a9d4e75854eb8da0012/include
libungif_EX_LIB := ungif
libungif_EX_USE := root_cxxdefaults zlib

ALL_TOOLS      += libunwind
libunwind_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libunwind/1.6.2-master-da573a24e8497382f776d9f95ef48336/include
libunwind_EX_LIB := unwind
libunwind_EX_USE := root_cxxdefaults

ALL_TOOLS      += libuuid
libuuid_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libuuid/2.34-0451b31e1b9a58c6aeefab41c18eea34/include
libuuid_EX_LIB := uuid
libuuid_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += libxml2
libxml2_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libxml2/2.9.10-8276cdd7c759f220df9a495a1e710246/include/libxml2
libxml2_EX_LIB := xml2
libxml2_EX_USE := root_cxxdefaults

ALL_TOOLS      += libxslt
libxslt_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libxslt/1.1.28-2147ef893b903a2b686c904d09f62ed1/include/libxslt
libxslt_EX_LIB := xslt

ALL_TOOLS      += libzmq
libzmq_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/libzmq/4.3.4-2a366197a57103a16875b0391d9093e6/include
libzmq_EX_LIB := zmq

ALL_TOOLS      += llvm-analyzer-ccompiler
llvm-analyzer-ccompiler_EX_USE := llvm-ccompiler

ALL_TOOLS      += llvm-analyzer-cxxcompiler
llvm-analyzer-cxxcompiler_EX_USE := llvm-cxxcompiler

ALL_TOOLS      += llvm-ccompiler
llvm-ccompiler_EX_USE := gcc-ccompiler

ALL_TOOLS      += llvm-cxxcompiler
llvm-cxxcompiler_EX_USE := gcc-cxxcompiler
llvm-cxxcompiler_EX_FLAGS_CXXFLAGS  := -Wno-c99-extensions -Wno-c++11-narrowing -D__STRICT_ANSI__ -Wno-unused-private-field -Wno-unknown-pragmas -Wno-unused-command-line-argument -Wno-unknown-warning-option -ftemplate-depth=512 -Wno-error=potentially-evaluated-expression -Wno-tautological-type-limit-compare -fsized-deallocation
llvm-cxxcompiler_EX_FLAGS_REM_CXXFLAGS  := -Wno-non-template-friend -Werror=format-contains-nul -Werror=maybe-uninitialized -Werror=unused-but-set-variable -Werror=return-local-addr -fipa-pta -frounding-math -mrecip -fno-crossjumping -fno-aggressive-loop-optimizations -funroll-all-loops
llvm-cxxcompiler_EX_FLAGS_REM_LTO_FLAGS  := -fipa-icf -flto-odr-type-merging

ALL_TOOLS      += llvm-f77compiler
llvm-f77compiler_EX_USE := gcc-f77compiler

ALL_TOOLS      += llvm
llvm_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/llvm/12.0.1-34ef258eea96cb5b111081adc51a827d/include
llvm_EX_LIB := clang
llvm_EX_FLAGS_CXXFLAGS  := -D_DEBUG -D_GNU_SOURCE -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -O3 -fomit-frame-pointer -fPIC -Wno-enum-compare -Wno-strict-aliasing -fno-rtti
llvm_EX_FLAGS_LDFLAGS  := -Wl,-undefined -Wl,suppress

ALL_TOOLS      += log4cplus
log4cplus_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/log4cplus/2.0.7-ae4d03ad227f093b8194f6d6c0e8322b/include
log4cplus_EX_LIB := log4cplusS

ALL_TOOLS      += lwtnn
lwtnn_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/lwtnn/2.13-c9325f37819a21f5dae660bb6c377604/include
lwtnn_EX_LIB := lwtnn
lwtnn_EX_USE := root_cxxdefaults eigen boost_system

ALL_TOOLS      += lz4
lz4_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/lz4/1.9.2-031da253be076d002e4d6af36bc64212/include
lz4_EX_LIB := lz4
lz4_EX_USE := root_cxxdefaults

ALL_TOOLS      += madgraph5amcatnlo
madgraph5amcatnlo_EX_USE := root_cxxdefaults gosamcontrib

ALL_TOOLS      += mctester
mctester_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/mctester/1.25.1-5e01e292375ccc378e0b5952479d64c9/include
mctester_EX_LIB := HEPEvent HepMCEvent MCTester
mctester_EX_USE := root_cxxdefaults root hepmc

ALL_TOOLS      += md5
md5_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/md5/1.0.0-e68283f2de2e2e709a0db99db3b53205/include
md5_EX_LIB := cms-md5

ALL_TOOLS      += meschach
meschach_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/meschach/1.2.pCMS1-a6c940d49300e01334c28ef7c2460c02/include
meschach_EX_LIB := meschach
meschach_EX_USE := root_cxxdefaults

ALL_TOOLS      += millepede
millepede_EX_USE := sockets pcre zlib

ALL_TOOLS      += mpfr
mpfr_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/mpfr-static/4.1.0-7b42d3d4fda2a86ff65c7f5831de4ca9/include
mpfr_EX_LIB := mpfr

ALL_TOOLS      += mxnet-predict
mxnet-predict_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/mxnet-predict/1.6.0-0bf665ed4554f7ec04f3af8e5fc5103e/include
mxnet-predict_EX_LIB := mxnet
mxnet-predict_EX_USE := openblas
mxnet-predict_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += numactl
numactl_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/numactl/2.0.14-9f463f2c31f4a09ab8785f9a5f851b80/include
numactl_EX_LIB := numa

ALL_TOOLS      += numpy-c-api
numpy-c-api_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-numpy/1.22.4-1e1240160543ee39175a125e9a5e9cae/c-api/core/include
numpy-c-api_EX_LIB := npymath
numpy-c-api_EX_USE := root_cxxdefaults

ALL_TOOLS      += nvperf
nvperf_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/cuda/11.5.2-91421d5831fe9d02c7a55c817b17889e/include
nvperf_EX_LIB := nvperf_host nvperf_target

ALL_TOOLS      += ofast-flag
ofast-flag_EX_FLAGS_CXXFLAGS  := -Ofast -fno-reciprocal-math -mrecip=none
ofast-flag_EX_FLAGS_NO_RECURSIVE_EXPORT  := 1

ALL_TOOLS      += onnxruntime
onnxruntime_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/onnxruntime/1.10.0-d4ea5c99a91518780757077276506fc0/include
onnxruntime_EX_LIB := onnxruntime
onnxruntime_EX_USE := protobuf cuda cudnn re2

ALL_TOOLS      += openblas
openblas_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/OpenBLAS/0.3.15-26c67b8b638762cfd2e2bcfc936e3ec7/include
openblas_EX_LIB := openblas

ALL_TOOLS      += opencv
opencv_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/opencv/4.5.5-def22b0c43b52ecc33443da0ced5a956/include
opencv_EX_LIB := opencv_core
opencv_EX_USE := libpng libjpeg-turbo zlib eigen openblas

ALL_TOOLS      += opengl
opengl_EX_LIB := GL GLU
opengl_EX_USE := x11

ALL_TOOLS      += openldap
openldap_EX_USE := db6

ALL_TOOLS      += openloops

ALL_TOOLS      += openmpi
openmpi_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/openmpi/4.1.5-2418ebaa5670c4f86488177011873558/include
openmpi_EX_LIB := mpi mpi_cxx

ALL_TOOLS      += openssl
openssl_EX_LIB := ssl crypto

ALL_TOOLS      += oracle
oracle_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/oracle/19.11.0.0.0dbru-092ba337a70a02c4734616a84f842897/include
oracle_EX_LIB := clntsh
oracle_EX_USE := root_cxxdefaults sockets

ALL_TOOLS      += oracleocci
oracleocci_EX_LIB := occi
oracleocci_EX_USE := oracle

ALL_TOOLS      += pacparser
pacparser_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pacparser/1.4.0-417b085cd92df49d866b76828d34095a/include
pacparser_EX_LIB := pacparser
pacparser_EX_USE := root_cxxdefaults

ALL_TOOLS      += pcre2
pcre2_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pcre2/10.36-88c59c64ebe54b36ea4626d5827d6026/include
pcre2_EX_LIB := pcre2
pcre2_EX_USE := root_cxxdefaults zlib bz2lib

ALL_TOOLS      += pcre
pcre_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pcre/8.43-5dcc901acc02f624b22dd9840b2357e8/include
pcre_EX_LIB := pcre
pcre_EX_USE := root_cxxdefaults zlib bz2lib

ALL_TOOLS      += photospp
photospp_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/photospp/3.64-4a9cdf1e1431c4ba7fa05e053796f38c/include
photospp_EX_LIB := Photospp PhotosppHepMC PhotosppHEPEVT
photospp_EX_USE := hepmc

ALL_TOOLS      += professor2
professor2_EX_USE := py3-numpy root yoda eigen

ALL_TOOLS      += protobuf
protobuf_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/protobuf/3.15.1-6884a7793e9b97ba2d2b356a32a0ab17/include
protobuf_EX_LIB := protobuf
protobuf_EX_USE := root_cxxdefaults

ALL_TOOLS      += py3-absl-py

ALL_TOOLS      += py3-anyio

ALL_TOOLS      += py3-appdirs

ALL_TOOLS      += py3-argon2-cffi-bindings

ALL_TOOLS      += py3-argon2-cffi

ALL_TOOLS      += py3-argparse

ALL_TOOLS      += py3-asn1crypto

ALL_TOOLS      += py3-astor

ALL_TOOLS      += py3-astroid

ALL_TOOLS      += py3-asttokens

ALL_TOOLS      += py3-astunparse

ALL_TOOLS      += py3-atomicwrites

ALL_TOOLS      += py3-attrs

ALL_TOOLS      += py3-autopep8

ALL_TOOLS      += py3-avro

ALL_TOOLS      += py3-awkward

ALL_TOOLS      += py3-backcall

ALL_TOOLS      += py3-backports-entry-points-selectable

ALL_TOOLS      += py3-beautifulsoup4

ALL_TOOLS      += py3-beniget

ALL_TOOLS      += py3-bleach

ALL_TOOLS      += py3-boost-histogram

ALL_TOOLS      += py3-bottleneck

ALL_TOOLS      += py3-cachecontrol

ALL_TOOLS      += py3-cachetools

ALL_TOOLS      += py3-cachy

ALL_TOOLS      += py3-certifi

ALL_TOOLS      += py3-cffi

ALL_TOOLS      += py3-chardet

ALL_TOOLS      += py3-charset-normalizer

ALL_TOOLS      += py3-cleo

ALL_TOOLS      += py3-click

ALL_TOOLS      += py3-clikit

ALL_TOOLS      += py3-cmsml

ALL_TOOLS      += py3-colorama

ALL_TOOLS      += py3-commonmark

ALL_TOOLS      += py3-contextlib2

ALL_TOOLS      += py3-correctionlib

ALL_TOOLS      += py3-coverage

ALL_TOOLS      += py3-cppy

ALL_TOOLS      += py3-crashtest

ALL_TOOLS      += py3-cryptography

ALL_TOOLS      += py3-cx-oracle

ALL_TOOLS      += py3-cycler

ALL_TOOLS      += py3-cython

ALL_TOOLS      += py3-debugpy

ALL_TOOLS      += py3-decorator

ALL_TOOLS      += py3-defusedxml

ALL_TOOLS      += py3-deprecated

ALL_TOOLS      += py3-deprecation

ALL_TOOLS      += py3-dill

ALL_TOOLS      += py3-distlib

ALL_TOOLS      += py3-docopt

ALL_TOOLS      += py3-docutils

ALL_TOOLS      += py3-downhill

ALL_TOOLS      += py3-dulwich

ALL_TOOLS      += py3-dxr

ALL_TOOLS      += py3-editables

ALL_TOOLS      += py3-entrypoints

ALL_TOOLS      += py3-executing

ALL_TOOLS      += py3-fastjsonschema

ALL_TOOLS      += py3-filelock

ALL_TOOLS      += py3-fire

ALL_TOOLS      += py3-flake8

ALL_TOOLS      += py3-flatbuffers

ALL_TOOLS      += py3-flawfinder

ALL_TOOLS      += py3-flit-core

ALL_TOOLS      += py3-flit

ALL_TOOLS      += py3-fonttools

ALL_TOOLS      += py3-funcsigs

ALL_TOOLS      += py3-future

ALL_TOOLS      += py3-gast

ALL_TOOLS      += py3-gitdb

ALL_TOOLS      += py3-gitpython

ALL_TOOLS      += py3-google-auth

ALL_TOOLS      += py3-google-auth-oauthlib

ALL_TOOLS      += py3-google-pasta

ALL_TOOLS      += py3-grpcio

ALL_TOOLS      += py3-grpcio-tools

ALL_TOOLS      += py3-h5py

ALL_TOOLS      += py3-hatchling

ALL_TOOLS      += py3-hepdata-lib

ALL_TOOLS      += py3-hepdata-validator

ALL_TOOLS      += py3-hep_ml

ALL_TOOLS      += py3-histbook

ALL_TOOLS      += py3-hist

ALL_TOOLS      += py3-histogrammar

ALL_TOOLS      += py3-histoprint

ALL_TOOLS      += py3-html5lib

ALL_TOOLS      += py3-idna

ALL_TOOLS      += py3-importlib-metadata

ALL_TOOLS      += py3-importlib-resources

ALL_TOOLS      += py3-iniconfig

ALL_TOOLS      += py3-ipaddress

ALL_TOOLS      += py3-ipykernel

ALL_TOOLS      += py3-ipython_genutils

ALL_TOOLS      += py3-ipython

ALL_TOOLS      += py3-ipywidgets

ALL_TOOLS      += py3-isort

ALL_TOOLS      += py3-jaraco-classes

ALL_TOOLS      += py3-jedi

ALL_TOOLS      += py3-jeepney

ALL_TOOLS      += py3-jinja2

ALL_TOOLS      += py3-joblib

ALL_TOOLS      += py3-jsonpickle

ALL_TOOLS      += py3-jsonschema

ALL_TOOLS      += py3-jupyter-client

ALL_TOOLS      += py3-jupyter-console

ALL_TOOLS      += py3-jupyter-core

ALL_TOOLS      += py3-jupyterlab-pygments

ALL_TOOLS      += py3-jupyterlab-widgets

ALL_TOOLS      += py3-jupyter-packaging

ALL_TOOLS      += py3-keras2onnx

ALL_TOOLS      += py3-keras-applications

ALL_TOOLS      += py3-keras

ALL_TOOLS      += py3-keras-preprocessing

ALL_TOOLS      += py3-keyring

ALL_TOOLS      += py3-kiwisolver

ALL_TOOLS      += py3-law

ALL_TOOLS      += py3-lazy-object-proxy

ALL_TOOLS      += py3-lizard

ALL_TOOLS      += py3-llvmlite

ALL_TOOLS      += py3-lockfile

ALL_TOOLS      += py3-luigi

ALL_TOOLS      += py3-lxml

ALL_TOOLS      += py3-lz4

ALL_TOOLS      += py3-mako

ALL_TOOLS      += py3-markdown

ALL_TOOLS      += py3-markupsafe

ALL_TOOLS      += py3-matplotlib-inline

ALL_TOOLS      += py3-matplotlib

ALL_TOOLS      += py3-mccabe

ALL_TOOLS      += py3-mistune

ALL_TOOLS      += py3-mock

ALL_TOOLS      += py3-more-itertools

ALL_TOOLS      += py3-mplhep-data

ALL_TOOLS      += py3-mplhep

ALL_TOOLS      += py3-mpmath

ALL_TOOLS      += py3-msgpack

ALL_TOOLS      += py3-nbclient

ALL_TOOLS      += py3-nbconvert

ALL_TOOLS      += py3-nbformat

ALL_TOOLS      += py3-nest-asyncio

ALL_TOOLS      += py3-networkx

ALL_TOOLS      += py3-notebook

ALL_TOOLS      += py3-numba

ALL_TOOLS      += py3-numexpr

ALL_TOOLS      += py3-numpy

ALL_TOOLS      += py3-oauthlib

ALL_TOOLS      += py3-onnxconverter-common

ALL_TOOLS      += py3-onnx

ALL_TOOLS      += py3-onnxmltools

ALL_TOOLS      += py3-opt-einsum

ALL_TOOLS      += py3-packaging

ALL_TOOLS      += py3-pandas

ALL_TOOLS      += py3-pandocfilters

ALL_TOOLS      += py3-parsimonious

ALL_TOOLS      += py3-parso

ALL_TOOLS      += py3-pastel

ALL_TOOLS      += py3-pathlib2

ALL_TOOLS      += py3-pathspec

ALL_TOOLS      += py3-pbr

ALL_TOOLS      += py3-pexpect

ALL_TOOLS      += py3-pickleshare

ALL_TOOLS      += py3-pillow

ALL_TOOLS      += py3-pip

ALL_TOOLS      += py3-pkgconfig

ALL_TOOLS      += py3-pkginfo

ALL_TOOLS      += py3-plac

ALL_TOOLS      += py3-platformdirs

ALL_TOOLS      += py3-pluggy

ALL_TOOLS      += py3-ply

ALL_TOOLS      += py3-poetry-core

ALL_TOOLS      += py3-poetry

ALL_TOOLS      += py3-poetry-plugin-export

ALL_TOOLS      += py3-prettytable

ALL_TOOLS      += py3-prometheus-client

ALL_TOOLS      += py3-prompt_toolkit

ALL_TOOLS      += py3-protobuf

ALL_TOOLS      += py3-prwlock

ALL_TOOLS      += py3-psutil

ALL_TOOLS      += py3-ptyprocess

ALL_TOOLS      += py3-pure-eval

ALL_TOOLS      += py3-pyasn1

ALL_TOOLS      += py3-pyasn1-modules

ALL_TOOLS      += py3-pybind11
py3-pybind11_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/py3-pybind11/2.10.0-f3c5316631d7b62dd5c5ad1242f1862a/lib/python3.9/site-packages/pybind11/include

ALL_TOOLS      += py3-pycodestyle

ALL_TOOLS      += py3-pycparser

ALL_TOOLS      += py3-pycuda

ALL_TOOLS      += py3-pycurl

ALL_TOOLS      += py3-pydantic

ALL_TOOLS      += py3-pydot

ALL_TOOLS      += py3-pyflakes

ALL_TOOLS      += py3-pygithub

ALL_TOOLS      += py3-pygments

ALL_TOOLS      += py3-pyjwt

ALL_TOOLS      += py3-pylev

ALL_TOOLS      += py3-pylint

ALL_TOOLS      += py3-py

ALL_TOOLS      += py3-pynacl

ALL_TOOLS      += py3-pyparsing

ALL_TOOLS      += py3-pyrsistent

ALL_TOOLS      += py3-pysocks

ALL_TOOLS      += py3-pysqlite3

ALL_TOOLS      += py3-pytest-cov

ALL_TOOLS      += py3-pytest

ALL_TOOLS      += py3-pytest-runner

ALL_TOOLS      += py3-python-daemon

ALL_TOOLS      += py3-python-dateutil

ALL_TOOLS      += py3-python-ldap

ALL_TOOLS      += py3-python-rapidjson

ALL_TOOLS      += py3-pythran

ALL_TOOLS      += py3-pytoml

ALL_TOOLS      += py3-pytools

ALL_TOOLS      += py3-pytz

ALL_TOOLS      += py3-pyyaml

ALL_TOOLS      += py3-pyzmq

ALL_TOOLS      += py3-regex

ALL_TOOLS      += py3-repoze-lru

ALL_TOOLS      += py3-requests

ALL_TOOLS      += py3-requests-oauthlib

ALL_TOOLS      += py3-requests-toolbelt

ALL_TOOLS      += py3-rich

ALL_TOOLS      += py3-rsa

ALL_TOOLS      += py3-scandir

ALL_TOOLS      += py3-schema

ALL_TOOLS      += py3-scikit-learn

ALL_TOOLS      += py3-scinum

ALL_TOOLS      += py3-scipy

ALL_TOOLS      += py3-seaborn

ALL_TOOLS      += py3-secretstorage

ALL_TOOLS      += py3-semantic-version

ALL_TOOLS      += py3-send2trash

ALL_TOOLS      += py3-setuptools

ALL_TOOLS      += py3-setuptools-rust

ALL_TOOLS      += py3-setuptools-scm

ALL_TOOLS      += py3-shellingham

ALL_TOOLS      += py3-simplegeneric

ALL_TOOLS      += py3-singledispatch

ALL_TOOLS      += py3-six

ALL_TOOLS      += py3-skl2onnx

ALL_TOOLS      += py3-smmap

ALL_TOOLS      += py3-sniffio

ALL_TOOLS      += py3-soupsieve

ALL_TOOLS      += py3-sqlalchemy

ALL_TOOLS      += py3-stack-data

ALL_TOOLS      += py3-stevedore

ALL_TOOLS      += py3-subprocess32

ALL_TOOLS      += py3-sympy

ALL_TOOLS      += py3-tables

ALL_TOOLS      += py3-tenacity

ALL_TOOLS      += py3-tensorboard-data-server

ALL_TOOLS      += py3-tensorboard

ALL_TOOLS      += py3-tensorboard-plugin-wit

ALL_TOOLS      += py3-tensorflow-estimator

ALL_TOOLS      += py3-tensorflow

ALL_TOOLS      += py3-termcolor

ALL_TOOLS      += py3-terminado

ALL_TOOLS      += py3-testpath

ALL_TOOLS      += py3-theano

ALL_TOOLS      += py3-threadpoolctl

ALL_TOOLS      += py3-tinycss2

ALL_TOOLS      += py3-tomli

ALL_TOOLS      += py3-tomli-w

ALL_TOOLS      += py3-tomlkit

ALL_TOOLS      += py3-toml

ALL_TOOLS      += py3-tornado

ALL_TOOLS      += py3-tqdm

ALL_TOOLS      += py3-traitlets

ALL_TOOLS      += py3-typed-ast

ALL_TOOLS      += py3-typing-extensions

ALL_TOOLS      += py3-uhi

ALL_TOOLS      += py3-uncertainties

ALL_TOOLS      += py3-uproot

ALL_TOOLS      += py3-urllib3

ALL_TOOLS      += py3-vector

ALL_TOOLS      += py3-virtualenv-clone

ALL_TOOLS      += py3-virtualenv

ALL_TOOLS      += py3-virtualenvwrapper

ALL_TOOLS      += py3-wcwidth

ALL_TOOLS      += py3-webencodings

ALL_TOOLS      += py3-werkzeug

ALL_TOOLS      += py3-wheel

ALL_TOOLS      += py3-widgetsnbextension

ALL_TOOLS      += py3-wrapt

ALL_TOOLS      += py3-xgboost

ALL_TOOLS      += py3-zipp

ALL_TOOLS      += pyclang

ALL_TOOLS      += pydata
pydata_EX_FLAGS_LDFLAGS  := $(PYDATA_BASE)/lib/pydata.o
pydata_EX_FLAGS_NO_RECURSIVE_EXPORT  := 1

ALL_TOOLS      += pyquen
pyquen_EX_LIB := pyquen
pyquen_EX_USE := pythia6 lhapdf

ALL_TOOLS      += pythia6_headers
pythia6_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pythia6/426-154c9fa9309a9a96c7e05f80622d33eb/include
pythia6_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += pythia6
pythia6_EX_LIB := pythia6 pythia6_dummy pythia6_pdfdummy
pythia6_EX_USE := pythia6_headers f77compiler

ALL_TOOLS      += pythia8
pythia8_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/pythia8/306-687274b5349576f1777e8bb247ec5a88/include
pythia8_EX_LIB := pythia8
pythia8_EX_USE := root_cxxdefaults cxxcompiler hepmc3 hepmc lhapdf
pythia8_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += python3
python3_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/python3/3.9.14-85e7e49ccb4e877d2ed5b0ae311c699f/include/python3.9
python3_EX_LIB := python3.9
python3_EX_USE := sockets

ALL_TOOLS      += python-paths

ALL_TOOLS      += python_tools

ALL_TOOLS      += qd_f_main
qd_f_main_EX_LIB := qd_f_main
qd_f_main_EX_USE := qd

ALL_TOOLS      += qd
qd_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/qd/2.3.13-3e887020a27b6bda66261a2b30e928f2/include
qd_EX_LIB := qdmod qd

ALL_TOOLS      += rdma-core
rdma-core_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/rdma-core/39.1-ffcb004e45628f4e0489cd8958dc427f/include

ALL_TOOLS      += re2
re2_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/re2/2021-06-01-bd9694d6b7214964900dc18f8d290757/include
re2_EX_LIB := re2

ALL_TOOLS      += rivet
rivet_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/rivet/3.1.7-18cc6082ff08237274a9407a2c215412/include
rivet_EX_LIB := Rivet
rivet_EX_USE := hepmc fastjet fastjet-contrib gsl yoda

ALL_TOOLS      += rocm
rocm_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/rocm/5.4.3-d077dcce910eb4e8562ad8848b8fcd3a/include
rocm_EX_LIB := amdhip64
rocm_EX_FLAGS_CPPDEFINES  := -D__HIP_PLATFORM_HCC__ -D__HIP_PLATFORM_AMD__
rocm_EX_FLAGS_ROCM_FLAGS  := -fgpu-rdc --offload-arch=gfx900 --offload-arch=gfx906 --offload-arch=gfx908 --offload-arch=gfx90a --target=x86_64-redhat-linux-gnu --gcc-toolchain=$(COMPILER_PATH)
rocm_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += rocm-rocrand
rocm-rocrand_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/rocm-rocrand/5.4.3-436b888cbc93d5eb5f4906eab77634c8/include
rocm-rocrand_EX_LIB := hiprand rocrand
rocm-rocrand_EX_USE := rocm
rocm-rocrand_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += roofitcore
roofitcore_EX_LIB := RooFitCore
roofitcore_EX_USE := rootcore roothistmatrix rootgpad rootminuit root_cxxdefaults

ALL_TOOLS      += roofit
roofit_EX_LIB := RooFit
roofit_EX_USE := roofitcore rootcore rootmath roothistmatrix

ALL_TOOLS      += roostats
roostats_EX_LIB := RooStats
roostats_EX_USE := roofitcore roofit rootcore roothistmatrix rootgpad

ALL_TOOLS      += rootcling
rootcling_EX_LIB := Core
rootcling_EX_USE := root_interface sockets pcre zlib xz
TOOLS_OVERRIDABLE_FLAGS  +=ROOTCLING_ARGS

ALL_TOOLS      += rootcore
rootcore_EX_LIB := Tree Net
rootcore_EX_USE := rootmathcore rootthread

ALL_TOOLS      += root_cxxdefaults

ALL_TOOLS      += rootdataframe
rootdataframe_EX_LIB := ROOTDataFrame
rootdataframe_EX_USE := rootcore rootgraphics roothistmatrix rootrio rootvecops

ALL_TOOLS      += rooteg
rooteg_EX_LIB := EG
rooteg_EX_USE := rootgraphics

ALL_TOOLS      += rooteve
rooteve_EX_LIB := Eve
rooteve_EX_USE := rootgeompainter rootrgl rootged

ALL_TOOLS      += rootfoam
rootfoam_EX_LIB := Foam
rootfoam_EX_USE := roothistmatrix

ALL_TOOLS      += rootged
rootged_EX_LIB := Ged
rootged_EX_USE := rootgui

ALL_TOOLS      += rootgeom
rootgeom_EX_LIB := Geom
rootgeom_EX_USE := rootrio rootmathcore

ALL_TOOLS      += rootgeompainter
rootgeompainter_EX_LIB := GeomPainter
rootgeompainter_EX_USE := rootgeom rootgraphics

ALL_TOOLS      += rootglew
rootglew_EX_LIB := GLEW

ALL_TOOLS      += rootgpad
rootgpad_EX_LIB := Gpad Graf
rootgpad_EX_USE := roothistmatrix

ALL_TOOLS      += rootgraphics
rootgraphics_EX_LIB := TreePlayer Graf3d Postscript
rootgraphics_EX_USE := rootgpad

ALL_TOOLS      += rootguihtml
rootguihtml_EX_LIB := GuiHtml
rootguihtml_EX_USE := rootgui rootinteractive

ALL_TOOLS      += rootgui
rootgui_EX_LIB := Gui
rootgui_EX_USE := rootgpad

ALL_TOOLS      += roothistmatrix
roothistmatrix_EX_LIB := Hist Matrix
roothistmatrix_EX_USE := rootcore

ALL_TOOLS      += roothistpainter
roothistpainter_EX_LIB := HistPainter
roothistpainter_EX_USE := roothistmatrix rootgpad rootimt

ALL_TOOLS      += roothtml
roothtml_EX_LIB := Html
roothtml_EX_USE := rootgpad

ALL_TOOLS      += rootimt
rootimt_EX_LIB := Imt
rootimt_EX_USE := rootthread tbb

ALL_TOOLS      += rootinteractive
rootinteractive_EX_LIB := Gui
rootinteractive_EX_USE := libjpeg-turbo libpng rootgpad rootrint

ALL_TOOLS      += root_interface
root_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/lcg/root/6.26.11-43b2b4148365ed22373c74503805094a/include
root_interface_EX_USE := root_cxxdefaults
root_interface_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += rootmathcore
rootmathcore_EX_LIB := MathCore
rootmathcore_EX_USE := rootcling

ALL_TOOLS      += rootmath
rootmath_EX_LIB := GenVector MathMore
rootmath_EX_USE := rootcore gsl

ALL_TOOLS      += rootminuit2
rootminuit2_EX_LIB := Minuit2
rootminuit2_EX_USE := rootgpad

ALL_TOOLS      += rootminuit
rootminuit_EX_LIB := Minuit
rootminuit_EX_USE := rootgpad

ALL_TOOLS      += root
root_EX_USE := rootphysics
root_EX_FLAGS_CXXMODULES  := 0
root_EX_FLAGS_GENREFLEX_FAILES_ON_WARNS  := --fail_on_warnings

ALL_TOOLS      += rootmlp
rootmlp_EX_LIB := MLP
rootmlp_EX_USE := rootgraphics

ALL_TOOLS      += rootntuple
rootntuple_EX_LIB := ROOTNTuple
rootntuple_EX_USE := rootvecops rootthread

ALL_TOOLS      += rootphysics
rootphysics_EX_LIB := Physics
rootphysics_EX_USE := roothistmatrix

ALL_TOOLS      += rootpy
rootpy_EX_USE := rootgraphics

ALL_TOOLS      += rootpymva
rootpymva_EX_LIB := PyMVA
rootpymva_EX_USE := roottmva numpy-c-api

ALL_TOOLS      += rootrflx
rootrflx_EX_USE := root_interface rootcling
rootrflx_EX_FLAGS_GENREFLEX_CPPFLAGS  := -DCMS_DICT_IMPL -D_REENTRANT -DGNUSOURCE -D__STRICT_ANSI__
rootrflx_EX_FLAGS_GENREFLEX_GCCXMLOPT  := -m64
TOOLS_OVERRIDABLE_FLAGS  +=GENREFLEX_CPPFLAGS

ALL_TOOLS      += rootrgl
rootrgl_EX_LIB := RGL
rootrgl_EX_USE := rootglew rootgui rootinteractive rootgraphics

ALL_TOOLS      += rootrint
rootrint_EX_LIB := Rint
rootrint_EX_USE := rootcling

ALL_TOOLS      += rootrio
rootrio_EX_LIB := RIO
rootrio_EX_USE := rootcling

ALL_TOOLS      += rootsmatrix
rootsmatrix_EX_LIB := Smatrix
rootsmatrix_EX_USE := rootcling

ALL_TOOLS      += rootspectrum
rootspectrum_EX_LIB := Spectrum
rootspectrum_EX_USE := roothistmatrix

ALL_TOOLS      += rootthread
rootthread_EX_LIB := Thread
rootthread_EX_USE := rootrio

ALL_TOOLS      += roottmva
roottmva_EX_LIB := TMVA
roottmva_EX_USE := rootmlp rootminuit

ALL_TOOLS      += rootvecops
rootvecops_EX_LIB := ROOTVecOps
rootvecops_EX_USE := rootcore

ALL_TOOLS      += rootx11
rootx11_EX_LIB := GX11
rootx11_EX_USE := rootcling

ALL_TOOLS      += rootxmlio
rootxmlio_EX_LIB := XMLIO
rootxmlio_EX_USE := rootrio

ALL_TOOLS      += rootxml
rootxml_EX_LIB := XMLParser
rootxml_EX_USE := rootcore libxml2

ALL_TOOLS      += sanitizer-flags-asan

ALL_TOOLS      += sanitizer-flags-ubsan

ALL_TOOLS      += scitokens-cpp
scitokens-cpp_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/scitokens-cpp/0.7.0-74812507075698569e96969a300f3f66/include
scitokens-cpp_EX_LIB := SciTokens
scitokens-cpp_EX_USE := sqlite libuuid curl

ALL_TOOLS      += self
self_EX_INCLUDE := /afs/cern.ch/user/m/mwitt/public/CPinHToTauTau/CMSSW_13_0_19/src /afs/cern.ch/user/m/mwitt/public/CPinHToTauTau/CMSSW_13_0_19/include/el9_amd64_gcc11/src /afs/cern.ch/user/m/mwitt/public/CPinHToTauTau/CMSSW_13_0_19/include/LCG /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/cmssw/CMSSW_13_0_19/src
self_EX_LIBDIR := /afs/cern.ch/user/m/mwitt/public/CPinHToTauTau/CMSSW_13_0_19/biglib/el9_amd64_gcc11 /afs/cern.ch/user/m/mwitt/public/CPinHToTauTau/CMSSW_13_0_19/lib/el9_amd64_gcc11 /afs/cern.ch/user/m/mwitt/public/CPinHToTauTau/CMSSW_13_0_19/external/el9_amd64_gcc11/lib /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/cmssw/CMSSW_13_0_19/biglib/el9_amd64_gcc11 /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/cmssw/CMSSW_13_0_19/lib/el9_amd64_gcc11 /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/cmssw/CMSSW_13_0_19/external/el9_amd64_gcc11/lib
self_EX_FLAGS_ALPAKA_BACKENDS  := cuda rocm serial
self_EX_FLAGS_CHECK_PRIVATE_HEADERS  := 1
self_EX_FLAGS_CODE_CHECK_ALPAKA_BACKEND  := serial
self_EX_FLAGS_DEFAULT_COMPILER  := gcc
self_EX_FLAGS_ENABLE_LTO  := 1
self_EX_FLAGS_ENABLE_PGO  := 0
self_EX_FLAGS_EXTERNAL_SYMLINK  := PATH LIBDIR CMSSW_SEARCH_PATH
self_EX_FLAGS_LLVM_ANALYZER  := llvm-analyzer
self_EX_FLAGS_NO_EXTERNAL_RUNTIME  := LD_LIBRARY_PATH PATH CMSSW_SEARCH_PATH
TOOLS_OVERRIDABLE_FLAGS  +=CPPDEFINES CXXFLAGS FFLAGS CFLAGS CPPFLAGS LDFLAGS CUDA_FLAGS CUDA_LDFLAGS LTO_FLAGS PGO_FLAGS ROCM_FLAGS ROCM_LDFLAGS
self_EX_FLAGS_SCRAM_TARGETS  := haswell sandybridge nehalem
self_EX_FLAGS_SKIP_TOOLS_SYMLINK  := cxxcompiler ccompiler f77compiler gcc-cxxcompiler gcc-ccompiler gcc-f77compiler llvm-cxxcompiler llvm-ccompiler llvm-f77compiler llvm-analyzer-cxxcompiler llvm-analyzer-ccompiler icc-cxxcompiler icc-ccompiler icc-f77compiler x11 dpm
self_EX_FLAGS_SYMLINK_DEPTH_CMSSW_SEARCH_PATH  := 2
self_ORDER := 20000
IS_PATCH:=

ALL_TOOLS      += sherpa
sherpa_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/sherpa/2.2.15-74d02d077f398b33bec77eb885db2b1d/include/SHERPA-MC
sherpa_EX_LIB := SherpaMain ToolsMath ToolsOrg
sherpa_EX_USE := root_cxxdefaults hepmc lhapdf qd blackhat fastjet sqlite openmpi openloops

ALL_TOOLS      += sigcpp
sigcpp_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/sigcpp/3.2.0-aa26ce8063258bda618dfc2c25b06c89/include/sigc++-3.0
sigcpp_EX_LIB := sigc-3.0
sigcpp_EX_USE := root_cxxdefaults

ALL_TOOLS      += sloccount

ALL_TOOLS      += sockets
sockets_EX_LIB := crypt dl rt

ALL_TOOLS      += sqlite
sqlite_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/sqlite/3.36.0-0f26675926fd468efdd431be2b62785e/include
sqlite_EX_LIB := sqlite3
sqlite_EX_USE := root_cxxdefaults

ALL_TOOLS      += starlight
starlight_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/starlight/r193-0682b464b6be741957f835e3dd810902/include
starlight_EX_LIB := Starlib
starlight_EX_USE := root_cxxdefaults clhep

ALL_TOOLS      += stdcxx-fs
stdcxx-fs_EX_LIB := stdc++fs

ALL_TOOLS      += tauolapp
tauolapp_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/tauolapp/1.1.8-f64b0c86b11d0a40b5f0141624b4a6ef/include
tauolapp_EX_LIB := TauolaCxxInterface TauolaFortran TauolaTauSpinner TauolaHepMC TauolaHEPEVT
tauolapp_EX_USE := root_cxxdefaults hepmc f77compiler pythia8 lhapdf

ALL_TOOLS      += tbbbind
tbbbind_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/tbb/v2021.8.0-ffea5a11aae93f99db95a8d25c58844d/include
tbbbind_EX_LIB := tbbbind_2_0
tbbbind_EX_USE := tbb
tbbbind_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += tbb
tbb_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/tbb/v2021.8.0-ffea5a11aae93f99db95a8d25c58844d/include
tbb_EX_LIB := tbb
tbb_EX_USE := root_cxxdefaults
tbb_EX_FLAGS_CPPDEFINES  := -DTBB_USE_GLIBCXX_VERSION=110201 -DTBB_SUPPRESS_DEPRECATED_MESSAGES -DTBB_PREVIEW_RESUMABLE_TASKS=1 -DTBB_PREVIEW_TASK_GROUP_EXTENSIONS=1
tbb_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += tcmalloc_minimal
tcmalloc_minimal_EX_LIB := tcmalloc_minimal

ALL_TOOLS      += tcmalloc
tcmalloc_EX_LIB := tcmalloc

ALL_TOOLS      += tensorflow-cc
tensorflow-cc_EX_LIB := tensorflow_cc
tensorflow-cc_EX_USE := tensorflow-framework eigen libpng sqlite

ALL_TOOLS      += tensorflow-c
tensorflow-c_EX_LIB := tensorflow
tensorflow-c_EX_USE := tensorflow-framework eigen libpng sqlite

ALL_TOOLS      += tensorflow-executable_run_options
tensorflow-executable_run_options_EX_LIB := executable_run_options
tensorflow-executable_run_options_EX_USE := tensorflow

ALL_TOOLS      += tensorflow-framework
tensorflow-framework_EX_LIB := tensorflow_framework
tensorflow-framework_EX_USE := tensorflow giflib zlib libjpeg-turbo protobuf

ALL_TOOLS      += tensorflow
tensorflow_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/tensorflow/2.6.4-d11d9da6622cff6ad896cf997c355738/include
tensorflow_EX_FLAGS_SYSTEM_INCLUDE  := 1

ALL_TOOLS      += tensorflow-runtime
tensorflow-runtime_EX_LIB := cpu_function_runtime
tensorflow-runtime_EX_USE := tensorflow

ALL_TOOLS      += tensorflow-tf2xla
tensorflow-tf2xla_EX_LIB := tf2xla
tensorflow-tf2xla_EX_USE := tensorflow

ALL_TOOLS      += tensorflow-xla_compiled_cpu_function
tensorflow-xla_compiled_cpu_function_EX_LIB := xla_compiled_cpu_function
tensorflow-xla_compiled_cpu_function_EX_USE := tensorflow

ALL_TOOLS      += thepeg
thepeg_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/thepeg/2.2.2-1a088b5f65da98c866d4cf0482fbfa9c/include
thepeg_EX_LIB := ThePEG LesHouches
thepeg_EX_USE := root_cxxdefaults lhapdf gsl

ALL_TOOLS      += tinyxml2
tinyxml2_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/tinyxml2/6.2.0-c2bad61e58f94d6db8f640afbd739be2/include
tinyxml2_EX_LIB := tinyxml2

ALL_TOOLS      += tkonlineswdb
tkonlineswdb_EX_LIB := DeviceDescriptions Fed9UDeviceFactory
tkonlineswdb_EX_USE := tkonlinesw oracle oracleocci

ALL_TOOLS      += tkonlinesw
tkonlinesw_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/tkonlinesw/4.2.0-1_gcc7-c67dc574ab14cf69b68b80d14c3f8413/include
tkonlinesw_EX_LIB := ICUtils Fed9UUtils
tkonlinesw_EX_USE := root_cxxdefaults xerces-c
tkonlinesw_EX_FLAGS_CXXFLAGS  := -DCMS_TK_64BITS

ALL_TOOLS      += triton-inference-client
triton-inference-client_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/triton-inference-client/2.11.0-44a844f581c7352f5bc91ff5f9e3e26a/include
triton-inference-client_EX_LIB := grpcclient
triton-inference-client_EX_USE := protobuf grpc cuda re2

ALL_TOOLS      += ucx
ucx_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/ucx/1.12.1-c9b53883e0f074b0369dd174026a7e8e/include
ucx_EX_LIB := ucp uct ucs ucm

ALL_TOOLS      += utm
utm_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/utm/utm_0.11.2-f53b43dd28a187aa608b8d058fcc755c/include
utm_EX_LIB := tmeventsetup tmtable tmxsd tmgrammar tmutil
utm_EX_USE := root_cxxdefaults

ALL_TOOLS      += valgrind
valgrind_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/valgrind/3.17.0-5d0361c2ec2c1fed0e1b1636c68aa3dd/include
valgrind_EX_USE := root_cxxdefaults

ALL_TOOLS      += vdt_headers
vdt_headers_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/cms/vdt/0.4.3-30d69e6d7a9d2a719cf771bf4b327768/include
vdt_headers_EX_USE := root_cxxdefaults

ALL_TOOLS      += vdt
vdt_EX_LIB := vdt
vdt_EX_USE := vdt_headers

ALL_TOOLS      += vecgeom_interface
vecgeom_interface_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/vecgeom/v1.1.17-10b2db97d6b1736305b7b352fc807529/include /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/vecgeom/v1.1.17-10b2db97d6b1736305b7b352fc807529/include/VecGeom
vecgeom_interface_EX_USE := root_cxxdefaults

ALL_TOOLS      += vecgeom
vecgeom_EX_LIB := vecgeom
vecgeom_EX_USE := vecgeom_interface

ALL_TOOLS      += x11
x11_EX_USE := sockets

ALL_TOOLS      += xerces-c
xerces-c_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xerces-c/3.1.3-96261f23c7d6fbfb7d59be544bd882f3/include
xerces-c_EX_LIB := xerces-c
xerces-c_EX_USE := root_cxxdefaults

ALL_TOOLS      += xgboost
xgboost_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xgboost/1.6.2-91ec2674a1cadbad486c39e4b005538c/include
xgboost_EX_LIB := xgboost

ALL_TOOLS      += xpmem
xpmem_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xpmem/v2.6.3-20220308-e4f1f864ac3800b7f784bbf7f543927d/include
xpmem_EX_LIB := xpmem

ALL_TOOLS      += xrdcl-record

ALL_TOOLS      += xrootd
xrootd_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xrootd/5.5.1-4928a053ea0ac3b7a02379c7b3391056/include/xrootd /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xrootd/5.5.1-4928a053ea0ac3b7a02379c7b3391056/include/xrootd/private
xrootd_EX_LIB := XrdUtils XrdCl
xrootd_EX_USE := root_cxxdefaults scitokens-cpp

ALL_TOOLS      += xtensor
xtensor_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xtensor/0.24.1-d8b77ecbecf86f2ed0762c65ec1e33c6/include
xtensor_EX_USE := xtl

ALL_TOOLS      += xtl
xtl_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xtl/0.7.4-7a266262fc00ce525fb46d766cb524cc/include

ALL_TOOLS      += xz
xz_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/xz/5.2.5-f2a869b3041e7ff7a6dcea149634b26b/include
xz_EX_LIB := lzma
xz_EX_USE := root_cxxdefaults

ALL_TOOLS      += yoda
yoda_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/yoda/1.9.7-c53cb3ac7f0a5881a7b8f2738e416aaa/include
yoda_EX_LIB := YODA
yoda_EX_USE := root_cxxdefaults

ALL_TOOLS      += zlib
zlib_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/zlib/1.2.11-3dfb2715f3608466b74431b80eb9d788/include
zlib_EX_LIB := z
zlib_EX_USE := root_cxxdefaults

ALL_TOOLS      += zstd
zstd_EX_INCLUDE := /cvmfs/cms.cern.ch/el9_amd64_gcc11/external/zstd/1.5.2-4341c217fdb8259db5ea57a7c59a8008/include
zstd_EX_LIB := zstd
zstd_EX_USE := root_cxxdefaults

