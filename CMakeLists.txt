# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

cmake_minimum_required(VERSION 3.27)

project(opm-flowdiagnostics-applications
  VERSION   2024.04
  LANGUAGES C CXX
)

enable_testing()

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_C_STANDARD          99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS        OFF)

set(CMAKE_DEBUG_POSTFIX "-d")

#----------------------------------------------------------------------------
# Prerequisites
#----------------------------------------------------------------------------

# Disable Git here in order to not run the internal Git-based version logic
# in LibECL's build system.  If we do not do this, we get a build failure in
# the sub-project due to the compiler not being able to use the string literal
#
#    "release/ResInsight/2017"
#
# as an integer in lib/util/ecl_version.cpp (ecl_version_get_major_version()).
#
# This is a hack.

set(CMAKE_DISABLE_FIND_PACKAGE_Git
  TRUE CACHE INTERNAL
  "Disable Git-based version logic in LibECL"
)

include(FetchContent)
FetchContent_Declare(LibECL
  GIT_REPOSITORY https://github.com/equinor/resdata.git
  GIT_TAG        539b836ac8fec7f1be9f3474419f189d54e10edf # 2.7.0
  SYSTEM
)

FetchContent_MakeAvailable(LibECL)

find_package(Boost
  1.44
  REQUIRED
  COMPONENTS "filesystem" "system" "unit_test_framework"
)

find_package(OpmFDEngine 2024.04 REQUIRED)

#----------------------------------------------------------------------------
# Main Library
#----------------------------------------------------------------------------

add_library(OpmFDApplications STATIC)

set_target_properties (OpmFDApplications PROPERTIES
  SOVERSION 1
  VERSION   2024.04
)

target_sources(OpmFDApplications
  PRIVATE
    opm/utility/ECLCaseUtilities.cpp
    opm/utility/ECLEndPointScaling.cpp
    opm/utility/ECLFluxCalc.cpp
    opm/utility/ECLGraph.cpp
    opm/utility/ECLPropertyUnitConversion.cpp
    opm/utility/ECLPropTable.cpp
    opm/utility/ECLPvtCommon.cpp
    opm/utility/ECLPvtCurveCollection.cpp
    opm/utility/ECLPvtGas.cpp
    opm/utility/ECLPvtOil.cpp
    opm/utility/ECLPvtWater.cpp
    opm/utility/ECLRegionMapping.cpp
    opm/utility/ECLResultData.cpp
    opm/utility/ECLSaturationFunc.cpp
    opm/utility/ECLTableInterpolation1D.cpp
    opm/utility/ECLUnitHandling.cpp
    opm/utility/ECLWellSolution.cpp

  PUBLIC FILE_SET applications_public_headers
    TYPE HEADERS
    FILES
      opm/utility/ECLCaseUtilities.hpp
      opm/utility/ECLEndPointScaling.hpp
      opm/utility/ECLFluxCalc.hpp
      opm/utility/ECLGraph.hpp
      opm/utility/ECLPhaseIndex.hpp
      opm/utility/ECLPiecewiseLinearInterpolant.hpp
      opm/utility/ECLPropertyUnitConversion.hpp
      opm/utility/ECLPropTable.hpp
      opm/utility/ECLPvtCommon.hpp
      opm/utility/ECLPvtCurveCollection.hpp
      opm/utility/ECLPvtGas.hpp
      opm/utility/ECLPvtOil.hpp
      opm/utility/ECLPvtWater.hpp
      opm/utility/ECLRegionMapping.hpp
      opm/utility/ECLResultData.hpp
      opm/utility/ECLSaturationFunc.hpp
      opm/utility/ECLTableInterpolation1D.hpp
      opm/utility/ECLUnitHandling.hpp
      opm/utility/ECLWellSolution.hpp

      opm/utility/imported/Units.hpp
)


target_link_libraries(OpmFDApplications
  PUBLIC
    OpmFDEngine::OpmFDEngine
    Boost::filesystem
    ecl
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

write_basic_package_version_file(OpmFDApplicationsConfigVersion.cmake
  VERSION ${OpmFDApplications_VERSION}
  COMPATIBILITY ExactVersion
)

install(TARGETS OpmFDApplications EXPORT OpmFDApplicationsTargets
  RUNTIME           # Following options apply to runtime artifacts.
    COMPONENT OpmFDApplications_Runtime
  LIBRARY           # Following options apply to library artifacts.
    COMPONENT OpmFDApplications_Runtime
    NAMELINK_COMPONENT OpmFDApplications_Development
  ARCHIVE           # Following options apply to archive artifacts.
    COMPONENT OpmFDApplications_Development
  FILE_SET applications_public_headers # Following options apply to file set HEADERS.
    COMPONENT OpmFDApplications_Development
)

install(EXPORT OpmFDApplicationsTargets
  FILE OpmFDApplicationsTargets.cmake
  NAMESPACE OpmFDApplications::
  DESTINATION share/cmake/OpmFDApplications
)

install(FILES
    "OpmFDApplicationsConfig.cmake"
    "${CMAKE_BINARY_DIR}/OpmFDApplicationsConfigVersion.cmake"
  DESTINATION
    share/cmake/OpmFDApplications
)

#----------------------------------------------------------------------------
# Example Applications
#----------------------------------------------------------------------------

add_library(FDExampleParameterHandling STATIC)
target_sources(FDExampleParameterHandling
  PRIVATE
    opm/utility/imported/ParameterGroup.cpp
    opm/utility/imported/ParameterRequirement.cpp
    opm/utility/imported/ParameterTools.cpp
    opm/utility/imported/Parameter.cpp

  PUBLIC FILE_SET parameter_public_headers
    TYPE HEADERS
    FILES
      opm/utility/imported/ParameterGroup.hpp
      opm/utility/imported/ParameterGroup_impl.hpp
      opm/utility/imported/ParameterMapItem.hpp
      opm/utility/imported/ParameterRequirement.hpp
      opm/utility/imported/ParameterStrings.hpp
      opm/utility/imported/ParameterTools.hpp
      opm/utility/imported/Parameter.hpp
)

add_executable(computeTracers)
target_sources(computeTracers
  PRIVATE
    examples/computeTracers.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(computeTracers
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
    OpmFDEngine::OpmFDEngine
    Boost::filesystem
)

#add_executable(dynamicCellProperty)
#target_sources(dynamicCellProperty
#  PRIVATE
#    examples/dynamicCellProperty.cpp
#    examples/exampleSetup.hpp
#)
#target_link_libraries(dynamicCellProperty
#  PUBLIC
#    OpmFDApplications
#    FDExampleParameterHandling
#    OpmFDEngine::OpmFDEngine
#    Boost::filesystem
#)

add_executable(extractFromRestart)
target_sources(extractFromRestart
  PRIVATE
    examples/extractFromRestart.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(extractFromRestart
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
)

add_executable(extractPropCurves)
target_sources(extractPropCurves
  PRIVATE
    examples/extractPropCurves.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(extractPropCurves
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
    OpmFDEngine::OpmFDEngine
    Boost::filesystem
)

install(TARGETS extractPropCurves
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

add_executable(computeFlowStorageCurve)
target_sources(computeFlowStorageCurve
  PRIVATE
    examples/computeFlowStorageCurve.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(computeFlowStorageCurve
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
    OpmFDEngine::OpmFDEngine
)

add_executable(computeLocalSolutions)
target_sources(computeLocalSolutions
  PRIVATE
    examples/computeLocalSolutions.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(computeLocalSolutions
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
    OpmFDEngine::OpmFDEngine
)

add_executable(computePhaseFluxes)
target_sources(computePhaseFluxes
  PRIVATE
    examples/computePhaseFluxes.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(computePhaseFluxes
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
    OpmFDEngine::OpmFDEngine
    Boost::filesystem
    Boost::system
)

add_executable(computeToFandTracers)
target_sources(computeToFandTracers
  PRIVATE
    examples/computeToFandTracers.cpp
    examples/exampleSetup.hpp
)
target_link_libraries(computeToFandTracers
  PUBLIC
    OpmFDApplications
    FDExampleParameterHandling
    OpmFDEngine::OpmFDEngine
)

#----------------------------------------------------------------------------
# Unit Tests
#----------------------------------------------------------------------------

add_executable(test_eclendpointscaling tests/test_eclendpointscaling.cpp)
target_link_libraries(test_eclendpointscaling
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_eclpropertyunitconversion tests/test_eclpropertyunitconversion.cpp)
target_link_libraries(test_eclpropertyunitconversion
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_eclproptable tests/test_eclproptable.cpp)
target_link_libraries(test_eclproptable
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_eclpvtcommon tests/test_eclpvtcommon.cpp)
target_link_libraries(test_eclpvtcommon
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_eclregionmapping tests/test_eclregionmapping.cpp)
target_link_libraries(test_eclregionmapping
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_eclsimple1dinterpolant tests/test_eclsimple1dinterpolant.cpp)
target_link_libraries(test_eclsimple1dinterpolant
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_eclunithandling tests/test_eclunithandling.cpp)
target_link_libraries(test_eclunithandling
  PRIVATE OpmFDApplications Boost::unit_test_framework
  )

add_executable(test_propcurves_integration tests/test_propcurves_integration.cpp)
target_link_libraries(test_propcurves_integration
  PUBLIC Boost::filesystem
  PRIVATE OpmFDApplications Boost::unit_test_framework
)

add_custom_target(copy_integration_test_input ALL
  COMMAND ${CMAKE_COMMAND} -E echo "Making integration test data available in output directory"
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/Integration"
    "${CMAKE_CURRENT_BINARY_DIR}/TestInput/Integration"
)

#----------------------------------------------------------------------------
# Individual Unit Tests
#----------------------------------------------------------------------------

add_test(Test.eclendpointscaling test_eclendpointscaling)
add_test(Test.eclproptable test_eclproptable)
add_test(Test.eclpvtcommon test_eclpvtcommon)
add_test(Test.eclregionmapping test_eclregionmapping)
add_test(Test.eclsimple1dinterpolant test_eclsimple1dinterpolant)
add_test(Test.eclunithandling test_eclunithandling)
add_test(Test.propcurves_integration test_propcurves_integration)
