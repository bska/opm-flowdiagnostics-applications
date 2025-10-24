include(CMakeFindDependencyMacro)

find_dependency(OpmFDEngine 2024.04)
find_dependency(Boost
  COMPONENTS "filesystem"
)
find_dependency(ecl)

include("${CMAKE_CURRENT_LIST_DIR}/OpmFDApplicationsTargets.cmake")
