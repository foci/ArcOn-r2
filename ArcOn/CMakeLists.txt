# Use make edit_cache for settings

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
cmake_policy(SET CMP0037 OLD)

FIND_PACKAGE(deal.II 8.3 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../trigger/ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#
# Are all dependencies fullfilled?
#
IF(NOT DEAL_II_WITH_PETSC OR NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_WITH_P4EST = ON
One or all of these are OFF in your installation but are required for this tutorial step."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

ADD_EXECUTABLE(./lib/$ENV{ArcOn_dimension}/ArcOn $ENV{ArcOn_DIR}/ArcOn/source/Main.cc)
DEAL_II_SETUP_TARGET(./lib/$ENV{ArcOn_dimension}/ArcOn)
INCLUDE_DIRECTORIES($ENV{ArcOn_DIR}/ArcOn/include)


