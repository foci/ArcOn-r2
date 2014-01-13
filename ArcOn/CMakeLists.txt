
    # Use make edit_cache for settings

    CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

    FIND_PACKAGE(deal.II 8.0 REQUIRED
      HINTS
        ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR} 
      )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

    DEAL_II_INITIALIZE_CACHED_VARIABLES()

    PROJECT(ArcOn)	
     
    ADD_EXECUTABLE(./lib/$ENV{ArcOn_dimension}/ArcOn $ENV{ArcOn_DIR}/ArcOn/source/Main.cc)

    DEAL_II_SETUP_TARGET(./lib/$ENV{ArcOn_dimension}/ArcOn)

    INCLUDE_DIRECTORIES($ENV{ArcOn_DIR}/ArcOn/include)

    #DEAL_II_INVOKE_AUTOPILOT() 
