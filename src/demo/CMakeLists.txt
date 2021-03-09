ADD_SUBDIRECTORY(libImage)

FIND_PACKAGE(OpenMP) 
if (OPENMP_FOUND) 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# Executable demo to estimate a homography matrix from SIFT matches and
# display the corresponding mosaic.
ADD_EXECUTABLE(demo_orsa_homography orsa_homography.cpp)
TARGET_LINK_LIBRARIES(demo_orsa_homography Orsa image sift)

IF(UNIX)
    SET_TARGET_PROPERTIES(demo_orsa_homography PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra")
ENDIF(UNIX)

INSTALL(FILES CMakeLists.txt
              cmdLine.h
              orsa_homography.cpp
              Rect.hpp
              siftMatch.hpp
              warping.hpp
        DESTINATION src/demo)

add_library(orsa SHARED
    liborsa.cpp 
    )

add_library(worsa SHARED
    wrapper_orsa.c
    )



TARGET_LINK_LIBRARIES(orsa Orsa image sift)
target_link_libraries(worsa orsa)
set_property(TARGET worsa PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(TARGET worsa POST_BUILD
   COMMAND gcc -I${CMAKE_CURRENT_SOURCE_DIR} -E ${CMAKE_CURRENT_SOURCE_DIR}/all.h >all.h
   COMMENT "Creating all.h"
)
