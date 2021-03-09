CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(mmm_orsa C CXX)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fPIC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -fPIC")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -march=native -fPIC")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")


# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else ()
    set (CMAKE_C_STANDARD 99)
endif ()

ENABLE_TESTING()

INCLUDE_DIRECTORIES(
  ./demo/
  ./third_party/
  ./
)

IF (WIN32)
  INCLUDE_DIRECTORIES(
    third_party/jpeg
    third_party/zlib
    third_party/png)
ENDIF (WIN32)

MACRO (UNIT_TEST NAME EXTRA_LIBS)
  ADD_EXECUTABLE(${NAME}_test ${NAME}_test.cpp)
  TARGET_LINK_LIBRARIES(${NAME}_test
                        ${EXTRA_LIBS} # Extra libs MUST be first.
                        CppUnitLite)
  ADD_TEST(${NAME}_test ${NAME}_test)
ENDMACRO (UNIT_TEST)

# Third parties libraries
ADD_SUBDIRECTORY(third_party)

# Orsa Core library
ADD_SUBDIRECTORY(libOrsa)

# External libraries required for the demo
ADD_SUBDIRECTORY(extras)

# Demo directory
ADD_SUBDIRECTORY(demo)

INSTALL(FILES CMakeLists.txt DESTINATION src)
INSTALL(FILES testing/testing.h DESTINATION src/testing)
