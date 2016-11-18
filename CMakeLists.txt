cmake_minimum_required(VERSION 2.8)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O2 -fopenmp")
# set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELEASE)

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
    set(CMAKE_CXX_COMPILER "g++")
endif()
if(WIN32)
    LINK_DIRECTORIES(C:\\lapack)
    set(CMAKE_CXX_COMPILER "g++-6")
endif()

project(alm)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(OpenMP)
if (OPENMP_FOUND) 
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Version numbers
file(READ ${PROJECT_SOURCE_DIR}/src/version.h version_file)
string(REGEX MATCH "ALAMODE_VERSION = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" alm_version ${version_file})
set(alm_version ${CMAKE_MATCH_1})
MESSAGE("ALM version: ${CMAKE_MATCH_1}")
set(serial "${alm_version}")
set(soserial "1")

# Source code
include_directories("${PROJECT_SOURCE_DIR}/src")
set(SOURCES ${PROJECT_SOURCE_DIR}/src/alm.cpp
            ${PROJECT_SOURCE_DIR}/src/alm_core.cpp
            ${PROJECT_SOURCE_DIR}/src/constraint.cpp
            ${PROJECT_SOURCE_DIR}/src/error.cpp
            ${PROJECT_SOURCE_DIR}/src/fcs.cpp
            ${PROJECT_SOURCE_DIR}/src/files.cpp
            ${PROJECT_SOURCE_DIR}/src/fitting.cpp
            ${PROJECT_SOURCE_DIR}/src/interaction.cpp
            ${PROJECT_SOURCE_DIR}/src/memory.cpp
            ${PROJECT_SOURCE_DIR}/src/patterndisp.cpp
            ${PROJECT_SOURCE_DIR}/src/symmetry.cpp
            ${PROJECT_SOURCE_DIR}/src/system.cpp
            ${PROJECT_SOURCE_DIR}/src/timer.cpp
            ${PROJECT_SOURCE_DIR}/src/writer.cpp)

# # Trick to use gcc to compile *.cpp: This doesn't work because of boost
# SET_SOURCE_FILES_PROPERTIES(${SOURCES} PROPERTIES LANGUAGE C)
# set(CMAKE_C_COMPILER "gcc")

# Executable
add_executable(alm ${PROJECT_SOURCE_DIR}/src/main.cpp
                   ${PROJECT_SOURCE_DIR}/src/alm_cui.cpp
		   ${PROJECT_SOURCE_DIR}/src/input_parser.cpp
		   ${PROJECT_SOURCE_DIR}/src/input_setter.cpp
		   ${SOURCES})
target_link_libraries(alm ${Boost_LIBRARIES} lapack)

# add_executable(alm ${PROJECT_SOURCE_DIR}/src/main.cpp)
# target_link_libraries(alm ${almcxx_static} ${Boost_LIBRARIES} lapack)

# Shared library
add_library(almcxx SHARED ${SOURCES})
target_link_libraries(almcxx ${Boost_LIBRARIES} lapack)
set_property(TARGET almcxx PROPERTY VERSION ${serial})
set_property(TARGET almcxx PROPERTY SOVERSION ${soserial})
#install(TARGETS almcxx LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib)
install(TARGETS almcxx DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Static link library
add_library(almcxx_static STATIC ${SOURCES})
set_property(TARGET almcxx_static PROPERTY VERSION ${serial})
set_property(TARGET almcxx_static PROPERTY SOVERSION ${soserial})
set_property(TARGET almcxx_static PROPERTY OUTPUT_NAME almcxx)
install(TARGETS almcxx_static ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Header file
install(FILES ${PROJECT_SOURCE_DIR}/src/alm.h DESTINATION ${PROJECT_SOURCE_DIR}/include)

