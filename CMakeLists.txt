cmake_minimum_required(VERSION 3.1)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O2")
# set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELEASE)

project(alm)
#set(CMAKE_CXX_COMPILER "g++")
option(WITH_SPARSE_SOLVER "Use sparse solver option" ON)

if(WIN32)
    LINK_DIRECTORIES(C:\\lapack)
	LINK_DIRECTORIES(C:\\spglib)
    include_directories(C:\\spglib)
	set(LAPACK_LIBRARIES lapack)
	set(spglib symspg)
endif()

if (UNIX)
    #set( ENV{BLA_VENDOR} "Intel10_64lp")
	include_directories("/Users/tadano/src/spglib/include")
    find_package(LAPACK REQUIRED)
	set(spglib "-L/Users/tadano/src/spglib/lib -lsymspg")
endif()

#
# Add openmp flag if available 
# 
find_package(OpenMP)
if (OPENMP_FOUND) 
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

if (WITH_SPARSE_SOLVER)
    add_definitions(-DWITH_SPARSE_SOLVER)
    find_package(Eigen3 REQUIRED)
    include_directories(${EIGEN3_INCLUDE_DIR})
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
            ${PROJECT_SOURCE_DIR}/src/constraint.cpp
            ${PROJECT_SOURCE_DIR}/src/fcs.cpp
            ${PROJECT_SOURCE_DIR}/src/files.cpp
            ${PROJECT_SOURCE_DIR}/src/fitting.cpp
            ${PROJECT_SOURCE_DIR}/src/input_parser.cpp
            ${PROJECT_SOURCE_DIR}/src/input_setter.cpp
            ${PROJECT_SOURCE_DIR}/src/interaction.cpp
            ${PROJECT_SOURCE_DIR}/src/lasso.cpp
            ${PROJECT_SOURCE_DIR}/src/patterndisp.cpp
            ${PROJECT_SOURCE_DIR}/src/rref.cpp
            ${PROJECT_SOURCE_DIR}/src/symmetry.cpp
            ${PROJECT_SOURCE_DIR}/src/system.cpp
            ${PROJECT_SOURCE_DIR}/src/timer.cpp
            ${PROJECT_SOURCE_DIR}/src/writer.cpp)

set(HEADERS ${PROJECT_SOURCE_DIR}/src/alm.h
            ${PROJECT_SOURCE_DIR}/src/constraint.h
            ${PROJECT_SOURCE_DIR}/src/fcs.h
            ${PROJECT_SOURCE_DIR}/src/files.h
            ${PROJECT_SOURCE_DIR}/src/fitting.h
            ${PROJECT_SOURCE_DIR}/src/input_parser.h
            ${PROJECT_SOURCE_DIR}/src/input_setter.h
            ${PROJECT_SOURCE_DIR}/src/interaction.h
            ${PROJECT_SOURCE_DIR}/src/lasso.h
            ${PROJECT_SOURCE_DIR}/src/patterndisp.h
            ${PROJECT_SOURCE_DIR}/src/rref.h
            ${PROJECT_SOURCE_DIR}/src/symmetry.h
            ${PROJECT_SOURCE_DIR}/src/system.h
            ${PROJECT_SOURCE_DIR}/src/timer.h
            ${PROJECT_SOURCE_DIR}/src/writer.h)

# # Trick to use gcc to compile *.cpp: This doesn't work because of boost
# SET_SOURCE_FILES_PROPERTIES(${SOURCES} PROPERTIES LANGUAGE C)
# set(CMAKE_C_COMPILER "gcc")

# Executable
add_executable(alm ${PROJECT_SOURCE_DIR}/src/main.cpp
                   ${PROJECT_SOURCE_DIR}/src/alm_cui.cpp
		   ${SOURCES} ${HEADERS})
target_link_libraries(alm ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${spglib})
target_link_libraries(alm ${LAPACK_LIBRARIES})
set_property(TARGET alm PROPERTY CXX_STANDARD 11)
set_property(TARGET alm PROPERTY CXX_STANDARD_REQUIRED ON)

# add_executable(alm ${PROJECT_SOURCE_DIR}/src/main.cpp)
# target_link_libraries(alm ${almcxx_static} ${Boost_LIBRARIES} lapack)

# Shared library
add_library(almcxx SHARED ${SOURCES})
#target_link_libraries(almcxx ${Boost_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(almcxx ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${spglib})
#target_link_libraries(almcxx ${LAPACK_LIBRARIES})
set_property(TARGET almcxx PROPERTY VERSION ${serial})
set_property(TARGET almcxx PROPERTY SOVERSION ${soserial})
set_property(TARGET almcxx PROPERTY CXX_STANDARD 11)
set_property(TARGET almcxx PROPERTY CXX_STANDARD_REQUIRED ON)
#install(TARGETS almcxx LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib)
install(TARGETS almcxx DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Static link library
add_library(almcxx_static STATIC ${SOURCES})
set_property(TARGET almcxx_static PROPERTY VERSION ${serial})
set_property(TARGET almcxx_static PROPERTY SOVERSION ${soserial})
set_property(TARGET almcxx_static PROPERTY OUTPUT_NAME almcxx)
set_property(TARGET almcxx_static PROPERTY CXX_STANDARD 11)
set_property(TARGET almcxx_static PROPERTY CXX_STANDARD_REQUIRED ON)
install(TARGETS almcxx_static ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# Header file
install(FILES ${PROJECT_SOURCE_DIR}/src/alm.h DESTINATION ${PROJECT_SOURCE_DIR}/include)

