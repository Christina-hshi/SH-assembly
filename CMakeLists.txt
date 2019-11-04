cmake_minimum_required(VERSION 3.0.2)

set(CMAKE_C_COMPILER "g++")
set(CMAKE_CXX_COMPILER "g++")

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_FLAGS "-std=gnu++11")
#SET(CMAKE_EXE_LINKER_FLAGS "-Ofast")
#set(COMPILER_FLAGS "-Wall -Ofast -m64 -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare")
set(COMPILER_FLAGS "-Wall -m64 -Wno-unused-result -Wno-strict-aliasing -Wno-unused-function -Wno-sign-compare")

#set(CMAKE_CXX_FLAGS "${COMPILER_FLAGS} -pg -pthread -g")
set(CMAKE_CXX_FLAGS_DEBUG "${COMPILER_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${COMPILER_FLAGS} -Ofast")

#set(CMAKE_C_FLAGS "-std=gnu11 ${COMPILER_FLAGS} -g -pg")
set(CMAKE_C_FLAGS_DEBUG "${COMPILER_FLAGS} -g")
set(CMAKE_C_FLAGS_RELEASE "${COMPILER_FLAGS} -Ofast")

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/Config.h.in"
#  "${PROJECT_BINARY_DIR}/Config.h"
#	)
# add the binary tree to the search path for include files
# so that we will find Config.h
include_directories(${PROJECT_BINARY_DIR})

# external libraries
include_directories($ENV{WD}Tools/miniconda2/include)
link_directories($ENV{WD}Tools/miniconda2/lib)
link_directories(/usr/lib64)

# check which compiler is used
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	# using Clang
	option(COMPILER_CLANG "using Clang to compile" ON)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	# using GCC
	option(COMPILER_GNU "using GNC to compile" ON)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	# using Intel C++
	option(COMPILER_INTEL "using INTEL complier" ON)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	# using Visual Studio C++
	option(COMPILER_MSVC "using MSVC to compile" ON)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

#add_subdirectory(gkm)
add_subdirectory(src)
add_subdirectory(base)
add_subdirectory(cqf)
#add_subdirectory(CQF-deNoise)
#add_subdirectory(contig)
#add_subdirectory(cqf)
#add_subdirectory(cqf_mt)
#add_subdirectory(cqf_mt_FC)
#add_subdirectory(scaffold)

#set(SH-objs $<TARGET_OBJECTS:gkm-objs> $<TARGET_OBJECTS:base-objs>)
#set(SH-objs $<TARGET_OBJECTS:base-objs>)
#add_library(SH-lib STATIC ${SH-objs})

set(SH-lib base-lib)

