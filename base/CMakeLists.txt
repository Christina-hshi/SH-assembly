include_directories(${PROJECT_SOURCE_DIR})

# add the binary tree to the search path for include files
# so that we will find Config.h
include_directories(${PROJECT_BINARY_DIR})

add_library(base-objs OBJECT
	Utility.cpp
	global.cpp
#	Hash.cpp
#	MurmurHash2.cpp
#	Fasta.cpp
#	Sequence.cpp	
#	Contig.cpp
#	Assemblier.cpp
	multithread_io.cpp
#	io.cpp
#	Tree_parallel.cpp
#	DNAStrInBits.cpp
	)
#	Tree.cpp

#set_properties(TARGET gkm-objs PROPERTIES CXX_STANDAR 11 COMPILE_FLAGS "-pthread")
#set_target_properties(gkm-objs PROPERTIES COMPILE_OPTIONS "-pthread")

# Visible to parent directory
set(base-objs $<TARGET_OBJECTS:base-objs> PARENT_SCOPE)

add_library(base-lib STATIC $<TARGET_OBJECTS:base-objs>)
target_include_directories(base-lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
