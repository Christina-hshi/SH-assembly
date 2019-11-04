include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

add_library(cqf-core-objs OBJECT
	gqf.c
	true2falseKmer_DP.cpp
)	

set(cqf-core-objs $<TARGET_OBJECTS:cqf-core-objs> PARENT_SCOPE)

add_library(cqf-core-lib STATIC $<TARGET_OBJECTS:cqf-core-objs>)

#message(${CMAKE_CURRENT_LIST_DIR}/../base)
#target_include_directories(cqf-core-lib PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)

