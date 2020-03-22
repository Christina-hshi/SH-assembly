include_directories(${PROJECT_SOURCE_DIR})

# add the binary tree to the search path for include files
# so that we will find Config.h
include_directories(${PROJECT_BINARY_DIR})

add_library(core-objs OBJECT
	unitig_graph.cpp
	)

#set_properties(TARGET core-objs PROPERTIES CXX_STANDAR 11 COMPILE_FLAGS "-pthread")
#set_target_properties(core-objs PROPERTIES COMPILE_OPTIONS "-pthread")

# Visible to parent directory
set(core-objs $<TARGET_OBJECTS:core-objs> PARENT_SCOPE)

add_library(core-lib STATIC $<TARGET_OBJECTS:core-objs>)
#target_include_directories(core-lib PUBLIC ${CMAKE_CURRENT_LIST_DIR})
