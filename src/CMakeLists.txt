include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

#gkm-lib
set(SH-lib core-lib base-lib)

add_executable(Contiger contig_assembly.cpp)
target_link_libraries(Contiger core-lib cqf-core-lib base-lib z bz2 boost_system boost_thread boost_iostreams boost_program_options pthread rt boost_timer boost_chrono tbb)

add_executable(CQF-deNoise CQF-deNoise.cpp)
target_link_libraries(CQF-deNoise cqf-core-lib base-lib z bz2 boost_system boost_thread boost_iostreams boost_program_options pthread rt boost_timer boost_chrono)

#add_executable(compare_contigs compare_contigs.cpp)
#target_link_libraries(compare_contigs cqf-core-lib base-lib boost_program_options)

#	#gkm-kernel
#	add_executable(gkm-kernel main_calKernel.cpp)
#	target_link_libraries(gkm-kernel ${SH-lib})
#
#	#test 
#	add_executable(BUSY test.cpp)
#	target_link_libraries(BUSY ${SH-lib})
#
#	#generate test data
#	add_executable(GEN generate_data.cpp)
#	target_link_libraries(GEN ${SH-lib})
#
#	#combine all libraries
#	add_library(combinedLib INTERFACE)
#	target_link_libraries(combinedLib INTERFACE contig-lib base-lib)
#
#	#perform contiging
#	add_executable(CONTIG contiging.cpp)
#	target_link_libraries(CONTIG combinedLib)
#	#target_link_libraries(CONTIG contig-lib)
#
#	#set_properties(TARGET gkm-kernel PROPERTIES CXX_STANDARD 11 COMPILE_FLAGS "-pthread")
#	#set_target_properties(gkm-kernel PROPERTIES COMPILE_OPTIONS "-pthread")
#
#	#test each sub-module
#	add_executable(DNAStrInBits_test DNAStrInBits_test.cpp)
#	target_link_libraries(DNAStrInBits_test base-lib)
#
#	add_executable(DeBruijnGraph_test DeBruijnGraph_test.cpp)
#	target_link_libraries(DeBruijnGraph_test combineLib)
#
#	add_executable(UnitigGraph_test UnitigGraph_test.cpp)
#	target_link_libraries(UnitigGraph_test contig-lib base-lib)
#
#	add_executable(splitHub_test splitHub_test.cpp)
#	target_link_libraries(splitHub_test contig-lib base-lib)
#
#	add_executable(CQF_test CQF_test.cpp)
#	target_link_libraries(CQF_test cqf-lib base-lib)
#
#	add_executable(CQF_mt_test CQF_mt_test.cpp)
#	target_link_libraries(CQF_mt_test cqf-mt-lib base-lib z bz2 boost_system-mt boost_thread-mt)
#
#	add_executable(CQF_mt_FC_test CQF_mt_FC_test.cpp)
#	target_link_libraries(CQF_mt_FC_test cqf-mt-FC-lib base-lib z bz2 boost_system-mt boost_thread-mt)
#
#	add_executable(CQF_contiger CQF_contiger.cpp)
#	target_link_libraries(CQF_contiger cqf-lib base-lib)
#
#	add_executable(CQF_mt_contiger CQF_mt_contiger.cpp)
#	target_link_libraries(CQF_mt_contiger cqf-mt-lib base-lib z bz2 boost_system-mt boost_thread-mt boost_iostreams)
#
#	add_executable(unordered_map_mt_test unordered_map_mt_test.cpp)
#	target_link_libraries(unordered_map_mt_test base-lib boost_thread-mt)
#
#	add_executable(vector_mt_test vector_mt_test.cpp)
#	target_link_libraries(vector_mt_test base-lib boost_thread-mt)
#
#	add_executable(SAM_test SAM_test.cpp)
#	target_link_libraries(SAM_test base-lib)
#
#	add_executable(try try.cpp)
#	target_link_libraries(try base-lib)
