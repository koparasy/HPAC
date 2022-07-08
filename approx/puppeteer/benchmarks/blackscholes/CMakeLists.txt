
include_directories(${CMAKE_CURRENT_BINARY_DIR})


set(SRC_FILES "blackscholes.cpp")
set(APP_NAME "blackscholes")

set(ACCURATE_BIN blck_accurate)
set(ANALYSIS_BIN blck_analysis)

add_executable(blck_analysis ${SRC_FILES})
target_link_libraries(blck_analysis ${LIBAPPROX_LIB_FILE})
target_include_directories(blck_analysis PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(blck_analysis PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fapprox )
set_target_properties(blck_analysis PROPERTIES LINK_FLAGS " -lapprox")

add_executable(blck_accurate ${SRC_FILES})
target_link_libraries(blck_accurate ${LIBAPPROX_LIB_FILE})
target_include_directories(blck_accurate PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(blck_accurate PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native )
set_target_properties(blck_accurate PROPERTIES LINK_FLAGS " -lapprox")

add_executable(blck_approx blackscholes_eval.cpp)
target_link_libraries(blck_approx ${LIBAPPROX_LIB_FILE})
target_include_directories(blck_approx PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(blck_approx PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native )
set_target_properties(blck_approx PROPERTIES LINK_FLAGS " -lapprox")

configure_file(blackscholes.yaml.in blackscholes.yaml)
configure_file(blackscholes.sh.in blackscholes_analysis.sh)

