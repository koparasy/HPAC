include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(SRC_FILES "hpccg_entire.cpp")

set(ACCURATE_BIN hpccg_entire_accurate)
set(ANALYSIS_BIN hpccg_entire_analysis)

add_executable(${ANALYSIS_BIN} ${SRC_FILES})
target_link_libraries(${ANALYSIS_BIN} ${LIBAPPROX_LIB_FILE})
target_include_directories(${ANALYSIS_BIN} PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(${ANALYSIS_BIN} PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fapprox -fopenmp -DWITH_QUANT)
set_target_properties(${ANALYSIS_BIN} PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

add_executable(${ACCURATE_BIN} ${SRC_FILES})
target_link_libraries(${ACCURATE_BIN} ${LIBAPPROX_LIB_FILE})
target_include_directories(${ACCURATE_BIN} PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(${ACCURATE_BIN} PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp -DWITH_QUANT)
set_target_properties(${ACCURATE_BIN} PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

configure_file(hpccg_entire.yaml.in hpccg_entire.yaml)

set(SRC_FILES "hpccg_slice.cpp")

set(ACCURATE_BIN hpccg_slice_accurate)
set(ANALYSIS_BIN hpccg_slice_analysis)

add_executable(${ANALYSIS_BIN} ${SRC_FILES})
target_link_libraries(${ANALYSIS_BIN} ${LIBAPPROX_LIB_FILE})
target_include_directories(${ANALYSIS_BIN} PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(${ANALYSIS_BIN} PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fapprox -fopenmp -DWITH_QUANT)
set_target_properties(${ANALYSIS_BIN} PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

add_executable(${ACCURATE_BIN} ${SRC_FILES})
target_link_libraries(${ACCURATE_BIN} ${LIBAPPROX_LIB_FILE})
target_include_directories(${ACCURATE_BIN} PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(${ACCURATE_BIN} PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp -DWITH_QUANT)
set_target_properties(${ACCURATE_BIN} PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

configure_file(hpccg_slice.yaml.in hpccg_slice.yaml)

set(SRC_FILES "hpccg_eval.cpp")

set(EVAL_BIN hpccg_eval)

add_executable(${EVAL_BIN} ${SRC_FILES})
target_link_libraries(${EVAL_BIN} ${LIBAPPROX_LIB_FILE})
target_include_directories(${EVAL_BIN} PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(${EVAL_BIN} PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fapprox -fopenmp -DWITH_QUANT)
set_target_properties(${EVAL_BIN} PROPERTIES LINK_FLAGS "-fopenmp -lapprox")
