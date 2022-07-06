
include_directories(${CMAKE_CURRENT_BINARY_DIR})


set(SRC_FILES "DCT_IDCT.cpp")

set(ACCURATE_BIN dct_accurate_wq)
set(ANALYSIS_BIN dct_analysis_wq)

add_executable(dct_analysis_wq ${SRC_FILES})
target_link_libraries(dct_analysis_wq ${LIBAPPROX_LIB_FILE})
target_include_directories(dct_analysis_wq PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(dct_analysis_wq PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fapprox -fopenmp -DWITH_QUANT)
set_target_properties(dct_analysis_wq PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

add_executable(dct_accurate_wq ${SRC_FILES})
target_link_libraries(dct_accurate_wq ${LIBAPPROX_LIB_FILE})
target_include_directories(dct_accurate_wq PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(dct_accurate_wq PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp -DWITH_QUANT)
set_target_properties(dct_accurate_wq PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

add_executable(dct_approx_wq DCT_IDCT_eval.cpp)
target_link_libraries(dct_approx_wq ${LIBAPPROX_LIB_FILE})
target_include_directories(dct_approx_wq PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(dct_approx_wq PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp -DWITH_QUANT)
set_target_properties(dct_approx_wq PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

configure_file(dct_wq.yaml.in dct_wq.yaml)

set(ACCURATE_BIN dct_accurate_nq)
set(ANALYSIS_BIN dct_analysis_nq)

add_executable(dct_analysis_nq ${SRC_FILES})
target_link_libraries(dct_analysis_nq ${LIBAPPROX_LIB_FILE})
target_include_directories(dct_analysis_nq PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(dct_analysis_nq PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fapprox -fopenmp)
set_target_properties(dct_analysis_nq PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

add_executable(dct_accurate_nq ${SRC_FILES})
target_link_libraries(dct_accurate_nq ${LIBAPPROX_LIB_FILE})
target_include_directories(dct_accurate_nq PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(dct_accurate_nq PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp)
set_target_properties(dct_accurate_nq PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

add_executable(dct_approx_nq DCT_IDCT_eval.cpp)
target_link_libraries(dct_approx_nq ${LIBAPPROX_LIB_FILE})
target_include_directories(dct_approx_nq PUBLIC ${APPROX_LIB_UTIL_DIR}/)
target_compile_options(dct_approx_nq PUBLIC ${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp)
set_target_properties(dct_approx_nq PROPERTIES LINK_FLAGS "-fopenmp -lapprox")

configure_file(dct_nq.yaml.in dct_nq.yaml)
