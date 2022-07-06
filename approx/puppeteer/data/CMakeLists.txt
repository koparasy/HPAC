
find_package(Python3 REQUIRED COMPONENTS Interpreter)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

set(OUT_NAME "figure_10")
file(GLOB LOGS perf_eval/adapt_vs_puppet/*.log)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/perf_eval)
message("${Magenta}Emitting ${CMAKE_CURRENT_BINARY_DIR}/perf_eval/${OUT_NAME} ${ColourReset}")
execute_process(COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/perf_eval/adapt_vs_puppet/generate_figure.py
  ${CMAKE_CURRENT_BINARY_DIR}/perf_eval/ ${OUT_NAME} ${LOGS})

set(OUT_NAME "figure_11")
file(GLOB LOGS perf_eval/thread_scaling_data/*.log)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/perf_eval)
message("${Magenta}Emitting ${CMAKE_CURRENT_BINARY_DIR}/perf_eval/${OUT_NAME} ${ColourReset}")
execute_process(COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/perf_eval/thread_scaling_data/scaling.py
  ${CMAKE_CURRENT_BINARY_DIR}/perf_eval/ ${OUT_NAME} ${LOGS})

set(OUT_NAME "figure_12")
file(GLOB LOGS perf_eval/hpccg_memory_overhead/*_memory.log)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/perf_eval)
message("${Magenta}Emitting ${CMAKE_CURRENT_BINARY_DIR}/perf_eval/${OUT_NAME} ${ColourReset}")
execute_process(COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/perf_eval/hpccg_memory_overhead/generate_figure.py
  ${CMAKE_CURRENT_BINARY_DIR}/perf_eval/ ${OUT_NAME} ${LOGS})


