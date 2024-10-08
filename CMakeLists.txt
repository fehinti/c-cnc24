# Minimal settings
cmake_minimum_required(VERSION 3.15)
set(VERSION "1.1.1")
project(C-CNC VERSION ${VERSION} LANGUAGES C)

# Build mode
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(NOTICE "Debug build, enabling all warnings")
  add_compile_options(-Wall)
  add_compile_definitions(DEBUG)
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "Release build")
else()
  message(FATAL_ERROR "Please set CMAKE_BUILD_TYPE to either Debug or Release!")
endif()

# Build settings
set(CMAKE_C_STANDARD 11)
set(SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
set(MAIN_DIR ${SOURCE_DIR}/main)
add_compile_definitions(_GNU_SOURCE)
file(GLOB LIB_SOURCES ${SOURCE_DIR}/*.c)
list(FILTER LIB_SOURCES EXCLUDE REGEX "block.c")
list(FILTER LIB_SOURCES EXCLUDE REGEX "program.c")
message(STATUS "Library source files: ${LIB_SOURCES}")
include_directories(/usr/local/include)
link_directories(/usr/local/lib)


# Configured files
configure_file(${SOURCE_DIR}/defines.h.in ${SOURCE_DIR}/defines.h)

# Target

# static library
add_library(ccnc_lib STATIC ${LIB_SOURCES})

# Test executables
add_executable(point_test ${SOURCE_DIR}/point.c)
target_compile_definitions(point_test PUBLIC POINT_MAIN)
target_link_libraries(point_test m mosquitto)

add_executable(machine_test ${LIB_SOURCES})
target_compile_definitions(machine_test PUBLIC MACHINE_MAIN)
target_link_libraries(machine_test m mosquitto)

add_executable(block_test ${LIB_SOURCES})
target_compile_definitions(block_test PUBLIC BLOCK_MAIN)
target_link_libraries(block_test m mosquitto)

add_executable(program_test ${LIB_SOURCES})
target_compile_definitions(program_test PUBLIC PROGRAM_MAIN)
target_link_libraries(program_test m mosquitto)

add_executable(ccnc ${MAIN_DIR}/ccnc.c)
target_link_libraries(ccnc ccnc_lib m mosquitto)