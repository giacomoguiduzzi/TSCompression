cmake_minimum_required(VERSION 3.30.0)
project(TSCompression)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type (default to Debug if not specified)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Include custom modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

# Project metadata
set(PROJECT_VENDOR "Giacomo Guiduzzi")
set(CMAKE_PROJECT_NAME "TSCompression")
set(PROJECT_CONTACT "giacomo.guiduzzi@unimore.it")

# SSE4.2 Check
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-msse4.2" SSE4_2_SUPPORTED)
if(SSE4_2_SUPPORTED)
    add_compile_options(-msse4.2)
    message(STATUS "CPU does support SSE4.2.")
else()
    message(STATUS "CPU does NOT support SSE4.2.")
endif()

# Add subdirectories
add_subdirectory(rlz/external/sdsl-lite-custom)

message(STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")

# Include directories
include_directories(
        ${CMAKE_SOURCE_DIR}/rlz
        ${CMAKE_SOURCE_DIR}/rlz/external/googletest/include
        ${CMAKE_SOURCE_DIR}/rlz/include
        ${CMAKE_BINARY_DIR}/rlz/external/sdsl-lite-custom/include
        ${CMAKE_BINARY_DIR}/rlz/external/sdsl-lite-custom/external/libdivsufsort/include
        ${CMAKE_SOURCE_DIR}/rlz/external/googletest/include
        ${CMAKE_SOURCE_DIR}/rlz/external/sdsl-lite-custom/external/FastPFor/headers
        ${STXXL_INCLUDE_DIRS}
)

# Add executable and link libraries
add_executable(TSCompression ts_compression.cpp
        ts_compression.h)
add_executable(compress_with_naive compress_with_naive.cpp ts_compression.h)
target_link_libraries(TSCompression sdsl divsufsort divsufsort64)
target_link_libraries(compress_with_naive sdsl divsufsort divsufsort64)