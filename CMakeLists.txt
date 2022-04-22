cmake_minimum_required(VERSION 3.16)

project(Epidemics-Mean-Field LANGUAGES CXX)

# *** Simulation library sources and headers
set(LIB_SOURCES
    epidemics/graph.cpp
    epidemics/random.cpp
    epidemics/simulation.cpp
    epidemics/NextReaction.cpp 
    epidemics/nMGA.cpp
    epidemics/utility.cpp)

set(LIB_HEADERS
    epidemics/types.h
    epidemics/graph.h
    epidemics/random.h
    epidemics/simulation.h
    epidemics/utility.h)

SET(LIB_PCH epidemics/stdafx.h)

# *** Command-line simulator sources and headers
set(SIMULATOR_SOURCES
    epidemics/analysis.cpp
    epidemics/main.cpp)

set(SIMULATOR_HEADERS
    epidemics/analysis.h)

# *** Unit tests sources and headers
set(TESTS_SOURCES
    epidemics/analysis.cpp
    epidemics/tests/main.cpp
    epidemics/tests/graph.cpp
    epidemics/tests/random.cpp
    epidemics/tests/nextreaction.cpp)

set(TESTS_HEADERS
    epidemics/catch/catch.hpp)

set(TESTS_PCH
    epidemics/tests/stdafx.h)

# *** Global build settings

# Includes are relative to the epidemics/ folder
include_directories(AFTER SYSTEM epidemics)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable all warnings
add_compile_options(-Wall)

# Dependencies
find_package(Boost REQUIRED)
link_libraries(Boost::boost)

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
link_libraries(Python3::Python Python3::NumPy Python3::Module)
add_compile_definitions(ENABLE_PLOTTING=1)

#include_directories(${Boost_INCLUDE_DIRS})

# *** Build command-line similator
add_executable(simulator
    ${SIMULATOR_SOURCES}
    ${LIB_SOURCES}
    ${SIMULATOR_HEADERS}
    ${LIB_HEADERS}
    ${LIB_PCH})
target_precompile_headers(simulator
    PUBLIC ${LIB_PCH})

# Enable address and undefined behaviour sanitizers for debug builds
target_compile_options(simulator PRIVATE
    "$<$<CONFIG:DEBUG>:-fsanitize=address>"
    "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
target_link_options(simulator PRIVATE
    "$<$<CONFIG:DEBUG>:-fsanitize=address>"
    "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")

# *** Build unit test runner
add_executable(tests
    ${TESTS_SOURCES}
    ${LIB_SOURCES}
    ${TESTS_HEADERS}
    ${TESTS_PCH}
    ${LIB_HEADERS}
    ${LIB_PCH})
target_precompile_headers(tests
    PUBLIC ${LIB_PCH}
    PUBLIC ${TESTS_PCH})
set_source_files_properties(epidemics/tests/main.cpp
    PROPERTIES SKIP_PRECOMPILE_HEADERS ON)

# Enable address and undefined behaviour sanitizers for all builds
target_compile_options(tests PRIVATE
    -fsanitize=address
    -fsanitize=undefined)
target_link_options(tests PRIVATE
    -fsanitize=address
    -fsanitize=undefined)
