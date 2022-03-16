cmake_minimum_required(VERSION 3.5)

project(Epidemics-Mean-Field LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIB_SOURCES
    epidemics/graph.cpp
    epidemics/random.cpp
    epidemics/simulation.cpp)

set(LIB_HEADERS
    epidemics/types.h
    epidemics/graph.h
    epidemics/random.h
    epidemics/simulation.h)

SET(LIB_PCH epidemics/stdafx.h)

set(SIMULATOR_SOURCES
    epidemics/analysis.cpp
    epidemics/main.cpp)

set(SIMULATOR_HEADERS
    epidemics/analysis.h)

set(TESTS_SOURCES
    epidemics/analysis.cpp
    epidemics/tests/main.cpp
    epidemics/tests/graph.cpp)

set(TESTS_HEADERS
    epidemics/catch/catch.hpp)

set(TESTS_PCH
    epidemics/tests/stdafx.h)

include_directories(AFTER SYSTEM epidemics)

add_executable(simulator
    ${SIMULATOR_SOURCES}
    ${LIB_SOURCES}
    ${SIMULATOR_HEADERS}
    ${LIB_HEADERS}
    ${LIB_PCH})
target_precompile_headers(simulator
    PUBLIC ${LIB_PCH})

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
