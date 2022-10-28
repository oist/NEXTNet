cmake_minimum_required(VERSION 3.16)

project(Epidemics-On-Networks LANGUAGES CXX)

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
    epidemics/algorithm.h
    epidemics/NextReaction.h
    epidemics/nMGA.h
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
    epidemics/tests/plot.cpp
    epidemics/tests/main.cpp
    epidemics/tests/graph.cpp
    epidemics/tests/random.cpp
    epidemics/tests/nextreaction.cpp
    epidemics/tests/nmga.cpp)

set(TESTS_HEADERS
    epidemics/tests/simulate.h
    epidemics/tests/analytical.h
    epidemics/tests/parallel.h
    epidemics/tests/plot.h
    epidemics/tests/gnuplot-iostream.h
    epidemics/catch/catch.hpp)

set(TESTS_PCH
    epidemics/tests/stdafx.h)

# *** Global build settings

# Includes are relative to the epidemics/ folder
include_directories(AFTER SYSTEM epidemics)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hide symbols
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN true)

# Enable all warnings
add_compile_options(-Wall)

# *** Dependencies

# Boost headers are mandatory, iostreamd and filesystem are optional
# and allow plots to be generated during test runs
find_package(Boost REQUIRED OPTIONAL_COMPONENTS iostreams filesystem)
link_libraries(Boost::boost)
#include_directories(${Boost_INCLUDE_DIRS})

# Intel Thread Building Blocks are optional, if present
# they speed up the tests by executing simulations in parallel
find_package(TBB)

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
if (TBB_FOUND)
    message("Thread Building Blocks: Enabled")
    target_link_libraries(tests TBB::tbb)
else()
    message("Thread Building Blocks: Disabled")
endif()
if (Boost_FILESYSTEM_FOUND AND Boost_IOSTREAMS_FOUND)
    message("Test plots: Enabled")
    target_link_libraries(tests Boost::boost Boost::iostreams Boost::filesystem)
    target_compile_definitions(tests PRIVATE ENABLE_PLOTTING=1)
else()
    message("Test plots: Disabled")
endif()

# Enable address and undefined behaviour sanitizers for all builds
#target_compile_options(tests PRIVATE
#    -fsanitize=address
#    -fsanitize=undefined)
#target_link_options(tests PRIVATE
#    -fsanitize=address
#    -fsanitize=undefined)
