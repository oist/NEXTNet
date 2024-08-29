cmake_minimum_required(VERSION 3.16)

project(Epidemics_On_Networks LANGUAGES CXX)

# *** Simulation library sources and headers
set(LIB_SOURCES
    epidemics/graph.cpp
    epidemics/dynamic_graph.cpp
    epidemics/brownian_proximity_graph.cpp
    epidemics/algorithm.cpp
    epidemics/random.cpp
    epidemics/NextReaction.cpp
    epidemics/NextReactionMeanField.cpp
    epidemics/nMGA.cpp
    epidemics/REGIR.cpp
    epidemics/utility.cpp)

set(LIB_HEADERS
    epidemics/types.h
    epidemics/graph.h
    epidemics/dynamic_graph.h
    epidemics/brownian_proximity_graph.h
    epidemics/random.h
    epidemics/algorithm.h
    epidemics/NextReaction.h
    epidemics/NextReactionMeanField.h
    epidemics/nMGA.h
    epidemics/REGIR.h
    epidemics/permutation.h
    epidemics/utility.h)

SET(LIB_PCH epidemics/stdafx.h)

# *** Command-line simulator sources and headers
set(SIMULATOR_SOURCES
    epidemics/analysis.cpp
    epidemics/main.cpp
    epidemics/programs/average.cpp
    epidemics/programs/scale_free.cpp
    epidemics/programs/benchmark.cpp
    epidemics/programs/trajectory.cpp
    epidemics/programs/customtrajectory.cpp
    epidemics/programs/sis_meanfield_gamma_gamma.cpp
    epidemics/programs/profile.cpp
    epidemics/programs/benchmark_priorityqueue.cpp
    epidemics/programs/custered_configuration_model_serrano.cpp)

set(SIMULATOR_HEADERS
    epidemics/analysis.h)

# *** Unit tests sources and headers
set(TESTS_SOURCES
    epidemics/analysis.cpp
    epidemics/tests/plot.cpp
    epidemics/tests/main.cpp
    epidemics/tests/utility.cpp
    epidemics/tests/graph.cpp
    epidemics/tests/dynamic_graph.cpp
    epidemics/tests/random.cpp
    epidemics/tests/nextreaction.cpp
    epidemics/tests/dynamic_nextreaction.cpp
    epidemics/tests/brownian_proximity_graph.cpp
    epidemics/tests/nmga.cpp
    epidemics/tests/REGIR.cpp)

set(TESTS_HEADERS
    epidemics/tests/statistics.h
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
include_directories(AFTER SYSTEM epidemics ext/dyndist)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hide symbols
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN true)

# Enable all warnings
add_compile_options(-Wall)

add_compile_definitions(${CMAKE_DEFINES})

# *** Dependencies

# Boost headers are mandatory, iostreamd and filesystem are optional
# and allow plots to be generated during test runs
find_package(Boost REQUIRED OPTIONAL_COMPONENTS iostreams filesystem)
link_libraries(Boost::boost)
include_directories(${Boost_INCLUDE_DIRS})
if (Boost_FILESYSTEM_FOUND)
    message("Boost filesystem: Found")
else()
    message("Boost filesystem: Not Found")
endif()
if (Boost_IOSTREAMS_FOUND)
    message("Boost iostreams: Found")
else()
    message("Boost iostreams: Not Found")
endif()

# Intel Thread Building Blocks are optional, if present
# they speed up the tests by executing simulations in parallel
if (NOT SERIALIZE)
    find_package(TBB)
    if (TBB_FOUND)
        message("Thread Building Blocks: Found")
        add_compile_definitions(TBB_SUPPRESS_DEPRECATED_MESSAGES)
    else()
        message("Thread Building Blocks: Not Found")
    endif()
endif()

# Gnuplot is optional, if present can generate plos
find_program(Gnuplot gnuplot)
if (Gnuplot)
    message("Gnuplot: Found (${Gnuplot})")
else()
    message("Gnuplot: Not Found")
endif()

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
if (SANITIZE STREQUAL "debug")
    message("Address and thread sanitizer during simulations: Debug builds")
    target_compile_options(simulator PRIVATE
        "$<$<CONFIG:DEBUG>:-fsanitize=address>"
        "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
    target_link_options(simulator PRIVATE
        "$<$<CONFIG:DEBUG>:-fsanitize=address>"
        "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
elseif (SANITIZE STREQUAL "always")
    message("Address and thread sanitizer during simulations: Always")
    target_compile_options(simulator PRIVATE
        "-fsanitize=address"
        "-fsanitize=undefined")
    target_link_options(simulator PRIVATE
        "-fsanitize=address"
        "-fsanitize=undefined")
else()
    message("Address and thread sanitizer during simulations: Disabled")
endif()

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
if ((NOT SERIALIZE) AND TBB_FOUND)
    message("Parallelization: Enabled")
    target_link_libraries(tests TBB::tbb)
    add_definitions(-DPARALLELIZE=1)
else()
    message("Parallelization: Disabled")
    add_definitions(-DPARALLELIZE=0)
endif()
if (Gnuplot AND Boost_FILESYSTEM_FOUND AND Boost_IOSTREAMS_FOUND)
    message("Plotting during tests: Enabled")
    target_link_libraries(tests Boost::boost Boost::iostreams Boost::filesystem)
    target_compile_definitions(tests PRIVATE
        ENABLE_PLOTTING=1
        GNUPLOT_DEFAULT_COMMAND="${Gnuplot} -persist"
     )
else()
    message("Plotting during tests: Disabled")
endif()

# Enable address and undefined behaviour sanitizers for all builds
if (SANITIZE STREQUAL "debug")
    message("Address and thread sanitizer during tests: Debug builds")
    target_compile_options(tests PRIVATE
        "$<$<CONFIG:DEBUG>:-fsanitize=address>"
        "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
    target_link_options(tests PRIVATE
        "$<$<CONFIG:DEBUG>:-fsanitize=address>"
        "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
elseif (SANITIZE STREQUAL "always")
    message("Address and thread sanitizer during tests: Always")
    target_compile_options(tests PRIVATE
        "-fsanitize=address"
        "-fsanitize=undefined")
    target_link_options(tests PRIVATE
        "-fsanitize=address"
        "-fsanitize=undefined")
else()
    message("Address and thread sanitizer during tests: Disabled")
endif()
