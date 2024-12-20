cmake_minimum_required(VERSION 3.16)

project(NEXTNet LANGUAGES CXX)

# *** Simulation library sources and headers
set(LIB_SOURCES
    nextnet/network.cpp
    nextnet/weighted_network.cpp
    nextnet/temporal_network.cpp
    nextnet/brownian_proximity_network.cpp
    nextnet/algorithm.cpp
    nextnet/random.cpp
    nextnet/NextReaction.cpp
    nextnet/NextReactionMeanField.cpp
    nextnet/nMGA.cpp
    nextnet/REGIR.cpp
    nextnet/utility.cpp)

set(LIB_HEADERS
    nextnet/types.h
    nextnet/network.h
    nextnet/weighted_network.h
    nextnet/temporal_network.h
    nextnet/brownian_proximity_network.h
    nextnet/random.h
    nextnet/algorithm.h
    nextnet/NextReaction.h
    nextnet/NextReactionMeanField.h
    nextnet/nMGA.h
    nextnet/REGIR.h
    nextnet/permutation.h
    nextnet/utility.h)

SET(LIB_PCH nextnet/stdafx.h)

# *** Command-line simulator sources and headers
set(SIMULATOR_SOURCES
    nextnet/analysis.cpp
    nextnet/main.cpp
    nextnet/programs/simulate.cpp)

set(SIMULATOR_HEADERS
    nextnet/analysis.h)

# *** Unit tests sources and headers
set(TESTS_SOURCES
    nextnet/analysis.cpp
    nextnet/tests/plot.cpp
    nextnet/tests/main.cpp
    nextnet/tests/utility.cpp
    nextnet/tests/network.cpp
    nextnet/tests/weighted_network.cpp
    nextnet/tests/temporal_network.cpp
    nextnet/tests/random.cpp
    nextnet/tests/nextreaction.cpp
    nextnet/tests/temporal_nextreaction.cpp
    nextnet/tests/brownian_proximity_network.cpp
    nextnet/tests/nmga.cpp
    nextnet/tests/REGIR.cpp)

set(TESTS_HEADERS
    nextnet/tests/statistics.h
    nextnet/tests/simulate.h
    nextnet/tests/analytical.h
    nextnet/tests/parallel.h
    nextnet/tests/plot.h
    nextnet/tests/gnuplot-iostream.h
    nextnet/catch/catch.hpp)

set(TESTS_PCH
    nextnet/tests/stdafx.h)

# *** Global build settings

# Includes are relative to the nextnet/ folder
include_directories(AFTER SYSTEM nextnet ext/dyndist ext/popl/include)

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
find_package(Boost CONFIG REQUIRED OPTIONAL_COMPONENTS iostreams filesystem)
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
else()
    # It seems that the mere presence of TBB on a system causes
    # the parallel STL backend that uses TBB to be selected, which
    # then causes linker errors if TBB is not actually linked. It
    # seems that defining this macro forces the serial backend and
    # prevents the linker errors. This is probably the wrong solution
    # but it works.
    add_compile_definitions(_PSTL_PAR_BACKEND_SERIAL)
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
set_source_files_properties(nextnet/tests/main.cpp
    PROPERTIES SKIP_PRECOMPILE_HEADERS ON)
# Test data
target_compile_definitions(tests PRIVATE TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/nextnet/tests/data")
message("Test data: ${PROJECT_SOURCE_DIR}/nextnet/tests/data")
if ((NOT SERIALIZE) AND TBB_FOUND)
    message("Parallelization: Enabled")
    target_link_libraries(tests TBB::tbb)
    target_link_libraries(simulator PRIVATE TBB::tbb)
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
