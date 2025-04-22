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
set(CLI_SOURCES
    nextnet/main.cpp
    nextnet/factories/factories.cpp
    nextnet/factories/network.cpp
    nextnet/factories/time.cpp
    nextnet/factories/algorithm.cpp
    nextnet/pstream/pstream.cpp)

set(CLI_HEADERS
    nextnet/factories/factories.h
    nextnet/factories/network.h
    nextnet/factories/time.h
    nextnet/factories/algorithm.h
    nextnet/pstream/pstream.h)

# *** Unit tests sources and headers
set(TESTS_SOURCES
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

# *** Bundled source that shouldn't be formatted
set(BUNDLED_HEADERS_AND_SOURCES
    nextnet/catch/catch.hpp
    nextnet/pstream/pstream.h
    nextnet/pstream/pstream.cpp)

# *** Variables

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Defaulting to CMAKE_BUILD_TYPE=RelWithDebInfo.")
	set(CMAKE_BUILD_TYPE Debug CACHE STRING "build type" FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
else()
	message(STATUS "CMAKE_BUILD_TYPE=" ${CMAKE_BUILD_TYPE})
endif()

if(NOT SANITIZE)
	set(SANITIZE "debug" CACHE STRING "whether to enable address and thread sanitizers" FORCE)
	set_property(CACHE SANITIZE PROPERTY STRINGS "always" "debug" "never")
endif()

if(NOT UNITTESTS)
	set(UNITTESTS NO CACHE BOOL "whether to enable target nextnet-unittests" FORCE)
	set_property(CACHE UNITTESTS PROPERTY TYPE BOOL)
endif()

if(NOT SANITIZE_UNITTESTS)
	set(SANITIZE_UNITTESTS "always" CACHE STRING "whether to enable address and thread sanitizer during tests" FORCE)
	set_property(CACHE SANITIZE_UNITTESTS PROPERTY STRINGS "always" "debug" "never")
endif()

if(NOT SINGLETHREADED)
	set(SINGLETHREADED NO CACHE BOOL "whether to parallelize unit tests" FORCE)
	set_property(CACHE SINGLETHREADED PROPERTY TYPE BOOL)
endif()

if(NOT CLANGFORMAT)
	set(CLANGFORMAT NO CACHE BOOL "whether to enable target clangformat" FORCE)
	set_property(CACHE CLANGFORMAT PROPERTY TYPE BOOL)
endif()

# *** Global build settings

# Includes are relative to the project folder
include_directories(AFTER SYSTEM . ext ext/dyndist ext/popl/include)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hide symbols
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN true)

# Enable all warnings, add CMAKE_DEFINES, and tell Boost that
# we don't have std::unary_function in C++17 mode
add_compile_options(-Wall)
add_compile_definitions(${CMAKE_DEFINES})

# Boost policies
cmake_policy(SET CMP0144 NEW)

# *** Dependencies

# Add cmake directory to CMake module search path
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Boost headers are mandatory, iostreamd and filesystem are optional
# and allow plots to be generated during test runs
find_package(Boost CONFIG REQUIRED OPTIONAL_COMPONENTS iostreams filesystem)
link_libraries(Boost::boost)
include_directories(${Boost_INCLUDE_DIRS})
# Some boost libraries attempt to use std::unary_function if this is
# not set, which causes problems since C++17 removed this class.
add_compile_definitions(NEXTNET_BOOST_NO_CXX98_FUNCTION_BASE)

# *** Build command-line similator
   message(STATUS "nextnet: enabling target")
add_executable(nextnet
    ${CLI_SOURCES}
    ${LIB_SOURCES}
    ${CLI_HEADERS}
    ${LIB_HEADERS}
    ${LIB_PCH})
target_precompile_headers(nextnet
    PUBLIC ${LIB_PCH})
# Keep nextnet single-threaded, avoid having to link TBB even if its present
target_compile_definitions(nextnet PRIVATE _PSTL_PAR_BACKEND_SERIAL)
target_compile_definitions(nextnet PRIVATE NEXTNET_BOOST_NO_CXX17_HDR_EXECUTION)

# On older versions of GCC, we have to link libstdc++fs
target_link_libraries(nextnet PRIVATE
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)

# Enable address and undefined behaviour sanitizers for debug builds
if (SANITIZE STREQUAL "debug")
    message(STATUS "nextnet: SANITIZE=${SANITIZE}, enabling addresses and thread sanitizer for debug builds")
    target_compile_options(nextnet PRIVATE
        "$<$<CONFIG:DEBUG>:-fsanitize=address>"
        "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
    target_link_options(nextnet PRIVATE
        "$<$<CONFIG:DEBUG>:-fsanitize=address>"
        "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
elseif (SANITIZE STREQUAL "always")
    message(STATUS "nextnet: SANITIZE=${SANITIZE}, enabling addresses and thread sanitizer for all builds")
    target_compile_options(nextnet PRIVATE
        "-fsanitize=address"
        "-fsanitize=undefined")
    target_link_options(nextnet PRIVATE
        "-fsanitize=address"
        "-fsanitize=undefined")
else()
    message(STATUS "nextnet: SANITIZE=${SANITIZE}, disabling addresses and thread sanitizer for all builds")
endif()

# *** Build unit test runner
if (UNITTESTS)
    message(STATUS "nextnet-unittests: UNITTESTS=${UNITTESTS}, enabling target")

    add_executable(nextnet-unittests
        ${TESTS_SOURCES}
        ${LIB_SOURCES}
        ${TESTS_HEADERS}
        ${TESTS_PCH}
        ${LIB_HEADERS}
        ${LIB_PCH})
    target_precompile_headers(nextnet-unittests
        PUBLIC ${LIB_PCH}
        PUBLIC ${TESTS_PCH})
    set_source_files_properties(nextnet/tests/main.cpp
        PROPERTIES SKIP_PRECOMPILE_HEADERS ON)
    # Test data
    target_compile_definitions(nextnet-unittests PRIVATE TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/nextnet/tests/data")
    message(STATUS "nextnet-unittests: test data in ${PROJECT_SOURCE_DIR}/nextnet/tests/data")

    # These boost libraries are optional
    if (Boost_FILESYSTEM_FOUND)
        message(STATUS "Boost filesystem: Found")
    else()
        message(STATUS "Boost filesystem: Not found")
    endif()
    if (Boost_IOSTREAMS_FOUND)
        message(STATUS "Boost iostreams: Found")
    else()
        message(STATUS "Boost iostreams: Not found")
    endif()

    # Gnuplot is optional, if present can generate plos
    find_program(Gnuplot gnuplot)
    if (Gnuplot)
        message(STATUS "Gnuplot: Found (${Gnuplot})")
    else()
        message(STATUS "Gnuplot: Not found")
    endif()

    if (Gnuplot AND Boost_FILESYSTEM_FOUND AND Boost_IOSTREAMS_FOUND)
        message(STATUS "nextnet-unittests: enabling plotting")
        target_link_libraries(nextnet-unittests Boost::boost Boost::iostreams Boost::filesystem)
        target_compile_definitions(nextnet-unittests PRIVATE ENABLE_PLOTTING=1)
	target_compile_definitions(nextnet-unittests PRIVATE GNUPLOT_DEFAULT_COMMAND="${Gnuplot} -persist")
    else()
        message(STATUS "nextnet-unittests: disabling plotting")
    endif()

    # Intel Thread Building Blocks are optional, if present
    # they speed up the tests by executing simulations in parallel
    if (NOT SINGLETHREADED)
        find_package(TBB)
		if (TBB_FOUND)
			message(STATUS "TBB: Found")
		else()
			message(STATUS "TBB: Not found")
		endif()
    endif()
    if (NOT SINGLETHREADED AND TBB_FOUND)
        message(STATUS "nextnet-unittests: SINGLETHREADED=${SINGLETHREADED}, enabling parallelization")
        target_link_libraries(nextnet-unittests TBB::tbb)
        target_compile_definitions(nextnet-unittests PRIVATE TBB_SUPPRESS_DEPRECATED_MESSAGES)
        target_compile_definitions(nextnet-unittests PRIVATE PARALLELIZE=1)
    else()
        message(STATUS "nextnet-unittests: SINGLETHREADED=${SINGLETHREADED}, disabling parallelization")
        # It seems that the mere presence of TBB on a system causes
        # the parallel STL backend that uses TBB to be selected, which
        # then causes linker errors if TBB is not actually linked. It
        # seems that defining this macro forces the serial backend and
        # prevents the linker errors. This is probably the wrong solution
        # but it works.
        target_compile_definitions(nextnet-unittests PRIVATE _PSTL_PAR_BACKEND_SERIAL)
        target_compile_definitions(nextnet-unittests PRIVATE NEXTNET_BOOST_NO_CXX17_HDR_EXECUTION)
        target_compile_definitions(nextnet-unittests PRIVATE PARALLELIZE=0)
    endif()

    # Enable address and undefined behaviour sanitizers for all builds
    if (SANITIZE_UNITTESTS STREQUAL "debug")
		message(STATUS "nextnet-unittests: SANITIZE_UNITTESTS=${SANITIZE_UNITTESTS}, enabling addresses and thread sanitizer for debug builds")
        target_compile_options(nextnet-unittests PRIVATE
            "$<$<CONFIG:DEBUG>:-fsanitize=address>"
            "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
        target_link_options(nextnet-unittests PRIVATE
            "$<$<CONFIG:DEBUG>:-fsanitize=address>"
            "$<$<CONFIG:DEBUG>:-fsanitize=undefined>")
    elseif (SANITIZE_UNITTESTS STREQUAL "always")
		message(STATUS "nextnet-unittests: SANITIZE_UNITTESTS=${SANITIZE_UNITTESTS}, enabling addresses and thread sanitizer for all builds")
        target_compile_options(nextnet-unittests PRIVATE
            "-fsanitize=address"
            "-fsanitize=undefined")
        target_link_options(nextnet-unittests PRIVATE
            "-fsanitize=address"
            "-fsanitize=undefined")
    else()
		message(STATUS "nextnet-unittests: SANITIZE_UNITTESTS=${SANITIZE_UNITTESTS}, disabling addresses and thread sanitizer for all builds")
    endif()
else()
    message(STATUS "nextnet-unittests: UNITTESTS=${UNITTESTS}, disabling target")
endif()

# *** Code formatter
if (CLANGFORMAT)
    message(STATUS "clangformat: CLANGFORMAT=${CLANGFORMAT}, enabling target")
    include(ClangFormat)

    # Format all sources except third-party files we just bundlet
    set(FORMAT_SOURCES
        ${CLI_SOURCES}
        ${LIB_SOURCES}
        ${CLI_HEADERS}
        ${LIB_HEADERS}
        ${TESTS_HEADERS}
        ${TESTS_SOURCES}
        ${TESTS_PCH}
        ${LIB_PCH})
    list(REMOVE_ITEM FORMAT_SOURCES ${BUNDLED_HEADERS_AND_SOURCES})

    # Create clangformat target, don't build by default
    clangformat_setup(${FORMAT_SOURCES})
    set_target_properties(clangformat PROPERTIES EXCLUDE_FROM_ALL TRUE)
else()
    message(STATUS "clangformat: CLANGFORMAT=${CLANGFORMAT}, disabling target")
endif()
