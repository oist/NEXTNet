cmake_minimum_required(VERSION 3.5)

project(Epidemics-Mean-Field LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    epidemics/main.cpp
    epidemics/graph.cpp
    epidemics/testing.cpp
    epidemics/savingData.cpp)
set(HEADERS
    epidemics/stdafx.h
    epidemics/types.h
    epidemics/graph.h
    epidemics/testing.h
    epidemics/savingData.h)
add_executable(simulator ${SOURCES} ${HEADERS})
