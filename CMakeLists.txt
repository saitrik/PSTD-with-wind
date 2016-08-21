cmake_minimum_required(VERSION 3.5)
project(untitled)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
INCLUDE_DIRECTORIES(C:/Users/strikootam/ClionProjects/untitled/fftw3)
set(SOURCE_FILES main.cpp wind.cpp wind.h spatderv.cpp spatderv.h PML.cpp PML.h)
add_executable(untitled ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(untitled C:/Users/strikootam/ClionProjects/untitled/fftw3/libfftw3-3.lib)

