cmake_minimum_required(VERSION 3.6)
project(class_matrix)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES Circuit.h parsing.h parsing.cpp main.cpp Matrix.h Circuit.cpp)
add_executable(class_matrix ${SOURCE_FILES})