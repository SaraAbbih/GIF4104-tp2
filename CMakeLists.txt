cmake_minimum_required(VERSION 3.16)
project(tp1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

add_executable(tp1 main.cpp miller-rabin-gmp.cpp miller-rabin-gmp.h)

target_link_libraries(tp1 gmp gmpxx)