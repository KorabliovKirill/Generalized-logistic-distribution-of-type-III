cmake_minimum_required(VERSION 3.31 FATAL_ERROR)

set (PROJECT_NAME)
project(${PROJECT_NAME} GLD3)

add_executable(${PROJECT_NAME}
   src/main.cpp
   src/distribution.cpp
   src/distribution.hpp
   src/mix_distribution.cpp
   src/mix_distribution.hpp
   src/emp_distribution.cpp
   src/emp_distribution.hpp
   src/spec_func.cpp
   src/spec_func.hpp
)