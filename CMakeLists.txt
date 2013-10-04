cmake_minimum_required(VERSION 2.8)

project( proto-eigen C CXX )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(Eigen3)
find_package(Boost 1.53)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})

add_executable(proto-expressions proto_expressions.cpp)
add_executable(proto-grammars proto_grammars.cpp)
add_executable(proto-transforms proto_transforms.cpp)
add_executable(eigen-calculator-problem eigen_calculator_problem.cpp)
add_executable(eigen-calculator-solution eigen_calculator_solution.cpp)

add_executable(user-proto-terminals user_proto_terminals.cpp)
set_target_properties(user-proto-terminals PROPERTIES COMPILE_FLAGS -std=c++11)

add_executable(user-proto-transpose user_proto_transpose.cpp)
add_executable(user-proto-matrix user_proto_matrix.cpp)
add_executable(fem-example fem_example.cpp)
add_executable(fem-expressions fem_expressions.cpp)
