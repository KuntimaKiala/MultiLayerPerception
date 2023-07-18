cmake_minimum_required(VERSION 3.6)
project(train) 
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (source_dir "${PROJECT_SOURCE_DIR}/src/")

file (GLOB source_files "${source_dir}/main.cpp" 
                        "${source_dir}/Matrix.cpp" 
                        "${source_dir}/Topology.cpp" 
                        "${source_dir}/SimpleNeuralNetWork.cpp" 
                        "${source_dir}/Neurons.cpp")

add_executable(train ${source_files})
