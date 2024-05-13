cmake_minimum_required(VERSION 3.17)
project(digitazionpp VERSION 0.0.1)

find_package(ROOT REQUIRED)

add_compile_options(-Wall -Wextra -O3)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

add_executable(digitazionpp "${PROJECT_SOURCE_DIR}/MC_data_new.cxx")
target_link_libraries(digitazionpp PUBLIC ${ROOT_LIBRARIES})