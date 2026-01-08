cmake_minimum_required(VERSION 3.10)

project(algebra)

### Include this package's and external library headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(/usr/local/include)

### Add subdirectories
# add_subdirectory(src)
add_subdirectory(apps)
