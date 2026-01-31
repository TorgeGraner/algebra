cmake_minimum_required(VERSION 3.15)
project(algebra)

include(CTest)
### Include this package's and external library headers
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(/usr/local/include)

### Add subdirectories
add_subdirectory(apps)

enable_testing()

add_subdirectory(tests)