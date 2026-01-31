project(decomposer)

add_executable(decomposer main.cpp)
target_include_directories(decomposer PUBLIC ${CMAKE_SOURCE_DIR}/include)