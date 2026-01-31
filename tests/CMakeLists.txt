project(tests)

set(LIB_TESTS
    testMatrix
    testJordanDecomp
    testPolynomial
)

foreach(TEST_NAME ${LIB_TESTS})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
