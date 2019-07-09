# cmake -D CMAKE_C_COMPILER=gcc-6 -D CMAKE_CXX_COMPILER=g++-6 ..

cmake_minimum_required(VERSION 3.7.2)

add_compile_options("-std=gnu++11")
add_compile_options("-O3")

Set(SourceDir "./src/")

Set(Headers "${SourceDir}first_lower_bound.h" "${SourceDir}/mathFunctions.h" "${SourceDir}/random_points.h" "${SourceDir}/reading.h")
Set(Sources "${SourceDir}/sampling.cpp" "${SourceDir}/mathFunctions.cpp" "${SourceDir}/random_points.cpp")

add_executable(main "${SourceDir}main.cpp" ${Sources} ${Headers})

target_link_libraries (main -pthread -lboost_system -lgsl -lgslcblas -lm -lboost_filesystem)

################################
# GTest
################################
ADD_SUBDIRECTORY (Tests/googletest/)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

################################
# Unit Tests
################################
# Add test cpp file
Set(TestSourceDir "./Tests/TestSrc/")
Set(TestSources "${TestSourceDir}testgtest.cpp")

add_executable(runUnitTests ${TestSources} ${Sources} ${Headers})

# Link test executable against gtest & gtest_main
target_link_libraries(runUnitTests gtest gtest_main -pthread -lboost_system)

set(UNIT_TEST runUnitTests)
add_custom_command(
     TARGET ${UNIT_TEST}
     COMMENT "Run tests"
     POST_BUILD 
     COMMAND ${UNIT_TEST}
)


