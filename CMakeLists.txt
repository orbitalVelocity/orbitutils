cmake_minimum_required (VERSION 2.8)

project (orbitutils)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(include)

file(GLOB src_files
	"src/*.cpp"
)

add_library(orbitutils	${src_files})
add_executable(test-kepler	${src_files}	test/test-kepler.cpp)
add_executable(test-oeconvert	${src_files}	test/test-oeconvert.cpp)
add_executable(test-lambert	${src_files}	test/test-lambert.cpp)
