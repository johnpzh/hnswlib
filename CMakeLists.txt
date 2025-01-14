cmake_minimum_required (VERSION 2.6)
project (hnsw_lib)

# Verbose Makefile
set(CMAKE_VERBOSE_MAKEFILE ON)

# C++ Standard
set(CMAKE_CXX_STANDARD 14)

# OpenMP
find_package(OpenMP REQUIRED)

# Thread
find_package(Threads)

include_directories("${PROJECT_BINARY_DIR}")

#set(SOURCE_EXE main.cpp)
#
#set(SOURCE_LIB sift_1b.cpp)

#add_library(sift_test STATIC ${SOURCE_LIB})
#target_link_libraries(sift_test PRIVATE OpenMP::OpenMP_CXX)


#add_executable(main ${SOURCE_EXE})
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET( CMAKE_CXX_FLAGS  "-Ofast -DNDEBUG -DHAVE_CXX0X -march=native -fpic -ftree-vectorize")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
SET( CMAKE_CXX_FLAGS  "-Ofast -lrt -DNDEBUG -DHAVE_CXX0X -march=native -fpic -w -ftree-vectorize -ftree-vectorizer-verbose=0" )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  SET( CMAKE_CXX_FLAGS  "-Ofast -lrt -DNDEBUG -DHAVE_CXX0X -march=native -fpic -w -ftree-vectorize" )
endif()

add_executable(test_updates examples/updates_test.cpp)
target_link_libraries (test_updates PRIVATE Threads::Threads)

#target_link_libraries(main sift_test)


add_executable(
        sift_1b
        sift_1b.cpp)
target_link_libraries(
        sift_1b
        PRIVATE OpenMP::OpenMP_CXX)

add_executable(
        sift_test
        sift_test.cpp)
target_link_libraries(
        sift_test
        PRIVATE OpenMP::OpenMP_CXX)

add_executable(
        hnsw_index
        app/hnsw_index.cpp)
target_link_libraries(
        hnsw_index
        PRIVATE OpenMP::OpenMP_CXX)

add_executable(
        hnsw_search
        app/hnsw_search.cpp)
target_link_libraries(
        hnsw_search
        PRIVATE OpenMP::OpenMP_CXX)

