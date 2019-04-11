cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)
enable_testing()

project (CGALAllegedBug C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

find_package(CGAL REQUIRED)

find_package(TBB REQUIRED COMPONENTS tbbmalloc)
find_package(Threads REQUIRED)

include_directories(
  SYSTEM
  ${CGAL_INCLUDE_DIRS}
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -std=c++17 -ggdb -pedantic -Wall -Wextra -Werror")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -std=c11 -ggdb -pedantic -Wall -Wextra -Werror")

add_executable (a.out
  main.cpp
)

target_link_libraries(a.out
  ${CMAKE_THREAD_LIBS_INIT}
  ${CGAL_LIBRARY}
  ${TBB_LIBRARIES}
  ${TBB_MALLOC_LIBRARIES}
  ${GMP_LIBRARIES}
)
