cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)
enable_testing()

project (CGALAllegedBug C CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

find_package(CGAL REQUIRED)
add_compile_definitions(CGAL_LINKED_WITH_TBB)

find_package(TBB REQUIRED COMPONENTS tbbmalloc)

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
  ${CGAL_LIBRARY}
  ${TBB_LIBRARIES}
  ${TBB_MALLOC_LIBRARIES}
  ${GMP_LIBRARIES}
)
