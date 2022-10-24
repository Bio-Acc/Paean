cmake_minimum_required(VERSION 3.10)

project(Paean LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
else()
    add_definitions(-DNDEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(OpenMP)
find_package(Thrust)
find_package(HTSlib)
find_package(Threads)
find_package(ZLIB)

include_directories(
    ${CMAKE_CURRENT_LIST_DIR}/include
    ${CMAKE_CURRENT_LIST_DIR}/third_party)

# gff
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/third_party/gff)

# for parse
add_library(parse STATIC
    ${CMAKE_SOURCE_DIR}/src/parse.cpp)
target_include_directories(parse PRIVATE
    ${HTSlib_INCLUDE_DIR}
    ${ZLIB_INCLUDE_DIRS})
target_link_libraries(parse
    ${HTSlib_LIBRARY} ${ZLIB_LIBRARIES}
    gff Threads::Threads)

# for prepare
add_executable(prepare ${CMAKE_SOURCE_DIR}/src/prepare.cpp)
target_link_libraries(prepare gff)

set(SRCS ${CMAKE_SOURCE_DIR}/src/sort.cpp)

# for paean
add_executable(paean ${CMAKE_SOURCE_DIR}/src/gene.cpp ${SRCS})
target_compile_definitions(paean PUBLIC THRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP
                                        THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP)
target_include_directories(paean PUBLIC ${THRUST_INCLUDE_DIR})
target_link_libraries(paean parse OpenMP::OpenMP_CXX)