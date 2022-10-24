cmake_minimum_required(VERSION 3.10)

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif()

if (POLICY CMP0104)
	cmake_policy(SET CMP0104 OLD)
endif()

project(Paean LANGUAGES CXX CUDA)

# standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O0 -g")
else()
    add_definitions(-DNDEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
endif()

set(CUDA_MIN_VERSION "9.2")
find_package(CUDA ${CUDA_MIN_VERSION})
if (NOT CUDA_FOUND)
	message(FATAL_ERROR "CUDA >= ${CUDA_MIN_VERSION} Required, Found ${CUDA_VERSION}.")
endif()

cuda_select_nvcc_arch_flags(ARCH_FLAGS "Auto")
string(REPLACE ";" " " ARCH_FLAGS "${ARCH_FLAGS}")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${ARCH_FLAGS}")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-declarations")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(CUB)
find_package(Thrust)
find_package(HTSlib)
find_package(Threads)
find_package(ZLIB)

set(CUDA_INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

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

# for cub sort
add_library(sort STATIC
    ${CMAKE_SOURCE_DIR}/src/sort.cu)
target_include_directories(sort PRIVATE
    ${CUDA_INCLUDE_DIRS}
    ${CUB_INCLUDE_DIR}
    ${THRUST_INCLUDE_DIR})
# if cuda version >= 11.0
if (CUDA_VERSION VERSION_GREATER_EQUAL 11.0)
    target_compile_definitions(sort PRIVATE THRUST_IGNORE_CUB_VERSION_CHECK)
    target_compile_definitions(sort PRIVATE THRUST_IGNORE_DEPRECATED_CPP_DIALECT)
endif()

# for prepare
add_executable(prepare ${CMAKE_SOURCE_DIR}/src/prepare.cpp)
target_link_libraries(prepare gff)

# for paean
add_executable(paean ${CMAKE_SOURCE_DIR}/src/gene.cu)
target_include_directories(paean PRIVATE
    ${CUDA_INCLUDE_DIRS})
target_link_libraries(paean parse sort)
