file(GLOB SRCS src/*.cpp)
add_library(gff STATIC ${SRCS})
target_include_directories(gff PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
