message(STATUS "CMake version found: ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(ggfftools VERSION 0.1.0 LANGUAGES CXX)
include(defaults)
include(external)
set(CMAKE_CXX_STANDARD 14)

set(GGFF_PROTO_FILES ggff.proto)
set(GGFF_PROTO_GENERATED ${CMAKE_BINARY_DIR}/ggff.pb.h ${CMAKE_BINARY_DIR}/ggff.pb.cc)
add_custom_command(OUTPUT ${GGFF_PROTO_GENERATED}
                   COMMENT "Generating C++ classes to encode/decode ggff protobuf messages"
                   COMMAND ${PROTOC_BINARY} -I${PROJECT_SOURCE_DIR} -I${PROTOBUF_INCLUDE_DIR} --cpp_out=${CMAKE_BINARY_DIR} ${GGFF_PROTO_FILES}
                   DEPENDS ${GGFF_PROTO_FILES} protobuf WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

add_library(ggff_proto STATIC ${GGFF_PROTO_GENERATED})
target_include_directories(ggff_proto SYSTEM PUBLIC ${CMAKE_BINARY_DIR} PRIVATE ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(ggff_proto PUBLIC ${LIBPROTOBUF})

add_library(ggff STATIC ggff.h ggff.cc)
target_link_libraries(ggff PUBLIC leveldb roaring_cpp absl::strings ggff_proto PRIVATE fmt absl::memory absl::flat_hash_map)
target_include_directories(ggff PUBLIC ${PROJECT_SOURCE_DIR} PRIVATE ${CMAKE_BINARY_DIR} ${croaring_SOURCE_DIR}/include)
set_target_properties(ggff PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE 1)

if(IPO_SUPPORTED AND "${CMAKE_BUILD_TYPE}" STREQUAL Release)
  set_property(TARGET ggff PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_executable(ggfftools main.cc)
target_link_libraries(ggfftools PRIVATE ggff)
target_include_directories(ggfftools PUBLIC ${PROJECT_SOURCE_DIR})
set_target_properties(ggfftools PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON POSITION_INDEPENDENT_CODE 1)

if(IPO_SUPPORTED AND "${CMAKE_BUILD_TYPE}" STREQUAL Release)
  set_property(TARGET ggfftools PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
