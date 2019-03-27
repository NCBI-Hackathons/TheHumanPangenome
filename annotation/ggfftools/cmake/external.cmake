# include modules
include(ProcessorCount)
include(ExternalProject)
include(FetchContent)

ProcessorCount(NumCores)
find_program(MAKE_EXE NAMES gmake nmake make)

# Overwrite the message() function to prevent subproject deps to write messages
function(message)
  if (NOT MESSAGE_QUIET)
    _message(${ARGN})
  endif()
endfunction()

set(MESSAGE_QUIET ON)

set(PROTOBUF_ROOT_DIR "${CMAKE_BINARY_DIR}/protobuf-3.7.0")
set(PROTOC_BINARY "${PROTOBUF_ROOT_DIR}/bin/protoc")
set(PROTOBUF_INCLUDE_DIR "${PROTOBUF_ROOT_DIR}/include")
set(LIBPROTOBUF "${PROTOBUF_ROOT_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}protobuf${CMAKE_STATIC_LIBRARY_SUFFIX}")
ExternalProject_Add(protobuf
                    URL https://github.com/protocolbuffers/protobuf/releases/download/v3.7.0/protobuf-cpp-3.7.0.tar.gz
                    URL_MD5 f1631a8e551e569273d78538f6ecf41c
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR} SOURCE_DIR ${PROTOBUF_ROOT_DIR} BUILD_IN_SOURCE 1
                    CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/scripts/configure.sh ${PROTOBUF_ROOT_DIR} ${CMAKE_C_COMPILER}
                    BUILD_COMMAND ${MAKE_EXE} -j${NumCores} INSTALL_COMMAND ${MAKE_EXE} install
                    BUILD_BYPRODUCTS ${LIBPROTOBUF} ${PROTOC_BINARY} LOG_DOWNLOAD ON LOG_CONFIGURE ON
                    LOG_BUILD ON USES_TERMINAL_DOWNLOAD OFF USER_TERMINAL_CONFIGURE OFF USES_TERMINAL_BUILD OFF)

FetchContent_Declare(abseil
                     URL https://github.com/abseil/abseil-cpp/archive/bf29470384a101b307873b26d358433138c857fc.tar.gz
                     URL_MD5 bb28901913baa9a8cdc51b9e0dce917a)
FetchContent_GetProperties(abseil)
if(NOT abseil_POPULATED)
  set(BUILD_TESTING OFF)
  FetchContent_Populate(abseil)
  add_subdirectory(${abseil_SOURCE_DIR} ${abseil_BINARY_DIR})
endif()

FetchContent_Declare(fmt
                     URL https://github.com/fmtlib/fmt/archive/5.3.0.tar.gz
                     URL_MD5 1015bf3ff2a140dfe03de50ee2469401)
FetchContent_GetProperties(fmt)
if(NOT fmt_POPULATED)
  FetchContent_Populate(fmt)
  add_subdirectory(${fmt_SOURCE_DIR} ${fmt_BINARY_DIR})
  target_include_directories(fmt SYSTEM INTERFACE ${fmt_SOURCE_DIR}/include)
endif()

FetchContent_Declare(level
                     URL https://github.com/google/leveldb/archive/7035af5fc36657447054617759854a726d31dbe0.tar.gz
                     URL_MD5 540885fc4a74d5cdd61b499de28c54bf)
FetchContent_GetProperties(level)
if(NOT level_POPULATED)
  FetchContent_Populate(level)
  set(LEVELDB_BUILD_TESTS OFF)
  set(LEVELDB_BUILD_BENCHMARKS OFF)
  add_subdirectory(${level_SOURCE_DIR} ${level_BINARY_DIR})
endif()

FetchContent_Declare(croaring
                     URL https://github.com/RoaringBitmap/CRoaring/archive/v0.2.60.tar.gz
                     URL_MD5 29602918e6890ffdeed84cb171857046)
FetchContent_GetProperties(croaring)
if(NOT croaring_POPULATED)
  FetchContent_Populate(croaring)
  set(ROARING_BUILD_STATIC ON)
  set(ENABLE_ROARING_TESTS OFF)
  add_subdirectory(${croaring_SOURCE_DIR} ${croaring_BINARY_DIR})

  add_library(roaring_cpp INTERFACE)
  target_compile_features(roaring_cpp INTERFACE cxx_std_14)
  target_include_directories(roaring_cpp SYSTEM INTERFACE ${croaring_SOURCE_DIR}/cpp)
  add_dependencies(roaring_cpp croaring)
endif()

unset(MESSAGE_QUIET)
