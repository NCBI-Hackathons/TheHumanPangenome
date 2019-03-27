set(CMAKE_POLICY_DEFAULT_CMP0012 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0074 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

# Check if IPO is supported
cmake_policy(SET CMP0069 NEW)
include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In source builds are not recommended. Please run cmake in a separate build directory")
endif()

if(NOT DEFINED CMAKE_BUILD_TYPE)
  message(STATUS "CMAKE_BUILD_TYPE not set. Setting Release build type as default")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
endif()

# Use ccache if found to cache previously built object files
find_program(CCACHE_EXE ccache)
if(CCACHE_EXE)
  message(STATUS "Found ccache in PATH. Using ccache to cache build artifacts")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_EXE})
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_EXE})
endif()
