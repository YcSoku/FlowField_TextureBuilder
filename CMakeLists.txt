cmake_minimum_required(VERSION 3.20)
project(FlowField_TextureBuilder)

set(CMAKE_CXX_STANDARD 20)

# Configure assets header file
configure_file(${CMAKE_SOURCE_DIR}/include/helpers/RootDir.h.in
        ${CMAKE_SOURCE_DIR}/include/helpers/RootDir.h)

# We need a CMAKE_DIR with some code to find external dependencies
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

message("system: ${CMAKE_SYSTEM_PREFIX_PATH}")

# stb
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/stb/include)
add_library(stb ${CMAKE_SOURCE_DIR}/thirdparty/stb/stb_image.cpp)
List(APPEND LIBS stb)

# quikGrid
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/quikgrid/include)
add_library(quikgrid
        ${CMAKE_SOURCE_DIR}/thirdparty/quikgrid/contour.cpp
        ${CMAKE_SOURCE_DIR}/thirdparty/quikgrid/gridxtyp.cpp
        ${CMAKE_SOURCE_DIR}/thirdparty/quikgrid/scatdata.cpp
        ${CMAKE_SOURCE_DIR}/thirdparty/quikgrid/surfgrid.cpp
        ${CMAKE_SOURCE_DIR}/thirdparty/quikgrid/xpand.cpp
        )
List(APPEND LIBS quikgrid)

# Nlohmann_json
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/json)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/json/include)
List(APPEND LIBS nlohmann_json::nlohmann_json)

# PROJ
find_package(PROJ REQUIRED CONFIG)
List(APPEND LIBS PROJ::proj)

# GDAL
find_package(GDAL REQUIRED)
List(APPEND LIBS GDAL::GDAL)

# Define the executable
aux_source_directory("source" SRC_FILES)
add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SRC_FILES} source/cesium_cartesian3.cpp source/FlowField.cpp)

# Define the include directories
target_include_directories( ${PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        )

# Define the link libraries
target_link_libraries(${PROJECT_NAME} ${LIBS})

