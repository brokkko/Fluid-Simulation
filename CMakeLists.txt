cmake_minimum_required(VERSION 3.21)
project(Fluid_Simulation)

set(CMAKE_CXX_STANDARD 20)
set(netcdf_DIR "D:/Program Files/netCDF 4.8.1/lib/cmake/netCDF")

include_directories(headers /usr/include D:/SFML/include)
include_directories(${netcdf_DIR}/../../../include)
link_directories(${netcdf_DIR}/../../../lib GL)

add_executable(Fluid_Simulation main.cpp)

set(OGLLibs  opengl32.lib glu32.lib)

set(SFML_ROOT D:/SFML)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake_modules")
find_package(SFML REQUIRED system window graphics network audio)

find_package(netcdf 4.8.1 CONFIG REQUIRED)

if (SFML_FOUND)
    include_directories(${SFML_INCLUDE_DIR})
    target_link_libraries(Fluid_Simulation PRIVATE netcdf ${SFML_LIBRARIES} ${SFML_DEPENDENCIES} ${OGLLibs})
endif()


# Copy dlls to build
if(WIN32)
    file(GLOB BINARY_DEP_DLLS "${SFML_INCLUDE_DIR}/../bin/*.dll")
    file(COPY ${BINARY_DEP_DLLS} DESTINATION ${CMAKE_BINARY_DIR})
    file(GLOB MINGW_DEP_DLLS "C:/mingw64/bin/*.dll")
    file(COPY ${MINGW_DEP_DLLS} DESTINATION ${CMAKE_BINARY_DIR})
endif()