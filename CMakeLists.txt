cmake_minimum_required(VERSION 2.8)
project(capdiff)

find_package(Boost COMPONENTS thread system filesystem)

# add OpenMP compiler flag
# enable optimization
# use c++11
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fopenmp -O3 -std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# add subdirectories containing header files to include path
include_directories(${CMAKE_SOURCE_DIR}/2D1C)
include_directories(${CMAKE_SOURCE_DIR}/FieldTool)
include_directories(${CMAKE_SOURCE_DIR}/Simulator)
include_directories(${CMAKE_SOURCE_DIR}/util)

# include the boost headers
include_directories(${BOOST_INCLUDE_DIR})

add_executable(2d1cp
        2D1C/2D1C.cpp # file containing main function
        2D1C/Plasma2D.cpp
        util/Matrix2D.cpp
        util/Matrix3D.cpp
        util/Progress.cpp
        util/RandomNumberGenerator.cpp)

add_executable(fieldtool
        FieldTool/FieldTool.cpp # file containing main function
        Simulator/util/Capillary2D.cpp
        Simulator/util/PrecalculatedField2D.cpp
        Simulator/util/CapillaryConfiguration.cpp
        Simulator/util/Compartement.cpp
        Simulator/util/MPICommunicator.cpp
        util/Progress.cpp
        util/Matrix2D.cpp
        util/Matrix3D.cpp)

add_executable(mcmc
        Main/Exec.cpp
        Simulator/SimulationEvaluator.cpp
        Simulator/util/IOManager.cpp
        Simulator/util/Capillary2D.cpp
        Simulator/util/PrecalculatedField2D.cpp
        Simulator/util/LiveField2D.cpp
        Simulator/util/CapillaryConfiguration.cpp
        Simulator/util/Compartement.cpp
        Simulator/util/MPICommunicator.cpp
        util/Progress.cpp
        util/Matrix2D.cpp
        util/Matrix3D.cpp
        util/RandomNumberGenerator.cpp)


target_link_libraries(2d1cp ${Boost_LIBRARIES})
target_link_libraries(fieldtool ${Boost_LIBRARIES})
target_link_libraries(mcmc ${Boost_LIBRARIES})