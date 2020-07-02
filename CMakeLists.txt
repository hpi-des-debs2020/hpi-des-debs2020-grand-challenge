cmake_minimum_required(VERSION 2.8.4)
project(grand-challenge)

#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 14)

find_package(Threads)
find_package(RapidJSON)
find_package(CURL REQUIRED)

include_directories(${RapidJSON})
include_directories(${CURL_INCLUDE_DIR})

add_executable(grand-challenge main.cpp Incremental_Detector.cpp Incremental_Detector.h Incremental_Detector_Disorder.cpp Incremental_Detector_Disorder.h query1_NonParallel.cpp query1_NonParallel.h query2_NonParallel.cpp query2_NonParallel.h)
target_link_libraries(grand-challenge ${RapidJSON})
target_link_libraries(grand-challenge ${CURL_LIBRARIES})
target_link_libraries(grand-challenge ${CMAKE_THREAD_LIBS_INIT})

