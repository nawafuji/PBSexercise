cmake_minimum_required(VERSION 2.8)

project(Exercise2)

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

# GLUT
find_package(GLUT)
if(NOT GLUT_FOUND)
	set(GLUT_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/GL)
	find_library(GLUT_LIBRARY glut32 PATHS ${CMAKE_CURRENT_LIST_DIR}/GL )
	message(${GLUT_LIBRARY})
	set(GLUT_LIBRARIES ${GLUT_LIBRARY})
endif(NOT GLUT_FOUND)
include_directories(${GLUT_INCLUDE_DIRS})

file(GLOB ex2_files
		${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/*.h
	)
	
add_executable(Exercise2 ${ex2_files})

target_link_libraries(Exercise2 ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
