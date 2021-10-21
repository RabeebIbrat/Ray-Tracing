cmake_minimum_required(VERSION 3.17)
project(Ray_Tracing)

set(CMAKE_CXX_STANDARD 14)

add_executable(Ray_Tracing 1605055_main.cpp 1605055_objects.h 1605055_glut_display.h 1605055_capture_lighting.h)
target_link_libraries(Ray_Tracing -lOpenGL32 -lfreeGlut -lglu32)