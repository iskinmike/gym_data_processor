# Copyright 2018, mike at myasnikov.mike@gmail.com
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required ( VERSION 2.8.12 )

project ( error_test_gym_data_processor CXX )

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -g -O0")

set(SOURCES
   src/main.cpp
    src/data_parser.cpp
    src/data_parser.hpp
    src/utils.cpp
    src/utils.hpp
)

set(LIBS
    # avformat
    # avcodec
    # swscale
    # avdevice
    # swresample
    # avutil
    pthread
    # X11
    # /home/mike/hobbies/learn_gl/glfw/build/src/libglfw.so
    # ${ARTOOL_SDK}/lib/libARX.so
    # SDL2
    # GLEW
    # GL
)


add_executable(${PROJECT_NAME} ${SOURCES})

# add_library(${PROJECT_NAME} SHARED ${SOURCES})

 add_custom_target ( dist ALL
     COMMAND ${CMAKE_COMMAND} -E remove_directory result
     COMMAND ${CMAKE_COMMAND} -E make_directory result
#     COMMAND ${CMAKE_COMMAND} -E make_directory dist/resources
     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_LIST_DIR}/resources/process_result.sh result/
 )


target_include_directories ( ${PROJECT_NAME} BEFORE PRIVATE
#        $ENV{WILTON_HOME}/devel/include
        # ${ARTOOL_SDK}/include
#        ${PLATFORM_INCLUDE}
        # /home/mike/hobbies/learn_gl/glfw/include
        # /home/mike/hobbies/learn_gl/glfw/src
#        ${STATICLIBS_HEADERS_DIRS}
)

target_link_libraries ( ${PROJECT_NAME} ${LIBS})

# target_compile_definitions(${PROJECT_NAME} PUBLIC -DARX_TARGET_PLATFORM_LINUX=1
                                           # PUBLIC -DDEBUG=1)

