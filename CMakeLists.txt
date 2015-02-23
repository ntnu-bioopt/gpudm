cmake_minimum_required(VERSION 2.8)
project(gpudm)
set(CMAKE_BUILD_TYPE Release)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(libchromophoreconfig)

#CUDA
find_package(CUDA)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_11 -Xptxas -v)
cuda_compile(CUDA_O src/inv.cu src/base.cu src/gpudm.cu src/misc.cu src/wlenhelper.cu)
cuda_add_library(gpudm ${CUDA_O})


add_executable(gpudm-bin src/main.cpp src/readimage.cpp)
TARGET_LINK_LIBRARIES(gpudm-bin gpudm chromophoreconfig boost_system)

#install
install (TARGETS gpudm DESTINATION lib)
install (FILES include/gpudm.h DESTINATION include)
