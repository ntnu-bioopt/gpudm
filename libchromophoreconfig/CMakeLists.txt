cmake_minimum_required(VERSION 2.6)

project(libchromophoreconfig)

#using Qt for xml and regex
FIND_PACKAGE(Qt4 REQUIRED)
SET(QT_USE_QTXML TRUE)
INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

set(CMAKE_BUILD_TYPE Release)

#create library
add_library(chromophoreconfig SHARED chromophores.cpp parseconfig.cpp absorption_props.cpp)
TARGET_LINK_LIBRARIES(chromophoreconfig ${QT_LIBRARIES})

#install
install(TARGETS chromophoreconfig DESTINATION lib)
install(FILES chromophores.h parseconfig.h absorption_props.h DESTINATION include)



