cmake_minimum_required(VERSION 2.8.11)

project(libchromophoreconfig)

#using Qt for xml and regex
set(CMAKE_AUTOMOC ON)
find_package(Qt5Core REQUIRED)
find_package(Qt5Xml REQUIRED)

set(CMAKE_BUILD_TYPE Release)

#create library
add_library(chromophoreconfig SHARED chromophores.cpp parseconfig.cpp absorption_props.cpp)
qt5_use_modules(chromophoreconfig Xml)

#install
install(TARGETS chromophoreconfig DESTINATION lib)
install(FILES chromophores.h parseconfig.h absorption_props.h DESTINATION include)



