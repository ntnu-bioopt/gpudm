#ifndef CONFIG_H_DEFINED
#define CONFIG_H_DEFINED

#include <string>
#include "chromophores.h"

//a static XML document tree resides in config.cpp and is read from and to

Chromophores parseChromophoreString(std::string chromString); //parse string of listed chromophores, water,melanin,cohb, etc, return properly set chromophore object. Blood is set by default.
void readConfigFile(std::string configFile); //read file into static XML document tree
void generateConfigFile(std::string configFile); //generate default config file FIXME: currently not implemented. 
std::string readStringSetting(std::string category, std::string property); //read string setting from specified property
float readNumberSetting(std::string category, std::string property); //read string setting and convert to number
float readNumberSetting(std::string category, std::string property, std::string type); //property type="input string" in addition to usual tags


const std::string CAT_SEP="/";






#endif
