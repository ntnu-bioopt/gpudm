//=======================================================================================================
//// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
//// Distributed under the MIT License.
//// (See accompanying file LICENSE or copy at
//// http://opensource.org/licenses/MIT)
////=======================================================================================================

#include "parseconfig.h"
#include <QString>
#include <QDomDocument>
#include <QDomElement>
#include <QFile>
#include <QIODevice>
#include <QDebug>
#include <iostream>
#include <QStringList>
using namespace std;

QDomDocument *configTree; //xml document representing config file

void readConfigFile(string configFile){
	configTree = new QDomDocument("configFile");
	QFile file(QString::fromStdString(configFile));
	if (!file.open(QIODevice::ReadOnly)){
		cerr << "Failed to open configuration file" << endl;
		exit(1);
	}	
	if (!configTree->setContent(&file)) {
		cerr << "Failed to parse XML in configuration file" << endl;
		file.close();
		exit(1);
	}
	file.close();
}

QDomElement createElem(const string tagname){
	return configTree->createElement(QString::fromStdString(tagname));
}

void generateConfigFile(string configFile){
}

QString readSetting(string categories, string property){
	//split categories
	QStringList splitted = QString::fromStdString(categories).split(QString::fromStdString(CAT_SEP));

	//loop through nesting of categories
	QDomElement categoryNode = configTree->documentElement().elementsByTagName(splitted[0]).at(0).toElement();
	for (int i=1; i < splitted.size(); i++){
		categoryNode = categoryNode.elementsByTagName(splitted[i]).at(0).toElement();
	}

	//get all properties
	QDomNodeList propertyNodes = categoryNode.elementsByTagName(QString::fromStdString(property));
	if (propertyNodes.size() == 0){
		cerr << "Error in config file, option missing: " << categories << "/" << property << endl;
		exit(1);
	}

	QString retStr;

	//return comma-separated list of properties at same nest level
	for (int i=0; i < propertyNodes.size(); i++){
		retStr += QString(propertyNodes.at(i).toElement().text());
		if ((propertyNodes.size() > 1) && (i < propertyNodes.size() - 1)){
			retStr += ",";
		}
	}
	return retStr;
}

QString readSetting(string categories, string property, string type){
	QStringList splitted = QString::fromStdString(categories).split(QString::fromStdString(CAT_SEP));

	QDomElement categoryNode = configTree->documentElement().elementsByTagName(splitted[0]).at(0).toElement();
	for (int i=1; i < splitted.size(); i++){
		categoryNode = categoryNode.elementsByTagName(splitted[i]).at(0).toElement();
	}
	QDomNodeList propertyNodes = categoryNode.elementsByTagName(QString::fromStdString(property));
	if (propertyNodes.size() == 0){
		cerr << "Error in config file, option missing: " << categories << "/" << property << endl;
		exit(1);
	}

	//find node with correct type attribute
	for (int i=0; i < propertyNodes.size(); i++){
		if (propertyNodes.at(i).toElement().attribute("type", "DEFAULT") == QString::fromStdString(type)){
			return propertyNodes.at(i).toElement().text();
		}
	}
	cerr << "Error in config file, option missing: " << categories << "/" << property << ", " << type << endl;
	exit(1);
}

string readStringSetting(string categories, string property){
	return readSetting(categories, property).toStdString();
}

float readNumberSetting(string categories, string property){
	return readSetting(categories, property).toFloat();
}

float readNumberSetting(std::string categories, std::string property, std::string type){
	return readSetting(categories, property, type).toFloat();
}

Chromophores parseChromophoreString(string chromString){
	Chromophores retChrom;
	QStringList splitted = QString::fromStdString(chromString).split(",");
	for (int i = 0; i < splitted.size(); i++){
		string str = splitted[i].toStdString();
		if (str == string("water")){
			retChrom.setWat();
		} else if (str == string("const")){
			retChrom.setKonst();
		} else if (str == string("melanin")){
			retChrom.setMel();
		} else if (str == string("bilirubin")){
			retChrom.setBil();
		} else if (str == string("betacarotene")){
			retChrom.setBet();
		} else if (str == string("methb")){
			retChrom.setMethb();
		} else if (str == string("cohb")){
			retChrom.setCOHb();
		} else if (str == string("fat")){
			retChrom.setFat();
		}
	}
	return retChrom;
}
