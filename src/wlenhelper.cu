//=======================================================================================================
// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//=======================================================================================================


#include "wlenhelper.h"
using namespace std;

//helper functions for common frame of reference for what indices should correspond to what wavelengths
int getStartInd(float *wlens, int bands, float startwlen){
	//find the wavelength indices corresponding to the specificed wavelength interval
	int startwlenInd;
	bool foundStart = false;
	for (int i=0; i < bands; i++){
		if ((wlens[i] > startwlen) && (!foundStart)){
			startwlenInd = i-1;
			foundStart = true;
			break;
		}
	}
	if (!foundStart){
		startwlenInd = 0;
	}
	return startwlenInd;
}
int getEndInd(float *wlens, int bands, float endwlen){
	//find the wavelength indices corresponding to the specificed wavelength interval
	int endwlenInd;
	bool foundEnd = false;
	for (int i=0; i < bands; i++){
		if ((wlens[i] > endwlen) && (!foundEnd)){
			endwlenInd = i-1;
			foundEnd = true;
			break;
		}
	}
	if (!foundEnd){
		endwlenInd = bands-1;
	}
	return endwlenInd;
}
