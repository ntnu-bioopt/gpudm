#include "chromophores.h"
#include "absorption_props.h"
#include <cmath>
#include "parseconfig.h"
#include <iostream>
using namespace std;


float melanin(float w){
	return pow((694.0f/w), 3.46);
}


Chromophores::Chromophores(){
	hasoxy = 1;
	hasdeoxy = 1;
	haswat = 0;
	hasbil = 0;
	hasbet = 0;
	hasmel = 0;
	haskonst = 0;
	hasmethb = 0;
	hascohb = 0;
	hasfat = 0;
	numEnd = 2;
		
	oxymax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "oxyfrac");
	deoxymax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "deoxyfrac");
	methbmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "methb");
	bilmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "bilirubin");
	betmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "betacarotene");
	watmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "water");
	cohbmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "cohb");
	melmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "melanin");
	konstmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "const");
	fatmax = readNumberSetting("visualization/chromophoreMaxVals", "chromophore", "fat");
}

std::string Chromophores::getName(int endmember){
	int k=0;
	float retval = 0;
	if (hasoxy){
		if (k == endmember){
			return string("oxyHb");
		}
		k++;
	} 
	if (hasdeoxy){
		if (k == endmember){
			return string("deoxyHb");
		}
		k++;
	}
	if (hasmethb){
		if (k == endmember){
			return string("methb");
		}
		k++;
	}
	if (hasbil){
		if (k == endmember){
			return string("bil");
		}
		k++;
	}
	if (hasbet){
		if (k == endmember){
			return string("bet");
		}
		k++;
	}
	if (haswat){
		if (k == endmember){
			return string("water");
		}
		k++;
	}
	if (hascohb){
		if (k == endmember){
			return string("hbco");
		}
		k++;
	}
	if (hasmel){
		if (k == endmember){
			return string("melanin");
		}
		k++;
	}
	if (haskonst){
		if (k == endmember){
			return string("const");
		}
		k++;
	}
	if (hasfat){
		if (k == endmember){
			return string("fat");
		}
		k++;
	}
	
}
		
float *Chromophores::getAbsArray(float wlen){
	float *arr = new float[numEnd];
	int k=0;
	if (hasoxy){
		arr[k] = muh_oxy_calc(wlen);
		k++;
	} 
	if (hasdeoxy){
		arr[k] = muh_deoxy_calc(wlen);
		k++;
	}
	if (hasmethb){
		arr[k] = methb(wlen);
		k++;
	}
	if (hasbil){
		arr[k] = bilirub(wlen);
		k++;
	}
	if (hasbet){
		arr[k] = betacarot(wlen);
		k++;
	}
	if (haswat){
		arr[k] = water(wlen);
		k++;
	}
	if (hascohb){
		arr[k] = HbCO(wlen);
		k++;
	}
	if (hasmel){
		arr[k] = melanin(wlen);
		k++;
	}
	if (haskonst){
		arr[k] = 1.0f;
		k++;
	}
	if (hasfat){
		arr[k] = fat(wlen);
		k++;
	}

	return arr;
}
int Chromophores::getMelInd(){
	int k=0;
	if (hasoxy){
		k++;
	}
	if (hasdeoxy){
		k++;
	}
	if (hasmethb){
		k++;
	}
	if (hasbil){
		k++;
	}
	if (hasbet){
		k++;
	}
	if (haswat){
		k++;
	}
	if (hascohb){
		k++;
	}
	if (hasmel){
		return k;
	}
	return -1;	
}

int Chromophores::getMethbInd(){
	int k=2;
	if (hasmethb){
		return k;
		k++;
	}
	return -1;	
	
}

int Chromophores::getKonstInd(){
	int k=0;
	if (hasoxy){
		k++;
	}
	if (hasdeoxy){
		k++;
	}
	if (hasmethb){
		k++;
	}
	if (hasbil){
		k++;
	}
	if (hasbet){
		k++;
	}
	if (haswat){
		k++;
	}
	if (hascohb){
		k++;
	}
	if (hasmel){
		k++;
	}
	if (haskonst){
		return k;
	}
	return -1;	
}

float Chromophores::getAbs(float w, int endmember){
	int k=0;
	float retval = 0;
	if (hasoxy){
		if (k == endmember){
			retval = muh_oxy_calc(w);
		}
		k++;
	} 
	if (hasdeoxy){
		if (k == endmember){
			retval = muh_deoxy_calc(w);
		}
		k++;
	}
	if (hasmethb){
		if (k == endmember){
			retval = methb(w);
		}
		k++;
	}
	if (hasbil){
		if (k == endmember){
			retval = bilirub(w);
		}
		k++;
	}
	if (hasbet){
		if (k == endmember){
			retval = betacarot(w);
		}
		k++;
	}
	if (haswat){
		if (k == endmember){
			retval = water(w);
		}
		k++;
	}
	if (hascohb){
		if (k == endmember){
			retval = HbCO(w);
		}
		k++;
	}
	if (hasmel){
		if (k == endmember){
			retval = melanin(w);
		}
		k++;
	}
	if (haskonst){
		if (k == endmember){
			retval = 1.0f;
		}
		k++;
	}
	if (hasfat){
		if (k == endmember){
			retval = fat(w);
		}
		k++;
	}
	return retval;
}

void Chromophores::checkAndSet(bool *val){
	if (!(*val)){
		*val = true;
		numEnd++;
	}
}

void Chromophores::checkAndUnset(bool *val){
	if ((*val)){
		*val = false;
		numEnd--;
	}
}

float Chromophores::getMaxVal(int endmember){
	int k=0;
	float retval = 0;
	if (hasoxy){
		if (k == endmember){
			retval = oxymax;
		}
		k++;
	} 
	if (hasdeoxy){
		if (k == endmember){
			retval = deoxymax;
		}
		k++;
	}
	if (hasmethb){
		if (k == endmember){
			retval = methbmax;
		}
		k++;
	}
	if (hasbil){
		if (k == endmember){
			retval = bilmax;
		}
		k++;
	}
	if (hasbet){
		if (k == endmember){
			retval = betmax;
		}
		k++;
	}
	if (haswat){
		if (k == endmember){
			retval = watmax;
		}
		k++;
	}
	if (hascohb){
		if (k == endmember){
			retval = cohbmax;
		}
		k++;
	}
	if (hasmel){
		if (k == endmember){
			retval = melmax;
		}
		k++;
	}
	if (haskonst){
		if (k == endmember){
			retval = konstmax;
		}
		k++;
	}
	if (hasfat){
		if (k == endmember){
			retval = fatmax;
		}
		k++;
	}
	return retval;
	
}
