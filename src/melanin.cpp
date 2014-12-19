#include "melanin.h"
#include <cmath>
MelaninType MELANIN_TYPE = SVAASAND;

float melanin(float w){
	switch(MELANIN_TYPE){
		case EUMELANIN:
			return exp(-kEu*(w - MELANIN_REFERENCE_WAVELENGTH)/MELANIN_REFERENCE_WAVELENGTH);
		break;
		case PHEOMELANIN:
			return exp(-kPheo*(w - MELANIN_REFERENCE_WAVELENGTH)/MELANIN_REFERENCE_WAVELENGTH);
		break;
		default:
			return pow(MELANIN_REFERENCE_WAVELENGTH*1.0/w, 3.46);
	}
}
