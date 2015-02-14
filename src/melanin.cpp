//=======================================================================================================
//// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
//// Distributed under the MIT License.
//// (See accompanying file LICENSE or copy at
//// http://opensource.org/licenses/MIT)
////=======================================================================================================

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
