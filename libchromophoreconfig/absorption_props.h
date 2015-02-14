//=======================================================================================================
//// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
//// Distributed under the MIT License.
//// (See accompanying file LICENSE or copy at
//// http://opensource.org/licenses/MIT)
////=======================================================================================================

#ifndef MUH_BLOOD_H_DEFINED
#define MUH_BLOOD_H_DEFINED


float muh_oxy_calc(float l);

float muh_deoxy_calc(float l);

float bilirub(float l);
float methb(float l);
float betacarot(float l);
float betacarot_schwieter(float l);
float water(float l);
float HbCO(float l);
float fat(float l);

#endif
