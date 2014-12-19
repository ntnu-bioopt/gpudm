#ifndef MELANIN_H_DEFINED
#define MELANIN_H_DEFINED
float melanin(float w);
enum MelaninType{PHEOMELANIN, EUMELANIN, SVAASAND};

const float MELANIN_REFERENCE_WAVELENGTH = 694.0f;

enum MelaninDetermination{USE_STRAIGHT_LINE, USE_MELANIN_CURVE};

#define PHEOMELANIN_GPU 0
#define SVAASAND_MELANIN_GPU 1
#define EUMELANIN_GPU 2
#define NONE_GPU 3

#define kPheo 4.780f
#define kEu 2.429f

const float SCALING_THRESHOLD_MELANINTYE_DETECTION = 0.90;

#endif
