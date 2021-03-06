//=======================================================================================================
//// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
//// Distributed under the MIT License.
//// (See accompanying file LICENSE or copy at
//// http://opensource.org/licenses/MIT)
////=======================================================================================================

__global__ void skindataDeinit(float *muam, float *bvf, float *oxy, float *melanin_type, float muamStd, float bvfStd, float oxyStd, size_t pitch);

__global__ void SetMelaninType(float *melanin_type, int inputType, size_t pitch);

__global__ void addUpAbsorption(float *mua, float *endvalues, float *AT, int startblockInd, int numBands, int endmembers, int ignoreInd, size_t pitch);

__global__ void calcOxyBvf(float *inputRes, float *oxyOut, float *bvfOut, size_t pitch);
