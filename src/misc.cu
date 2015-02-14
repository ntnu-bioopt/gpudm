//=======================================================================================================
// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//=======================================================================================================

#include "melanin.h"
__global__ void skindataDeinit(float *muam, float *bvf, float *oxy, float *melanin_type, float muamStd, float bvfStd, float oxyStd, size_t pitch){
	int ind = blockIdx.x*pitch+ threadIdx.x;
	muam[ind] = muamStd;
	bvf[ind] = bvfStd;
	oxy[ind] = oxyStd;
	melanin_type[ind] = SVAASAND_MELANIN_GPU;

}

__global__ void SetMelaninType(float *melanin_type, int inputType, size_t pitch){
	melanin_type[blockIdx.x*pitch + threadIdx.x] = inputType;	
}




__global__ void calcOxyBvf(float *inputRes, float *oxyOut, float *bvfOut, size_t pitch){
	int ind = blockIdx.x*pitch + threadIdx.x;
	float foxy = inputRes[ind];
	float fdeoxy = inputRes[ind + gridDim.x*pitch];
	float bvf = foxy + fdeoxy;
	float oxy = fdividef(foxy, bvf);
	bvf = (bvf >= 0)*(bvf <= 1)*bvf;
	oxy = (oxy >= 0)*(oxy <= 1)*oxy;
	oxyOut[ind] = oxy;
	bvfOut[ind] = bvf;
}


__global__ void addUpAbsorption(float *mua, float *endvalues, float *AT, int startblockInd, int numBands, int endmembers, int ignoreInd, size_t pitch){
	int ind = blockIdx.x*pitch + threadIdx.x;
	int pixind = ind + gridDim.x*(blockIdx.y+startblockInd)*pitch;
	float muaval = 0;
	for (int i=0; i < endmembers; i++){
		__shared__ float temp;
		if (threadIdx.x == 0){
			temp = AT[i*numBands + blockIdx.y];
		}
		__syncthreads();
		muaval += temp*endvalues[ind+i*gridDim.x*pitch]*(i != ignoreInd);
	}
	mua[pixind] = muaval;
}
