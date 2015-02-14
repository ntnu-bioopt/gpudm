//=======================================================================================================
// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//=======================================================================================================


///SCA KIS
#define NUM_ITERATIONS 600
#define METHB_TRESH 1

//implementation of SCA using shared memory
template<int MAX_ENDMEMBERS, int MAX_BLOCKDIM>
__global__ void SCA(float *AT, float *H, float *inputmua, float *x, float *mu, int startind, int numBands, int numEndmembers, size_t pitch){
	__shared__ float sh_x[MAX_ENDMEMBERS][MAX_BLOCKDIM];
	__shared__ float sh_mu[MAX_ENDMEMBERS][MAX_BLOCKDIM];

	//if (blockIdx.x == 5){
	//calculate initial values for mu and x
	int pixind = pitch*blockIdx.x + threadIdx.x;
	
	for (int i=0; i < numEndmembers; i++){
		float temp = 0;
		for (int j=0; j < numBands; j++){
			__shared__ float currAT;
			if (threadIdx.x == 0){
				currAT = AT[i*numBands + j];
			}
			__syncthreads();
			float inp = inputmua[pixind + gridDim.x*pitch*(j+startind)];
			temp += currAT*inp;
		}
		float xstart = 0;
		float mustart = -1*temp;
		sh_x[i][threadIdx.x] = xstart;
		sh_mu[i][threadIdx.x] = mustart;
	}
	__shared__ float sh_H[MAX_ENDMEMBERS][MAX_ENDMEMBERS];
	if (threadIdx.x == 0){
		for (int i=0; i < numEndmembers; i++){
			for (int j=0; j < numEndmembers; j++){
				sh_H[i][j] = H[i*numEndmembers + j];
			}
		}	
	}
	__syncthreads();


	//iterate
	for (int i=0; i < NUM_ITERATIONS; i++){
		for (int j=0; j < numEndmembers; j++){
			float xprev;
			float muprev;
			//extract previous values
			xprev = sh_x[j][threadIdx.x];
			muprev = sh_mu[j][threadIdx.x];

			//calculate new value
			float temp = xprev - fdividef(muprev, sh_H[j][j]);
			temp = (temp > 0)*temp;
			
			//save new value to x
			sh_x[j][threadIdx.x] = temp;

			//calculate new mu
			for (int k=0; k < numEndmembers; k++){
				float addition = (temp - xprev)*sh_H[k][j];
				sh_mu[k][threadIdx.x] += addition;
			}
		}
	}
	for (int i=0; i < numEndmembers; i++){
		float val;
		val = sh_x[i][threadIdx.x];
		x[pixind + i*gridDim.x*pitch] = val;
	}
}


//implementation using registers
//YOLO
#define ATmultMua(Arg) i=Arg;\
		temp = 0;\
		for (int j=0; j < numBands; j++){\
			__shared__ float currAT;\
			if (threadIdx.x == 0){\
				currAT = AT[i*numBands + j];\
			}\
			__syncthreads();\
			float inp = inputmua[pixind + gridDim.x*pitch*(j+startind)];\
			temp += currAT*inp;\
		}\
		x##Arg = 0;\
		/*temp = ((i != currMethbind) || (hasMethb))*(temp);*/\
		mu##Arg = -1*temp
		/*mu_sh[Arg][threadIdx.x] = -1*temp*/

#define CalcX(Arg) j=Arg;\
		xprev = x##Arg;\
		muprev = mu##Arg;\
		/*muprev = mu_sh[Arg][threadIdx.x];*/\
		temp = xprev - fdividef(muprev, sh_H[j][j]);\
		temp = (temp > 0)*temp;\
		x##Arg = temp

#define CalcMu(Arg) k=Arg;\
		addition = (temp - xprev)*sh_H[k][j];\
		mu##Arg += addition
		/*mu_sh[Arg][threadIdx.x] += addition*/

#define updateMu CalcMu(0);\
		CalcMu(1);\
		if (MAX_ENDMEMBERS >= 3)\
			CalcMu(2);\
		if (MAX_ENDMEMBERS >= 4)\
			CalcMu(3);\
		if (MAX_ENDMEMBERS >= 5)\
			CalcMu(4);\
		if (MAX_ENDMEMBERS >= 6)\
			CalcMu(5);\
		if (MAX_ENDMEMBERS >= 7)\
			CalcMu(6)
		
//se opp for makrohelvete, dette ble skrevet for å se om det ble noe raskere om jeg brukte av registrene istedet for shared memory men DESSVERRE var det ikke mulig å opprette arrays i registrene uten like mye kuk som dette. ISRA-artikkelen lyver.  
template<int MAX_ENDMEMBERS>
__global__ void SCAFast(float *AT, float *H, float *inputmua, float *x, float *mu, int startind, int numBands, int numEndmembers, size_t pitch){
	//calculate initial values for mu and x
	int pixind = pitch*blockIdx.x + threadIdx.x;
	
	float temp;
	float x0, x1, x2, x3, x4, x5, x6;
	float mu0, mu1, mu2, mu3, mu4, mu5, mu6;
	int i;

	//calculate initial coordinates
	ATmultMua(0);
	ATmultMua(1);
	if (MAX_ENDMEMBERS >= 3)
		ATmultMua(2);
	if (MAX_ENDMEMBERS >= 4)
		ATmultMua(3);
	if (MAX_ENDMEMBERS >= 5)
		ATmultMua(4);
	if (MAX_ENDMEMBERS >= 6)
		ATmultMua(5);
	if (MAX_ENDMEMBERS >= 7)
		ATmultMua(6);
	

	__shared__ float sh_H[MAX_ENDMEMBERS][MAX_ENDMEMBERS];
	if (threadIdx.x == 0){
		for (int i=0; i < numEndmembers; i++){
			for (int j=0; j < numEndmembers; j++){
				sh_H[i][j] = H[i*numEndmembers + j];
			}
		}	
	}
	__syncthreads();


	//iterate
	for (int i=0; i < NUM_ITERATIONS; i++){
		float xprev, muprev, temp, addition;
		int k,j;
		CalcX(0);
		updateMu;

		CalcX(1);
		updateMu;

		if (MAX_ENDMEMBERS >= 3){
			CalcX(2);
			updateMu;
		}

		if (MAX_ENDMEMBERS >= 4){
			CalcX(3);
			updateMu;
		}
		
		if (MAX_ENDMEMBERS >= 5){
			CalcX(4);
			updateMu;
		}
		
		if (MAX_ENDMEMBERS >= 6){
			CalcX(5);
			updateMu;
		}
	
		if (MAX_ENDMEMBERS >= 7){
			CalcX(6);
			updateMu;
		}	
	}
	
	x[pixind + 0*gridDim.x*pitch] = x0;
	x[pixind + 1*gridDim.x*pitch] = x1;
	if (MAX_ENDMEMBERS >= 3)
		x[pixind + 2*gridDim.x*pitch] = x2;
	if (MAX_ENDMEMBERS >= 4)
		x[pixind + 3*gridDim.x*pitch] = x3;
	if (MAX_ENDMEMBERS >= 5)
		x[pixind + 4*gridDim.x*pitch] = x4;
	if (MAX_ENDMEMBERS >= 6)
		x[pixind + 5*gridDim.x*pitch] = x5;
	if (MAX_ENDMEMBERS >= 7)
		x[pixind + 6*gridDim.x*pitch] = x6;
}
