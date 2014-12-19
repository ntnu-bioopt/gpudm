#include "melanin.h"

#define div13 1.0f/3.0f
#include "inv.h"
#define A 0.14386f
//just be aware that this A is actually incorrect. The diffuse reflectance should also be calculated using Rd = gamma*(1-Rsp), but be also aware that the boundary conditions are just an approximation. This choice of A gives us a reasonably good approximation to monte carlo simulations. 

//#define A 0.17f
#define d1 0.0001f
#define de 0.0001f
#define NUM_ITERATIONS 15
__global__ void ReflIsoL2InvertMuae(float *muae, float *muse, float *muad, float *musd, float *gcol, float *lineData, size_t pitch, int startblockInd){
	#ifdef INV_MUAE_BLOCKDIM_64_QUICKFIX
	//NB FIXME: indeksuthentingen er bygget for en mindre blokkdimensjon enn 160
	int egBlockInd = (blockIdx.x*64 + threadIdx.x)/160; //indeksen til blokkene slik de egentlig var ment å være i originalgriddet
	int ind = ((gridDim.x*(blockIdx.y+startblockInd))*pitch*2)/5 + egBlockInd*pitch + blockIdx.x*64-egBlockInd*32*5 + threadIdx.x; //the pitch is made for the largest possible block size, but this uses too many registers for that block size to be valid and needs something less. Therefore some index hackin'
	#else
	int ind = (gridDim.x*(blockIdx.y+startblockInd) + blockIdx.x)*pitch+ threadIdx.x;
	#endif

	float musr2 = gcol[ind]; //ps dette er ikke musr, det er g FIXME

	//reduced scattering coefficients
	float musr1 = muse[ind]*(1.0f-musr2);
	musr2 = musd[ind]*(1.0f-musr2);

	//move mua into shared memory
	float mua1 = muae[ind];
	float mua2 = muad[ind];
	
	//diffusion constant
		
	float D2 = fdividef(1.0f,3.0f*(musr2 + mua2));
	float del2 = sqrtf(fdividef(D2,mua2));

	float res;
	float derivmuae;
	float currLineData = lineData[ind];
	float f2 = 1.0f + fdividef(del2 ,(D2 * 3.0f));
	float f6 = D2-del2*fdividef(del2, (D2*9.0f));

	for (int i=0; i < NUM_ITERATIONS; ++i){
		float D1 = fdividef(1.0f,3.0f*(musr1 + mua1));
		float del1 = sqrtf(fdividef(D1,mua1));
		float div1D1D1 = fdividef(1.0f, D1*D1);

		float coshval = coshf(fdividef(d1, del1));
		float sinhval = sinhf(fdividef(d1, del1));

		//result
		float f1 = (del1*del1 * fdividef(del2 , 3.0f) - del1*del1 * D2) * coshval + (powf(del1, 3.0f) * fdividef(D2 , D1*3.0f) - del1 * del2 * D1) * sinhval;
		float f3 = fdividef(musr2, musr1) * del2 * del2 * D1 * (del1*del1 * fdividef(div1D1D1,9.0f) - 1.0f) + del1*del1 * f6;
		float f4 = expf(-fdividef(d1 , (D1 *3.0f)));
		float f5 = del1*del1 * fdividef(div1D1D1,9.0f) - 1.0f;
		float f7 = D1 * del1 * (D2 + del2 * A) * coshval + (del2 * D1*D1 + D2 * del1*del1 * A) * sinhval;

		float fact = fdividef((f1 * f2 + f3 * f4) , (f5 * f2 * f7));

		res = del1 * musr1 * A * fact;
		
		
		float dD1d1 = -3.0f*D1*D1;
		float ddel1d1 = (dD1d1*mua1-D1)*fdividef(1.0f, mua1*mua1)*fdividef(1.0f, 2.0f*del1);

		
		float df1d1 = (fdividef(2.0f , 3.0f) * del1 * del2 * ddel1d1 - 2.0f * del1 * D2 * ddel1d1) * coshval - (del1*del1 * fdividef(del2 , 3.0f) - del1*del1 * D2) * sinhval * d1 * fdividef(1.0f, del1*del1) * ddel1d1 + (del1*del1 * fdividef(D2 , D1) * ddel1d1 - powf(del1, 3.0f) * D2 * (div1D1D1) * fdividef(dD1d1 , 3.0f) - ddel1d1 * del2 * D1 - del1 * del2 * dD1d1) * sinhval - (powf(del1, 3.0f) * fdividef(D2 , (D1 * 3.0f)) - del1 * del2 * D1) * coshval * d1 * fdividef(1.0f, del1*del1) * ddel1d1;
		float df3d1 = fdividef(musr2 , musr1) * del2 * del2 * dD1d1 * (del1*del1 * fdividef(div1D1D1,9.0f) - 1.0f) + fdividef(musr2 , musr1) * del2 * del2 * D1 * (fdividef(2.0f , 9.0f) * del1 * (div1D1D1) * ddel1d1 - fdividef(2.0f , 9.0f) * del1*del1 * powf(D1, -3.0f) * dD1d1) + 2.0f * del1 * (f6) * ddel1d1;
		float df4d1 = d1 * (div1D1D1) * dD1d1 * f4 *fdividef(1.0f, 3.0f);
		float df5d1 = fdividef(2.0f , 9.0f) * del1 * (div1D1D1) * ddel1d1 - fdividef(2.0f , 9.0f) * del1*del1 * powf(D1, -3.0f) * dD1d1;
		float df7d1 = dD1d1 * del1 * (D2 + del2 * A) * coshval + D1 * ddel1d1 * (D2 + del2 * A) * coshval - fdividef(D1 , del1) * (D2 + del2 * A) * sinhval * d1 * ddel1d1 + (2.0f * del2 * D1 * dD1d1 + 2.0f * D2 * del1 * A * ddel1d1) * sinhval - (del2 * D1*D1 + D2 * del1*del1 * A) * coshval * d1 * fdividef(1.0f, del1*del1) * ddel1d1;

		derivmuae = ddel1d1 * musr1 * A * fact + del1 * musr1 * A * (df1d1 * f2 + df3d1 * f4 + f3 * df4d1) *fdividef(1.0f, (f5 * f2 * f7)) - del1 * musr1 * A * fact * fdividef(1.0f, f5) * df5d1 - del1 * musr1 * A * fact * fdividef(1.0f, f7) * df7d1;
		//newton's method
		mua1 = mua1 - fdividef(res-currLineData, derivmuae);

		//correction in case muad wants to be negative, which we seriously don't want
		mua1 = mua1*(1-signbit(mua1)) + signbit(mua1);
	}
	muae[ind] = mua1;
}

//takes in pre-allocated arrays and the arrays containing the bases of the different absorption coefficients, fills the skin data arrays with the optical properties
__global__ void calcSkinData(float *wlens, float *oxy_arr, float *Bd_arr, float *muam694_arr, float *melanintype_arr, float *muae, float *muse, float *muad, float *musd, 
				float *muh_oxy, float *muh_deoxy, float *melanin_base, float *musm, float *musr, float *musb_base, size_t pitch, int startblockind){
	//walk down the lines, walk along the blocks, walk along the threads inside the block
	int index = (gridDim.x*(blockIdx.y+startblockind) + blockIdx.x)*pitch + threadIdx.x;

	//absorption properties
	float H = 0.41;
	float H0 = 0.45;
	float Be = 0.002;


	int chromInd = blockIdx.x*pitch + threadIdx.x;
	float oxy = oxy_arr[chromInd];
	float Bd = Bd_arr[chromInd];
	float muam694 = muam694_arr[chromInd];
	int melanintype = melanintype_arr[chromInd];	


	float mua_other = 25; //FIXME
	float muab_blood = (muh_oxy[index]*oxy + muh_deoxy[index]*(1-oxy))*fdividef(H,H0);

	__shared__ float wlen;
	if (threadIdx.x == 0){
		wlen = wlens[blockIdx.y+startblockind];
	}
	__syncthreads();

	float mua_melanin = muam694*((melanintype == SVAASAND_MELANIN_GPU)*powf(fdividef(694.0f,wlen), 3.46) + (melanintype == EUMELANIN_GPU)*expf(-kEu*fdividef(wlen-694.0f,694.0f)) + (melanintype == PHEOMELANIN_GPU)*expf(-kPheo*fdividef(wlen-694.0f,694.0f)));
	muae[index] = mua_melanin + muab_blood*Be + mua_other*(1-Be);
	muad[index] = muab_blood*Bd + mua_other*(1-Bd);

	//scattering properties
	float c_ray = 1.05e12;
	float c_mie = 105;
	float must = musm[index]*c_mie*100 + musr[index]*c_ray*100;
	
	//float f = 0.64;
	float aMie = 18.780, bMie = 0.22, aRay = 17.6;
	must = 100*(aMie*pow(wlen/500.0f, -bMie) + aRay*pow(wlen/500, -4));
	float gcol = 0.62 + wlen*29e-5;
	must /= (1-gcol);

	//float b = 0.91;
	//must = 3000*(f*powf(wlen/500.0f,-4) + (1-f)*powf(wlen/500.0f, -b))/(1-(0.62 + wlen*29e-5));
	
	float musb685 = 55.09e-12;
	float ve = 1.25e-16;
	float musb = musb685*H*(1-H)*(1.4-H)*fdividef(1.0f,ve)*musb_base[index];
	muse[index] = must;//*(1-Be);//+musb*Be;
	musd[index] = must;//*(1-Bd);//+musb*Bd;
}

__global__ void test(float *muae, float *muse, float *muad, float *musd, float *gcol, float *lineData, size_t pitch, int startblockind){
	int ind = (gridDim.x*(blockIdx.y+startblockind) + blockIdx.x)*pitch+ threadIdx.x;
	float musr2 = gcol[ind]; //ps dette er ikke musr, det er g FIXME

	//reduced scattering coefficients
	float musr1 = muse[ind]*(1.0f-musr2);
	musr2 = musd[ind]*(1.0f-musr2);

	//move mua into local memory
	float mua1 = muae[ind];
	float mua2 = muad[ind];

	float currLineData = lineData[ind];
	float res;
	float derivmuad;
	
	//diffusion constant
	float D1 = fdividef(1.0f,3.0f*(musr1 + mua1));

	float musr2dmusr1 = fdividef(musr2,musr1); //musr2 divided by musr1, keep for derivative calc


	//optical penetration depth
	float del1 = sqrtf(fdividef(D1,mua1));
	float sinhval = sinhf(fdividef(de, del1));
	float coshval = coshf(fdividef(de, del1));
	float expval = expf(-fdividef(de,D1)*div13);

	float D2 = fdividef(1.0f,3.0f*(musr2 + mua2));
	float del2 = sqrtf(fdividef(D2,mua2));

	//from Svaasand 1995	
	//calculate the reflectance value
	float f1 = (del1*del1*del2*div13-del1*del1*D2)*coshval+(del1*del1*del1*fdividef(D2,D1)*div13 - del1*del2*D1)*sinhval; //keeping for derivative calc
	float f2 = 1.0f+fdividef(del2, D2)*div13; //keeping for derivative calc
	float f3 = (musr2dmusr1*del2*del2*D1*(del1*fdividef(del1,D1*D1)*div13*div13-1.0f)+del1*del1*(D2-del2*fdividef(del2,D2)*div13*div13))*expval; //keeping for derivative calc
	float f4 = del1*fdividef(del1,D1*D1)*div13*div13 - 1.0f; //keep for derivative calc, f4
	float f5 = fdividef(del2,D2)*div13+1.0f; //keep for derivative calc, f5
	float f6 = D1*del1*(D2+del2*A)*coshval+(D1*D1*del2 + D2*del1*del1*A)*sinhval; //keep for derivative calc, f6
	float num = del1*musr1*A*(f1*f2+f3);
	float denom = f4*f5*f6;
	res = fdividef(num, denom);

	//calculate the derivative with respect to muad
	float dD2dmuad = -3.0f*D2*D2;
	float ddel2dmuad = (dD2dmuad*mua2-D2)*fdividef(1.0f, mua2*mua2)*fdividef(1.0f, 2.0f*del2);
	float df2dmuad = div13*(ddel2dmuad*D2 - del2*dD2dmuad)*fdividef(1.0f, D2*D2);
	derivmuad = (del1*musr1*A*((coshval*(del1*del1*div13*ddel2dmuad - del1*del1*dD2dmuad) + sinhval*(del1*del1*del1*fdividef(1.0f, D1)*div13*dD2dmuad - del1*D1*ddel2dmuad))*f2 + f1*df2dmuad + expval*(musr2dmusr1*2.0f*del2*ddel2dmuad*D1*(fdividef(del1*del1,9.0f*D1*D1)-1.0f) + del1*del1*(dD2dmuad + fdividef(1.0f, 9.0f*mua2*mua2))))*denom - (f4*(df2dmuad*f6 + D1*del1*(dD2dmuad + ddel2dmuad*A)*coshval + (D1*D1*ddel2dmuad + dD2dmuad*del1*del1*A)*sinhval*f5)*num))*fdividef(1.0f, denom*denom);

	//calculate next mua2
	//newton's method
	mua2 = mua2 - fdividef(res-currLineData, derivmuad);

	//correction in case muad wants to be negative, which we seriously don't want
	mua2 = mua2*(1-signbit(mua2)) + signbit(mua2);
	
}


#define BLOCK_DIM 160
__global__ void ReflIsoL2InvertMuad(float *muae, float *muse, float *muad, float *musd, float *gcol, float *lineData, size_t pitch, int startblockind){
	int ind = (gridDim.x*(blockIdx.y+startblockind) + blockIdx.x)*pitch+ threadIdx.x;
	float musr2 = gcol[ind]; //ps dette er ikke musr, det er g FIXME

	//reduced scattering coefficients
	float musr1 = muse[ind]*(1.0f-musr2);
	musr2 = musd[ind]*(1.0f-musr2);

	//move mua into local memory
	float mua1 = muae[ind];
	float mua2 = muad[ind];

	float currLineData = lineData[ind];
	float res;
	float derivmuad;
	
	//diffusion constant
	float D1 = fdividef(1.0f,3.0f*(musr1 + mua1));

	float musr2dmusr1 = fdividef(musr2,musr1); //musr2 divided by musr1, keep for derivative calc


	//optical penetration depth
	float del1 = sqrtf(fdividef(D1,mua1));
	float sinhval = sinhf(fdividef(de, del1));
	float coshval = coshf(fdividef(de, del1));
	float expval = expf(-fdividef(de,D1)*div13);

	for (int i=0; i < NUM_ITERATIONS; ++i){
		//result
		//res = (musr1 * sqrtf(1.0f / (musr1 + mua1) / mua1) * sqrtf(3.0f) * ((-musr1 / 2.0f + mua1) * (3 * musr2 + mua2) * mua2 * sqrtf((1 / (3 * musr2 + mua2) / mua2)) + (3.0f / 2.0f * musr2 - mua2) * mua1 - 3.0f / 2.0f * musr1 * mua2) * sinh(de * sqrtf(3.0f) * powf(1.0f / (musr1 + mua1) / mua1, -1.0f / 2.0f)) + 4.0f * (-3.0f / 8.0f * musr2 + mua2) * musr1 * cosh(de * sqrtf(3.0f) * powf(1.0f / (musr1 + mua1) / mua1, -1.0f / 2.0f)) + expf(-de * (musr1 + mua1)) * (3.0f * musr2 * mua1 - 4.0f * musr1 * mua2)) * sqrtf(1.0f / (musr1 + mua1) / mua1) * mua1 * sqrtf(3.0f) * (musr1 + mua1) * A / (-musr1 / 2.0f + mua1) / (3.0f + sqrtf((1 / (3 * musr2 + mua2) / mua2)) * (3 * musr2 + mua2)) / ((mua1 * (3 * musr2 + mua2) * sqrtf((1 / (3 * musr2 + mua2) / mua2)) + 3.0f * A * (musr1 + mua1)) * sinh(de * sqrtf(3.0f) * powf(1.0f / (musr1 + mua1) / mua1, -1.0f / 2.0f)) + sqrtf(1.0f / (musr1 + mua1) / mua1) * mua1 * sqrtf(3.0f) * (musr1 + mua1) * cosh(de * sqrtf(3.0f) * powf(1.0f / (musr1 + mua1) / mua1, -1.0f / 2.0f)) * (1.0f + A * (3 * musr2 + mua2) * sqrtf((1 / (3 * musr2 + mua2) / mua2)))) / mua2;

		//derivative
		//derivmuad = -powf((3 * musr2 + mua2), -1.0f / 2.0f) * A * musr2 * (expf(-de * (musr1 + mua1)) * (3.0f * sqrtf(mua2) * (A * (musr1 + mua1) + mua1 / 9.0f + 4.0f / 9.0f * musr1) * mua1 * sqrtf((3 * musr2 + mua2)) + (A * (musr1 + mua1) + mua1) * (2.0f * musr1 * mua2 + (3.0f / 2.0f * musr2 + mua2) * mua1)) * sinh(de * sqrtf(3.0f) * sqrtf(musr1 + mua1) * sqrtf(mua1)) + 2.0f * (expf(-de * (musr1 + mua1)) * (2.0f / 3.0f * ((A / 4.0f + 3.0f / 4.0f) * powf(mua1, 3.0f / 2.0f) + musr1 * sqrtf(mua1) * A) * sqrtf(mua2) * sqrtf((3 * musr2 + mua2)) + (A + 1.0f / 3.0f) * ((3.0f / 4.0f * musr2 + mua2 / 2.0f) * powf(mua1, 3.0f / 2.0f) + mua2 * musr1 * sqrtf(mua1))) * cosh(de * sqrtf(3.0f) * sqrtf(musr1 + mua1) * sqrtf(mua1)) - 5.0f / 4.0f * (mua2 + 3.0f / 5.0f * sqrtf((3 * musr2 + mua2)) * sqrtf(mua2) + 3.0f / 0.10e2 * musr2) * (A + 1.0f / 3.0f) * musr1 * sqrtf(mua1)) * sqrtf(3.0f) * sqrtf(musr1 + mua1)) * sqrtf(3.0f) * sqrtf(musr1 + mua1) * powf(mua2, 3.0f / 2.0f) * sqrtf(mua1) * powf(mua2 + sqrtf((3 * musr2 + mua2)) * sqrtf(mua2) / 3.0f, -2.0f) / (-musr1 / 2.0f + mua1) * powf((mua1 * sqrtf((3 * musr2 + mua2)) * sqrtf(mua2) / 3.0f + A * mua2 * (musr1 + mua1)) * sinh(de * sqrtf(3.0f) * sqrtf(musr1 + mua1) * sqrtf(mua1)) + sqrtf(3.0f) * sqrtf(musr1 + mua1) * sqrtf(mua1) * cosh(de * sqrtf(3.0f) * sqrtf(musr1 + mua1) * sqrtf(mua1)) * (mua2 + sqrtf((3 * musr2 + mua2)) * sqrtf(mua2) * A) / 3.0f, -2.0f) / 9.0f;
		//res += currLineData*(i+1);
		//derivmuad += currLineData*2*(i+1);
		//res += gcol[ind];
		float D2 = fdividef(1.0f,3.0f*(musr2 + mua2));
		float del2 = sqrtf(fdividef(D2,mua2));

		//from Svaasand 1995	
		//calculate the reflectance value
		float f1 = (del1*del1*del2*div13-del1*del1*D2)*coshval+(del1*del1*del1*fdividef(D2,D1)*div13 - del1*del2*D1)*sinhval; //keeping for derivative calc
		float f2 = 1.0f+fdividef(del2, D2)*div13; //keeping for derivative calc
		float f3 = (musr2dmusr1*del2*del2*D1*(del1*fdividef(del1,D1*D1)*div13*div13-1.0f)+del1*del1*(D2-del2*fdividef(del2,D2)*div13*div13))*expval; //keeping for derivative calc
		float f4 = del1*fdividef(del1,D1*D1)*div13*div13 - 1.0f; //keep for derivative calc, f4
		float f5 = fdividef(del2,D2)*div13+1.0f; //keep for derivative calc, f5
		float f6 = D1*del1*(D2+del2*A)*coshval+(D1*D1*del2 + D2*del1*del1*A)*sinhval; //keep for derivative calc, f6
		float num = del1*musr1*A*(f1*f2+f3);
		float denom = f4*f5*f6;
		res = fdividef(num, denom);

		//calculate the derivative with respect to muad
		float dD2dmuad = -3.0f*D2*D2;
		float ddel2dmuad = (dD2dmuad*mua2-D2)*fdividef(1.0f, mua2*mua2)*fdividef(1.0f, 2.0f*del2);
		float df2dmuad = div13*(ddel2dmuad*D2 - del2*dD2dmuad)*fdividef(1.0f, D2*D2);
		derivmuad = (del1*musr1*A*((coshval*(del1*del1*div13*ddel2dmuad - del1*del1*dD2dmuad) + sinhval*(del1*del1*del1*fdividef(1.0f, D1)*div13*dD2dmuad - del1*D1*ddel2dmuad))*f2 + f1*df2dmuad + expval*(musr2dmusr1*2.0f*del2*ddel2dmuad*D1*(fdividef(del1*del1,9.0f*D1*D1)-1.0f) + del1*del1*(dD2dmuad + fdividef(1.0f, 9.0f*mua2*mua2))))*denom - (f4*(df2dmuad*f6 + D1*del1*(dD2dmuad + ddel2dmuad*A)*coshval + (D1*D1*ddel2dmuad + dD2dmuad*del1*del1*A)*sinhval*f5)*num))*fdividef(1.0f, denom*denom);

		//calculate next mua2
		//newton's method
		mua2 = mua2 - fdividef(res-currLineData, derivmuad);

		//correction in case muad wants to be negative, which we seriously don't want
		mua2 = mua2*(1-signbit(mua2)) + signbit(mua2);
	}
	muad[ind] = mua2;
}

//find straight line through input mua, return as the value at wavelength w 
__global__ void StraightLine(float *wavelengths, float *mua, float wlen, float *res, int startwlenind, int endwlenind, size_t pitch){
	float xbar = 0;
	float ybar = 0;
	float xybar = 0;
	float x2bar = 0;
	int num = 0;
	for (int i=startwlenind; i < endwlenind; i++){
		num++;
		int ind = (gridDim.x*i + blockIdx.x)*pitch + threadIdx.x;
		__shared__ float w;
		if (threadIdx.x == 0){
			w = wavelengths[i];
		}
		__syncthreads();
		//float w = wavelengths[i];
		float abs = mua[ind];
		xbar += w;
		ybar += abs;
		xybar += w*abs;
		x2bar += w*w;
	}
	xbar /= num;
	ybar /= num;
	xybar /= num;
	x2bar /= num;
	float a = (xybar - xbar*ybar)/(x2bar - xbar*xbar);
	float b = ybar - a*xbar;
	int ind = blockIdx.x*pitch + threadIdx.x;
	res[ind] = a*wlen + b;
}

__global__ void MultVector(float *multVec, float *mua, float *res, float factor, int startwlenind, int endwlenind, size_t pitch){
	int i=0;
	float restemp = 0;
	for (i=startwlenind; i < endwlenind; i++){
		int ind = (gridDim.x*i + blockIdx.x)*pitch + threadIdx.x;
		restemp += mua[ind]*multVec[ind];
	}
	restemp *= factor;
	int ind = blockIdx.x*pitch + threadIdx.x;
	res[ind] = restemp;
}


__global__ void ReflIsoL2(float *muae, float *muse, float *muad, float *musd, float *gcol, float *res, size_t pitch, int startblockind){
	int ind = (gridDim.x*(blockIdx.y+startblockind) + blockIdx.x)*pitch+ threadIdx.x;
	float musr2 = gcol[ind]; //ps dette er ikke musr, det er g FIXME

	//reduced scattering coefficients
	float musr1 = muse[ind]*(1.0f-musr2);
	musr2 = musd[ind]*(1.0f-musr2);

	//move mua into local memory
	float mua1 = muae[ind];
	float mua2 = muad[ind];

	res[ind] = ReflIsoL2Calc(mua1, musr1, mua2, musr2);
}


__device__ float ReflIsoL2Calc(float mua1, float musr1, float mua2, float musr2){
	//diffusion constant
	float D1 = fdividef(1.0f,3.0f*(musr1 + mua1));

	float musr2dmusr1 = fdividef(musr2,musr1); //musr2 divided by musr1, keep for derivative calc


	//optical penetration depth
	float del1 = sqrtf(fdividef(D1,mua1));
	float sinhval = sinhf(fdividef(de, del1));
	float coshval = coshf(fdividef(de, del1));
	float fact1 = del1*musr1*A;
	float expval = expf(-fdividef(de,D1)*div13);

	float D2 = fdividef(1.0f,3.0f*(musr2 + mua2));
	float del2 = sqrtf(fdividef(D2,mua2));

	//from Svaasand 1995	
	//calculate the reflectance value
	float f1 = (del1*del1*del2*div13-del1*del1*D2)*coshval+(del1*del1*del1*fdividef(D2,D1)*div13 - del1*del2*D1)*sinhval; //keeping for derivative calc
	float f2 = 1.0f+fdividef(del2, D2)*div13; //keeping for derivative calc
	float f3 = (musr2dmusr1*del2*del2*D1*(del1*fdividef(del1,D1*D1)*div13*div13-1.0f)+del1*del1*(D2-del2*fdividef(del2,D2)*div13*div13))*expval; //keeping for derivative calc
	float f4 = del1*fdividef(del1,D1*D1)*div13*div13 - 1.0f; //keep for derivative calc, f4
	float f5 = fdividef(del2,D2)*div13+1.0f; //keep for derivative calc, f5
	float f6 = D1*del1*(D2+del2*A)*coshval+(D1*D1*del2 + D2*del1*del1*A)*sinhval; //keep for derivative calc, f6
	float num = fact1*(f1*f2+f3);
	float denom = f4*f5*f6;
	return fdividef(num,denom);
}

__global__ void ReflIsoL2ErrorCheck(float *muae, float *muse, float *musd, float *gcol, float *AT, float *x, int endmembers, int numbands, float *inputres, float *outputres, size_t pitch, int startblockind, int diff){
	float error = 0;
	for (int i=0; i < numbands; i++){
		int ind = (gridDim.x*(i+startblockind) + blockIdx.x)*pitch+ threadIdx.x;
		float musr2 = gcol[ind]; //ps dette er ikke musr, det er g FIXME

		//reduced scattering coefficients
		float musr1 = muse[ind]*(1.0f-musr2);
		musr2 = musd[ind]*(1.0f-musr2);

		//calculate mua
		float mua1 = muae[ind];
		float mua2 = 0;
		float temp;
		for (int j=0; j < endmembers; j++){
			__shared__ float Aval;
			if (threadIdx.x == 0){
				Aval = AT[j*numbands + i + diff];
			}
			__syncthreads();
			mua2 += Aval*x[(gridDim.x*j + blockIdx.x)*pitch + threadIdx.x];
		}
		float res = ReflIsoL2Calc(mua1, musr1, mua2, musr2);
		float measval = inputres[ind];
		error += (res-measval)*(res-measval);
	}
	outputres[blockIdx.x*pitch + threadIdx.x] = sqrtf(error);
	
}
