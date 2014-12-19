
//invert muae from input reflectance
__global__ void ReflIsoL2InvertMuae(float *muae, float *muse, float *muad, float *musd, float *gcol, float *lineData, size_t pitch, int startblockInd);
__global__ void calcSkinData(float *wlens, float *oxy_arr, float *Bd_arr, float *muam694_arr, float *melanintype_arr, float *muae, float *muse, float *muad, float *musd, 
				float *muh_oxy, float *muh_deoxy, float *melanin_base, float *musm, float *musr, float *musb_base, size_t pitch, int startblockind);
__global__ void MultVector(float *multVec, float *mua, float *res, float factor, int startwlenind, int endwlenind, size_t pitch);
__global__ void ReflIsoL2InvertMuad(float *muae, float *muse, float *muad, float *musd, float *gcol, float *lineData, size_t pitch, int startblockind);
__global__ void ReflIsoL2(float *muae, float *muse, float *muad, float *musd, float *gcol, float *res, size_t pitch, int startblockind);
__device__ float ReflIsoL2Calc(float mua1, float musr1, float mua2, float musr2);
__global__ void ReflIsoL2ErrorCheck(float *muae, float *muse, float *musd, float *gcol, float *AT, float *x, int endmembers, int numbands, float *inputres, float *outputres, size_t pitch, int startblockind, int diff);
__global__ void StraightLine(float *wavelengths, float *mua, float wlen, float *res, int startwlenInd, int endwlenInd, size_t inputPitch);
