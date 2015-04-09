//=======================================================================================================
//// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
//// Distributed under the MIT License.
//// (See accompanying file LICENSE or copy at
//// http://opensource.org/licenses/MIT)
////=======================================================================================================

#ifndef GPUDM_H_DEFINED
#define GPUDM_H_DEFINED

#include <chromophores.h>
#include "melanin.h"

#include <stdlib.h>

//base spectral dependencies
typedef struct{
	float *muh_oxy;
	float *muh_deoxy;
	float *melanin_base; //contains (695/lambda)^3.46
	float *gcol; //anisotropy factor for collagen
	float *musm; //unreduced mie scattering
	float *musr; //unreduced rayleigh scattering coefficient
	float *musb_base; //unreduced blood scattering coefficient
	float *wavelengths; //wavelengths
} GPUSkinBase;

//FIXME: simplify parameter list
void gpuskinbase_initialize(GPUSkinBase *base, int samples, int bands, size_t byteWidthPerBlock, int numberOfBlocks, float *wavelengths);
void gpuskinbase_free(GPUSkinBase *base);

//skin data
typedef struct{
	float *muam694;
	float *bvf;
	float *oxy;
	float *melanin_type; //SVAASAND_MELANIN_GPU, EUMELANIN_GPU or PHEOMELANIN_GPU
	float *download(int samples); //download all arrays to an host array
	size_t pitch;
	int height;
	size_t byteWidth;
} GPUSkinData;

//FIXME: simplify parameter list (bytewidthperblock wtf)
void gpuskindata_initialize(GPUSkinData *data, size_t byteWidthPerBlock, int numberOfBlocks, int thrPerBl);
void gpuskindata_download_from_gpu(GPUSkinData *data, int samples, float *outputHostArray);
void gpuskindata_free(GPUSkinData *data);


//optical properties
typedef struct{
	float *muae;
	float *muse;
	float *muad;
	float *musd;
} GPUOpticalProps;

void gpuopticalprops_initialize(GPUOpticalProps *props, size_t byteWidthPerBlock, int numberOfBlocks);
void gpuopticalprops_free(GPUOpticalProps *props);

//parameters relevant for non-negative least squares fitting on the GPU
typedef struct{
	//SCA parameters
	float *SCA_GPUATA; //H matrix on the GPU
	float *SCA_GPUA; //transpose of chromophore matrix on GPU
	float *SCA_mu; //den mu-greia i algoritmen anae rikke hva den heter
	
	//fitting parameters	
	int fitting_startWlenInd; //start index for fitting
	int fitting_endWlenInd;
	int fitting_numEndmembers;
	int fitting_numWlens;
	float *fitting_res;
	Chromophores fitting_chrom;

	//image properties
	int image_samples;
	float *image_wlens;
	int image_bands;

	//GPU properties
	int gpu_threadsPerBlock;
	size_t gpu_pitch;
	size_t gpu_arrayByteWidth;
	int gpu_arrayHeight;
} GNLSQParams;


//contains parameters and GPU arrays relevant for GPU-DM (estimation of skin optical parameters using the GPU and diffusion model)
typedef struct{
	GPUSkinData *skinData; //skin data (muam694, oxy, bvf, ...)
	GPUOpticalProps *optProps; //optical properties (mus, mua, ...), either fitted or calculated
	GPUSkinBase *skinBase; //base wavelength dependencies of skin properties
	
	float *currRefl; //gpu allocated array for containing reflectance
	int absfit_numIterations; //number of iterations for newton's method, estimation of absorption

	//image properties
	int image_bands;
	int image_samples;
	float *image_wlens;

	//GPU array properties
	//this is for convenience in the handling of the full samples x wavelengths-arrays, and will not necessarily be correct for samples x 1 or samples x endmembers-arrays
	int gpu_threadsPerBlock;
	size_t gpu_arrayByteWidth;
	int gpu_arrayHeight;
	size_t gpu_pitch;

	//melanin fitting variables
	float melMeth_factor; //parameter for melanin curve fitting
	float *melMeth_melcurve; //parameter for melanin curve fitting
	float *melMeth_gpumeltype; //melanin type. 0: svaasand, 1: pheomelanin, 2: eumelanin
	int melMeth_numIterations; //number of iterations of the entire melanin method, with increasing melanin values as initial values
	int melMeth_useCurveFromIterNumber; //from which iteration of the method to use the melanin absorption curve for the fitting of the melanin
	
	//default reinitialization variables
	float default_bvf, default_muam, default_oxy;

	//wavelength range for melanin fitting	
	int melMeth_startWlenInd;
	int melMeth_endWlenInd;

	//non-negative least squares fitting parameters for the separate wavelength intervals
	//will also contain the resulting fitting values
	GNLSQParams *lsq_w450;
	GNLSQParams *lsq_w530;
	GNLSQParams *lsq_w700;

	//nnls fitting parameters for dermal absorption in the melanin fitting routines
	GNLSQParams *melMeth_lsq;
} GPUDMParams;

//initialize GPUDMParams
//reads chromophore choices etc from config file initialized elsewhere
void gpudm_initialize(GPUDMParams *params, int samples, int bands, float *wlens);

//doesn't read parameters from config file, uses constant parameters instead
void gpudm_initialize_useconstants(GPUDMParams *params, int samples, int bands, float *wlens);

//free gpu memory
void gpudm_free(GPUDMParams *params);
	
//upload reflectance and reinitialize skindata arrays to standard values
void gpudm_reinitialize(GPUDMParams *params, float *host_reflectance);

//entrance function: given an array of host allocated reflectance data, estimate skin optical parameters, all the way from epidermal melanin to dermal properties
void gpudm_fit_reflectance(GPUDMParams *params, float *host_reflectance);

//estimate melanin given parameters and reflectance set in params
void gpudm_estimate_melanin(GPUDMParams *params);

//estimate dermal absorption coefficient from current reflectance data, using current muae
void gpudm_estimate_muad(GPUDMParams *params, int startWlenInd, int endWlenInd);

//estimate epidermal absorption coefficient from current reflectance data, using current muad
void gpudm_estimate_muae(GPUDMParams *params, int startWlenInd, int endWlenInd);

//calculate optical properties (muae, muad, muse, musd, ...) given the current skin data
void gpudm_calc_opticalprops(GPUDMParams *params, int startWlenInd, int endWlenInd);

//download bandarray (i.e. reflectance array, muad array, etc) to output
void gpudm_download_bandarray(GPUDMParams *params, float *gpu_array, float *host_output);

//run forward simulation, allocate outputRefl_host and return in this array
void gpudm_forward_simulation(GPUDMParams *params, float **outputRefl_host);
void gpudm_forward_simulation_singleres(GPUDMParams *params, GNLSQParams *gnlsq, float *refl);

//functions for downloading results to host arrays
void gpudm_download_melanin(GPUDMParams *params, float *host_muam);
void gpudm_download_muad(GPUDMParams *params, float *host_muad);
void gpudm_download_450res(GPUDMParams *params, float *host_res);
void gpudm_download_530res(GPUDMParams *params, float *host_res);
void gpudm_download_700res(GPUDMParams *params, float *host_res);



//initialize parameters
void gnlsq_initialize(GNLSQParams *params, int threadsPerBlock, int samples, int bands, Chromophores chrom, float *wlens, float startWlen, float endWlen);

//free gpu arrays
void gnlsq_free(GNLSQParams *params);

//fit defined chromophores to input absorption using SCA
//output will be contained within params, use other functions to extract
void gnlsq_fit_sca(GNLSQParams *params, float *absorption);

//reconstruct absorption into input array from results contained in params. Ignore parameter number ignoreIndex if set
void gnlsq_reconstruct_absorption(GNLSQParams *params, float *gpu_absorption, int ignoreIndex = -1);

//update input skinData based on lsqfit results
void gnlsq_update_skindata(GNLSQParams *params, GPUSkinData *skinData);

//download lsqfit results to array
void gnlsq_download_res(GNLSQParams *params, float *host_res);









#endif
