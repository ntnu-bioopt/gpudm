#include "gpudm.h"
#include <iostream>
#include <absorption_props.h>
#include <fstream>
#include <iostream>
#include "inv.h"
#include <ctime>
#include <parseconfig.h>
#include "wlenhelper.h"
#include "sca.cu"

#include "misc.h"
using namespace std;


void gpudm_free(GPUDMParams *params){
	gpuskindata_free(params->skinData);
	gpuopticalprops_free(params->optProps);
	gpuskinbase_free(params->skinBase);
	
	cudaFree(params->currRefl);
	cudaFree(params->melMeth_melcurve);
	cudaFree(params->melMeth_gpumeltype);

	gnlsq_free(params->lsq_w450);
	gnlsq_free(params->lsq_w530);
	gnlsq_free(params->lsq_w700);
	gnlsq_free(params->melMeth_lsq);
}

void gpudm_initialize_useconstants(GPUDMParams *params, int samples, int bands, float *wlens){
	//int numIntervals = 4;
	//_forwardRefls = new float[bands*numIntervals];
	params->image_samples = samples;
	params->image_bands = bands;

	//melanin estimation properties
	float startwlen = 730;
	float endwlen = 820;
	
	params->melMeth_numIterations = 2;
	params->melMeth_useCurveFromIterNumber = 3;
	params->melMeth_startWlenInd = getStartInd(wlens, bands, startwlen);
	params->melMeth_endWlenInd = getEndInd(wlens, bands, endwlen);
	
	//standard values for skin properties, for deinitialization
	params->default_muam = 100;
	params->default_bvf = 0.01;
	params->default_oxy = 0.8;

	//device array properties
	params->gpu_threadsPerBlock = 160;
	params->gpu_arrayByteWidth = params->gpu_threadsPerBlock*sizeof(float);
	params->gpu_arrayHeight = params->image_samples*params->image_bands/params->gpu_threadsPerBlock;

	//create buffers
	cudaMallocPitch(&params->currRefl, &(params->gpu_pitch), params->gpu_arrayByteWidth, params->gpu_arrayHeight);

	//skindata arrays
	params->skinBase = new GPUSkinBase;
	gpuskinbase_initialize(params->skinBase, params->image_samples, params->image_bands, params->gpu_arrayByteWidth, params->gpu_arrayHeight, wlens);

	params->skinData = new GPUSkinData;
	gpuskindata_initialize(params->skinData, params->gpu_arrayByteWidth, samples/params->gpu_threadsPerBlock, params->gpu_threadsPerBlock);

	params->optProps = new GPUOpticalProps;
	gpuopticalprops_initialize(params->optProps, params->gpu_arrayByteWidth, params->gpu_arrayHeight);

	//absorption coefficient estimation properties
	params->absfit_numIterations = 15;
	
	//melanin unmixing properties
	float *melcurveHost = new float[params->image_samples*params->image_bands];
	params->melMeth_factor = 0;
	for (int i=params->melMeth_startWlenInd; i < params->melMeth_endWlenInd; i++){
		float temp = pow(694.0f/wlens[i], 3.46);
		params->melMeth_factor += temp*temp;
		for (int j=0; j < params->image_samples; j++){
			melcurveHost[i*params->image_samples + j] = temp;
		}
	}
	params->melMeth_factor = 1.0f/params->melMeth_factor;
	cudaMallocPitch(&(params->melMeth_melcurve), &(params->gpu_pitch), params->gpu_arrayByteWidth, params->gpu_arrayHeight);
	cudaMemcpy2D(params->melMeth_melcurve, params->gpu_pitch, melcurveHost, params->gpu_arrayByteWidth, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyHostToDevice);
	delete [] melcurveHost;


	//chromophores to be used in the various fitting intervals and methods
	Chromophores melchrom;
	melchrom.setMel();
	melchrom.setWat();
	
	Chromophores chrom450;
	chrom450.setMel();
	chrom450.setBil();
	chrom450.setBet();
	chrom450.setKonst();
	
	Chromophores chrom530;
	chrom530.setMel();
	chrom530.setKonst();

	Chromophores chrom700;
	chrom700.setMel();
	chrom700.setWat();
	chrom700.setFat();
	chrom700.setKonst();


	//lsqfitting parameters for each wavelength interval
	params->lsq_w450 = new GNLSQParams;
	params->lsq_w530 = new GNLSQParams;
	params->lsq_w700 = new GNLSQParams;


	gnlsq_initialize(params->lsq_w450, params->gpu_threadsPerBlock, samples, bands, chrom450, wlens, 460, 530);
	
	gnlsq_initialize(params->lsq_w530, params->gpu_threadsPerBlock, samples, bands, chrom530, wlens, 510, 590);
	
	gnlsq_initialize(params->lsq_w700, params->gpu_threadsPerBlock, samples, bands, chrom700, wlens, 690, 820);

	//lsqfitting parameters for the melanin method
	params->melMeth_lsq = new GNLSQParams;
	gnlsq_initialize(params->melMeth_lsq, params->gpu_threadsPerBlock, samples, bands, melchrom, wlens, startwlen, endwlen);
}


void gpudm_initialize(GPUDMParams *params, int samples, int bands, float *wlens){
	//int numIntervals = 4;
	//_forwardRefls = new float[bands*numIntervals];
	params->image_samples = samples;
	params->image_bands = bands;

	//melanin estimation properties
	float startwlen = readNumberSetting("inversion/melanin/wavelengthInterval", "start");
	float endwlen = readNumberSetting("inversion/melanin/wavelengthInterval", "end");
	
	params->melMeth_numIterations = readNumberSetting("inversion/melanin", "iterations");
	params->melMeth_useCurveFromIterNumber = readNumberSetting("inversion/melanin", "useMelaninCurveFromIterNum");
	params->melMeth_startWlenInd = getStartInd(wlens, bands, startwlen);
	params->melMeth_endWlenInd = getEndInd(wlens, bands, endwlen);
	
	//standard values for skin properties, for deinitialization
	params->default_muam = readNumberSetting("inversion/melanin/initialValues", "muam694");
	params->default_bvf = readNumberSetting("inversion/melanin/initialValues", "bvf");
	params->default_oxy = readNumberSetting("inversion/melanin/initialValues", "oxy");

	//device array properties
	params->gpu_threadsPerBlock = readNumberSetting("inversion", "threadsPerBlock");
	params->gpu_arrayByteWidth = params->gpu_threadsPerBlock*sizeof(float);
	params->gpu_arrayHeight = params->image_samples*params->image_bands/params->gpu_threadsPerBlock;

	//create buffers
	cudaMallocPitch(&params->currRefl, &(params->gpu_pitch), params->gpu_arrayByteWidth, params->gpu_arrayHeight);

	//skindata arrays
	params->skinBase = new GPUSkinBase;
	gpuskinbase_initialize(params->skinBase, params->image_samples, params->image_bands, params->gpu_arrayByteWidth, params->gpu_arrayHeight, wlens);

	params->skinData = new GPUSkinData;
	gpuskindata_initialize(params->skinData, params->gpu_arrayByteWidth, samples/params->gpu_threadsPerBlock, params->gpu_threadsPerBlock);

	params->optProps = new GPUOpticalProps;
	gpuopticalprops_initialize(params->optProps, params->gpu_arrayByteWidth, params->gpu_arrayHeight);

	//absorption coefficient estimation properties
	params->absfit_numIterations = 15;
	
	//melanin unmixing properties
	float *melcurveHost = new float[params->image_samples*params->image_bands];
	params->melMeth_factor = 0;
	for (int i=params->melMeth_startWlenInd; i < params->melMeth_endWlenInd; i++){
		float temp = pow(694.0f/wlens[i], 3.46);
		params->melMeth_factor += temp*temp;
		for (int j=0; j < params->image_samples; j++){
			melcurveHost[i*params->image_samples + j] = temp;
		}
	}
	params->melMeth_factor = 1.0f/params->melMeth_factor;
	cudaMallocPitch(&(params->melMeth_melcurve), &(params->gpu_pitch), params->gpu_arrayByteWidth, params->gpu_arrayHeight);
	cudaMemcpy2D(params->melMeth_melcurve, params->gpu_pitch, melcurveHost, params->gpu_arrayByteWidth, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyHostToDevice);
	delete [] melcurveHost;


	//chromophores to be used in the various fitting intervals and methods
	Chromophores melchrom = parseChromophoreString(readStringSetting("inversion/melanin/chromophores", "chromophore"));
	Chromophores chrom450 = parseChromophoreString(readStringSetting("inversion/dermis/interval450/chromophores", "chromophore"));
	Chromophores chrom530 = parseChromophoreString(readStringSetting("inversion/dermis/interval530/chromophores", "chromophore"));	
	Chromophores chrom700 = parseChromophoreString(readStringSetting("inversion/dermis/interval700/chromophores", "chromophore"));	


	//lsqfitting parameters for each wavelength interval
	params->lsq_w450 = new GNLSQParams;
	params->lsq_w530 = new GNLSQParams;
	params->lsq_w700 = new GNLSQParams;


	gnlsq_initialize(params->lsq_w450, params->gpu_threadsPerBlock, samples, bands, chrom450, wlens, readNumberSetting("inversion/dermis/interval450/wavelengthInterval", "start"), readNumberSetting("inversion/dermis/interval450/wavelengthInterval", "end"));
	
	gnlsq_initialize(params->lsq_w530, params->gpu_threadsPerBlock, samples, bands, chrom530, wlens, readNumberSetting("inversion/dermis/interval530/wavelengthInterval", "start"), readNumberSetting("inversion/dermis/interval530/wavelengthInterval", "end"));
	
	gnlsq_initialize(params->lsq_w700, params->gpu_threadsPerBlock, samples, bands, chrom700, wlens, readNumberSetting("inversion/dermis/interval700/wavelengthInterval", "start"), readNumberSetting("inversion/dermis/interval700/wavelengthInterval", "end"));

	//lsqfitting parameters for the melanin method
	params->melMeth_lsq = new GNLSQParams;
	gnlsq_initialize(params->melMeth_lsq, params->gpu_threadsPerBlock, samples, bands, melchrom, wlens, startwlen, endwlen);
}


void gpudm_reinitialize(GPUDMParams *params, float *reflectance){
	//reinitialize skin data arrays before parameter estimation
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = 1;
	dim3 dimGrid = dim3(block_number_x, block_number_y);
	dim3 dimBlock = dim3(params->gpu_threadsPerBlock);
	size_t inputPitch = params->gpu_pitch/sizeof(float);
	skindataDeinit<<<dimGrid, dimBlock>>>(params->skinData->muam694, params->skinData->bvf, params->skinData->oxy, params->skinData->melanin_type, params->default_muam, params->default_bvf, params->default_oxy, inputPitch);

	//upload input reflectance to GPU
	cudaMemcpy2D(params->currRefl, params->gpu_pitch, reflectance, params->gpu_arrayByteWidth, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyHostToDevice);
}

void gpudm_fit_reflectance(GPUDMParams *params, float *reflectance){
	//upload reflectance and deinitialize
	gpudm_reinitialize(params, reflectance);
	
	//find required melanin
	gpudm_estimate_melanin(params);
	
	//estimate muad
	gpudm_calc_opticalprops(params, 0, params->image_bands-1);
	gpudm_estimate_muad(params, 0, params->image_bands-1);

	//fit muad in the defined wavelength intervals
	gnlsq_fit_sca(params->lsq_w450, params->optProps->muad);
	gnlsq_fit_sca(params->lsq_w530, params->optProps->muad);
	gnlsq_fit_sca(params->lsq_w700, params->optProps->muad);
}

void gpudm_estimate_melanin(GPUDMParams *params){
	float *muam694 = params->skinData->muam694;
	float *muae = params->optProps->muae;
	float *muad = params->optProps->muad;

	int startWlenInd = params->melMeth_startWlenInd;
	int endWlenInd = params->melMeth_endWlenInd;

	for (int i=0; i < params->melMeth_numIterations; i++){
		//one run of the melanin estimation method

		//estimate dermal absorption coefficient muad
		gpudm_calc_opticalprops(params, startWlenInd, endWlenInd);
		gpudm_estimate_muad(params, startWlenInd, endWlenInd);
	
		//fit muad
		gnlsq_fit_sca(params->melMeth_lsq, muad);

		//update skindata to get new bvf, oxy
		gnlsq_update_skindata(params->melMeth_lsq, params->skinData);

		//estimate muad from fitted coefficients, ignore melanin
		gnlsq_reconstruct_absorption(params->melMeth_lsq, muad, params->melMeth_lsq->fitting_chrom.getMelInd());

		//estimate epidermal absorption coefficient muae
		gpudm_estimate_muae(params, startWlenInd, endWlenInd);

		//gpu parameters for unmixing of melanin
		int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
		int block_number_y = 1;
		size_t inputPitch = params->gpu_pitch/sizeof(float);
		dim3 dimGrid(block_number_x, block_number_y);
		dim3 dimBlock(params->gpu_threadsPerBlock);

		//fit melanin to epidermal absorption coefficient using:
		if (i < params->melMeth_useCurveFromIterNumber){
			//straight line
			StraightLine<<<dimGrid, dimBlock>>>(params->skinBase->wavelengths, muae, MELANIN_REFERENCE_WAVELENGTH, muam694, params->melMeth_startWlenInd, params->melMeth_endWlenInd, inputPitch);
		} else {
			//or the melanin absorption curve
			MultVector<<<dimGrid, dimBlock>>>(params->melMeth_melcurve, muae, muam694, params->melMeth_factor, params->melMeth_startWlenInd, params->melMeth_endWlenInd, inputPitch);
		}
	}
}


void gpudm_download_bandarray(GPUDMParams *params, float *gpuarray, float *output){
	cudaMemcpy2D(output, params->gpu_arrayByteWidth, gpuarray, params->gpu_pitch, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyDeviceToHost);
}


void gpudm_estimate_muad(GPUDMParams *params, int startWlenInd, int endWlenInd){
	float *inputRefl = params->currRefl;
	float *muad = params->optProps->muad;

	//gpu grid 
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = endWlenInd - startWlenInd;
	dim3 dimGrid = dim3(block_number_x, block_number_y);
	dim3 dimBlock = dim3(params->gpu_threadsPerBlock);
	size_t inputPitch = params->gpu_pitch/sizeof(float);

	//estimate muad
	ReflIsoL2InvertMuad<<<dimGrid, dimBlock>>>(params->optProps->muae, params->optProps->muse, muad, params->optProps->musd, params->skinBase->gcol, inputRefl, inputPitch, startWlenInd);
}


void gpudm_calc_opticalprops(GPUDMParams *params, int startWlenInd, int endWlenInd){
	//gpu grid 
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = endWlenInd - startWlenInd;
	dim3 dimGrid = dim3(block_number_x, block_number_y);
	dim3 dimBlock = dim3(params->gpu_threadsPerBlock);
	size_t inputPitch = params->gpu_pitch/sizeof(float);

	//calculate initial optical properties given current skin data
	calcSkinData<<<dimGrid, dimBlock>>>(params->skinBase->wavelengths, params->skinData->oxy, params->skinData->bvf, params->skinData->muam694, params->skinData->melanin_type, params->optProps->muae, params->optProps->muse, params->optProps->muad, params->optProps->musd, params->skinBase->muh_oxy, params->skinBase->muh_deoxy, params->skinBase->melanin_base, params->skinBase->musm, params->skinBase->musr, params->skinBase->musb_base, inputPitch, startWlenInd);
	
}


void gpudm_estimate_muae(GPUDMParams *params, int startWlenInd, int endWlenInd){
	float *inputRefl = params->currRefl;
	float *muae = params->optProps->muae;

	#ifdef INV_MUAE_BLOCKDIM_64_QUICKFIX
	//the invmuae kernel has a lot of computations: will not be able to run on older gpus. Quickfix for decreasing number of threads per block (see also inside the kernel)
	int threadsPerBlock = 64;
	int block_number_x = params->image_samples/threadsPerBlock;
	int block_number_y = endWlenInd - startWlenInd;
	
	dim3 dimBlock = dim3(threadsPerBlock);
	dim3 dimGrid = dim3(block_number_x, block_number_y);
	#else 
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = endWlenInd - startWlenInd;
	dim3 dimGrid = dim3(block_number_x, block_number_y);
	dim3 dimBlock = dim3(params->gpu_threadsPerBlock);
	#endif

	size_t inputPitch = params->gpu_pitch/sizeof(float);

	ReflIsoL2InvertMuae<<<dimGrid, dimBlock>>>(muae, params->optProps->muse, params->optProps->muad, params->optProps->musd, params->skinBase->gcol, inputRefl, inputPitch, startWlenInd);
}


void gpudm_forward_simulation_singleres(GPUDMParams *params, GNLSQParams *gnlsq, float *refl){
	//create new lsqfitting object for all available wavelengths
	float startwlen = gnlsq->image_wlens[0];
	float endwlen = gnlsq->image_wlens[gnlsq->image_bands-1];
	GNLSQParams lsq_wall;
	gnlsq_initialize(&lsq_wall, params->gpu_threadsPerBlock, params->image_samples, params->image_bands, gnlsq->fitting_chrom, gnlsq->image_wlens, startwlen, endwlen);

	//copy fitting result from input lsqfitting object to object containing all wavelengths
	cudaMemcpy2D(lsq_wall.fitting_res, lsq_wall.gpu_pitch, gnlsq->fitting_res, gnlsq->gpu_pitch, gnlsq->gpu_arrayByteWidth, gnlsq->gpu_arrayHeight, cudaMemcpyDeviceToDevice);

	//generate muad
	gnlsq_reconstruct_absorption(&lsq_wall, params->optProps->muad);

	//run forward simulation (with previously set mus, muam, etc)
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = params->image_bands;
	dim3 dimGrid = dim3(block_number_x, block_number_y);
	dim3 dimBlock = dim3(params->gpu_threadsPerBlock);
	size_t inputPitch = params->gpu_pitch/sizeof(float);
	ReflIsoL2<<<dimGrid, dimBlock>>>(params->optProps->muae, params->optProps->muse, params->optProps->muad, params->optProps->musd, params->skinBase->gcol, refl, inputPitch, 0);

	gnlsq_free(&lsq_wall);
}

#define NUM_INTERVALS 3
void gpudm_forward_simulation(GPUDMParams *params, float **outputRefl_host){
	float *refl = params->currRefl;
	*outputRefl_host = NULL;
	*outputRefl_host = new float[params->image_samples*params->image_bands*NUM_INTERVALS];
	
	//forward simulation 450nm  parameters
	gpudm_forward_simulation_singleres(params, params->lsq_w450, refl);
	cudaMemcpy2D(*outputRefl_host, params->gpu_arrayByteWidth, refl, params->gpu_pitch, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyDeviceToHost);

	//forward simulation using current 530 nm parameters
	gpudm_forward_simulation_singleres(params, params->lsq_w530, refl);
	cudaMemcpy2D(*outputRefl_host + params->image_samples*params->image_bands, params->gpu_arrayByteWidth, refl, params->gpu_pitch, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyDeviceToHost);

	//forward simulation using current 700 nm parameters
	gpudm_forward_simulation_singleres(params, params->lsq_w700, refl);
	cudaMemcpy2D(*outputRefl_host + 2*params->image_samples*params->image_bands, params->gpu_arrayByteWidth, refl, params->gpu_pitch, params->gpu_arrayByteWidth, params->gpu_arrayHeight, cudaMemcpyDeviceToHost);

}

void gpudm_download_melanin(GPUDMParams *params, float *host_muam){
	cudaMemcpy2D(host_muam, params->skinData->byteWidth, params->skinData->muam694, params->skinData->pitch, params->skinData->byteWidth, params->skinData->height, cudaMemcpyDeviceToHost);
}

void gpudm_download_muad(GPUDMParams *params, float *host_muad){
	gpudm_download_bandarray(params, params->optProps->muad, host_muad);
}

void gpudm_download_450res(GPUDMParams *params, float *host_res){
	gnlsq_download_res(params->lsq_w450, host_res);
}

void gpudm_download_530res(GPUDMParams *params, float *host_res){
	gnlsq_download_res(params->lsq_w530, host_res);
}

void gpudm_download_700res(GPUDMParams *params, float *host_res){
	gnlsq_download_res(params->lsq_w700, host_res);
}

void gnlsq_free(GNLSQParams *params){
	cudaFree(params->SCA_GPUATA);
	cudaFree(params->SCA_GPUA);
	cudaFree(params->SCA_mu);
	cudaFree(params->fitting_res);
}




		


void gnlsq_initialize(GNLSQParams *params, int threadsPerBlock, int samples, int bands, Chromophores chrom, float *wlens, float startWlen, float endWlen){
	//image properties
	params->image_wlens = wlens;
	params->image_samples = samples;
	params->image_bands = bands;

	//gpu properties
	params->gpu_threadsPerBlock = threadsPerBlock;
	params->gpu_arrayByteWidth = sizeof(float)*threadsPerBlock;
	params->gpu_arrayHeight = samples*chrom.getNumEndmembers()/threadsPerBlock;


	//fitting parameters	
	params->fitting_numEndmembers = chrom.getNumEndmembers();
	params->fitting_startWlenInd = getStartInd(wlens, bands, startWlen);
	params->fitting_endWlenInd = getEndInd(wlens, bands, endWlen);
	params->fitting_numWlens = params->fitting_endWlenInd - params->fitting_startWlenInd;
	params->fitting_chrom = chrom; 


	//SCA specific parameters

	//absorption matrix
	int numWlens = params->fitting_numWlens;
	int numEndmembers = params->fitting_numEndmembers;
	float **chromMat = new float*[numWlens];
	for (int i=0; i < numWlens; i++){
		float wlen = wlens[params->fitting_startWlenInd + i];
		chromMat[i] = params->fitting_chrom.getAbsArray(wlen);
	}

	//H matrix in SCA fit
	float *ATA = new float[numEndmembers*numEndmembers];
	for (int i=0; i < numEndmembers; i++){
		for (int j=0; j < numEndmembers; j++){
			ATA[i*numEndmembers + j] = 0;
		}
	}
	for (int i=0; i < numEndmembers; i++){
		for (int k=0; k < numEndmembers; k++){
			float temp = 0;
			for (int j=0; j < numWlens; j++){
				temp += chromMat[j][i]*chromMat[j][k];
			}
			ATA[k*numEndmembers + i] = temp;
		}	
	}

	//transfer ATA to GPUATA on gpu
	cudaMalloc(&(params->SCA_GPUATA), sizeof(float)*numEndmembers*numEndmembers);
	cudaMemcpy(params->SCA_GPUATA, ATA, sizeof(float)*numEndmembers*numEndmembers, cudaMemcpyHostToDevice);
	delete [] ATA;

	//transfer absorption matrix A to GPU
	float *chromHost = new float[numEndmembers*numWlens];
	for (int i=0; i < numEndmembers; i++){
		for (int j=0; j < numWlens; j++){
			chromHost[i*numWlens + j] = chromMat[j][i];
		}
	}


	cudaMalloc(&(params->SCA_GPUA), sizeof(float)*numEndmembers*numWlens);
	cudaMemcpy(params->SCA_GPUA, chromHost, sizeof(float)*numEndmembers*numWlens, cudaMemcpyHostToDevice);


	//allocate mu variable for SCA
	size_t pitch;
	cudaMallocPitch(&(params->SCA_mu), &pitch, sizeof(float)*threadsPerBlock, numEndmembers*samples/threadsPerBlock);
	
	//allocate variable for storing the fitting results
	cudaMallocPitch(&(params->fitting_res), &pitch, sizeof(float)*threadsPerBlock, numEndmembers*samples/threadsPerBlock);
	params->gpu_pitch = pitch;

	//cleanup
	delete [] chromHost;
	for (int i=0; i < numWlens; i++){
		delete [] chromMat[i];
	}
	delete [] chromMat;
}



//do fitting using SCA
void gnlsq_fit_sca(GNLSQParams *params, float *absorption){
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = 1;
	dim3 dimBlock(params->gpu_threadsPerBlock);
	dim3 dimGrid(block_number_x, block_number_y);
	size_t inputPitch = params->gpu_pitch/sizeof(float);
	switch(params->fitting_numEndmembers){
		case 1:
		cerr << "ERROR! Only one endmember in SCA unmixing, you are doing something wrong." << endl;
		exit(1);
		case 2:
		SCAFast<2><<<dimGrid, dimBlock>>>(params->SCA_GPUA, params->SCA_GPUATA, absorption, params->fitting_res, params->SCA_mu, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, inputPitch);
		break;
		case 3:
		SCAFast<3><<<dimGrid, dimBlock>>>(params->SCA_GPUA, params->SCA_GPUATA, absorption, params->fitting_res, params->SCA_mu, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, inputPitch);
		break;
		case 4:
		SCAFast<4><<<dimGrid, dimBlock>>>(params->SCA_GPUA, params->SCA_GPUATA, absorption, params->fitting_res, params->SCA_mu, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, inputPitch);
		break;
		case 5:
		SCAFast<5><<<dimGrid, dimBlock>>>(params->SCA_GPUA, params->SCA_GPUATA, absorption, params->fitting_res, params->SCA_mu, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, inputPitch);
		break;
		case 6:
		SCAFast<6><<<dimGrid, dimBlock>>>(params->SCA_GPUA, params->SCA_GPUATA, absorption, params->fitting_res, params->SCA_mu, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, inputPitch);
		break;
		case 7:
		SCAFast<7><<<dimGrid, dimBlock>>>(params->SCA_GPUA, params->SCA_GPUATA, absorption, params->fitting_res, params->SCA_mu, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, inputPitch);
		break;
		default:
		cerr << "ERROR! This program was created in an infernally idiotic way. You need to add more cases in the switch statement in SCAFitting::doFitting, you tried to fit more parameters than was defined by the switch statement. (This is needed because of just in time compilation and the need for constants in the initialization of the CUDA kernel)" << endl;
		exit(1);
	}
}


void gnlsq_reconstruct_absorption(GNLSQParams *params, float *absorption, int ignoreIndex){
	//grid boilerplate
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = params->fitting_numWlens;
	dim3 dimBlock(params->gpu_threadsPerBlock);
	dim3 dimGrid(block_number_x, block_number_y);
	size_t inputPitch = params->gpu_pitch/sizeof(float);
		
	addUpAbsorption<<<dimGrid, dimBlock>>>(absorption, params->fitting_res, params->SCA_GPUA, params->fitting_startWlenInd, params->fitting_numWlens, params->fitting_numEndmembers, ignoreIndex, inputPitch);
}


void gnlsq_update_skindata(GNLSQParams *params, GPUSkinData *skinData){
	//grid boilerplate
	int block_number_x = params->image_samples/params->gpu_threadsPerBlock;
	int block_number_y = 1;
	dim3 dimBlock(params->gpu_threadsPerBlock);
	dim3 dimGrid(block_number_x, block_number_y);
	size_t inputPitch = params->gpu_pitch/sizeof(float);
		
		
	calcOxyBvf<<<dimGrid, dimBlock>>>(params->fitting_res, skinData->oxy, skinData->bvf, inputPitch);
	
}


void gnlsq_download_res(GNLSQParams *params, float *res){
	cudaMemcpy2D(res, sizeof(float)*params->gpu_threadsPerBlock, params->fitting_res, params->gpu_pitch, params->gpu_threadsPerBlock*sizeof(float), params->fitting_numEndmembers*params->image_samples/params->gpu_threadsPerBlock, cudaMemcpyDeviceToHost);
}

