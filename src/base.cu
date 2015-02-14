//=======================================================================================================
// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//=======================================================================================================



#include <iostream>
#include <absorption_props.h>
#include <melanin.h>
#include "gpudm.h"
using namespace std;




void gpuskinbase_initialize(GPUSkinBase *base, int samples, int bands, size_t byteWidth, int height, float *wavelengths){
	int arrayByteSize = sizeof(float)*samples*bands;
	size_t pitch;
	
	//wavelength dependencies of scattering and absorption coefficients
	cudaMallocPitch(&(base->muh_oxy), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->muh_deoxy), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->melanin_base), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->musm), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->musr), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->musb_base), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->gcol), &pitch, byteWidth, height);
	cudaMallocPitch(&(base->wavelengths), &pitch, byteWidth, height);


	float *muh_oxy = base->muh_oxy;
	float *muh_deoxy = base->muh_deoxy;
	float *melanin_base = base->melanin_base;
	float *musm = base->musm;
	float *musr = base->musr;
	float *musb_base = base->musb_base;
	float *gcol = base->gcol;

	//calculate the spectrum arrays on the host, try to transfer concurrently by swapping between initiating the cudaMemCpy and doing calculations on the CPUi
	float *muh_oxy_host = (float*)malloc(arrayByteSize);
	float *muh_deoxy_host = (float*)malloc(arrayByteSize);
	float *melanin_base_host = (float*)malloc(arrayByteSize);
	float *gcol_host = (float*)malloc(arrayByteSize);
	float *musm_host = (float*)malloc(arrayByteSize);
	float *musr_host = (float*)malloc(arrayByteSize);
	float *musb_base_host = (float*)malloc(arrayByteSize);
	float *wavelengths_host = (float*)malloc(arrayByteSize);

	//prepare the wavelength dependent parts of the absorption and scattering coefficients
	//it might seem inefficient to assign the exact same value to the whole range of spatial positions for each wavelength, but even if it is all a memory waste, it will pay off in speed later on because of memory coalescing in the GPU. Also tried to load it into shared memory across all the thread blocks, but had to make the threads diverge in order to do that and would not work well with multiple shared variables
	//Would be the exact same latency issues. whether 32 threads load 32 different variables in one go and in parallell, or one of the thread loads the whole thing into the cache and picks one variable and broadcasts it to the rest will give the exact same latency, plus it will have to stall threads. 
	float lambda;
	int position;
	float muh_oxy_temp, muh_deoxy_temp, gcol_temp, musm_temp, musr_temp, musb_temp, melanin_temp;
	int spatial_size = samples;
	for (int i=0; i < bands; i++){
		lambda = wavelengths[i];
		muh_oxy_temp = muh_oxy_calc(lambda);
		muh_deoxy_temp = muh_deoxy_calc(lambda);
		gcol_temp = 0.62 + lambda*29e-5;
		musm_temp = (1-1.745e-3*lambda + 9.843e-7*lambda*lambda)/(1-gcol_temp);
		musr_temp = pow(lambda, -4);
		musb_temp = pow((685/(lambda*1.0)), 0.37);
		melanin_temp = melanin(lambda);
		for (int j=0; j < samples; j++){
			position = i*spatial_size + j;
			muh_oxy_host[position] = muh_oxy_temp;
			muh_deoxy_host[position] = muh_deoxy_temp;
			melanin_base_host[position] = melanin_temp;
			gcol_host[position] = gcol_temp;
			musm_host[position] = musm_temp;
			musr_host[position] = musr_temp;
			musb_base_host[position] = musb_temp/(1-gcol_temp)*(1-0.996);
		}
		wavelengths_host[i] = lambda;

	}

	//copy the wavelength dependencies to the GPU device
	cerr << cudaMemcpy2D(muh_oxy, pitch, muh_oxy_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(muh_deoxy, pitch, muh_deoxy_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(melanin_base, pitch, melanin_base_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(musm, pitch, musm_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(musr, pitch, musr_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(musb_base, pitch, musb_base_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(gcol, pitch, gcol_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);
	cerr << cudaMemcpy2D(base->wavelengths, pitch, wavelengths_host, byteWidth, byteWidth, height, cudaMemcpyHostToDevice);

	//free arrays that are no longer needed
	free(muh_oxy_host);
	free(muh_deoxy_host);
	free(melanin_base_host);
	free(musm_host);
	free(musr_host);
	free(musb_base_host);
	free(gcol_host);
	free(wavelengths_host);
}

void gpuskindata_initialize(GPUSkinData *data, size_t byteWidth, int height, int thrPerB){
	size_t *pitch = &(data->pitch);

	cerr << cudaMallocPitch(&(data->oxy), pitch, byteWidth, height);
	cerr << cudaMallocPitch(&(data->bvf), pitch, byteWidth, height);
	cerr << cudaMallocPitch(&(data->muam694), pitch, byteWidth, height);
	cerr << cudaMallocPitch(&(data->melanin_type), pitch, byteWidth, height);
	data->byteWidth = byteWidth;
	data->height = height;
}

void gpuopticalprops_initialize(GPUOpticalProps *props, size_t byteWidth, int height){
	size_t pitch;
	cerr << cudaMallocPitch(&(props->muae), &pitch, byteWidth, height);
	cerr << cudaMallocPitch(&(props->muse), &pitch, byteWidth, height);
	cerr << cudaMallocPitch(&(props->muad), &pitch, byteWidth, height);
	cerr << cudaMallocPitch(&(props->musd), &pitch, byteWidth, height);
}

void gpuskindata_download_from_gpu(GPUSkinData *data, int samples, float *outputHostArray){
	float *muamhost = (float*)malloc(data->byteWidth*data->height);
	float *bvfhost = (float*)malloc(data->byteWidth*data->height);
	float *oxyhost = (float*)malloc(data->byteWidth*data->height);
	cerr << cudaMemcpy2D(muamhost, data->byteWidth, data->muam694, data->pitch, data->byteWidth, data->height, cudaMemcpyDeviceToHost);
	cerr << cudaMemcpy2D(bvfhost, data->byteWidth, data->bvf, data->pitch, data->byteWidth, data->height, cudaMemcpyDeviceToHost);
	cerr << cudaMemcpy2D(oxyhost, data->byteWidth, data->oxy, data->pitch, data->byteWidth, data->height, cudaMemcpyDeviceToHost);
	
	for (int i=0; i < samples; i++){
		cout << muamhost[i] << " ";
	}
	cout << endl;

	memcpy(outputHostArray, bvfhost, samples*sizeof(float));
	memcpy(outputHostArray+samples, oxyhost, samples*sizeof(float));

	free(muamhost);
	free(bvfhost);
	free(oxyhost);
}

void gpuskinbase_free(GPUSkinBase *base){
	cudaFree(base->muh_oxy);
	cudaFree(base->muh_deoxy);
	cudaFree(base->melanin_base);
	cudaFree(base->gcol);
	cudaFree(base->musm);
	cudaFree(base->musr);
	cudaFree(base->musb_base);
	cudaFree(base->wavelengths);
}

void gpuskindata_free(GPUSkinData *data){
	cudaFree(data->muam694);
	cudaFree(data->bvf);
	cudaFree(data->oxy);
	cudaFree(data->melanin_type);
}

void gpuopticalprops_free(GPUOpticalProps *props){
	cudaFree(props->muae);
	cudaFree(props->muad);
	cudaFree(props->musd);
	cudaFree(props->muse);
}
