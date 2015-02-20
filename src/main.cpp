#include "readimage.h"
#include "gpudm.h"
#include "parseconfig.h"
#include "hyperspectral.h"

#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
	if (argc < 3){
		cout << "Usage: " << argv[0] << " infilename outfilename" << endl;
		exit(1);
	}

	readConfigFile("configfile.xml");


	//read image to BIL-interleaved float-array
	char *filename = argv[1];
	size_t offset;
	HyspexHeader header;
	readHeader(filename, &header);
	
	ImageSubset subset;
	subset.startSamp = 0;
	subset.endSamp = header.samples;
	subset.startLine = 0;
	subset.endLine = header.lines;
	
	float *data = new float[header.lines*header.samples*header.bands];
	readImage(filename, &header, subset, data);

	float *wlens = new float[header.bands];
	for (int i=0; i < header.bands; i++){
		wlens[i] = header.wlens[i];
	}


	int lines = header.lines;
	int bands = header.bands;
	int samples = 1600;


	//prepare gpudm arrays
	GPUDMParams params;
	gpudm_initialize(&params, samples, bands, wlens);

	//prepare output arrays
	float *res_530 = new float[samples*lines*params.lsq_w530->fitting_numEndmembers];
	float *res_700 = new float[samples*lines*params.lsq_w700->fitting_numEndmembers];

	for (int i=0; i < lines; i++){
		float *lineReflOrig = data + i*header.samples*bands;
		float *lineRefl = new float[samples*bands]();
		for (int j=0; j < samples; j++){
			for (int k=0; k < bands; k++){
				lineRefl[k*samples + j] = lineReflOrig[k*header.samples + j];
			}
		}
	

		gpudm_fit_reflectance(&params, lineRefl);

		//result from interval around 500 nm
		float *res_line_530 = res_530 + samples*params.lsq_w530->fitting_numEndmembers*i;
		gpudm_download_530res(&params, res_line_530);

		//result from interval around 700 nm
		float *res_line_700 = res_700 + samples*params.lsq_w700->fitting_numEndmembers*i;
		gpudm_download_700res(&params, res_line_700);
		delete [] lineRefl;
	}

	//write results to hyperspectral files
	Hyperspectral *hyp530res = new Hyperspectral(res_530, lines, samples, params.lsq_w530->fitting_numEndmembers, BIL);
	hyp530res->writeToFile(string(argv[2]) + "_530res");

	Hyperspectral *hyp700res = new Hyperspectral(res_700, lines, samples, params.lsq_w700->fitting_numEndmembers, BIL);
	hyp700res->writeToFile(string(argv[2]) + "_700res");



	delete [] res_530;
	delete [] res_700;
	delete [] data;

	delete hyp530res;
	delete hyp700res;

	gpudm_free(&params);
}
