//=======================================================================================================
//// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
//// Distributed under the MIT License.
//// (See accompanying file LICENSE or copy at
//// http://opensource.org/licenses/MIT)
////=======================================================================================================

#include "readimage.h"
#include "gpudm.h"
#include "parseconfig.h"
#include <sstream>

#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
	if (argc < 3){
		cout << "Usage: " << argv[0] << " infilename outfilename [path to config file (optional, configfile.xml in current folder otherwise assumed)]" << endl;
		exit(1);
	}

	std::string configpath = "configfile.xml";
	if (argc > 3){
		configpath = argv[3];
	}

	readConfigFile(configpath);


	//read image to BIL-interleaved float-array
	char *filename = argv[1];
	size_t offset;
	HyspexHeader header;
	hyperspectral_read_header(filename, &header);

	ImageSubset subset;
	subset.startSamp = 0;
	subset.endSamp = header.samples;
	subset.startLine = 0;
	subset.endLine = header.lines;

	float *data = new float[header.lines*header.samples*header.bands];
	hyperspectral_read_image(filename, &header, subset, data);

	float *wlens = new float[header.bands];
	for (int i=0; i < header.bands; i++){
		wlens[i] = header.wlens[i];
	}

	int lines = header.lines;
	int bands = header.bands;
	int samples = 1600; //due to some stupid decisions in threadsPerBlock and subsequent arrays, this is fixed regardless of image width

	//prepare gpudm arrays
	GPUDMParams params;
	gpudm_initialize(&params, samples, bands, wlens);

	//prepare output arrays
	float *res_530 = new float[samples*lines*params.lsq_w530->fitting_numEndmembers];
	float *res_700 = new float[samples*lines*params.lsq_w700->fitting_numEndmembers];
	float *res_mel = new float[samples*lines];
	float *res_muad = new float[samples*lines*bands];

	for (int i=0; i < lines; i++){
		float *lineReflEg = data + i*header.samples*bands;
		float *lineRefl = new float[samples*bands]();
		for (int j=0; j < bands; j++){
			for (int k=0; k < header.samples; k++){
				lineRefl[j*samples + k] = lineReflEg[j*header.samples + k];
			}
		}

		gpudm_fit_reflectance(&params, lineRefl);

		//result from interval around 500 nm
		float *res_line_530 = res_530 + samples*params.lsq_w530->fitting_numEndmembers*i;
		gpudm_download_530res(&params, res_line_530);

		//result from interval around 700 nm
		float *res_line_700 = res_700 + samples*params.lsq_w700->fitting_numEndmembers*i;
		gpudm_download_700res(&params, res_line_700);

		//melanin result
		float *res_line_mel = res_mel + samples*1*i;
		gpudm_download_melanin(&params, res_line_mel);

		//muad result
		float *res_line_muad = res_muad + samples*bands*i;
		gpudm_download_muad(&params, res_line_muad);

		delete [] lineRefl;
	}

	//write results to hyperspectral files
	//parameters are in the sequence defined by the Chromophores class (see libchromophoreconfig)

	//530nm interval
	string outfilename = string(argv[2]) + "_530res";
	vector<string> bandnames;
	vector<float> bandnums;
	ostringstream *description = new ostringstream;
	*description << "Dermal chromophore values extracted from " << wlens[params.lsq_w530->fitting_startWlenInd] << "-" << wlens[params.lsq_w530->fitting_endWlenInd] << "nm";
	for (int i=0; i < params.lsq_w530->fitting_numEndmembers; i++){
		bandnums.push_back(i);
		bandnames.push_back(params.lsq_w530->fitting_chrom.getName(i));
	}
	hyperspectral_write_header(outfilename.c_str(), params.lsq_w530->fitting_numEndmembers, samples, lines, bandnums, description->str(), bandnames);
	hyperspectral_write_image(outfilename.c_str(), params.lsq_w530->fitting_numEndmembers, samples, lines, res_530);

	//700nm interval
	outfilename = string(argv[2]) + "_700res";
	bandnames.clear();
	bandnums.clear();
	delete description;
	description = new ostringstream;
	*description << "Dermal chromophore values extracted from " << wlens[params.lsq_w700->fitting_startWlenInd] << "-" << wlens[params.lsq_w700->fitting_endWlenInd] << "nm";
	for (int i=0; i < params.lsq_w700->fitting_numEndmembers; i++){
		bandnums.push_back(i);
		bandnames.push_back(params.lsq_w700->fitting_chrom.getName(i));
	}
	hyperspectral_write_header(outfilename.c_str(), params.lsq_w700->fitting_numEndmembers, samples, lines, bandnums, description->str(), bandnames);
	hyperspectral_write_image(outfilename.c_str(), params.lsq_w700->fitting_numEndmembers, samples, lines, res_700);


	//melanin values
	outfilename = string(argv[2]) + "_melres";
	bandnames.clear();
	bandnums.clear();
	delete description;
	description = new ostringstream;
	*description << "Epidermal melanin amount in terms of absorption in m-1 at 694 nm.";
	bandnames.push_back("melanin");
	bandnums.push_back(0);
	hyperspectral_write_header(outfilename.c_str(), 1, samples, lines, bandnums, description->str(), bandnames);
	hyperspectral_write_image(outfilename.c_str(), 1, samples, lines, res_700);

	//muad
	outfilename = string(argv[2]) + "_muadres";
	bandnames.clear();
	bandnums.clear();
	delete description;
	description = new ostringstream;
	*description << "Derived dermal absorption coefficients.";
	for (int i=0; i < bands; i++){
		bandnums.push_back(wlens[i]);
		bandnames.push_back("wavelength");
	}
	hyperspectral_write_header(outfilename.c_str(), bands, samples, lines, bandnums, description->str(), bandnames);
	hyperspectral_write_image(outfilename.c_str(), bands, samples, lines, res_muad);



	delete [] res_530;
	delete [] res_700;
	delete [] res_mel;
	delete [] data;

	gpudm_free(&params);
}	
