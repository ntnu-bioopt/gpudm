//=======================================================================================================
// Copyright 2015 Asgeir Bjorgan, Lise Lyngsnes Randeberg, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//=======================================================================================================

#ifndef READIMAGE_H_DEFINED
#define READIMAGE_H_DEFINED
#include <vector>
#include <string>

typedef struct {
	int samples;
	int bands;
	int lines;
	int offset;
	std::vector<float> wlens;
	int datatype;
} HyspexHeader;

//for specifying subset of hyperspectral image. (samp -> cols, line -> rows)  
typedef struct {
	int startSamp;
	int endSamp;
	int startLine;
	int endLine;
} ImageSubset;

//set contents of supplied HyspexHeader to information contained in [filename, extension removed].hdr
void hyperspectral_read_header(char *filename, HyspexHeader *header);

//read image in filename according to information supplied in header, subset, into float array supplied in data. 
void hyperspectral_read_image(char *filename, HyspexHeader *header, ImageSubset subset, float *data);

//write header and image. 
void hyperspectral_write_header(const char *filename, int numBands, int numPixels, int numLines, std::vector<float> wlens, std::string description, std::vector<std::string> bandnames);
void hyperspectral_write_image(const char *filename, int bands, int samples, int lines, float *data);


#endif
