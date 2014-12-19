
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <hyPipe/DataObjectBase.h>
#include <hyPipe/DataObjectFactory.h>
#include <hyPipe/ProcessingStage.h>
#include <hyPipe/TypeChecker.h>
#include <hyPipe/InputPort.h>
#include <hyPipe/OutputPort.h>
#include <hyLight/ImageMetadata.h>
#include <string>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <vector>

class RGBVisualize;
class CalibrationStage : public hypipe::ProcessingStage, public hypipe::TypeChecker, public hypipe::DataObjectFactory{
	public:
		CalibrationStage(int bands, int samples, int linesPerBlock, std::string specificationFile, RGBVisualize *rgbStage, std::vector<float> wlens);
		std::string checkType(hypipe::DataObjectBase* d);
		void addFloatImageConnection(hypipe::InputPort *p){_outputPort->addConnection(p);};
		hypipe::InputPort* getImageInputPort(){return _inputPort;};
		void createPorts();
		virtual hypipe::DataObjectBase* constructData();
		
		void detectan(float *inputImage, int height, int width, int *x1, int *x2, int *y1, int *y2);
		float* sumCalArray_with99slab(float *totImg, int startxCalSlab, int endxCalSlab, int startyCalSlab, int endyCalSlab);

		void stats(float *line, float *mean, float *stdev);

		void divideCalArrayBySpecs(float *calArray);
		void updateCalArray(float *lineData, float *calArray, int num, float *diff);
	protected:
		virtual void execute();
	private:
		std::string _specificationFile;
		bool _shouldFindSlab;	
		hypipe::InputPort *_inputPort; //input port, float image
		hypipe::OutputPort *_outputPort; //output port, float image

		int _samples; //samples in the output image
		int _imageInputSamples; //samples in the input image
		int _linesPerBlock;
		int _bands;

		int _stdBand;

		float* _calarray;
		std::vector<float> _wlens;
		float _calfactor;
		int _threadsPerBlock;
		float *_buffer;
		RGBVisualize *rgbVis; //used for some inter-stage communication: will get messages of read lines when reflectance standard has been detected and summed
};
