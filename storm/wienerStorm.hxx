/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*      Copyright 2010-2011 by Joachim Schleicher and Ullrich Koethe    */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*    joachim.schleicher@iwr.uni-heidelberg.de                          */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#include <vigra/stdconvolution.hxx>
#include <vigra/convolution.hxx>
#include <vigra/recursiveconvolution.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/fftw3.hxx> 
#include <vigra/localminmax.hxx>
#include <vigra/splineimageview.hxx>
#include <vigra/multi_fft.hxx>

#include <ctime>
//#include <set>
#include <fstream>
//#include <iomanip>
#include <vector>
//#include <algorithm>
#ifdef OPENMP_FOUND
    #include <omp.h>
#endif //OPENMP_FOUND
#include "util.h"
#include "fftfilter.hxx"
#include "myimportinfo.h"
//#include <iostream>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>
//#include <opengm/operations/adder.hxx>
#include <opengm/inference/graphcut.hxx>
//#include <opengm/operations/minimizer.hxx>
#include <opengm/inference/auxiliary/minstcutboost.hxx>

#include <opengm/operations/adder.hxx>

#include <opengm/functions/truncated_squared_difference.hxx>
#include <opengm/inference/auxiliary/minstcutkolmogorov.hxx>

#include <Rmath.h>

using namespace vigra;
using namespace vigra::functor;

/**
 * @file dSTORM data processing for localization microscopy
 * 
 * This file contains functions to localizes single molecule 
 * point-spread-functions by interpolation of each image in an
 * image stack after prefiltering. A Wiener filter can be learned
 * from good-quality input data and afterwards lied to
 * low SNR-measurements.
 * 
 * For algorithmic details and performance measurements please refer to
 * the diploma thesis available at
 * http://hci.iwr.uni-heidelberg.de/Staff/jschleic/
 * 
 * Some helper functions can be found in util.h, and the filter in Fourier 
 * domain is applied using fftw-wrappers in fftfilter.h.
 * 
 * @date 2010-2011 Diploma thesis J. Schleicher
 */

//--------------------------------------------------------------------------
// helper classes and functions
//--------------------------------------------------------------------------

/**
 * BSpline coefficients but no prefiltering
 */

double factorial_log(double in){
	double erg = 0;
	for(int i = 1; i<=in;i++){
		erg +=log((double)i);
	}
	return erg;
}

template <class T>
void poissonDistribution_log(std::vector<T>& xin, std::vector<T>& yout, T lamb){
	if(xin.size()!=yout.size()){std::cout<<"Warning x and y vector must have the same length!"<<std::endl;}
	int numberElements = xin.size();
	int minVal = (int) xin[0];
	for(int i = 0 ; i<numberElements;i++){
		double p1 = xin[i] * log(lamb);//log(pow(lamb, xin[i]));
		//std::cout<<"p1: "<<p1<<" ";
		double p2 = factorial_log(xin[i]);
		//std::cout<<"p2: "<<p2<<" ";
		double p3 = -lamb;//exp(-lamb);
		//std::cout<<"p3: "<<p3<<" ";
		double p4 = p1 - p2;
		//std::cout<<"p4: "<<p4<<" ";
		//double p5 = p4;//exp(p4);
		//std::cout<<"p5: "<<p5<<" ";
		yout[i] = p4 + p3;//(p5+p3);
		//std::cout<<"p5 * p3: "<<p5 * p3<<" ";
	}
}

template <class T>
void poissonlastValue_log(T xmax, T& y_lastout, T lamb){
	xmax = (int)xmax;
	double p1 = xmax * log(lamb);//log(pow(lamb, xin[i]));
	//std::cout<<"p1: "<<p1<<" ";
	double p2 = factorial_log(xmax);
	//std::cout<<"p2: "<<p2<<" ";
	double p3 = -lamb;//exp(-lamb);
	//std::cout<<"p3: "<<p3<<" ";
	double p4 = p1 - p2;
	//std::cout<<"p4: "<<p4<<" ";
	//double p5 = p4;//exp(p4);
	//std::cout<<"p5: "<<p5<<" ";
	y_lastout = p4 + p3;//(p5+p3);
	//std::cout<<"p5 * p3: "<<p5 * p3<<" ";

}

double factorial(double in){
	double erg = 1;
	for(int i = 1; i<=in;i++){erg *=(double)i;}
	return erg;
}

size_t variableIndex(const size_t x, const size_t y, const size_t nx){return x + nx *y;}

template <class T>
void performGraphcut(MultiArray<2, T>& array, T lambda){
	using namespace opengm;
	using namespace std;

	const size_t nx = array.size(0);
	const size_t ny = array.size(1);
	const size_t numberOfLabels = 2;
	//double lambda = 0.1;


	typedef SimpleDiscreteSpace<size_t, size_t> Space;
	Space space(nx * ny, numberOfLabels);
	typedef GraphicalModel<double, Adder, OPENGM_TYPELIST_2(
			ExplicitFunction<double> ,PottsFunction<double>), Space> Model;
	double maximum = 1.05;//max2DArr(array);
	Model gm(space);
	for(size_t y = 0; y < ny; ++y){
		for(size_t x = 0; x < nx; ++x){
			const size_t shape[] = {numberOfLabels};
			ExplicitFunction<double> f(shape, shape + 1);

			//if(array(x,y)>0){f(0) = array(x,y); f(1) = 0.0;}
			//if(array(x,y)<=0){f(1) = -array(x,y); f(0) = 0.0;}
			f(0) = array(x,y);
			f(1) = (maximum - array(x,y));

			Model::FunctionIdentifier fid = gm.addFunction(f);
			size_t variableIndices[] = {variableIndex(x,y,nx)};
			gm.addFactor(fid, variableIndices, variableIndices +1);
		}
	}

	PottsFunction<double> f(numberOfLabels, numberOfLabels, 0.0, lambda);
	Model::FunctionIdentifier fid = gm.addFunction(f);


	for(size_t y = 0; y < ny; y++){
		for(size_t x = 0; x < nx; ++x){
			if(x+1 < nx){
				size_t variableIndices[] = {variableIndex(x,y,nx), variableIndex(x+1,y,nx)};
				sort(variableIndices, variableIndices + 2);
				gm.addFactor(fid, variableIndices, variableIndices + 2);
			}
			if(y + 1 < ny){
				size_t variableIndices[] = {variableIndex(x,y,nx), variableIndex(x,y+1,nx)};
				sort(variableIndices, variableIndices + 2);
				gm.addFactor(fid, variableIndices, variableIndices + 2);
			}
		}
	}

	typedef opengm::external::MinSTCutKolmogorov<size_t, double> MinStCutType;

	typedef opengm::GraphCut<Model, opengm::Minimizer, MinStCutType> MinCut;

	std::vector<Model::LabelType> v;
	{
	MinCut mincut(gm);
	mincut.infer();
	mincut.arg(v);
	}


	int counter = 0;
	//std::vector<Model::LabelType>::iterator it;
	//for(it = v.begin(); it != v.end(); it++){
	//	std::cout<<it<<" ";
	//}

	for(int i = 0; i< ny;i++){
		for(int j = 0; j< nx;j++){
			array(j,i) = v[counter];
			counter += 1;
		}
	}
}

template <int ORDER, class T>
class BSplineWOPrefilter
 : public BSpline<ORDER, T> {

    public:
    /** Prefilter coefficients
        (array has zero length, since image is already prefiltered).
    */
    ArrayVector<double> const & prefilterCoefficients() const
    {
        static ArrayVector<double> b;
        return b;
    }
};

/**
 * Class to keep an image coordinate with corresponding pixel value
 * 
 * This corresponds to a vigra::Point2D with an additional value at that coordinate.
 */
template <class VALUETYPE>
class Coord{
    public:
        typedef VALUETYPE value_type;
        Coord(const int x_,const int y_,const VALUETYPE val_,const VALUETYPE asym_=1., const VALUETYPE snr_ = 0)
            : x(x_), y(y_), val(val_), asymmetry(asym_), signalNoiseRatio(snr_) {  }
        int x;
        int y;
        VALUETYPE val;
        VALUETYPE asymmetry;
        VALUETYPE signalNoiseRatio;

        bool operator<(const Coord<VALUETYPE>& c2) const {
            return ((this->y==c2.y)&&(this->x < c2.x)) || (this->y < c2.y);
        }

        // Coords are equal, if they're at exactly the same position and have the same value
        bool operator==(const Coord<VALUETYPE>& c2) const {
            bool ret = (this->x == c2.x) && (this->y == c2.y) && (this->val == c2.val);
            return ret;
        }
};


// hack to push the coordinates into an array instead of marking them in 
// a target image.
// This is used as an accessor although it doesn't access the pixel values ;-)
// To work on ROIs, a global offset can be set with setOffset().
template <class T, class ITERATOR>
class VectorPushAccessor{
    public:
        typedef typename T::value_type value_type;
        VectorPushAccessor(std::set<T>& arr, ITERATOR it_start)
            : m_arr(arr), m_it_start(it_start), m_offset() {        }

        T const &   operator() (ITERATOR const &i) const {
            return NumericTraits<T>::zero();
        }
        template<class V>
        void    set (V const & /*value*/, ITERATOR const &i) {
            int x = i.x+m_offset.x;
            int y = i.y-m_it_start.y+m_offset.y;
            typename T::value_type val = *i;
            T c (x,y,val);
            m_arr.insert(c);
        }
        void setOffset(Diff2D offset) {
            m_offset = offset;
        }

    private:
        std::set<T>& m_arr;
        ITERATOR m_it_start;
        Diff2D m_offset;
};

/**
 * Draw coordinates from all frames into the result image
 */
template <class C, class Image>
void drawCoordsToImage(const std::vector<std::set<C> >& coords, Image& res) {
    res = 0;
    typename std::vector<std::set<C> >::const_iterator it;
    //  loop over the images
    for(it = coords.begin(); it != coords.end(); ++it) {
        drawCoordsToImage( *it, res);
    }
}

/**
 *  Draw coordinates detected in one frame into the resulting image
 */
template <class C, class Image>
void drawCoordsToImage(const std::set<C>& coords, Image& res) {
    //  loop over the coordinates
    typename std::set<C>::iterator it2;

    for(it2 = coords.begin(); it2 != coords.end(); it2++) {
        const C& c = *it2;
        res(c.x, c.y) += c.val;
    }
}

template <class C>
int saveCoordsFile(const std::string& filename, const std::vector<std::set<C> >& coords, 
            const MultiArrayShape<3>::type & shape, const int factor) {
    int numSpots = 0;
    std::set<Coord<float> >::iterator it2;
    std::ofstream cfile (filename.c_str());
    cfile << shape[0] << " " << shape[1] << " " << shape[2] << std::endl;
    cfile << std::fixed; // fixed instead of scientific format
    for(unsigned int j = 0; j < coords.size(); j++) {
        for(it2=coords[j].begin(); it2 != coords[j].end(); it2++) {
            numSpots++;
            const Coord<float>& c = *it2;
            cfile << std::setprecision(3) << (float)c.x/factor << " " << (float)c.y/factor << " "
                << j << " " << std::setprecision(1) << c.val << " " << std::setprecision(3) << c.asymmetry << " "
                << c.signalNoiseRatio << std::endl;
        }
    }
    cfile.close();
    return numSpots;
}

/** 
 * finds the value, so that the given percentage of pixels is above / below that value.
 */
template <class Image>
void findMinMaxPercentile(Image& im, double minPerc, double& minVal, double maxPerc, double& maxVal) {
    std::vector<typename Image::value_type> v;
    for(int y=0; y<im.height(); y++) {
        for(int x=0; x<im.width(); x++) {
            v.push_back(im[y][x]);
        }
    }
    std::sort(v.begin(),v.end());
    minVal=v[(int)(v.size()*minPerc)];
    maxVal=v[(int)(v.size()*maxPerc)];
}

/**
 * Add asymmetry to the coordinates list
 */
template <class SrcIterator, class SrcAccessor, class T>
inline void determineAsymmetry(triple<SrcIterator, SrcIterator, SrcAccessor> s,
        std::set<Coord<T> >& coords,
        const int factor) {
    determineAsymmetry(s.first, s.second, s.third, coords, factor);
}

template <class SrcIterator, class SrcAccessor, class T>
void determineAsymmetry(SrcIterator srcUpperLeft,
        SrcIterator srcLowerRight,
        SrcAccessor acc,
        std::set<Coord<T> >& coords,
        const int factor) {
    vigra::SplineImageView<3,T> sview(srcUpperLeft, srcLowerRight, acc, true);
    std::set<Coord<float> > newcoords;
    std::set<Coord<float> >::iterator it2;
    for(it2 = coords.begin(); it2 != coords.end(); it2++) {
        const Coord<float>& c = *it2;
        T sxx = sview.dxx((float)(c.x)/factor, (float)(c.y)/factor);
        T syy = sview.dyy((float)(c.x)/factor, (float)(c.y)/factor);
        T sxy = sview.dxy((float)(c.x)/factor, (float)(c.y)/factor);
        // calculate the eigenvalues
        T ev1 = (sxx+syy)/2. - sqrt((sxx+syy)*(sxx+syy)/4. + sxy*sxy - sxx*syy);
        T ev2 = (sxx+syy)/2. + sqrt((sxx+syy)*(sxx+syy)/4. + sxy*sxy - sxx*syy);
        Coord<float> cc (c.x, c.y, c.val, ev1/ev2);
        newcoords.insert(cc); // copy for now. Hack hack hack...
    }
    coords=newcoords;
}

template <class SrcIterator, class SrcAccessor, class T>
inline void determineSNR(triple<SrcIterator, SrcIterator, SrcAccessor> s,
        std::set<Coord<T> >& coords,
        const int factor ) {
    determineSNR(s.first, s.second, s.third, coords, factor );
}

template <class SrcIterator, class SrcAccessor, class T>
void determineSNR(SrcIterator srcUpperLeft,
        SrcIterator srcLowerRight,
        SrcAccessor acc,
        std::set<Coord<T> >& coords,
        const int factor) {
    vigra::SplineImageView<3,T> sview(srcUpperLeft, srcLowerRight, acc, true);
    std::set<Coord<float> > newcoords;
    std::set<Coord<float> >::iterator it2;
    for(it2 = coords.begin(); it2 != coords.end(); it2++) {
        const Coord<float>& c = *it2;
        T SNR = c.val; //SNR = amplitude_signal/standardDeviation_noise, stddev_noise = 1 => SNR = intensity maximum
        Coord<float> cc (c.x, c.y, c.val, c.asymmetry, SNR);
        newcoords.insert(cc); // copy for now. Hack hack hack...
    }
    coords=newcoords;
}

//template <class T>
struct transformationFunctor{
	transformationFunctor(float a, float intercept, float minInt = 0): a(a), intercept(intercept), C(-2/a*std::sqrt(a*minInt+ intercept)){}
	float operator()(float poissInt){return 2/a*std::sqrt(a*poissInt + intercept)+ C;}

private:
	float a;
	float intercept;
	float C;
};


template <class T>
void findCorrectionCoefficients(const MyImportInfo& info, std::vector<T>& parameterTrafo, bool second = false) {

    unsigned int stacksize = info.shape(2);
    unsigned int w = info.shape(0);
    unsigned int h = info.shape(1);

    int numberPoints = 2000;
    int listPixelCoordinates[numberPoints];
    T b;

    findGoodPixelsForScellamApproach(info, b, numberPoints, listPixelCoordinates, parameterTrafo);
    //std::cout<<"skellam candidates found"<<std::endl;
    //for(int n = 0; n < numberPoints;n++){
    //	std::cout<<listPixelCoordinates[n]<<'\n';
    //}
    T meanValues[numberPoints];
    T skellamParameters[numberPoints];
    //std::cout<<"before skellam found"<<std::endl;

    calculateSkellamParameters(info, listPixelCoordinates, meanValues, skellamParameters, numberPoints,parameterTrafo);

    std::ofstream selectedPoints;

	char filename[1000];
	sprintf(filename, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/selectedPoints.txt");
	selectedPoints.open(filename);

	std::cout<<"[";
	for(int n = 0; n < numberPoints;n++){
		selectedPoints<<meanValues[n]<<" "<<skellamParameters[n]<<std::endl;
	}
	selectedPoints.close();


    findBestFit(info, meanValues, skellamParameters, numberPoints, parameterTrafo);
    //std::cout<<"bestfit found"<<std::endl;
    MultiArray<3, T> tempMA(Shape3(w,h,1));
	T minVal = 999999999999;
	for(int k = 0; k<stacksize; k++){
		readBlock(info, Shape3(0,0,k),Shape3(w,h,1), tempMA);
		FindMinMax<T> minmax;
	    inspectMultiArray(srcMultiArrayRange(tempMA), minmax);
	    if(minVal > minmax.min){minVal = minmax.min;}
	}

	if(minVal < parameterTrafo[1]){parameterTrafo[1] = minVal;}
	//parameterTrafo[1] = minVal;   // sets offset to the minimum of all pixels, to avoid negative values during calculation of Poisson distributions
	std::cout<<"minimum: "<< minVal<<std::endl;
}

template <class T>
void findBestFit(const MyImportInfo& info,T meanValues[],T skellamParameters[],int numberPoints, std::vector<T>& parameterTrafo){
	int mode = 0; //0: take only lowest point from each interval and fit line through this points
	switch(mode){
		case 0:{
			//std::cout<<"begin bestfit"<<std::endl;
			int numberIntervals = std::min(10, numberPoints);
			//std::cout<<numberIntervals;
			T meanMin = 9999999, meanMax = 0;
			for(int i = 0; i< numberPoints; i++){
				if(meanValues[i]<meanMin){meanMin = meanValues[i];}
				if(meanValues[i]>meanMax){meanMax = meanValues[i];}
			}
			T intervalSize = (meanMax - meanMin)/ numberIntervals;
			std::vector<T> temp_vec(meanValues, meanValues + numberPoints);
			//std::cout<<"middle best fit"<<std::endl;
			std::vector<int> iter(numberPoints); //contains information about permutation from sort
			linearSequence(iter.begin(), iter.end());
			indexSort(temp_vec.begin(), temp_vec.end(), iter.begin());

			std::vector<T> minMeanValuesIntervals;
			std::vector<T> minSkellamValuesIntervals;
			T highValueForMinDetection = 99999999;
			T localMinMean;
			T currentMeanValue;
			//std::cout<<"middle2 best fit"<<std::endl;
			//from each interval only the point with the lowest skellam parameter is selected
			for(int j = 0; j < numberIntervals; j++){
				localMinMean = highValueForMinDetection;
				for(int i = 0; i < numberPoints; i++){
					if(meanValues[i] >= meanMin + j * intervalSize and meanValues[i] < meanMin + (j+1) * intervalSize){
						if(skellamParameters[i] < localMinMean){localMinMean = skellamParameters[i]; currentMeanValue = meanValues[i];}
					}
				}
				if(localMinMean != highValueForMinDetection){
					minMeanValuesIntervals.push_back(currentMeanValue);
					minSkellamValuesIntervals.push_back(localMinMean);
				}
			}
			std::cout<<"from find best fit"<< std::endl;
			std::cout<<std::endl<<"[";
			for(int i = 0; i< minMeanValuesIntervals.size(); i++){std::cout<<minMeanValuesIntervals[i]<<",";}
			std::cout<<"]"<<std::endl;
			std::cout<<"[";
			for(int i = 0; i< minSkellamValuesIntervals.size(); i++){std::cout<<minSkellamValuesIntervals[i]<<",";}
			std::cout<<"]"<<std::endl;
			//std::cin.get();
			T slope, gof, intercept;
			//std::cout<<"before ransac"<<std::endl;
			doRansac(minMeanValuesIntervals, minSkellamValuesIntervals, 400, slope, intercept, gof);
			//std::cout<<"after ransac"<<std::endl;
			std::cout<<"slope: "<< slope<<" x0: "<<-intercept/slope<<" intercept: "<<intercept<<" lowest gof:"<<gof<<std::endl;

			parameterTrafo[0] = slope; parameterTrafo[1] = -intercept/slope; parameterTrafo[2] = intercept;
		}
	}

}

//the best fit through the lowes points of each interval is found here
template <class T>
void doRansac(std::vector<T>& MeanValues, std::vector<T>& SkellamValues,
		int numberIterations, T& bestSlope, T& bestIntercept, T& lowestGof){
	//std::cout<<"begin ransac"<<std::endl;
	int numberPoints = MeanValues.size();
	int minNumberPoints = 5, maxNumberPoints = numberPoints;
	if(minNumberPoints > maxNumberPoints){std::cout<<"not enough points found by findGoodPixelsForScellamApproach!"<<std::endl;}
	T slope, intercept, gof;
	std::vector<int> indicesVector(numberPoints);
	linearSequence(indicesVector.begin(), indicesVector.end());
	//for(int i = 0; i< numberPoints;i++){std::cout<<indicesVector[i]<<std::endl;}
	int nbrPointsChosen;

	std::vector<T> collectionSlopes;
	std::vector<T> collectionIntercept;
	std::vector<T> collectionGof;
	//std::cout<<"middle ransac"<<std::endl;
	for(int niter = 0; niter<numberIterations; niter++){
		//std::cout<<niter<<" ";
		nbrPointsChosen = rand()%(maxNumberPoints - minNumberPoints + 1) + minNumberPoints;
		random_shuffle(indicesVector.begin(), indicesVector.begin() + nbrPointsChosen);

		fitLine(MeanValues, SkellamValues, indicesVector, nbrPointsChosen, slope, intercept, gof);

		collectionSlopes.push_back(slope);
		collectionIntercept.push_back(intercept);
		collectionGof.push_back(gof);
	}
	//std::cout<<"middle2 ransac"<<std::endl;
	T minGof = collectionGof[0];
	int indexMinGof = 0;
	for(int i = 0; i < numberIterations; i++){
		if(collectionGof[i] < minGof){minGof = collectionGof[i]; indexMinGof = i;}
	}
	//std::cout<<"end Ransac"<<std::endl;
	bestSlope = collectionSlopes[indexMinGof];
	bestIntercept = collectionIntercept[indexMinGof];
	lowestGof = collectionGof[indexMinGof];
}

//fits a line taking only the chosen points into account
template <class T>
void fitLine(std::vector<T>& MeanValues, std::vector<T>& SkellamValues,
		std::vector<int>& indicesVector,int& nbrPointsChosen,
		T &slope, T &intercept, T &gof){
	int numberPoints = MeanValues.size();
	typename std::vector<T>::iterator it;
	//std::cout<<"begin fit line"<<std::endl;
	std::vector<T> tempVx;
	std::vector<T> tempVy;
	//tempVx.reserve(100);
	//std::cout<<"nbrPointsChosen: "<<nbrPointsChosen<<std::endl;
	for(int i = 0; i<nbrPointsChosen; i++){
		tempVx.push_back(MeanValues[indicesVector[i]]);
		tempVy.push_back(SkellamValues[indicesVector[i]]);
	}
	//std::cout<<"middle fit line"<<std::endl;
	T mXY = 0, mX = 0,mY = 0, mXX = 0, mXmX = 0;// mXY: mean of all products x*y, mX: mean x, mY: meanY, mXX: mean of pow(x,2), mXmX: pow(mX,2)
	for(int i = 0; i < nbrPointsChosen; i++){
		mXY += (tempVx[i] * tempVy[i]);
		mX += tempVx[i];
		mY += tempVy[i];
		mXX += (tempVx[i]*tempVx[i]);
	}
	mXY /= nbrPointsChosen;
	mX /= nbrPointsChosen;
	mY /= nbrPointsChosen;
	mXX /= nbrPointsChosen;
	slope = (mXY - mX*mY)/(mXX - mX*mX); // calculates slope based on total least squares
	//std::cout<<mXY<<" "<<mX<<" "<<mY<<" "<< mXX  << std::endl;
	intercept = mY - slope * mX;
	gof = 0;
	for(int i = 0; i< nbrPointsChosen;i++){
		gof += pow((tempVy[i] - (intercept + slope * tempVx[i])),2);
	}
	gof /= numberPoints;
	//std::cout<<"slope: "<< slope<<" intercept:"<<intercept<<" goodness: "<<gof<<std::endl;
}

template <class T>
T median(MultiArray<3, T> array){
	int size = 0;
	typename MultiArray<3,T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){size += 1;}
	std::vector<int> iter(size); //contains information permutation from sort
	linearSequence(iter.begin(), iter.end());
	indexSort(array.begin(), array.end(), iter.begin());
	if(size == 1){return array[0];}
	if(size%2 == 0){return (array[iter[size/2]] + array[iter[size/2-1]])/2;}
	if(size%2 == 1){return array[iter[size/2]];}
}
template <class T>
T median(MultiArray<2, T> array){
	int size = 0;
	typename MultiArray<2,T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){size += 1;}
	std::vector<int> iter(size); //contains information permutation from sort
	linearSequence(iter.begin(), iter.end());
	indexSort(array.begin(), array.end(), iter.begin());
	if(size == 1){return array[0];}
	if(size%2 == 0){return (array[iter[size/2]] + array[iter[size/2-1]])/2;}
	if(size%2 == 1){return array[iter[size/2]];}
}
template <class T>
T meanArr(MultiArray<3, T> array){
	int counter = 0;
	T meanVal = 0;
	typename MultiArray<3,T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){
		meanVal += array[counter];
		counter += 1;
	}
	return meanVal/counter;
}

template <class T>
T meanArr(MultiArray<2, T> array){
	int counter = 0;
	T meanVal = 0;
	typename MultiArray<2,T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){
		meanVal += array[counter];
		counter += 1;
	}
	return meanVal/counter;
}

template <class T>
T meanArr(BasicImage<T> array){
	int counter = 0;
	T meanVal = 0;
	typename BasicImage<T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){
		meanVal += *it;//array[counter];
		counter += 1;
	}
	return meanVal/counter;
}

template <class T>
T meanArr(BasicImageView<T> array){
	int counter = 0;
	T meanVal = 0;
	typename BasicImage<T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){
		meanVal += *it;//array[counter];
		counter += 1;
	}
	return meanVal/counter;
}


template <class T>
T minArr(MultiArray<3, T> array){
	int counter = 0;
	T minVal = 999999999;
	typename MultiArray<3,T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){
		if(array[counter]<minVal){minVal = array[counter];};
		counter += 1;
	}
	return minVal;
}

template <class T>
T max2DArr(MultiArray<2, T> array){
	int counter = 0;
	T maxVal = -999999999;
	typename MultiArray<2,T>::iterator it;
	for(it = array.begin(); it != array.end(); it++){
		if(array[counter]>maxVal){maxVal = array[counter];};
		counter += 1;
	}
	return maxVal;
}

void printIntensities(const MyImportInfo& info, int* vecw, int* vech, int nbrPoints, double a, double b){
	unsigned int stacksize2 = info.shape(2);
	unsigned int w = info.shape(0);
	unsigned int h = info.shape(1);
	//std::cout<<"w: "<<w<<" h: "<<h<<std::endl;
	//std::cin.get();
	MultiArray<3, float> mean_im_temp(Shape3(w,h,1));
	MultiArray<3, float> mean_im_temp2(Shape3(w,h,1));


	std::ofstream origimg[nbrPoints];

	char temp[1000];
	for(int i = 0;i < nbrPoints; i++){
		sprintf(temp, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/outputPoisTest/pos0_%d_pos1%d.txt", vecw[i], vech[i]);
		origimg[i].open(temp);
		std::cout<<"vecw["<<i<<"]="<<vecw[i]<<" vech["<<i<<"]="<<vech[i]<<std::endl;
		std::cout<<"a:"<< a<<" b" <<b<<std::endl;
	}



	for(int f = 0; f< stacksize2;f++){
		readBlock(info, Shape3(0,0,f),Shape3(w,h,1), mean_im_temp2);
		for(int i = 0; i <nbrPoints; i++){
			origimg[i]<< (mean_im_temp2(vecw[i], vech[i])-b)/a<<std::endl;
		}
	}

	for(int i=0; i<nbrPoints; i++){
		origimg[i].close();
	}
	std::cout<<"coordinates saved"<<std::endl;
}

template <class T>
void getMask2(BasicImage<T>& data, unsigned int w,unsigned int h, MultiArray<2, T>& mask,int framenumber){
	float alpha = 0.05;

	std::ofstream origimg, filter;

	char filtername[1000], origname[1000];
	sprintf(filtername, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/filterGauss%d.txt", framenumber);
	sprintf(origname, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/origBild%d.txt", framenumber);
	filter.open (filtername);
	origimg.open(origname);

	T meanData = meanArr(data);

	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			mask(i,j) = isSignalGauss(data(i,j), alpha, meanData);
			//filter << i<<"  "<< j<<"  "<< mask(i,j)<<std::endl;
			//origimg << i<<"  "<< j<<" "<<data(i,j)<<std::endl;
		}
	}
	origimg.close();
	filter.close();
	//std::cin.get();
}

template <class T>
void getMask2(BasicImageView<T>& data, unsigned int w,unsigned int h, MultiArray<2, T>& mask,int framenumber){
	float alpha = 0.05;

	std::ofstream origimg, filter;

	char filtername[1000], origname[1000];
	sprintf(filtername, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/filterGauss2%d.txt", framenumber);
	sprintf(origname, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/origBild2%d.txt", framenumber);
	filter.open (filtername);
	origimg.open(origname);

	T meanData = meanArr(data);

	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			mask(i,j) = isSignalGauss(data(i,j), alpha, meanData);
			//filter << i<<"  "<< j<<"  "<< mask(i,j)<<std::endl;
			//origimg << i<<"  "<< j<<" "<<data(i,j)<<std::endl;
		}
	}
	origimg.close();
	filter.close();
	//std::cin.get();
}


template <class T>
T isSignalGauss(T corrInt, T alpha, T mu){
	T sig = 1;
	T rest = corrInt - (int)corrInt;
	T factor = 10.0;
	int numberValues;
	std::vector<T> x;
	if(corrInt < mu){
		numberValues = (int)(3*sig *factor);
		for(int i=0;i<numberValues;i++){
			x.push_back(corrInt - 3*sig + i/factor);
		}
	}
	else{
		numberValues = (int)3*sig*factor + (corrInt - mu)*factor;
		for(int i=0;i<numberValues;i++){
			x.push_back( mu - 3*sig + i/factor);
		}
	}
	//std::cout<<"corrInt: "<<corrInt<<" numberValues:"<<numberValues<<std::endl;
	std::vector<T> y(numberValues);
	gaussDistribution(x,y,sig, mu);
	/*std::cout<<"mu: "<<mu<<std::endl;
	std::cout<<"[";
	for(int i =0; i< numberValues;i++){
		std::cout<<x[i]<<", ";
	}
	std::cout<<"]"<<std::endl;

	std::cout<<"[";
	for(int i =0; i< numberValues;i++){
		std::cout<<y[i]<<", ";
	}
	std::cout<<"]"<<std::endl;
	int kk;
	std::cin>>kk;*/
	T cdf = 0;
	for(int i = 0; i<y.size();i++){cdf += y[i];}
	//cdf = cdf + rest*y[y.size()-1]; //if corrInt is something like 13.9 this adds 0.9 times y[last]
	cdf /= factor;
	//std::cout<<"cdf: "<< cdf<<std::endl<<std::endl;

	if(cdf < 1-alpha){return 0;}
	if(cdf > 1-alpha and cdf <= 1){return (cdf - (1-alpha))/alpha;}
	else{return 1;} //case if very bright spot

}

template <class T>
void gaussDistribution(std::vector<T>& xin, std::vector<T>& yout, T sig,T mu){
	if(xin.size()!=yout.size()){std::cout<<"Warning x and y vector must have the same length!"<<std::endl;}
	int numberElements = xin.size();
	int minVal = (int) xin[0];
	for(int i = 0 ; i<numberElements;i++){
		yout[i] = 1/(std::sqrt(2*3.14159265)*sig)*std::exp(-0.5*std::pow((xin[i]-mu),2)/(std::pow(sig,2)));
		//std::cout<<"std::sqrt(2*pi)*sig: "<<(std::sqrt(2*3.14159265)*sig)<<" exp(-0.5*(x-mu)**2/sig**2): "<<std::exp(-0.5*std::pow((xin[i]-mu),2)/(std::pow(sig,2)))<<std::endl;
		//std::cout<<"-0.5*std::pow((xin[i]-mu),2): "<<-0.5*std::pow((xin[i]-mu),2)<<" (std::pow(sig,2)): "<<(std::pow(sig,2))<<" -0.5*std::pow((xin[i]-mu),2)/(std::pow(sig,2)): "<<-0.5*std::pow((xin[i]-mu),2)/(std::pow(sig,2))<<std::endl;
		//std::cout<<"xin[i]"<< xin[i]<<" mu: "<<mu<<std::endl;
	}

}

template<class T>
void showPoisson(const MyImportInfo& info, std::vector<T>& parameterTrafo){
	unsigned int stacksize = info.shape(2);
	unsigned int w = info.shape(0);
	unsigned int h = info.shape(1);

	std::ofstream f1,f2,f3,f4,f5;
	std::ofstream bf1,bf2,bf3,bf4,bf5;
	char nf1[1000],nf2[1000],nf3[1000],nf4[1000],nf5[1000];
	char bnf1[1000],bnf2[1000],bnf3[1000],bnf4[1000],bnf5[1000];
	sprintf(nf1, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/f1.txt");
	sprintf(nf2, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/f2.txt");
	sprintf(nf3, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/f3.txt");
	sprintf(nf4, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/f4.txt");
	sprintf(nf5, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/f5.txt");

	sprintf(bnf1, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/beforef1.txt");
	sprintf(bnf2, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/beforef2.txt");
	sprintf(bnf3, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/beforef3.txt");
	sprintf(bnf4, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/beforef4.txt");
	sprintf(bnf5, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/beforef5.txt");

	f1.open(nf1);
	f2.open(nf2);
	f3.open(nf3);
	f4.open(nf4);
	f5.open(nf5);

	bf1.open(bnf1);
	bf2.open(bnf2);
	bf3.open(bnf3);
	bf4.open(bnf4);
	bf5.open(bnf5);

	MultiArray<3, T> tempMA(Shape3(w,h,1));
	for(int k = 0; k<stacksize; k++){
		readBlock(info, Shape3(0,0,k),Shape3(w,h,1), tempMA);
		if(true){
			for(int i=0; i< w; i++){
				for(int j = 0; j< h;j++){
					if(i == 1 & j == 1){
						bf1 << tempMA(i,j,0);
					}
					if(i == 10 & j == 10){
						bf2 << tempMA(i,j,0)<<", ";
					}
					if(i == 100 & j == 100){
						bf3 << tempMA(i,j,0)<<", ";
					}
					if(i == 20 & j == 20){
						bf4 << tempMA(i,j,0)<<", ";
					}
					if(i == 5 & j == 5){
						bf5 << tempMA(i,j,0)<<", ";
					}

					tempMA(i,j,0) = (tempMA(i,j,0) - parameterTrafo[1])/parameterTrafo[0];
					if(i == 1 & j == 1){
						f1 << tempMA(i,j,0)<<", ";
					}
					if(i == 10 & j == 10){
						f2 << tempMA(i,j,0)<<", ";
					}
					if(i == 100 & j == 100){
						f3 << tempMA(i,j,0)<<", ";
					}
					if(i == 20 & j == 20){
						f4 << tempMA(i,j,0)<<", ";
					}
					if(i == 5 & j == 5){
						f5 << tempMA(i,j,0)<<", ";
					}

				}
			}
		}
	}
	f1.close();
	f2.close();
	f3.close();
	f4.close();
	f5.close();

	bf1.close();
	bf2.close();
	bf3.close();
	bf4.close();
	bf5.close();

}

//get Mask calculates the probabilities for each pixel of the current frame to be foreground,
//there are two different methods available: logliklihood, based on cumulative distribution function
template <class T>
void getMask(MultiArrayView<2,T >& array, const MyImportInfo& info, MultiArray<3, T>& PoissonMeans ,int framenumber,MultiArray<2,T>& mask, std::vector<T> parameterTrafo){
	unsigned int stacksize =  info.shape(2);
	unsigned int w = info.shape(0);
	unsigned int h = info.shape(1);
	//MultiArray<2,T> array(Shape2(w,h));
	//MultiArray<3, T> tempMA(Shape3(w,h,1));
	//readBlock(info, Shape3(0,0,framenr),Shape3(w,h,1), tempMA);
	T at = 4.38, bt = 380, it = 460.6, mt =20;
	/*for(int i=0; i<999999999; i++){
		std::cout<<factorial(100)<<std::endl;
		std::cin.get();
		isSignal(it,at,bt,mt);
	}*/

	//T muMedian =  median(PoissonMeans);

	//std::cout<<"median:"<<muMedian<<std::endl;
	//T muMean = meanArr(PoissonMeans);
	T muMin = minArr(PoissonMeans);
	//std::cout<<"mean:"<<muMean<<std::endl;

	vigra::FindMinMax<T> arrMinmax;
	inspectImage(srcImageRange(array), arrMinmax);

	T alpha = 0.01;
	//std::cout<<"max: "<< arrMinmax.max<<" poisMeanMin: "<<muMin;
	//T factor2ndDist = 1.3; // should be dependent of snr
	T factor2ndDist = 1+arrMinmax.max/muMin/10; // should be dependent of snr
	//std::cout<<" factor2ndDist: "<<factor2ndDist<<std::endl;
	std::ofstream origimg, maski;


	bool writeMatrices = false;
	if(writeMatrices){
	char origimgname[1000];
	char maskname[1000];
	sprintf(origimgname, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/origimgbeforefilter%d.txt", framenumber);
	sprintf(maskname, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/maskPoiss%d.txt", framenumber);
	origimg.open (origimgname);
	maski.open (maskname);
	}
	double cdf = 0;
	for(int i = 0; i< w; i++){

		for(int j = 0; j < h; j++){

			//mask(i,j) = isSignal(array(i,j), muMin, factor2ndDist);
			//mask(i,j) = isSignal_fast(array(i,j), muMin, factor2ndDist);
			//std::cin.get();
			//mask(i,j) = isSignalCDF(array(i,j), muMin, alpha);
			//mask(i,j) = isSignalCDF(array(i,j), PoissonMeans(i,j), alpha);
			cdf = ppois(array(i,j), PoissonMeans(i,j), 1,0);

			if (cdf > 1- alpha and cdf < 1){
				mask(i,j) = (cdf - (1-alpha))/alpha;
			}
			else if (cdf == 1){mask(i,j) = 1;}
			else if (cdf < 1- alpha){
				mask(i,j) = 0;
			}

			//printf("%f     ",ppois(132, 140.1, 1,0));
			//printf("%f     ",ppois(array(i,j), PoissonMeans(i,j), 1,0));
			//std::cout<<i<<"  "<<j<<"  "<<"ppois("<<array(i,j)<<", "<<PoissonMeans(i,j)<<", 1,0): "<<ppois(array(i,j), PoissonMeans(i,j), 1,0)<<" mask(i,j): "<<mask(i,j)<<std::endl;
			if (writeMatrices){
			maski << i<<"  "<< j<<"  "<< mask(i,j)<<std::endl;
			origimg << i<<"  "<< j<<"  "<< array(i,j)<<std::endl;}
			//std::cout<<i<<"  "<<j<<"  "<<PoissonMeans(i,j)<<" mask(i,j): "<<mask(i,j)<<std::endl;
		}
	}
	if (writeMatrices){
	origimg.close();
	maski.close();}
	T beta = 0.5;

	//performGraphcut(mask,beta);
	gaussianSmoothing(srcImageRange(mask), destImage(mask), parameterTrafo[3]/2.);
	transformImage(srcImageRange(mask),destImage(mask),
	        		ifThenElse(Arg1()>Param(0.2), Param(1.), Param(0.)));


	std::ofstream maskafter;
	bool writeMatrices2 = false;
	if(writeMatrices2){
		char maskaftername[1000];
		sprintf(maskaftername, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/maskafterGraphcutPoiss%d.txt", framenumber);
		maskafter.open (maskaftername);

	for(int i = 0; i< w; i++){
		for(int j = 0; j < h; j++){
			if (writeMatrices2){
			maskafter << i<<"  "<< j<<"  "<< mask(i,j)<<std::endl;}
		}
	}
	//std::cin.get();
	//std::cout<<"stop before"<<std::endl;
	maskafter.close();
	}
}

//assuming a poission distribution, each intensity can be compared with the probability
//of the cdf, if the value is above 1-alpha a non zero value is returned
template <class T>
T isSignalCDF(T corrInt, T lamb, T alpha){
	//T corrInt = (intensity - b)/ a;
	if(corrInt < 1){return 0;}
	T rest = corrInt - (int)corrInt;
	std::vector<T> x((int)corrInt);
	std::vector<T> y((int)corrInt);
	linearSequence(x.begin(), x.end());

	poissonDistribution(x,y,lamb);
	std::cout<<"corrInt: "<< corrInt<<" lamb: "<< lamb<<std::endl;
	for(int i = 0; i<corrInt; i++){
		std::cout<<y[i]<< " , ";
	}
	std::cout<<std::endl;
	T cdf = 0;
	for(int i = 0; i<y.size();i++){cdf += y[i];}
	cdf = cdf + rest*y[y.size()-1]; //if corrInt is something like 13.9 this adds 0.9 times y[last]
	//std::cout<<"cdf: "<< cdf<<std::endl;
	if(cdf < 1-alpha){return 0;}
	if(cdf > 1-alpha and cdf <= 1){
		//std::cout<<"alpha: "<<alpha<<" cdf: "<<cdf<<" (cdf - (1-alpha))/alpha: "<<(cdf - (1-alpha))/alpha<<std::endl;
		return (cdf - (1-alpha))/alpha;}
	else{std::cout<<"corrInt :"<<corrInt<<" alpha: "<<alpha<<" cdf: "<<cdf<<" (cdf - (1-alpha))/alpha: "<<(cdf - (1-alpha))/alpha<<std::endl;
	return 1;} //case if very bright spot

}

//implementation of poissonDistribution
template <class T>
void poissonDistribution(std::vector<T>& xin, std::vector<T>& yout, T lamb){
	if(xin.size()!=yout.size()){std::cout<<"Warning x and y vector must have the same length!"<<std::endl;}
	int numberElements = xin.size();
	int minVal = (int) xin[0];
	for(int i = 0 ; i<numberElements;i++){
		double p1 = log(pow(lamb, xin[i]));
		//std::cout<<"p1: "<<p1<<" ";
		double p2 = log(factorial(xin[i]));
		//std::cout<<"p2: "<<p2<<" ";
		double p3 = exp(-lamb);
		//std::cout<<"p3: "<<p3<<" ";
		double p4 = p1 - p2;
		//std::cout<<"p4: "<<p4<<" ";
		double p5 = exp(p4);
		//std::cout<<"p5: "<<p5<<" ";
		yout[i] = (p5*p3);
		//std::cout<<"p5 * p3: "<<p5 * p3<<" ";
	}

}

//two poisson distributions are compared, the two distributions are, the distribution around
//the estimated mean for the background pixels and a second distribution with a mean scaled with the
//factor factor2ndDist, each intensity  is explained better by either the background or the
//"foreground" distribution.
template <class T>
T isSignal(T corrInt, T lamb_,T factor2ndDist){
	if(int(corrInt)== 0){return 0;}
	T lamb = lamb_;
	T prob = 0;
	//T factor2ndDist = 1.5; //determines the shift of the 2nd distribution
	T limit = .5;		//determines the value in case signal/no signal for sure
	//T corrInt = (intensity - b)/ a;
	std::vector<T> x((int)corrInt);
	linearSequence(x.begin(), x.end());
	std::vector<T> y1(x.size());
	std::vector<T> y2(x.size());
	//poissonDistribution(x,y1,lamb);
	//poissonDistribution(x,y2,lamb * factor2ndDist);
	poissonDistribution_log(x,y1,lamb);
	poissonDistribution_log(x,y2,lamb * factor2ndDist);

	int lasty = x.size()-1;

	//if(corrInt < 1){prob =  -limit;} // is done by if(int(corrInt) == 0
	if(y1[lasty] == 0 && y2[lasty] > 0){prob = limit;}
	else if(y2[lasty] == 0 && y1[lasty] > 0){prob =  -limit;}
	else if(y2[lasty] == 0 && y1[lasty] == 0){prob =  limit;}
	//else{prob = -log(y1[lasty]/y2[lasty]);}
	else{prob = y2[lasty] - y1[lasty];}

	if(prob > limit or std::isnan(prob)){prob =  limit;}
	if(prob < -limit or std::isnan(-prob)){prob = -limit;}

	prob = (prob + limit)/(2*limit);
	//if(prob < 0.7){prob = 0;}
	return prob;
}

template <class T>
T isSignal_fast(T corrInt, T lamb_,T factor2ndDist){
	if(int(corrInt)== 0){return 0;}
	T lamb = lamb_;
	T prob = 0;
	T factor = 5;
	T limit = .5*factor;		//determines the value in case signal/no signal for sure
	T y1,y2;

	poissonlastValue_log(corrInt, y1, lamb);
	poissonlastValue_log(corrInt, y2, lamb * factor2ndDist);

	//if(corrInt < 1){prob =  -limit;} // is done by if(int(corrInt) == 0
	if(y1 == 0 && y2 > 0){prob = limit;}
	else if(y2 == 0 && y1 > 0){prob =  -limit;}
	else if(y2 == 0 && y1 == 0){prob =  limit;}
	//else{prob = -log(y1[lasty]/y2[lasty]);}
	else{prob = y2 - y1;}

	if(prob > limit or std::isnan(prob)){prob =  limit;}
	if(prob < -limit or std::isnan(-prob)){prob = -limit;}

	prob = (prob + limit)/(2*limit);
	//if(prob < 0.7){prob = 0;}
	return prob;
}

//estimates the mean of the poisson distributions
template <class T>
void getPoissonDistributions(const MyImportInfo& info,T a,T offset,MultiArray<3,T> &PoissonMeans){
	unsigned int stacksize = info.shape(2); // 100
	unsigned int w = info.shape(0);
	unsigned int h = info.shape(1);
	int g;
	MultiArray<3,T> tempMA(Shape3(w,h,1));
	for(int k = 0; k< stacksize; k++){
		//std::cout<<"framenr: "<<k<<std::endl;
		readBlock(info, Shape3(0,0,k),Shape3(w,h,1), tempMA);
		for(int i = 0; i<w;i++){
			for(int j=0; j<h; j++){
				tempMA(i,j,0) = (tempMA(i,j,0) - offset)/a;
			}
		}
		combineTwoMultiArrays(srcMultiArrayRange(PoissonMeans), srcMultiArray(tempMA), destMultiArray(PoissonMeans),
                std::plus<T>());
	}
	for(int i = 0; i<w;i++){
		for(int j=0; j<h; j++){
			PoissonMeans(i,j,0) = PoissonMeans(i,j,0)/ stacksize;
			//std::cout<<PoissonMeans(i,j,0)<<" ";
		}
	}
	//std::cin>>a;
}

//calcuates the skellam parameters (in principle variances) for the pixels chosen previously.
template <class T>
void calculateSkellamParameters(const MyImportInfo& info, int listPixelCoordinates[],T meanValues[],T skellamParameters[],
		int numberPoints, std::vector<T> parameterTrafo){

	unsigned int stacksize = info.shape(2);
	unsigned int w = info.shape(0);
	unsigned int h = info.shape(1);
	std::vector< std::vector <T> > intensities(numberPoints, std::vector<T>(stacksize));
	MultiArray<3, T> temp_arr(Shape3(w,h,1));
	for(int f = 0; f < stacksize; f++){
		readBlock(info, Shape3(0,0,f),Shape3(w,h,1), temp_arr);
		for(int i = 0; i < numberPoints; i++){
			intensities[i][f] = temp_arr[listPixelCoordinates[i]];
		}
	}
	for(int i = 0; i < numberPoints; i++){
		meanValues[i] = 0;
		for(int f = 0; f < stacksize; f++){meanValues[i] += intensities[i][f];}
		meanValues[i] /= stacksize;
	}
	T sigma[numberPoints];
	for(int i = 0; i < numberPoints; i++){
		skellamParameters[i] = (intensities[i][stacksize-1] - intensities[i][0])/stacksize;
		sigma[i] = 0;
		//std::cout<<"before "<<skellamParameters[i]<<" "<< intensities[i][0]<<"  "<< intensities[i][1]<<'\n';
		for(int f = 0; f < stacksize-1; f++){
			sigma[i] += std::pow((skellamParameters[i] - intensities[i][f] + intensities[i][f+1]),2);
		}
		skellamParameters[i]+= sigma[i];
		skellamParameters[i] /= 2*(stacksize - 1);
	}

}

//To estimate the gain factor points with different mean intensities are needed. This functions searches for
//good candidates, it tries to find pixels with as much different mean values as possible.
template <class T>
void findGoodPixelsForScellamApproach(const MyImportInfo& info, T & tmpType,int numberPoints, int listPixelCoordinates[],
		std::vector<T> parameterVector){
	unsigned int stacksize = info.shape(2);
	unsigned int w = info.shape(0);
	unsigned int h = info.shape(1);

	int numberFramesForMean = stacksize;

	MultiArray<3, T> mean_im_temp(Shape3(w,h,1));


	MultiArray<3, T> mean_im_temp2(Shape3(w,h,1));

	std::vector<T> pixelValues;
	int pixelPerFrame = w*h;
	int minMean = 999999999, maxMean = 0;
	std::vector<int> iter2;

    for(int f = 0; f< numberFramesForMean;f++){
    	readBlock(info, Shape3(0,0,f),Shape3(w,h,1), mean_im_temp2);
    	vigra::combineTwoMultiArrays(srcMultiArrayRange(mean_im_temp2),
    						srcMultiArray(mean_im_temp),
    						destMultiArray(mean_im_temp),
    		                std::plus<int>());
    }
    for(int i = 0; i< pixelPerFrame;i++){mean_im_temp[i] = mean_im_temp[i]/numberFramesForMean;}
    FindMinMax<T> minmax;

    inspectMultiArray(srcMultiArrayRange(mean_im_temp), minmax);
    std::cout<<"min: "<<minmax.min<<" max: "<<minmax.max<<std::endl;

    std::vector<int> iter(pixelPerFrame); //contains information permutation from sort
	linearSequence(iter.begin(), iter.end());
	indexSort(mean_im_temp.begin(), mean_im_temp.end(), iter.begin());

	std::cout << pixelPerFrame<<"  "<<mean_im_temp[iter[pixelPerFrame-20]]<<"  "<<mean_im_temp[iter[pixelPerFrame-30]]<<'\n';

	int mode = 1; 	//0: chose points with lower and higher mean intensities
					//1: take points from interval

	switch(mode){
		case 0:{
			int lowerPerc, higherPerc;
			float percentage_low = 10/100.;
			float percentage_high = 1/100.;

			std::cout<<std::rand()%10<<'\n';

			//findIndicesBackgroundAndBeads(lowerPerc, higherPerc, percentage, mean_im_temp, iter, minmax.min, minmax.max);
			lowerPerc = (int)(percentage_low * pixelPerFrame);
			higherPerc = (int)(pixelPerFrame * (1-percentage_high));
			std::cout<<"minVal: "<<minmax.min<<" lowerPerc: "<<lowerPerc<<" lowerVal: "<<mean_im_temp[iter[lowerPerc]]<<'\n';
			std::cout<<"maxVal: "<<minmax.max<<" higherPerc: "<<higherPerc<<" higherVal: "<<mean_im_temp[iter[higherPerc]]<<'\n';

			for(int n = 0;n < numberPoints;n++){
				if(std::rand()%2 == 0){listPixelCoordinates[n] = iter[std::rand()%lowerPerc];}
				else{listPixelCoordinates[n] = iter[std::rand()%(pixelPerFrame - higherPerc) + higherPerc];}
				//std::cout<<listPixelCoordinates[n]<<'\n';
			}
		}
		case 1:{
			int intervalCounter = 0;
			T intervalDistance = (minmax.max - minmax.min)/ numberPoints;
			std::cout<<"intervalDist: "<<intervalDistance<<'\n';
			for(int i = 0; i< pixelPerFrame; i++){
				if(mean_im_temp[iter[i]] > minmax.min + intervalDistance * intervalCounter){
					listPixelCoordinates[intervalCounter] = iter[i];
					//std::cout<<"listat counter: "<<listPixelCoordinates[intervalCounter]<<'\n';
					//std::cout<<"intervalCounter: "<<intervalCounter<<'\n';
					intervalCounter += 1;
					//std::cout<<"intervalCounter: "<<intervalCounter<<'\n';
					if (intervalCounter == numberPoints){break;}
				}
			}
			std::cout<<"intervalCounter ende: "<<intervalCounter<<'\n';
			if (intervalCounter < numberPoints){	//this is the case if the last for loop finished before enough points were found
				intervalCounter -= 1;
				for(int ic = intervalCounter; ic < numberPoints; ic++){listPixelCoordinates[ic] = iter[std::rand()%pixelPerFrame];}
			}
		}
	}
	//for(int n = 0; n < numberPoints;n++){
	//	std::cout<<listPixelCoordinates[n]<<"  intensities: "<<mean_im_temp[listPixelCoordinates[n]]<<'\n';
	//}
}

template <class T>
void applyTransformation(MultiArray<2, T>& array,int w, int h,T a,T offset){
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			array(i,j) = (array(i,j) - offset)/a;
		}
	}
}

template <class T>
void applyMask(BasicImage<T>& img, MultiArrayView<2, T>& mask, int frnr){
	int w = mask.shape(0);
	int h = mask.shape(1);

	std::ofstream imgAfterMaskFile, imgBeforeMaskFile;
	bool writeMatrices = false;
	if(writeMatrices){

		char imgAfterMask[1000];
		char imgBeforeMask[1000];
		sprintf(imgAfterMask, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/filterdImgAfterMask%d.txt", frnr);
		sprintf(imgBeforeMask, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/filterdImgBeforeMask%d.txt", frnr);
		imgAfterMaskFile.open(imgAfterMask);
		imgBeforeMaskFile.open(imgBeforeMask);}

	for(int i = 0; i< w; i++){
		for(int j = 0; j<h; j++){
			imgBeforeMaskFile << i<<"  "<< j<<"  "<< img(i,j)<<std::endl;
			if(mask(i,j)== 0.0){img(i,j)= 0;}
			imgAfterMaskFile << i<<"  "<< j<<"  "<< img(i,j)<<std::endl;

		}
	}
	if (writeMatrices){
		imgAfterMaskFile.close();
		imgBeforeMaskFile.close();}
	//std::cin.get();

}

//--------------------------------------------------------------------------
// GENERATE WIENER FILTER
//--------------------------------------------------------------------------
// Most of the following functions are available twice: Once taking a
// MultiArrayView (data) as input and once with MyImportInfo (filepointer).
// Since the algorithms work on single-frames only, there is no need to
// put the complete dataset into RAM but every frame can be read from disk
// when needed. The class MyImportInfo transparently handles hdf5 and sif
// input file pointers.

/**
 * Calculate Power-Spektrum
 */
template <class T, class DestIterator, class DestAccessor>
void powerSpectrum(const MultiArrayView<3, T>& array, 
                   DestIterator res_ul, DestAccessor res_acc,
                   std::vector<T> parameterTrafo) {
    unsigned int stacksize = array.size(2);
    unsigned int w = array.size(0);
    unsigned int h = array.size(1); 
    vigra::DImage ps(w, h);
    vigra::DImage ps_center(w, h);
    ps = 0;
    typename MultiArray<2,T>::iterator it;
    T a = parameterTrafo[0];
    T b = parameterTrafo[1];

    transformationFunctor tF(parameterTrafo[3], parameterTrafo[5], parameterTrafo[4]);
    for(unsigned int i = 0; i < stacksize; i++) {
        MultiArrayView <2, T> array2 = array.bindOuter(i); // select current image
        for(it = array2.begin(); it != array2.end(); it++){
        	if((*it-b)/a > 0){
        		*it = (*it-b)/a;
        		*it = (2*std::sqrt(*it + 3.0/8.0));
        	}
        	else{*it = 2*std::sqrt(0+3./8);}
        }

        //        for(int x0 = 0; x0 < w; x0++){
//			for(int x1 = 0; x1 < h; x1++){
//				//std::cout<<array2(x0,x1,0)<< "<- before ";
//				if(((array2(x0,x1) - parameterTrafo[1])/parameterTrafo[0])>0){
//				array2(x0,x1) = (array2(x0,x1) - parameterTrafo[1])/parameterTrafo[0];}
//				else{array2(x0,x1) = 0;}
//				//std::cout<<"a: "<< parameterTrafo[0]<< "offset: "<< parameterTrafo[1]<< "intersect: "<< parameterTrafo[2]<< "a2: "<< parameterTrafo[3]<< "offset2: "<< parameterTrafo[4]<< "intercept2: "<< parameterTrafo[5]<<" poissInt: "<<poissInt<<std::endl;
//				array2(x0,x1) = (2*std::sqrt(array2(x0,x1) + 3.0/8.0));
//				//std::cout<<array2(x0,x1,0)<< "<- after ";
//			}
//		}
        BasicImageView<T> input = makeBasicImageView(array2);  // access data as BasicImage     

        vigra::FFTWComplexImage fourier(w, h);
        fourierTransform(srcImageRange(input), destImage(fourier));
        
        // there is no squared magnitude accessor, so we use the magnitude here
        vigra::combineTwoImages(srcImageRange(ps), 
                srcImage(fourier, FFTWSquaredMagnitudeAccessor<double>()), 
                destImage(ps), Arg1()+Arg2());
        
        helper::progress(i, stacksize); // report progress
    }

    std::cout <<"info"<< std::endl;

    moveDCToCenter(srcImageRange(ps), destImage(ps_center));
    vigra::transformImage(
            srcImageRange(ps_center), 
            destIter(res_ul, res_acc), 
            Arg1() / Param(stacksize));
}

/**
 * Calculate Power-Spektrum
 */
template <class DestIterator, class DestAccessor, class T>
void powerSpectrum(const MyImportInfo& info, 
                   DestIterator res_ul, DestAccessor res_acc, std::vector<T> parameterTrafo) {
    unsigned int stacksize = info.shapeOfDimension(2);
    unsigned int w = info.shapeOfDimension(0);
    unsigned int h = info.shapeOfDimension(1); 
   // typedef float T;
    MultiArray<3, T> im(Shape3(w,h,1));
    vigra::DImage ps(w, h);
    vigra::DImage ps_center(w, h);
    ps = 0;
    //std::default_random_engine generator;
    //std::normal_distribution<double> distribution(0.0,1.0);

    transformationFunctor tF(1, 3./8., 0);
    for(unsigned int i = 0; i < stacksize; i++) {
        readBlock(info, Shape3(0,0,i), Shape3(w,h,1), im);
        for(int x0 = 0; x0 < w; x0++){
        	for(int x1 = 0; x1 < h; x1++){
        		//std::cout<<im(x0,x1,0)<< "<- before ";
        		if(((im(x0,x1,0) - parameterTrafo[1])/parameterTrafo[0])>0){
        		im(x0,x1,0) = (im(x0,x1,0) - parameterTrafo[1])/parameterTrafo[0];}
        		else{im(x0,x1,0) = 0;}
        		//std::cout<<"a: "<< parameterTrafo[0]<< "offset: "<< parameterTrafo[1]<< "intersect: "<< parameterTrafo[2]<< "a2: "<< parameterTrafo[3]<< "offset2: "<< parameterTrafo[4]<< "intercept2: "<< parameterTrafo[5]<<" poissInt: "<<poissInt<<std::endl;
        		im(x0,x1,0) = 2*std::sqrt(im(x0,x1,0) + 3.0/8.0);//tF(poissInt);
        		//im(x0,x1,0) = distribution(generator);
        		//std::cout<<im(x0,x1,0)<< "<- after ";
        	}
        }

        MultiArrayView <2, T> array2 = im.bindOuter(0); // select current image
        BasicImageView<T> input = makeBasicImageView(array2);  // access data as BasicImage     

        vigra::FFTWComplexImage fourier(w, h);
        fourierTransform(srcImageRange(input), destImage(fourier));
        
        // there is no squared magnitude accessor, so we use the magnitude here
        vigra::combineTwoImages(srcImageRange(ps), 
                srcImage(fourier, FFTWSquaredMagnitudeAccessor<double>()), 
                destImage(ps), Arg1()+Arg2());

        helper::progress(i, stacksize); // report progress
    }
    std::cout <<std::endl;
    moveDCToCenter(srcImageRange(ps), destImage(ps_center));

    vigra::transformImage(
            srcImageRange(ps_center),
            destIter(res_ul, res_acc),
            Arg1() / Param(stacksize));
}


template <class DestIterator, class DestAccessor, class StormDataSet,class T>
inline
void powerSpectrum(
                   const StormDataSet& im,
                   pair<DestIterator, DestAccessor> ps, std::vector<T> parameterTrafo)
{
    powerSpectrum(im, ps.first, ps.second, parameterTrafo);
}

/**
 * Estimate Noise Power
 */
template <class SrcIterator, class SrcAccessor>
typename SrcIterator::value_type estimateNoisePower(int w, int h,
        SrcIterator is, SrcIterator end, SrcAccessor as)
{
    typedef double sum_type; // use double here since the sum can get larger than float range

    vigra::FindSum<sum_type> sum;   // init functor
    vigra::FindSum<sum_type> sumROI;   // init functor
    vigra::BImage mask(w,h);
    mask = 0;
    // TODO: this size should depend on image dimensions!
    int borderWidth = 10;
    for(int y = borderWidth; y < h-borderWidth; y++) {  // select center of the fft image
        for(int x = borderWidth; x < w-borderWidth; x++) {
            mask(x,y) = 1;
        }
    }
    vigra::inspectImage(is, end, as, sum);
    vigra::inspectImageIf(is, end, as, mask.upperLeft(), mask.accessor(), sumROI);

    sum_type s = sum() - sumROI();
    return s / (w*h - (w-2*borderWidth)*(h-2*borderWidth));
}


template <class SrcIterator, class SrcAccessor>
inline
typename SrcIterator::value_type estimateNoisePower(int w, int h,
                   triple<SrcIterator, SrcIterator, SrcAccessor> ps)
{
    return estimateNoisePower(w, h, ps.first, ps.second, ps.third);
}

/**
 * Construct Wiener Filter using noise power estimated 
 * at high frequencies.
 */
// Wiener filter is defined as 
// H(f) = (|X(f)|)^2/[(|X(f)|)^2 + (|N(f)|)^2]
// where X(f) is the power of the signal and 
// N(f) is the power of the noise
// (e.g., see http://cnx.org/content/m12522/latest/)
template <class T, class DestImage, class StormData>
void constructWienerFilter(StormData& im, 
                DestImage& dest, std::vector<T>& parameterVector) {

    int w = im.shape(0);
    int h = im.shape(1);
    std::cout<<w<< "  "<<h<<std::endl;
    BasicImage<T> ps(w,h);
    powerSpectrum(im, destImage(ps), parameterVector);

    std::ofstream ostreamPS;
	bool writeMatrices = false;
	if(writeMatrices){
		char fnamePS[1000];
		sprintf(fnamePS, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/meanPS.txt");
		ostreamPS.open(fnamePS);}
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			//std::cout<<i<<" 2 "<<j<<std::endl;
			ostreamPS << i<<"  "<< j<<"  "<< ps(i,j)<<std::endl;
		}
	}
	ostreamPS.close();

	std::cout<<std::endl;
	vigra::exportImage(srcImageRange(ps), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/meanPS.png");


    T noise = estimateNoisePower(w,h,srcImageRange(ps));
    std::cout<<"noise: "<<noise<<std::endl;


    transformImage(srcImageRange(ps),destImage(ps),
        		ifThenElse(Arg1()-Param(noise)>Param(0.), Arg1()-Param(noise), Param(0.)));

    double totvar = 0;
	double dx;
	double dy;

	double sum = 0;
	std::cout<<std::endl;

	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			sum += ps(i,j);
			dx = (i - w/2.)/(w/2.);
			dy = (j - h/2.)/(h/2.);
			//std::cout<<"dx: "<<dx<<" dy: "<<dy<<std::endl;
			totvar += ps(i,j) * (dx * dx + dy*dy);
			if(i == w/2){std::cout<<ps(i,j)<<" ,";}
		}
	}
	std::cout<<std::endl;

	double varx = 0;
	double sum2 = 0;
	for(int i = 0; i< w; i++){
		sum2 += ps(i, h/2);
		dx = ((i-w/2.)/(w/2.));
		varx += ps(i,h/2) * dx *dx;
	}

	varx = varx / (sum2-ps(w/2,h/2));
	std::cout<<"Variance: "<< varx<<" Sigma spatial domain 1D: "<<1/(2*3.14*std::sqrt(varx/2.)) <<std::endl;

	double vary = 0;
	double sum3 = 0;
	for(int i = 0; i< h; i++){
		sum3 += ps(w/2, i);
		dy = ((i-h/2.)/(h/2.));
		vary += ps(w/2,i) * dy *dy;
	}

	vary = vary / (sum3-ps(w/2,h/2));
	std::cout<<"Variance: "<< vary<<" Sigma spatial domain 1Dy: "<<1/(2*3.14*std::sqrt(vary/2.)) <<std::endl;


	totvar = totvar / ((sum-ps(w/2,h/2)) *2.); //totvar is the sum of varx and vary in case of symmetric signal
	//double korr2 = std::sqrt(w/50.), korr3 = std::sqrt(h/50.);
	//double korr4 = 1/(korr2*korr3);
	double sigmaSpatial = 1/(2*3.14*std::sqrt(totvar/2.));
	parameterVector[3] = 1/(2*3.14*std::sqrt(vary/2.));//--sigmaSpatial;

	std::cout<<"Variance: "<< totvar<<" Sigma spatial domain: "<<sigmaSpatial <<" sum:"<<sum<<"ps(w,h/2): "<<ps(w/2,h/2)<<std::endl;



    transformImage(srcImageRange(ps),
            destImage(ps),
            ifThenElse(Arg1()-Param(noise)>Param(0.), (Arg1()-Param(noise))/Arg1(), Param(0.)));
    moveDCToUpperLeft(srcImageRange(ps), destImage(dest));



/*
	vigra::FFTWComplexImage spatial(w, h), ps2(w,h), ps3(w,h), fourier(w,h), fourier2(w,h), spatialSq(w, h), spatialSq2(w,h),fourier_center(w,h);
	vigra::FFTWComplexImage::iterator it2;
	vigra::FFTWRealAccessor< double > acc;


	MultiArray<2, FFTWComplex<double> >::iterator it3;
	//MultiArray<2, FFTWComplex<double> > wienerFourierDomain(Shape2(w,h));
	MultiArray<2, FFTWComplex<double> > wienerSpatialDomain(Shape2(w,h));
	MultiArray<2, FFTWComplex<double> > wienerFourierDomain(Shape2(w,h));
	typename DestImage::iterator itdest= dest.begin();
	int i=0;

	for (it3 = wienerFourierDomain.begin(); it3 != wienerFourierDomain.end(); it3++){
		FFTWComplex<double> temp22(*itdest,0);
		*it3 = temp22;
		//std::cout<<*it3<<"  "<<i<<std::endl;
		itdest +=1;
		i++;

	}
	moveDCToUpperLeft(wienerFourierDomain);
	fourierTransformInverse(wienerFourierDomain, wienerSpatialDomain);
	moveDCToCenter(wienerSpatialDomain);
//	moveDCToCenter(srcImageRange(wienerSpatialDomain), destImage(wienerSpatialDomain));
//	FFTWComplex<double> temp22;
//		i = 0;
//		itdest= dest.begin();
//		for (it3 = wienerSpatialDomain.begin(); it3 != wienerSpatialDomain.end(); it3++){
//		//for (it3 = wienerFourierDomain.begin(); i < 10000; it3++){
//			temp22 = *it3;
//			//std::cout<<temp22<<" #$# "<< temp22[0]<<std::endl;
//			//std::cout<<i<<std::endl;
//			*itdest = temp22[0];
//			itdest += 1;
//			i++;
//		}
//	std::cout<<wienerSpatialDomain(10,10)<<std::endl;
//	//transformImage(srcImageRange(wienerSpatialDomain),
//	//					destImage(wienerSpatialDomain), Arg1()*Arg1());
//	std::cout<<wienerSpatialDomain(10,10)<<std::endl;
//
	for(it3 = wienerSpatialDomain.begin(); it3 != wienerSpatialDomain.end(); it3++){
		FFTWComplex<double> temp = *it3;
		//std::cout<<*it3<<temp[0]<<" "<<temp[0]*temp[0]<<"  ";
		FFTWComplex<double> temp2(temp[0]*temp[0],0); *it3 = temp2;
	}
//	std::cout<<wienerSpatialDomain(10,10)<<std::endl;

	fourierTransform(wienerSpatialDomain, wienerFourierDomain);
//
	//moveDCToCenter(srcImageRange(wienerFourierDomain), destImage(wienerFourierDomain));
	//for (itdest = dest.begin(); itdest != dest.end(); itdest++){*itdest = 0;}

	FFTWComplex<double> temp22;
	i = 0;
	itdest= dest.begin();
	double maximum = 0;
	for (it3 = wienerFourierDomain.begin(); it3 != wienerFourierDomain.end(); it3++){
	//for (it3 = wienerFourierDomain.begin(); i < 10000; it3++){
		temp22 = *it3;
		//std::cout<<temp22<<" #$# "<< temp22[0]<<std::endl;
		//std::cout<<i<<std::endl;
		if(temp22[0]>=0){*itdest = 50*temp22[0];}
		if(temp22[0]<0){*itdest = -50*temp22[0];}
		if(*itdest > maximum){maximum = *itdest;}
		itdest += 1;
		i++;
	}
	for(itdest = dest.begin(); itdest != dest.end(); itdest++){
		*itdest = *itdest / maximum;
	}

*/


//	int counter = 0;
//	it3 = wienerFourierDomain.begin();
//	for(int i = 0; i< w; i++){
//		for(int j = 0; j< h/2; j++){
//			temp22 = *it3;
//			dest(i,j) = temp22[0];
//			it3++;
//			std::cout<<counter<<std::endl;
//			counter++;
//		}
//	}
//	it3 = wienerFourierDomain.end();
//	for(int i = 0; i< w; i++){
//		for(int j = h/2; j< h; j++){
//			temp22 = *it3;
//			dest(i,j) = temp22[0];
//			it3--;
//			std::cout<<counter<<std::endl;
//			counter++;
//		}
//	}
//	for (it3 = wienerFourierDomain.end(); it3 != wienerFourierDomain.begin(); it3--){
//		temp22 = *it3;
//		std::cout<<temp22<<" #$# "<< temp22[0]<<std::endl;
//		std::cout<<i<<std::endl;
//		*itdest = temp22[0];
//		itdest += 1;
//		i++;
//	}


//	i = 0;
//	itdest= dest.begin();
//	for (it2 = ps2.begin(); i <w*h; it2++){
//		acc.set(*itdest,it2);
//		itdest +=1;
//		i++;
//		//std::cout<<*itdest<<"  ";
//	}
//	moveDCToUpperLeft(srcImageRange(ps2), destImage(ps3));
//
//	vigra::DImage temp(w, h);
//	temp = 0;
//	transformImage(srcImageRange(ps3, FFTWMagnitudeAccessor<double>()),
//					destImage(temp), Arg1());
//	vigra::exportImage(srcImageRange(temp), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/Wienerps2.png");
//	temp = 0;
//
//	i=0;
//
//
//
//	fourierTransformInverse(srcImageRange(ps3), destImage(spatial));
//	for(it2 = spatial.begin(); i<w*h; it2++){
//		std::cout<<*it2<<" ";
//		i++;
//	}
//	vigra::DImage temp2(w, h);
//	temp = 0;
//	transformImage(srcImageRange(spatial, FFTWMagnitudeAccessor<double>()),
//						destImage(temp2), Arg1());
//	vigra::exportImage(srcImageRange(temp2), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/wienerFilterSpatial.png");
//	temp = 0;
//
//    //vigra::FFTWComplexImage realSpace(w, h);
//    //vigra::fourierTransformInverse(srcImageRange(dest).first,srcImageRange(dest).second, srcImageRange(dest).third, destImage(realSpace).first, destImage(realSpace).second);
//    //vigra::exportImage(srcImageRange(realSpace), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/imgRealWienerFilter.png");
//	double norm_factor = 1./(w*h);
//
//
//    transformImage(srcImageRange(spatial, FFTWMagnitudeAccessor<double>()), destImage(spatialSq), Arg1()*Arg1());
//    transformImage(srcImageRange(spatialSq, FFTWMagnitudeAccessor<double>()),
//    					destImage(temp), Arg1());
//	vigra::exportImage(srcImageRange(temp), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/wienerFilterSpatialSq.png");
//
//    i=0;
////    for(it2 = spatialSq.begin(); i<w*h; it2++){
////		std::cout<<*it2<<" ";
////		i++;
////	}
//    //vigra::DImage fourier_center(w,h);
//
//    moveDCToUpperLeft(srcImageRange(spatialSq), destImage(spatialSq2));
//    fourierTransform(srcImageRange(spatialSq2), destImage(fourier));
//
//	vigra::FindMinMax<T> spatialSqMinmax;
//
//	vigra::inspectImage(srcImageRange(fourier, FFTWMagnitudeAccessor<double>()), spatialSqMinmax);
//	norm_factor = spatialSqMinmax.max;
//
//	std::cout<<"min: "<<spatialSqMinmax.min<<" max: "<<spatialSqMinmax.max<<std::endl;
//	transformImage(srcImageRange(fourier, FFTWMagnitudeAccessor<double>()), destImage(fourier2), Arg1()/Param(norm_factor));
//
//	vigra::inspectImage(srcImageRange(fourier2, FFTWMagnitudeAccessor<double>()), spatialSqMinmax);
//	norm_factor = spatialSqMinmax.max;
//
//	std::cout<<"min: "<<spatialSqMinmax.min<<" max: "<<spatialSqMinmax.max<<std::endl;
//
//    transformImage(srcImageRange(fourier2, FFTWMagnitudeAccessor<double>()),
//           					destImage(dest), Arg1());
////    vigra::exportImage(srcImageRange(temp), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/wienerFilterFourier.png");
////
////    //moveDCToCenter(srcImageRange(fourier), destImage(fourier_center));
//////    //moveDCToUpperLeft(srcImageRange(fourier), destImage(fourier_center));
//////    i=0;
//////
////	temp=0;
////    transformImage(srcImageRange(fourier,  <double>()), destImage(dest), Arg1());
////    transformImage(srcImageRange(fourier, FFTWSquaredMagnitudeAccessor<double>()),
////        					destImage(temp), Arg1());
////    vigra::exportImage(srcImageRange(temp), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/wienerFilterFourierFinal.png");


	std::ofstream ostreamWF;

	if(writeMatrices){
		char fnameWF[1000];
		sprintf(fnameWF, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/Wienerfilter.txt");
		ostreamWF.open(fnameWF);}
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			//std::cout<<i<<" 2 "<<j<<std::endl;
			ostreamWF << i<<"  "<< j<<"  "<< dest(i,j)<<std::endl;
		}
	}
	ostreamWF.close();

}


/** 
 Generate a filter for enhancing the image quality in fourier space.
 Either using constructWienerFilter() or by loading the given file.
  
 @param filter if this file exists, load it. Otherwise create a filter
        from the data and save it to file 'filter'
 @param in 3-dimensional measurement as MultiArrayView<3,float> or MyImageInfo
*/
template <class T, class StormDataSet>
void generateFilter(StormDataSet& in, BasicImage<T>& filter, const std::string& filterfile,
		std::vector<T>& parameterVector) {
    bool constructNewFilter = true;
    if(filterfile != "" && helper::fileExists(filterfile)) {
        vigra::ImageImportInfo filterinfo(filterfile.c_str());

        if(filterinfo.isGrayscale())
        {
            std::cout << "using filter from file " << filterfile << std::endl;
            vigra::BasicImage<T> filterIn(filterinfo.width(), filterinfo.height());
            vigra::importImage(filterinfo, destImage(filterIn)); // read the image
            vigra::resizeImageSplineInterpolation(srcImageRange(filterIn), destImageRange(filter));
            constructNewFilter = false;
        }
        else
        {
            std::cout << "filter image must be grayscale" << std::endl;
        }
    }
    if(constructNewFilter) {
        std::cout << "generating wiener filter from the data" << std::endl;
        constructWienerFilter<T>(in, filter, parameterVector);
        std::cout << "wiener filter constructed"<<parameterVector[3]<<std::endl;
        vigra::exportImage(srcImageRange(filter), filterfile.c_str()); // save to disk
    }
    
}

//--------------------------------------------------------------------------
// STORM DATA PROCESSING
//--------------------------------------------------------------------------

/** 
 * Estimate Background level and subtract it from the image
 */
template <class Image>
void subtractBackground(Image& im, Image& bg) {
    float sigma = 10.; // todo: estimate from data

    vigra::recursiveSmoothX(srcImageRange(im), destImage(bg), sigma);
    vigra::recursiveSmoothY(srcImageRange(bg), destImage(bg), sigma);
    vigra::combineTwoImages(srcImageRange(im), srcImage(bg), destImage(im), Arg1()-Arg2());
}

/**
 * Prefilter for BSpline-Interpolation
 */
template <int ORDER, class Image>
void prefilterBSpline(Image& im) {
    ArrayVector<double> const & b = BSpline<ORDER, double>().prefilterCoefficients();

    for(unsigned int i=0; i<b.size(); ++i)
    {
        recursiveFilterX(srcImageRange(im), destImage(im), b[i], BORDER_TREATMENT_REFLECT);
        recursiveFilterY(srcImageRange(im), destImage(im), b[i], BORDER_TREATMENT_REFLECT);
    }
}

/**
 * Localize Maxima of the spots and return a list with coordinates
 * 
 * This is the actual loop over a microscopic image stack to 
 * reconstruct a super-resolution image out of single molecule detections.
 * 
 * The localization is done on per-frame basis in wienerStormSingleFrame()
 * 
 * @param im MultiArrayView on the actual image data
 */
template <class T>
void wienerStorm(const MultiArrayView<3, T>& im, const BasicImage<T>& filter, 
            std::vector<std::set<Coord<T> > >& maxima_coords, std::vector<T> parameterTrafo,
            const T threshold=800, const int factor=8, const int mylen=9,
            const std::string &frames="", const char verbose=0) {
    
    unsigned int stacksize = im.size(2);
    unsigned int w = im.size(0);
    unsigned int h = im.size(1);
    unsigned int w_xxl = factor*(w-1)+1;
    unsigned int h_xxl = factor*(h-1)+1;
    unsigned int i_stride=1;
    int i_beg=0, i_end=stacksize;
    if(frames!="") {
        helper::rangeSplit(frames, i_beg, i_end, i_stride);
        if(i_beg < 0) i_end = stacksize+i_beg; // allow counting backwards from the end
        if(i_end < 0) i_end = stacksize+i_end; // allow counting backwards from the end
        if(verbose) std::cout << "processing frames [" << i_beg << ":" 
            << i_end << ":" << i_stride << "]" << std::endl;
    }
    


    // TODO: Precondition: res must have size (factor*(w-1)+1, factor*(h-1)+1)
    // filter must have the size of input


    // initialize fftw-wrapper; create plans
    BasicImageView<T> sampleinput = makeBasicImageView(im.bindOuter(0));  // access first frame as BasicImage
    FFTFilter<T> fftwWrapper(sampleinput);

    std::cout << "Finding the maximum spots in the images...2" << std::endl;
    helper::progress(-1,-1); // reset progress

    //over all images in stack
    transformationFunctor tF(parameterTrafo[3],parameterTrafo[5],parameterTrafo[4]);
    #pragma omp parallel for schedule(static, CHUNKSIZE)
    for(int i = i_beg; i < i_end; i+=i_stride) {
        MultiArrayView <2, T> array = im.bindOuter(i); // select current image
        for(int x0 = 0; x0 < w; x0++){
			for(int x1 = 0; x1 < h; x1++){
				//T poissInt = (array(x0,x1,0) - parameterTrafo[1])/parameterTrafo[0];
				//array(x0,x1,0) = tF(poissInt);
			}
		}
        wienerStormSingleFrame(array, filter, maxima_coords[i],
                fftwWrapper, // TODO (this is no real function argument but should be global)
                threshold, factor, mylen, verbose);

        #ifdef OPENMP_FOUND
        if(omp_get_thread_num()==0) { // master thread
            helper::progress(i+1, i_end); // update progress bar
        }
        #else
            helper::progress(i+1, i_end); // update progress bar
        #endif //OPENMP_FOUND
    }
    std::cout << std::endl;
}

/**
 * Localize Maxima of the spots and return a list with coordinates
 *
 * This is the actual loop over a microscopic image stack to
 * reconstruct a super-resolution image out of single molecule detections.
 *
 * The localization is done on per-frame basis in wienerStormSingleFrame()
 *
 * @param info MyImportInfo file info containing the image stack
 */

template <class T>
void wienerStorm(const MyImportInfo& info, const BasicImage<T>& filter,
            std::vector<std::set<Coord<T> > >& maxima_coords, std::vector<T> parameterTrafo,
            MultiArray<3, T>& PoissonMeans,
            const T threshold=800, const int factor=8, const int mylen=9,
            const std::string &frames="", const char verbose=0) {

    unsigned int stacksize = info.shape(2);
    unsigned int w = info.shape(0);
    unsigned int h = info.shape(1);
    unsigned int w_xxl = factor*(w-1)+1;
    unsigned int h_xxl = factor*(h-1)+1;
    unsigned int i_stride=1;
    int i_beg=0, i_end=stacksize;
    if(frames!="") {
        helper::rangeSplit(frames, i_beg, i_end, i_stride);
        if(i_beg < 0) i_end = stacksize+i_beg; // allow counting backwards from the end
        if(i_end < 0) i_end = stacksize+i_end; // allow counting backwards from the end
        if(verbose) std::cout << "processing frames [" << i_beg << ":"
            << i_end << ":" << i_stride << "]" << std::endl;
    }

    // TODO: Precondition: res must have size (factor*(w-1)+1, factor*(h-1)+1)
    // filter must have the size of input
    MultiArray<3, T> im(Shape3(w,h,1));

    // initialize fftw-wrapper; create plans
    readBlock(info, Shape3(0,0,0), Shape3(w,h,1), im);
    BasicImageView<T> sampleinput = makeBasicImageView(im.bindOuter(0));  // access first frame as BasicImage
    FFTFilter<T> fftwWrapper(srcImageRange(sampleinput));

    #ifndef STORM_QT // silence stdout
    std::cout << "Finding the maximum spots in the images..." << std::endl;
    #endif // STORM_QT
    helper::progress(-1,-1); // reset progress

    transformationFunctor tF(1, 3./8,0);
    //over all images in stack
    #pragma omp parallel for schedule(static, CHUNKSIZE) firstprivate(im)
    for(int i = i_beg; i < i_end; i+=i_stride) {
        readBlock(info, Shape3(0,0,i), Shape3(w,h,1), im);
        MultiArrayView <2, T> array = im.bindOuter(0); // select current image

        for(int x0 = 0; x0 < w; x0++){
        	for(int x1 = 0; x1 < h; x1++){array(x0,x1) = (array(x0,x1) - parameterTrafo[1])/parameterTrafo[0];}
        }

        MultiArray<2, T> mask(Shape2(info.shape()[0],info.shape()[1]));
        //applyTransformation(array, w, h, a, offset);
        getMask(array, info, PoissonMeans, i, mask, parameterTrafo);

        for(int x0 = 0; x0 < w; x0++){
			for(int x1 = 0; x1 < h; x1++){array(x0,x1) = tF(array(x0,x1));} // this does Anscombe transformation
		}

        wienerStormSingleFrame(array, filter, maxima_coords[i],
                fftwWrapper, // TODO (this is no real function argument but should be global)
                mask, i, parameterTrafo,
                threshold, factor, mylen, verbose);
        #ifdef OPENMP_FOUND
        if(omp_get_thread_num()==0) { // master thread
            helper::progress(i+1, i_end); // update progress bar
        }
        #else
            helper::progress(i+1, i_end); // update progress bar
        #endif //OPENMP_FOUND
    }
    #ifndef STORM_QT // silence stdout
    std::cout << std::endl;
    #endif // STORM_QT
}

template <class T>
void wienerStormSingleFrame(const MultiArrayView<2, T>& in, const BasicImage<T>& filter,
            std::set<Coord<T> >& maxima_coords,
            FFTFilter<T> & fftwWrapper, MultiArray<2, T>& mask, int framenumber, std::vector<T> parameterTrafo,
            const T threshold=800, const int factor=8, const int mylen=9,
            const char verbose=0) {

    unsigned int w = in.shape(0); // width
    unsigned int h = in.shape(1); // height

    BasicImage<T> filtered(w,h);
    BasicImage<T> bg(w,h), bg2(w,h);        // background
    const int mylen2 = mylen/2;
    unsigned int w_roi = factor*(mylen-1)+1;
    unsigned int h_roi = factor*(mylen-1)+1;
    BasicImage<T> im_xxl(w_roi, h_roi);

    BasicImageView<T> input = makeBasicImageView(in);  // access data as BasicImage

    //fft, filter with Wiener filter in frequency domain, inverse fft, take real part
    BasicImageView<T> filteredView(filtered.data(), filtered.size());

    BasicImage<T> unfiltered(w,h);

    vigra::copyImage(srcImageRange(input), destImage(unfiltered));
    //std::cout<<"using a sigma of:"<<parameterTrafo[3]<<std::endl;

//    vigra::FindMinMax<T> Minmax;
//    vigra::inspectImage(srcImageRange(input), Minmax);




    //gaussianSmoothing(srcImageRange(input), destImage(filteredView), parameterTrafo[3]);
    gaussianSmoothing(srcImageRange(input), destImage(filteredView), 1.4);


//    vigra::FindMinMax<T> Minmax2;
//    vigra::inspectImage(srcImageRange(filteredView), Minmax2);
//    std::cout<<"before max:"<<Minmax.max<<" after max:"<<Minmax2.max<<std::endl;
//    std::cout<<"before min:"<<Minmax.min<<" after min:"<<Minmax2.min<<std::endl;

    //fftwWrapper.applyFourierFilter(srcImageRange(input), srcImage(filter), destImage(filteredView));
    //vigra::exportImage(srcImageRange(filtered), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/imgAfterWienerFilter.png");
    //~ vigra::gaussianSmoothing(srcImageRange(input), destImage(filtered), 1.2);
    subtractBackground(filtered, bg);
    subtractBackground(unfiltered, bg2);

    //vigra::FindMinMax<T> bgMinmax;
    //vigra::inspectImage(srcImageRange(bg), bgMinmax);


    //T baseline = bgMinmax.min;
    //std::cout<<"baseline: "<< baseline<<'\n';

	applyMask(filtered, mask, framenumber);

	vigra::FindMinMax<T> filteredMinMax;
	inspectImage(srcImageRange(filtered), filteredMinMax);
	//T thresh = (filteredMinMax.max - filteredMinMax.min)*0.10+ filteredMinMax.min;

    std::set<Coord<T> > maxima_candidates_vect;  // we use a set for the coordinates to automatically squeeze duplicates
                                                 // (from overlapping ROIs)
    VectorPushAccessor<Coord<T>, typename BasicImage<T>::const_traverser> maxima_candidates(maxima_candidates_vect, filtered.upperLeft());
    vigra::localMaxima(srcImageRange(filtered), destImage(filtered, maxima_candidates), vigra::LocalMinmaxOptions().threshold(threshold));

    VectorPushAccessor<Coord<T>, typename BasicImage<T>::const_traverser> maxima_acc(maxima_coords, im_xxl.upperLeft());
    //upscale filtered image regions with spline interpolation
    std::set<Coord<float> >::iterator it2;

	//vigra::exportImage(srcImageRange(filtered), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/filtered.png");
	//vigra::exportImage(srcImageRange(unfiltered), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/unfiltered.png");
//	vigra::FindMinMax<T> Minmax3;
//	vigra::inspectImage(srcImageRange(unfiltered), Minmax3);
//	std::cout<<"before max:"<<Minmax3.max<<" after max:"<<Minmax3.max<<std::endl;
//	std::cout<<"before min:"<<Minmax3.min<<" after min:"<<Minmax3.min<<std::endl;
    for(it2=maxima_candidates_vect.begin(); it2 != maxima_candidates_vect.end(); it2++) {
            Coord<float> c = *it2;

            if(unfiltered(c.x,c.y)<1) { // skip very low signals with SNR lower 3
            	std::cout<<"value skipped: "<<unfiltered(c.x,c.y)<<" bg(x,y): "<<bg(c.x,c.y)<<std::endl;
                continue;
            }
            //std::cout<<"value take: "<<filtered(c.x,c.y)<<" bg(x,y): "<<bg(c.x,c.y)<<std::endl;
            Diff2D roi_ul (c.x-mylen2, c.y-mylen2);
            Diff2D roi_lr (c.x-mylen2+mylen, c.y-mylen2+mylen);

            Diff2D xxl_ul (0, 0);  // offset in xxl image
            Diff2D xxl_lr (0, 0);

            // Maxima-Candidates near the border
            if(c.x-mylen2<0 || c.y-mylen2<0 || c.x-mylen2+mylen>(int)w || c.y-mylen2+mylen>(int)h) {
                Diff2D _roi_ul (
                    ((c.x-mylen2)<0) ? 0 : (c.x-mylen2),
                    ((c.y-mylen2)<0) ? 0 : (c.y-mylen2) );
                Diff2D _roi_lr (
                    ((c.x-mylen2+mylen)>(int)w) ? w : (c.x-mylen2+mylen),
                    ((c.y-mylen2+mylen)>(int)h) ? h : (c.y-mylen2+mylen) );

                xxl_ul += (_roi_ul-roi_ul)*factor; // offset in xxl image
                xxl_lr += (_roi_lr-roi_lr)*factor;
                roi_ul = _roi_ul;
                roi_lr = _roi_lr;

            }

            vigra::resizeImageSplineInterpolation(
                    srcIterRange(filtered.upperLeft()+roi_ul, filtered.upperLeft()+roi_lr),
                    destIterRange(im_xxl.upperLeft()+xxl_ul, im_xxl.lowerRight()+xxl_lr),
                    BSplineWOPrefilter<3,double>());
            // find local maxima that are above a given threshold
            // at least the values should be above background+baseline
            // here we include only internal pixels, no border
            // to get every maximum only once, the maxima are pushed into a std::set
            maxima_acc.setOffset(Diff2D(factor*(c.x-mylen2), factor*(c.y-mylen2)));
            //std::cout<<roi_ul<<"  "<<roi_lr<<"  "<<xxl_ul<<"  "<<xxl_lr<<"  "<< framenumber<<std::endl;
            vigra::localMaxima(srcIterRange(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), im_xxl.lowerRight()+xxl_lr-Diff2D(factor,factor)),
                    destIter(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), maxima_acc), vigra::LocalMinmaxOptions().threshold(threshold));
    }
//    typename std::set<Coord<T> >::iterator iter3;
//    for(iter3 = maxima_coords.begin(); iter3 != maxima_coords.end(); iter3++){
//    	Coord<float> c = *iter3;
//    	std::cout<<c.x<<" "<<c.y<<" "<<c.val<<" "<<framenumber<<std::endl;
//    }

    determineAsymmetry(srcImageRange(unfiltered), maxima_coords, factor);
    determineSNR(srcImageRange(unfiltered), maxima_coords, factor);
}
