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
#include <vigra/accumulator.hxx>
#include <vigra/multi_resize.hxx>

#include <ctime>
//#include <set>
#include <fstream>
//#include <iomanip>
#include <vector>
#include <valarray>

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

#include <Rinternals.h>
#include <Rinterface.h>
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
            const MultiArrayShape<3>::type & shape, int factor, float pxSize) {
    int numSpots = 0;
    std::set<Coord<float> >::iterator it2;
    std::ofstream cfile (filename.c_str());
    cfile << shape[0] << " " << shape[1] << " " << shape[2] << std::endl;
    cfile << std::fixed; // fixed instead of scientific format
    for(unsigned int j = 0; j < coords.size(); j++) {
        for(it2=coords[j].begin(); it2 != coords[j].end(); it2++) {
            numSpots++;
            const Coord<float>& c = *it2;
            cfile << std::setprecision(3) << (float)c.x/factor * pxSize << " " << (float)c.y/factor * pxSize << " "
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
	sprintf(filename, "/home/herrmannsdoerfer/master/workspace/output/selectedPoints.txt");
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
    int nbins = 10;

    std::string rScript(info.executableDir());
    rScript.append("/").append(STORM_RSCRIPT);
    if (!helper::fileExists(rScript)) {
        rScript.clear();
        rScript.append(STORM_RSCRIPT_DIR).append(STORM_RSCRIPT);
    }

    SEXP fun, t, tmp, tmp2;
    PROTECT(tmp = Rf_ScalarInteger(42));

    PROTECT(fun = t = Rf_allocList(2));
    SET_TYPEOF(fun, LANGSXP);
    SETCAR(t, Rf_install("set.seed"));
    t = CDR(t);
    SETCAR(t, tmp);
    Rf_eval(fun, R_GlobalEnv);
    UNPROTECT(2);

    PROTECT(tmp = Rf_mkString(rScript.c_str()));
    PROTECT(fun = t = Rf_allocList(2));
    SET_TYPEOF(fun, LANGSXP);
    SETCAR(t, Rf_install("parse"));
    t = CDR(t);
    SETCAR(t, tmp);
    SET_TAG(t, Rf_install("file"));
    PROTECT(tmp2 = Rf_eval(fun, R_GlobalEnv));
    for (R_len_t i = 0; i < Rf_length(tmp2); ++i) {
        Rf_eval(VECTOR_ELT(tmp2, i), R_GlobalEnv);
    }
    UNPROTECT(3);

    PROTECT(tmp = Rf_allocMatrix(REALSXP, numberPoints, 2));
    double *mat = REAL(tmp);
    for (int row = 0; row < numberPoints; ++row) {
        mat[row] = meanValues[row];
    }
    for (int row = 0; row < numberPoints; ++row) {
        mat[numberPoints + row] = skellamParameters[row];
    }

    SEXP bins;
    PROTECT(bins = Rf_ScalarInteger(nbins));
    PROTECT(fun = t = Rf_allocList(3));
    SET_TYPEOF(fun, LANGSXP);
    SETCAR(t, Rf_install("fit.storm.points"));
    t = CDR(t);
    SETCAR(t, tmp);
    t = CDR(t);
    SETCAR(t, bins);
    PROTECT(tmp = Rf_eval(fun, R_GlobalEnv));

    if (tmp == R_NilValue) {
        std::cerr << "Fit could not be performed, exiting..." << std::endl;
        std::exit(1);
    } else {
        double *coefs = REAL(tmp);
        parameterTrafo[0] = coefs[1];
        parameterTrafo[1] = -coefs[0] / coefs[1];
        parameterTrafo[2] = coefs[0];
        std::cout<<"slope: "<< parameterTrafo[0]<<" x0: "<<parameterTrafo[1] << " intercept: "<<parameterTrafo[2] << std::endl;
    }
    UNPROTECT(4);
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
		sprintf(temp, "/home/herrmannsdoerfer/master/workspace/output/intensities/pos0_%d_pos1_%d.txt", vecw[i], vech[i]);
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

//get Mask calculates the probabilities for each pixel of the current frame to be foreground,
//there are two different methods available: logliklihood, based on cumulative distribution function
template <class T>
void getMask(MultiArrayView<2,T >& array, int w, int h, int stacksize, MultiArrayView<2, T>& PoissonMeans ,int framenumber,MultiArray<2,T>& mask, std::vector<T> parameterTrafo){

	T alpha = 0.05;

	T factor2ndDist; //= 1+arrMinmax.max/muMin/10; // should be dependent of snr

	std::ofstream origimg, maski;

	bool writeMatrices = (framenumber == 1);
	if(writeMatrices){
		char origimgname[1000], maskname[1000];
		sprintf(origimgname, "/home/herrmannsdoerfer/master/workspace/output/outputMyStorm/origimgbeforefilter%d.txt", framenumber);
		sprintf(maskname, "/home/herrmannsdoerfer/master/workspace/output/outputMyStorm/maskPoiss%d.txt", framenumber);

		origimg.open (origimgname);
		maski.open (maskname);}
	double cdf = 0;
	for(int i = 0; i< w; i++){

		for(int j = 0; j < h; j++){
//			factor2ndDist = 1+arrMinmax.max/PoissonMeans(i,j)/10;
//			mask(i,j) = isSignal_fast(array(i,j), muMin, factor2ndDist);
			cdf = ppois(array(i,j), PoissonMeans(i,j), 1,0);

			if (cdf > 1- alpha and cdf < 1){
				mask(i,j) = (cdf - (1-alpha))/alpha;
			}
			else if (cdf == 1){mask(i,j) = 1;}
			else if (cdf < 1- alpha){mask(i,j) = 0;}

			if (writeMatrices){
				maski << i<<"  "<< j<<"  "<< mask(i,j)<<std::endl;}
		}
	}
	if (writeMatrices){
	origimg.close();
	maski.close();}
	T beta = 0.5;

	//performGraphcut(mask,beta);
	gaussianSmoothing(srcImageRange(mask), destImage(mask), parameterTrafo[3]);
	transformImage(srcImageRange(mask),destImage(mask),
	        		ifThenElse(Arg1()>Param(2/(2*3.14*parameterTrafo[3])), Param(1.), Param(0.)));

	std::ofstream maskafter;
	bool writeMatrices2 = (framenumber == 1);
	if(writeMatrices2){
		char maskaftername[1000];

		sprintf(maskaftername, "/home/herrmannsdoerfer/master/workspace/output/outputMyStorm/maskafterGraphcutPoiss%d.txt", framenumber);
		maskafter.open (maskaftername);

	for(int i = 0; i< w; i++){
		for(int j = 0; j < h; j++){
			if (writeMatrices2){
			maskafter << i<<"  "<< j<<"  "<< mask(i,j)<<std::endl;}
		}
	}
	maskafter.close();
	}
}

//two poisson distributions are compared, the two distributions are, the distribution around
//the estimated mean for the background pixels and a second distribution with a mean scaled with the
//factor factor2ndDist, each intensity  is explained better by either the background or the
//"foreground" distribution.
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
void getSmoothedPixelMeans(const MyImportInfo& info,T a,T offset,MultiArray<3,T> &regionMeans){
    unsigned int stacksize = info.shape(2);
    unsigned int w = info.shape(0);
    unsigned int h = info.shape(1);

    unsigned int blocks = 30;

    MultiArray<3, T> temp_arr(Shape3(w,h,blocks));
    auto *test = new MultiArray<3, T>(Shape3(w, h, blocks));
    delete test;
    MultiArray<3, T> tempRegionMeans(Shape3(std::ceil(w / (float)blocks), std::ceil(h / (float)blocks), std::ceil(stacksize / (float)blocks)));
    MultiArray<3, int> labels(Shape3(w, h, blocks));
    regionMeans.reshape(Shape3(w, h, stacksize));

    for (int y = 0, f = 0; y < std::ceil(h / (float)blocks); ++y) {
        for (int x = 0; x < std::ceil(w / (float)blocks); ++x, ++f) {
            vigra::Shape3 index(std::min((x + 1) * blocks, w), std::min((y + 1) * blocks, h), blocks);
            auto roi = labels.subarray(vigra::Shape3(x * blocks, y * blocks, 0), index);
            vigra::initMultiArray(destMultiArrayRange(roi), f);
        }
    }

    for(int f = 0; f < std::ceil(stacksize / (float)blocks); f++) {
        int zShape;
        MultiArrayView<3, T> temp_arrView;
        MultiArrayView<3, int> labelsView;
        if ((f + 1) * blocks <= stacksize) {
            zShape = blocks;
            temp_arrView = temp_arr; // should be assignment in MultiArrayView => no copying
            labelsView = labels;
        } else {
            zShape = stacksize % blocks;
            temp_arrView = temp_arr.subarray(Shape3(0, 0, 0), Shape3(w, h, zShape));
            labelsView = labels.subarray(Shape3(0, 0, 0), Shape3(w, h, zShape));
        }
        readBlock(info, Shape3(0,0,f * blocks),Shape3(w,h, zShape), temp_arrView);
        vigra::transformMultiArray(srcMultiArrayRange(temp_arrView), destMultiArrayRange(temp_arrView), [&a, &offset](T p){return (p - offset) / a;});
        vigra::acc::AccumulatorChainArray<typename CoupledIteratorType<3, T, int>::type::value_type, vigra::acc::Select<vigra::acc::DataArg<1>, vigra::acc::LabelArg<2>, vigra::acc::StandardQuantiles<vigra::acc::AutoRangeHistogram<0>>>> accChain;
        auto iter = vigra::createCoupledIterator(temp_arrView, labelsView);
        auto iterEnd = iter.getEndIterator();
        vigra::acc::extractFeatures(iter, iterEnd, accChain);
        for (int y = 0, i = 0; y < std::ceil(h / (float)blocks); ++y) {
            for (int x = 0; x < std::ceil(w / (float)blocks); ++x, ++i) {
                tempRegionMeans(x, y, f) = 0.75 * vigra::acc::get<vigra::acc::StandardQuantiles<vigra::acc::AutoRangeHistogram<0>>>(accChain, i)[3];
                auto test = tempRegionMeans(x, y, f);

            }
        }
    }
    vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(tempRegionMeans), destMultiArrayRange(regionMeans), vigra::BSpline<3>());
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

	for(int f = 0; f < stacksize; f++) {
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

	MultiArray<3, T> mean_im_temp(Shape3(w,h,1)), mean_im_temp2(Shape3(w,h,1));

	std::vector<T> pixelValues;
	int pixelPerFrame = w*h;
	int minMean = INT_MAX, maxMean = 0;
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
}

template <class T>
void applyMask(BasicImage<T>& img, MultiArrayView<2, T>& mask, int frnr){
	int w = mask.shape(0);
	int h = mask.shape(1);

	std::ofstream imgAfterMaskFile, imgBeforeMaskFile;
	bool writeMatrices = true;
	if(writeMatrices){

		char imgAfterMask[1000];
		char imgBeforeMask[1000];

		sprintf(imgAfterMask, "/home/herrmannsdoerfer/master/workspace/output/outputMyStorm/filterdImgAfterMask%d.txt", frnr);
		sprintf(imgBeforeMask, "/home/herrmannsdoerfer/master/workspace/output/outputMyStorm/filterdImgBeforeMask%d.txt", frnr);

		imgAfterMaskFile.open(imgAfterMask);
		imgBeforeMaskFile.open(imgBeforeMask);}

	for(int i = 0; i< w; i++){
		for(int j = 0; j<h; j++){
			if(writeMatrices){imgBeforeMaskFile << i<<"  "<< j<<"  "<< img(i,j)<<std::endl;}
			if(mask(i,j)== 0.0){img(i,j)= 0;}
			if(writeMatrices){imgAfterMaskFile << i<<"  "<< j<<"  "<< img(i,j)<<std::endl;}

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

    for(unsigned int i = 0; i < stacksize; i++) {
        MultiArrayView <2, T> array2 = array.bindOuter(i); // select current image
        for(it = array2.begin(); it != array2.end(); it++){
        	if((*it-b)/a > 0){
        		*it = (*it-b)/a;
        		*it = (2*std::sqrt(*it + 3.0/8.0));
        	}
        	else{*it = 2*std::sqrt(0+3./8);}
        }

        BasicImage<T> bg(w,h);        // background
		vigra::exportImage(srcImageRange(array2), "/home/herrmannsdoerfer/master/workspace/output/vorBGSubstr.png");
		subtractBackground(array2, bg);
		vigra::exportImage(srcImageRange(array2), "/home/herrmannsdoerfer/master/workspace/output/nachBGSubstr.png");

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

    transformationFunctor tF(1, 3./8., 0);
    for(unsigned int i = 0; i < stacksize; i++) {
        readBlock(info, Shape3(0,0,i), Shape3(w,h,1), im);
        for(int x0 = 0; x0 < w; x0++){
        	for(int x1 = 0; x1 < h; x1++){
        		if(((im(x0,x1,0) - parameterTrafo[1])/parameterTrafo[0])>0){
        		im(x0,x1,0) = (im(x0,x1,0) - parameterTrafo[1])/parameterTrafo[0];}
        		else{im(x0,x1,0) = 0;}
        		im(x0,x1,0) = 2*std::sqrt(im(x0,x1,0) + 3.0/8.0);//tF(poissInt);
        	}
        }

        BasicImage<T> bg(w,h);        // background
        BasicImageView<T> bgV(bg.data(), w,h);
        BasicImageView<T> tmp = makeBasicImageView(im.bindOuter(0));

		vigra::exportImage(srcImageRange(tmp), "/home/herrmannsdoerfer/master/workspace/output/vorBGSubstr.png");
		subtractBackground(tmp, bgV);
		vigra::exportImage(srcImageRange(tmp), "/home/herrmannsdoerfer/master/workspace/output/nachBGSubstr.png");

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
template <class T, class StormData>
void constructWienerFilter(StormData& im, std::vector<T>& parameterVector) {

    int w = im.shape(0);
    int h = im.shape(1);
    std::cout<<w<< "  "<<h<<std::endl;
    BasicImage<T> ps(w,h);
    powerSpectrum(im, destImage(ps), parameterVector);

    std::ofstream ostreamPS;
	bool writeMatrices = false;
	if(writeMatrices){
		char fnamePS[1000];

		//sprintf(fnamePS, "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/meanPS.txt");

		ostreamPS.open(fnamePS);}
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			//std::cout<<i<<" 2 "<<j<<std::endl;
			ostreamPS << i<<"  "<< j<<"  "<< ps(i,j)<<std::endl;
		}
	}
	ostreamPS.close();

	std::cout<<std::endl;

	//vigra::exportImage(srcImageRange(ps), "/home/herrmannsdoerfer/master/workspace/myStorm/storm/build/output/meanPS.png");

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
	double vary = 0;
	double sum3 = 0;

	for(int i = 0; i< w; i++){
		sum2 += ps(i, h/2);
		dx = ((i-w/2.)/(w/2.));
		varx += ps(i,h/2) * dx *dx;

		sum3 += ps(w/2, i);
		dy = ((i-h/2.)/(h/2.));
		vary += ps(w/2,i) * dy *dy;
	}
	varx = varx / (sum2-ps(w/2,h/2));
	vary = vary / (sum3-ps(w/2,h/2));

	std::cout<<"Variance: "<< varx<<" Sigma spatial domain 1D: "<<1/(2*3.14*std::sqrt(varx/2.)) <<std::endl;
	std::cout<<"Variance: "<< vary<<" Sigma spatial domain 1Dy: "<<1/(2*3.14*std::sqrt(vary/2.)) <<std::endl;

	totvar = totvar / ((sum-ps(w/2,h/2)) *2.); //totvar is the sum of varx and vary in case of symmetric signal
	double sigmaSpatial = 1/(2*3.14*std::sqrt(totvar/2.));
	parameterVector[3] = sigmaSpatial;

	std::cout<<"Variance: "<< totvar<<" Sigma spatial domain: "<<sigmaSpatial <<" sum:"<<sum<<"ps(w,h/2): "<<ps(w/2,h/2)<<std::endl;
}


/**
 Generate a filter for enhancing the image quality in fourier space.
 Either using constructWienerFilter() or by loading the given file.

 @param filter if this file exists, load it. Otherwise create a filter
        from the data and save it to file 'filter'
 @param in 3-dimensional measurement as MultiArrayView<3,float> or MyImageInfo
*/
template <class T, class StormDataSet>
void generateFilter(StormDataSet& in, std::vector<T>& parameterVector) {
    bool constructNewFilter = true;
    if(!in.getFilterfile().empty() && helper::fileExists(in.getFilterfile())) {
        std::ifstream filter(in.getFilterfile());
        double filterWidth;
        filter >> filterWidth;
        parameterVector[3] = filterWidth;
        if (!filter.fail() && !filter.bad()) {
            constructNewFilter = false;
            std::cout << "using filter from file " << in.getFilterfile() << std::endl;
        }
        filter.close();
    }
    if(constructNewFilter) {
        std::cout << "generating wiener filter from the data" << std::endl;
        constructWienerFilter(in, parameterVector);
        std::cout << "wiener filter constructed"<<parameterVector[3]<<std::endl;
        std::ofstream filter(in.getFilterfile());
        filter << parameterVector[3] << std::endl;
        filter.close();
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
void wienerStorm(const MyImportInfo& info, std::vector<std::set<Coord<T> > >& maxima_coords,
			std::vector<T> parameterTrafo, MultiArray<3, T>& PoissonMeans) {

    unsigned int stacksize = info.shape(2);
    unsigned int w = info.shape(0);
    unsigned int h = info.shape(1);
    unsigned int w_xxl = info.getFactor()*(w-1)+1;
    unsigned int h_xxl = info.getFactor()*(h-1)+1;
    unsigned int i_stride=1;
    int i_beg=0, i_end=stacksize;
    if(!info.getFrameRange().empty()) {
        helper::rangeSplit(info.getFrameRange(), i_beg, i_end, i_stride);
        if(i_beg < 0) i_end = stacksize+i_beg; // allow counting backwards from the end
        if(i_end < 0) i_end = stacksize+i_end; // allow counting backwards from the end
        if(info.verbose) std::cout << "processing frames [" << i_beg << ":"
            << i_end << ":" << i_stride << "]" << std::endl;
    }

    // TODO: Precondition: res must have size (info.getFactor()*(w-1)+1, info.getFactor()*(h-1)+1)
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

        MultiArray<2, T> mask(Shape2(w,h));
        MultiArrayView<2, T> poissView = PoissonMeans.bindOuter(i);
        getMask(array, w,h,stacksize, poissView, i, mask, parameterTrafo);

        for(int x0 = 0; x0 < w; x0++){
			for(int x1 = 0; x1 < h; x1++){array(x0,x1) = tF(array(x0,x1));} // this does Anscombe transformation
		}

        wienerStormSingleFrame(info, array, maxima_coords[i],
                fftwWrapper, // TODO (this is no real function argument but should be global)
                mask, i, parameterTrafo);
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
void wienerStormSingleFrame(const MyImportInfo& info, const MultiArrayView<2, T>& in, std::set<Coord<T> >& maxima_coords,
            FFTFilter<T> & fftwWrapper, MultiArray<2, T>& mask, int framenumber, std::vector<T> parameterTrafo)
          	{

    unsigned int w = in.shape(0); // width
    unsigned int h = in.shape(1); // height

    BasicImage<T> filtered(w,h);
    BasicImage<T> bg(w,h), bg2(w,h);        // background

    int factor = info.getFactor();
    int mylen = info.getRoilen();

    const int mylen2 = mylen/2;
    unsigned int w_roi = factor*(mylen-1)+1;
    unsigned int h_roi = factor*(mylen-1)+1;
    BasicImage<T> im_xxl(w_roi, h_roi);

    BasicImageView<T> input = makeBasicImageView(in);  // access data as BasicImage

    //fft, filter with Wiener filter in frequency domain, inverse fft, take real part
    BasicImageView<T> filteredView(filtered.data(), filtered.size());

    BasicImage<T> unfiltered(w,h);

    vigra::copyImage(srcImageRange(input), destImage(unfiltered));

    gaussianSmoothing(srcImageRange(input), destImage(filteredView), parameterTrafo[3]/2.);
    //gaussianSmoothing(srcImageRange(input), destImage(filteredView), 1.4);

    //fftwWrapper.applyFourierFilter(srcImageRange(input), srcImage(filter), destImage(filteredView));
    subtractBackground(filtered, bg);
    subtractBackground(unfiltered, bg2);

	applyMask(filtered, mask, framenumber);

	vigra::FindMinMax<T> filteredMinMax;
	inspectImage(srcImageRange(filtered), filteredMinMax);
	//T thresh = (filteredMinMax.max - filteredMinMax.min)*0.10+ filteredMinMax.min;

    std::set<Coord<T> > maxima_candidates_vect;  // we use a set for the coordinates to automatically squeeze duplicates
                                                 // (from overlapping ROIs)
    VectorPushAccessor<Coord<T>, typename BasicImage<T>::const_traverser> maxima_candidates(maxima_candidates_vect, filtered.upperLeft());
    vigra::localMaxima(srcImageRange(filtered), destImage(filtered, maxima_candidates), vigra::LocalMinmaxOptions().threshold(0));

    VectorPushAccessor<Coord<T>, typename BasicImage<T>::const_traverser> maxima_acc(maxima_coords, im_xxl.upperLeft());
    //upscale filtered image regions with spline interpolation
    std::set<Coord<float> >::iterator it2;


//	vigra::exportImage(srcImageRange(filtered), "/home/herrmannsdoerfer/master/workspace/output/filtered.png");
//	vigra::exportImage(srcImageRange(unfiltered), "/home/herrmannsdoerfer/master/workspace/output/unfiltered.png");

    for(it2=maxima_candidates_vect.begin(); it2 != maxima_candidates_vect.end(); it2++) {
            Coord<float> c = *it2;
            //std::cout<<"value not skipped: "<<unfiltered(c.x,c.y)<<" bg(x,y): "<<bg(c.x,c.y)<<std::endl;
            if(unfiltered(c.x,c.y)<3) { // skip very low signals with SNR lower 3
            	//std::cout<<"value skipped: "<<unfiltered(c.x,c.y)<<" bg(x,y): "<<bg(c.x,c.y)<<std::endl;
                continue;
            }
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
                    destIter(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), maxima_acc), vigra::LocalMinmaxOptions().threshold(0));
    }
    determineAsymmetry(srcImageRange(unfiltered), maxima_coords, factor);
    determineSNR(srcImageRange(unfiltered), maxima_coords, factor);
}
