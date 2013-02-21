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
#include <vigra/multi_impex.hxx>

#include <ctime>
//#include <set>
#include <fstream>
//#include <iomanip>
#include <vector>
#include <valarray>
#include <limits>

//#include <algorithm>
#ifdef OPENMP_FOUND
    #include <omp.h>
#endif //OPENMP_FOUND
#include "util.h"
#include "fftfilter.hxx"
#include "dataparams.h"
//#include <iostream>

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
template <class T, class S, class ITERATOR>
class VectorPushAccessor{
    public:
        typedef typename T::value_type value_type;
        VectorPushAccessor(std::set<T>& arr, ITERATOR it_start, int factor, const MultiArrayView<2, S> &mask)
            : m_arr(arr), m_it_start(it_start), m_offset(), m_factor(factor), m_mask(mask) {}

        T const &   operator() (ITERATOR const &i) const {
            return NumericTraits<T>::zero();
        }
        template<class V>
        void    set (V const & /*value*/, ITERATOR const &i) {
            int x = i.x+m_offset.x;
            int y = i.y-m_it_start.y+m_offset.y;
            if (m_mask(std::round(x / m_factor), std::round(y / m_factor)) > 0) {
                m_arr.insert(T(x, y, *i));
            }
        }
        void setOffset(Diff2D offset) {
            m_offset = offset;
        }

    private:
        std::set<T>& m_arr;
        ITERATOR m_it_start;
        Diff2D m_offset;
        float m_factor; //needs to be float for division, avoids casting
        const MultiArrayView<2, S> &m_mask;
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
int saveCoordsFile(const DataParams &params, std::ofstream &cfile, const std::vector<std::set<C> >& coords) {
    int numSpots = 0;
    std::set<Coord<float> >::iterator it2;
    cfile << params.shape(0) << " " << params.shape(1) << " " << params.shape(2) << std::endl;
    cfile << std::fixed; // fixed instead of scientific format
    for(unsigned int j = 0; j < coords.size(); j++) {
        for(it2=coords[j].begin(); it2 != coords[j].end(); it2++) {
            numSpots++;
            const Coord<float>& c = *it2;
            cfile << std::setprecision(3) << (float)c.x/params.getFactor() * params.getPixelSize() << " " << (float)c.y/params.getFactor() * params.getPixelSize() << " "
                << j << " " << std::setprecision(1) << c.val << " " << std::setprecision(3) << c.asymmetry << " "
                << c.signalNoiseRatio << std::endl;
        }
    }
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
class transformationFunctor {
public:
	transformationFunctor(float a, float intercept, float minInt = 0): a(a), intercept(intercept), C(-2/a*std::sqrt(a*minInt+ intercept)){}
	float operator()(float poissInt) const {return 2/a*std::sqrt(a*poissInt + intercept)+ C;}

private:
	float a;
	float intercept;
	float C;
};

template <class T>
void findBestFit(DataParams &params,T meanValues[],T skellamParameters[],int numberPoints){
    int nbins = 10;

    std::string rScript(params.executableDir());
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
        params.setSlope(coefs[1]);
        params.setIntercept(-coefs[0] / coefs[1]);
        std::cout<<"slope: "<< params.getSlope()<<" x0: "<<params.getIntercept() << std::endl;
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

void printIntensities(const DataParams &params, int* vecw, int* vech, int nbrPoints){
    unsigned int stacksize2 = params.shape(2);
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);
    float a = params.getSlope(), b = params.getIntercept();
	//std::cout<<"w: "<<w<<" h: "<<h<<std::endl;
	//std::cin.get();
	MultiArray<3, float> mean_im_temp(Shape3(w,h,1));
	MultiArray<3, float> mean_im_temp2(Shape3(w,h,1));


	std::ofstream origimg[nbrPoints];

	char temp[1000];
	for(int i = 0;i < nbrPoints; i++){
		sprintf(temp, "/home/herrmannsdoerfer/tmpOutput/pos0_%d_pos1_%d.txt", vecw[i], vech[i]);
		origimg[i].open(temp);
		std::cout<<"vecw["<<i<<"]="<<vecw[i]<<" vech["<<i<<"]="<<vech[i]<<std::endl;
		std::cout<<"a:"<< a<<" b" <<b<<std::endl;
	}



	for(int f = 0; f< stacksize2;f++){
        params.readBlock(Shape3(0,0,f),Shape3(w,h,1), mean_im_temp2);
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
void getMask(const DataParams &params, const BasicImage<T>& array, int framenumber, MultiArray<2,T>& mask){
	double cdf = params.getMaskThreshold();
    vigra::transformImage(srcImageRange(array), destImage(mask), [&cdf](T p) {return p >= cdf ? 1 : 0;});

    vigra::IImage labels(array.width(), array.height());
    unsigned int nbrCC = vigra::labelImageWithBackground(srcImageRange(mask), destImage(labels), false, 0);
    std::valarray<int> bins(0, nbrCC + 1);
    auto fun = [&bins](int32_t p){++bins[p];};
    vigra::inspectImage(srcImageRange(labels), fun);
    vigra::transformImage(srcImageRange(labels), destImage(mask), [&params, &bins](T p) {if(!p || bins[p] < 3.14 * std::pow(params.getSigma(), 2)) return 0; else return 1;});
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

/**
 * Estimate Background level and subtract it from the image
 */
template <class Image>
void subtractBackground(Image& im) {
    float sigma = 10.; // todo: estimate from data
    BasicImage<typename Image::value_type> bg(im.size());
    vigra::recursiveSmoothX(srcImageRange(im), destImage(bg), sigma);
    vigra::recursiveSmoothY(srcImageRange(bg), destImage(bg), sigma);
    vigra::combineTwoImages(srcImageRange(im), srcImage(bg), destImage(im), Arg1()-Arg2());
}



//--------------------------------------------------------------------------
// GENERATE WIENER FILTER
//--------------------------------------------------------------------------
// Since the algorithms work on single-frames only, there is no need to
// put the complete dataset into RAM but every frame can be read from disk
// when needed. The class MyImportInfo transparently handles hdf5 and sif
// input file pointers.

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
template <class T>
void constructWienerFilter(DataParams &params, BasicImage<T> &ps) {

    int w = params.shape(0);
    int h = params.shape(1);
    std::cout<<w<< "  "<<h<<std::endl;

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
    params.setSigma(sigmaSpatial);

    std::cout<<"Variance: "<< totvar<<" Sigma spatial domain: "<<sigmaSpatial <<" sum:"<<sum<<"ps(w,h/2): "<<ps(w/2,h/2)<<std::endl;
}


//To estimate the gain factor points with different mean intensities are needed. This functions searches for
//good candidates, it tries to find pixels with as much different mean values as possible.
template <class T>
void estimateParameters(DataParams &params) {
    bool needFilter = !(params.getSkellamFramesSaved() && params.getSigmaSaved());
    bool needSkellam = !(params.getSkellamFramesSaved() && params.getSlopeSaved() && params.getInterceptSaved());
    if (!needFilter && !needSkellam)
        return;
    unsigned int stacksize = params.getSkellamFrames();
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);

    int maxNumberPoints = 2000;

    T minVal;
    MultiArray<3, T> meanArr, *img = new MultiArray<3, T>(Shape3(w, h, 1)), *lastVal, anscombeTransf;
    MultiArray<2, vigra::acc::AccumulatorChain<T, vigra::acc::Select<vigra::acc::Sum, vigra::acc::Variance>>> skellamArr(Shape2(w, h));
    unsigned int passes;
    int *listPixelCoordinates;
    T *meanValues;
    T *skellamParameters;
    if (needSkellam) {
        meanArr.reshape(Shape3(w, h, 1));
        lastVal = new MultiArray<3, T>(Shape3(w, h, 1));
        skellamArr.reshape(Shape2(w, h));
        passes = skellamArr(0, 0).passesRequired();
        minVal = std::numeric_limits<T>::max();
        listPixelCoordinates = (int*)std::malloc(maxNumberPoints * sizeof(int));
        meanValues = (T*)std::malloc(maxNumberPoints * sizeof(T));
        skellamParameters = (T*)std::malloc(maxNumberPoints * sizeof(T));
    }
    if (needFilter)
        anscombeTransf.reshape(Shape3(w, h, 1));

    transformationFunctor tF(1, 3./8., 0);
    vigra::DImage ps;
    vigra::DImage ps_center;
    if (needFilter) {
        ps.resize(w, h, 0.0);
        ps_center.resize(w, h);
    }

    for(int f = 0; f< stacksize;f++){
        params.readBlock(Shape3(0,0,f),Shape3(w,h,1), *img);
        vigra::combineTwoMultiArrays(srcMultiArrayRange(*img), srcMultiArray(meanArr), destMultiArray(meanArr), std::plus<T>());
        if (needFilter) {
            // calculate PowerSpectrum
            vigra::transformMultiArray(srcMultiArrayRange(*img), destMultiArrayRange(anscombeTransf), tF);
            MultiArrayView <2, T> array = anscombeTransf.bindOuter(0);
            auto arrayView = makeBasicImageView(array);
            subtractBackground(arrayView);
            vigra::FFTWComplexImage fourier(w, h);
            fourierTransform(srcImageRange(array), destImage(fourier));
            vigra::combineTwoImages(srcImageRange(ps),
                                    srcImage(fourier, FFTWSquaredMagnitudeAccessor<double>()),
                                    destImage(ps), Arg1()+Arg2());
        }
        if (needSkellam) {
            if (f > 0) {
                FindMinMax<T> minmax;
                inspectMultiArray(srcMultiArrayRange(*img), minmax);
                if(minVal > minmax.min)
                minVal = minmax.min;
                vigra::combineTwoMultiArrays(srcMultiArrayRange(*lastVal), srcMultiArrayRange(*img), destMultiArrayRange(*lastVal), std::minus<T>());
                auto imgIter = lastVal->begin(), imgEnd = lastVal->end();
                auto skellamIter = skellamArr.begin(), skellamEnd = skellamArr.end();
                for (; imgIter != imgEnd && skellamIter != skellamEnd; ++imgIter, ++skellamIter) {
                    for (int n = 1; n <= passes; ++n) {
                        skellamIter->updatePassN(*imgIter, n);
                    }
                }
            }
            auto *tmp = lastVal;
            lastVal = img;
            img = tmp;
        }
    }
    if (needFilter) {
        moveDCToCenter(srcImageRange(ps), destImage(ps_center));
        vigra::transformImage(srcImageRange(ps_center), destImage(ps),
                            Arg1() / Param(stacksize));
        constructWienerFilter(params, ps);
    }
    if (needSkellam) {
        vigra::transformMultiArray(srcMultiArrayRange(meanArr), destMultiArrayRange(meanArr), [&stacksize](T p){return p / stacksize;});

        FindMinMax<T> minmax;
        inspectMultiArray(srcMultiArrayRange(meanArr), minmax);
        if(params.getVerbose()){
        std::cout<<"min: "<<minmax.min<<" max: "<<minmax.max<<std::endl;}

        std::vector<int> iter(w * h); //contains information permutation from sort
        linearSequence(iter.begin(), iter.end());
        indexSort(meanArr.begin(), meanArr.end(), iter.begin());

        int intervalCounter = 0;
        T intervalDistance = (minmax.max - minmax.min)/ maxNumberPoints;

        for(int i = 0; i< w * h && intervalCounter < maxNumberPoints; i++) {
            if(meanArr[iter[i]] > minmax.min + intervalDistance * intervalCounter) {
                listPixelCoordinates[intervalCounter] = iter[i];
                meanValues[intervalCounter] = meanArr[iter[i]];
                skellamParameters[intervalCounter] = 0.5 * (vigra::acc::get<vigra::acc::Sum>(skellamArr[iter[i]]) / stacksize + vigra::acc::get<vigra::acc::Variance>(skellamArr[iter[i]]));
                ++intervalCounter;
            }
        }
        findBestFit(params, meanValues, skellamParameters, intervalCounter);

        if(params.getIntercept() > 0)
            params.setIntercept(std::min(minVal, params.getIntercept()));
        else
            params.setIntercept(minVal);
        delete lastVal;
    }
    delete img;
}

template <class T>
void applyMask(BasicImage<T>& img, MultiArrayView<2, T>& mask, int frnr){
	int w = mask.shape(0);
	int h = mask.shape(1);

	std::ofstream imgAfterMaskFile, imgBeforeMaskFile;
	bool writeMatrices = frnr == 1;
	if(writeMatrices){

		char imgAfterMask[1000];
		char imgBeforeMask[1000];

		sprintf(imgAfterMask, "/home/herrmannsdoerfer/tmpOutput/filterdImgAfterMask%d.txt", frnr);
		sprintf(imgBeforeMask, "/home/herrmannsdoerfer/tmpOutput/filterdImgBeforeMask%d.txt", frnr);

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
// STORM DATA PROCESSING
//--------------------------------------------------------------------------

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

template <class T>
void getPoissonLabelsArray(const DataParams &params, MultiArray<3, T> &labels) {
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);
    unsigned int tChunkSize = params.getTChunkSize();
    unsigned int xyChunkSize = params.getXYChunkSize();
    unsigned int xChunks = std::ceil(w / (float)xyChunkSize);
    unsigned int yChunks = std::ceil(h / (float)xyChunkSize);
    labels.reshape(Shape3(w, h, tChunkSize));
    for (int x = 0, l = 0; x < xChunks; ++x) {
        for (int y = 0; y < yChunks; ++y, ++l) {
            vigra::Shape3 index(std::min((x + 1) * xyChunkSize, w), std::min((y + 1) * xyChunkSize, h), tChunkSize);
            auto roi = labels.subarray(vigra::Shape3(x * xyChunkSize, y * xyChunkSize, 0), index);
            vigra::initMultiArray(destMultiArrayRange(roi), l);
        }
    }
}

template <class T, class L>
void getPoissonMeansForChunk(const DataParams &params, const MultiArrayView<3, L> &labels, const MultiArrayView<3, T> &img, MultiArrayView<2, T> &regionMeans) {
    vigra::acc::AccumulatorChainArray<typename CoupledIteratorType<3, T, L>::type::value_type, vigra::acc::Select<vigra::acc::DataArg<1>, vigra::acc::LabelArg<2>, vigra::acc::StandardQuantiles<vigra::acc::AutoRangeHistogram<0>>>> accChain;
    auto iter = vigra::createCoupledIterator(img, labels);
    auto iterEnd = iter.getEndIterator();
    vigra::acc::extractFeatures(iter, iterEnd, accChain);
    for (int x = 0, i = 0; x < regionMeans.shape()[0]; ++x) {
        for (int y = 0; y < regionMeans.shape()[1]; ++y, ++i) {
            regionMeans(x, y) = vigra::acc::get<vigra::acc::StandardQuantiles<vigra::acc::AutoRangeHistogram<0>>>(accChain, i)[3] - 0.3333;
        }
    }
}

template <class T, class F>
void processChunk(const DataParams &params, MultiArray<3, T> &srcImage,
                  MultiArrayView<3, T> &poissonMeans, F tF, int &currframe, int middleChunk,
                  std::vector<std::set<Coord<T>>>& maxima_coords) {
    unsigned int middleChunkFrame = middleChunk * params.getTChunkSize();
    #pragma omp parallel for schedule(runtime) shared(srcImage, poissonMeans, maxima_coords)
    for (int f = 0; f < srcImage.shape()[2]; ++f) {
        auto currSrc = srcImage.bindOuter(f);
        auto currPoisson = poissonMeans.bindOuter(middleChunkFrame + f);
        vigra::combineTwoMultiArrays(srcMultiArrayRange(currSrc), srcMultiArray(currPoisson), destMultiArray(currSrc),
                                     [&tF](T srcPixel, T poissonPixel){return tF(srcPixel) - tF(poissonPixel);});
        wienerStormSingleFrame(params, currSrc, maxima_coords[currframe + f], currframe + f);
    }
    currframe += srcImage.shape()[2];
}

template <class T, class L>
void readChunk(const DataParams &params, MultiArray<3, T>** srcImage,
               MultiArrayView<3, T> &poissonMeansRaw, MultiArrayView<3, T> &poissonMeans,
               const MultiArrayView<3, L> &poissonLabels, int chunk) {
    unsigned int xChunks = std::ceil(params.shape(0) / (float)params.getXYChunkSize());
    unsigned int yChunks = std::ceil(params.shape(1) / (float)params.getXYChunkSize());
    unsigned int chunksInMemory = params.getChunksInMemory();
    unsigned int middleChunk = std::floor(0.5 * chunksInMemory);
    auto *tmp = srcImage[0];
    for (int i = 0; i < middleChunk; ++i) {
        srcImage[i] = srcImage[i + 1];
    }
    srcImage[middleChunk] = tmp;
    params.readBlock(Shape3(0, 0, chunk * params.getTChunkSize()), tmp->shape(), *tmp);
    vigra::transformMultiArray(srcMultiArrayRange(*tmp), destMultiArrayRange(*tmp), [&params](T p){return (p - params.getIntercept()) / params.getSlope();});
    for (int z = 0; z < chunksInMemory - 1; ++z) {
        for (int x = 0; x < xChunks; ++x) {
            for (int y = 0; y < yChunks; ++y) {
                poissonMeansRaw(x, y, z) = poissonMeansRaw(x, y, z + 1);
            }
        }
    }
    auto currRawMean = poissonMeansRaw.bindOuter(chunksInMemory - 1);
    getPoissonMeansForChunk(params, poissonLabels, *tmp, currRawMean);
    vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(poissonMeansRaw), destMultiArrayRange(poissonMeans), vigra::BSpline<3>());
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
void wienerStorm(DataParams &params, std::vector<std::set<Coord<T> > >& maxima_coords) {

    unsigned int stacksize = params.shape(2);
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);
    unsigned int i_stride=1;
    int i_beg=0, i_end=stacksize;
    unsigned int tChunks = std::ceil((i_end - i_beg) / (float)params.getTChunkSize());
    unsigned int xChunks = std::ceil(w / (float)params.getXYChunkSize());
    unsigned int yChunks = std::ceil(h / (float)params.getXYChunkSize());
    unsigned int chunksInMemory = params.getChunksInMemory();
    unsigned int lastTChunkIndex = params.getTChunkSize() - 1;
    unsigned int middleChunk = std::floor(0.5 * chunksInMemory);
    unsigned int middleChunkFrame = middleChunk * params.getTChunkSize();
    /*if(!params.getFrameRange().empty()) {
     h elper::rangeSplit(params*.getFrameRange(), i_beg, i_end, i_stride);
        if(i_beg < 0) i_end = stacksize+i_beg; // allow counting backwards from the end
        if(i_end < 0) i_end = stacksize+i_end; // allow counting backwards from the end
        if(params.verbose) std::cout << "processing frames [" << i_beg << ":"
            << i_end << ":" << i_stride << "]" << std::endl;
    }*/

    // TODO: Precondition: res must have size (params.getFactor()*(w-1)+1, params.getFactor()*(h-1)+1)
    // filter must have the size of input

    #ifndef STORM_QT // silence stdout
    std::cout << "Finding the maximum spots in the images..." << std::endl;
    #endif // STORM_QT
    helper::progress(-1,-1); // reset progress
#ifdef OPENMP_FOUND
    omp_set_schedule(omp_sched_dynamic, omp_get_num_threads / params.getTChunkSize());
#endif

    transformationFunctor tF(1, 3./8,0);

    MultiArray<3, T>** srcImage = (MultiArray<3, T>**)std::malloc(chunksInMemory * sizeof(MultiArray<3, T>*));
    for (int i = 0; i < chunksInMemory; ++i) {
        srcImage[i] = new MultiArray<3, T>(Shape3(w, h, params.getTChunkSize()));
    }
    MultiArray<3, int> poissonLabels;
    MultiArray<3, T> poissonMeansRaw(Shape3(xChunks, yChunks, chunksInMemory));
    MultiArray<3, T> poissonMeans(Shape3(w, h, params.getTChunkSize() * chunksInMemory));
    getPoissonLabelsArray(params, poissonLabels);

    int currframe = 0, chunk;

    for (chunk = 0; chunk < chunksInMemory; ++chunk) {
        params.readBlock(Shape3(0, 0, chunk * params.getTChunkSize()), srcImage[chunk]->shape(), *srcImage[chunk]);
        auto currRawMean = poissonMeansRaw.bindOuter(chunk);
        vigra::transformMultiArray(srcMultiArrayRange(*srcImage[chunk]), destMultiArrayRange(*srcImage[chunk]), [&params](T p){return (p - params.getIntercept()) / params.getSlope();});
        getPoissonMeansForChunk(params, poissonLabels, *srcImage[chunk], currRawMean);

    }
    vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(poissonMeansRaw), destMultiArrayRange(poissonMeans), vigra::BSpline<3>());
    for (int c = 0; c <= middleChunk; ++c) {
        processChunk(params, *srcImage[c], poissonMeans, tF, currframe, c, maxima_coords);
        helper::progress(currframe, i_end);
    }
    if (chunksInMemory % 2) {
        for (int c = 0; c < middleChunk; ++c) {
            delete srcImage[c];
            srcImage[c] = srcImage[c + middleChunk];
        }
        srcImage[middleChunk] = srcImage[chunksInMemory - 1];
    } else {
        for (int c = 0; c < middleChunk - 1; ++c) {
            srcImage[c] = srcImage[c + middleChunk - 1];
        }
        srcImage[middleChunk - 1] = srcImage[chunksInMemory - 2];
        srcImage[middleChunk] = srcImage[chunksInMemory - 1];
    }

    int lastChunkSize = stacksize % params.getTChunkSize();
    for (; chunk < (lastChunkSize ? tChunks - 1 : tChunks); ++chunk) {
        readChunk(params, srcImage, poissonMeansRaw, poissonMeans, poissonLabels, chunk);
        processChunk(params, *srcImage[0], poissonMeans, tF, currframe, middleChunk, maxima_coords);
        helper::progress(currframe, i_end);
    }
    if (lastChunkSize) {
        srcImage[0]->reshape(Shape3(w, h, lastChunkSize));
        Shape3 labelsShape = poissonLabels.shape();
        labelsShape[2] = lastChunkSize;
        auto lastPoissonLabelsView = poissonLabels.subarray(Shape3(0, 0, 0), labelsShape);
        readChunk(params, srcImage, poissonMeansRaw, poissonMeans, lastPoissonLabelsView, chunk);
        processChunk(params, *srcImage[0], poissonMeans, tF, currframe, middleChunk, maxima_coords);
        helper::progress(currframe, i_end);
    }
    delete srcImage[0];
    for (int c = middleChunk + 1; c < chunksInMemory; ++c) {
        int cIndex = c - middleChunk;
        processChunk(params, *srcImage[cIndex], poissonMeans, tF, currframe, c, maxima_coords);
        helper::progress(currframe, i_end);
        delete srcImage[cIndex];
    }
    std::free(srcImage);
}

template <class T>
void wienerStormSingleFrame(const DataParams &params, const MultiArrayView<2, T>& in, std::set<Coord<T> >& maxima_coords, int framenumber)
{
    int w = params.shape(0); // width
    int h = params.shape(1); // height

    //std::cout<<"frame: "<<framenumber<<std::endl;

    BasicImage<T> filtered(w,h);

    int factor = params.getFactor();
    int mylen = params.getRoilen();

    const int mylen2 = mylen/2;
    unsigned int w_roi = factor*(mylen-1)+1;
    unsigned int h_roi = factor*(mylen-1)+1;
    BasicImage<T> im_xxl(w_roi, h_roi);

    BasicImageView<T> input = makeBasicImageView(in);  // access data as BasicImage

    //fft, filter with Wiener filter in frequency domain, inverse fft, take real part
    BasicImageView<T> filteredView(filtered.data(), filtered.size());

    BasicImage<T> unfiltered(w,h);

    vigra::copyImage(srcImageRange(input), destImage(unfiltered));

    std::ofstream beforeimg, afterimg;
    char beforefilter[1000], afterfilter[1000];
    sprintf(beforefilter, "/home/herrmannsdoerfer/tmpOutput/frameData/beforefilter%d.txt", framenumber);
    sprintf(afterfilter, "/home/herrmannsdoerfer/tmpOutput/frameData/afterfilter%d.txt", framenumber);


//     beforeimg.open (beforefilter);
//     for (int i = 0; i < w; i++) {
//         for( int j = 0; j< h; j++) {
//             beforeimg <<i<<" "<<j<<" "<< input(i,j)<<std::endl;
//         }
//     }
//     beforeimg.close();
    gaussianSmoothing(srcImageRange(input), destImage(filteredView), std::pow(std::pow(params.getSigma(), 2)-std::pow(0.85, 2), .5));

//     afterimg.open (afterfilter);
//     for (int i = 0; i < w; i++) {
//         for( int j = 0; j< h; j++) {
//             afterimg <<i<<"  "<<j<<"  "<< filteredView(i,j)<<std::endl;
//         }
//     }
//     afterimg.close();

    vigra::FindMinMax<T> filteredMinMax;
    inspectImage(srcImageRange(filtered), filteredMinMax);
    MultiArray<2, T> mask(Shape2(w, h));
    getMask(params, filtered, framenumber, mask);
    std::set<Coord<T> > maxima_candidates_vect;  // we use a set for the coordinates to automatically squeeze duplicates
                                                 // (from overlapping ROIs)
    VectorPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_candidates(maxima_candidates_vect, filtered.upperLeft(), 1, mask);
    vigra::localMaxima(srcImageRange(filtered), destImage(filtered, maxima_candidates), vigra::LocalMinmaxOptions().neighborhood(4));

    VectorPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_acc(maxima_coords, im_xxl.upperLeft(), factor, mask);
    //upscale filtered image regions with spline interpolation
    std::set<Coord<float> >::iterator it2;


//	vigra::exportImage(srcImageRange(filtered), "/home/herrmannsdoerfer/master/workspace/output/filtered.png");
//	vigra::exportImage(srcImageRange(unfiltered), "/home/herrmannsdoerfer/master/workspace/output/unfiltered.png");

    for(it2=maxima_candidates_vect.begin(); it2 != maxima_candidates_vect.end(); it2++) {
            Coord<float> c = *it2;
            if(unfiltered(c.x,c.y)<0.5 or mask(c.x,c.y) == 0.0) { // skip very low signals with SNR lower 3 or maxima not covered by the mask
                //std::cout<<"value skipped: "<<unfiltered(c.x,c.y)<<std::endl;
                continue;
            }
            //std::cout<<"value not skipped: "<<unfiltered(c.x,c.y)<<std::endl;
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
                               destIter(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), maxima_acc), vigra::LocalMinmaxOptions().neighborhood(4));
    }
    determineAsymmetry(srcImageRange(unfiltered), maxima_coords, factor);
    determineSNR(srcImageRange(unfiltered), maxima_coords, factor);
}
