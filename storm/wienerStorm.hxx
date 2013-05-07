/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*      Copyright 2010-2013 by Joachim Schleicher, Ullrich Koethe,	    */
/*					Frank Herrmannsdoerfer and Ilia Kats				*/
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

#ifndef WIENERSTORM_HXX
#define WIENERSTORM_HXX

#include <vigra/convolution.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/splineimageview.hxx>
#include <vigra/multi_fft.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/multi_resize.hxx>
#include <vigra/multi_impex.hxx>

#include <set>
#include <fstream>
#include <iomanip>
#include <vector>
#include <valarray>
#include <limits>
#include <atomic>
#include <thread>
#include <mutex>
#include <algorithm>
#include <numeric>
#ifdef OPENMP_FOUND
    #include <omp.h>
#endif //OPENMP_FOUND
#include "util.h"
#include "dataparams.h"

#define R_NO_REMAP
#include <Rembedded.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>

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

extern std::mutex wienerStorm_R_mutex;


enum WienerStormStage{CameraParameters, ParameterCheck, PSFWidth, Localization};
class ProgressFunctor
{
public:
    ProgressFunctor() : m_abort(false), m_finished(false) {};
    virtual ~ProgressFunctor(){};
    virtual void setStage(WienerStormStage) = 0;
    virtual void setStackSize(int) = 0;
    virtual void frameFinished(int) = 0;
    virtual void abort() {m_abort = true;};
    bool getAbort() {return m_abort;};
    void setFinished() {m_finished = true;};
    bool isFinished() {return m_finished;};

protected:
    std::atomic<bool> m_abort;
    std::atomic<bool> m_finished; // to indicate when abortion is complete
};

//--------------------------------------------------------------------------
// helper classes and functions
//--------------------------------------------------------------------------

/**
 * BSpline coefficients but no prefiltering
 */

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
class SetPushAccessor{
    public:
        typedef typename T::value_type value_type;
        SetPushAccessor(std::set<T>& arr, ITERATOR it_start, int factor, const MultiArrayView<2, S> &mask)
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

template <class T, class ITERATOR>
class VectorPushAccessor{
public:
    typedef typename T::value_type value_type;
    VectorPushAccessor(std::vector<T>& arr, ITERATOR it_start)
    : m_arr(arr), m_it_start(it_start) {}

    T const &   operator() (ITERATOR const &i) const {
        return NumericTraits<T>::zero();
    }
    template<class V>
    void    set (V const & /*value*/, ITERATOR const &i) {
        int x = i.x-m_it_start.x;
        int y = i.y-m_it_start.y;
        m_arr.push_back(T(x, y, *i));
    }

private:
    std::vector<T>& m_arr;
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
int saveCoordsFile(const DataParams &params, std::ofstream &cfile, const std::vector<std::set<C> >& coords) {
    int numSpots = 0;
    std::set<Coord<float> >::iterator it2;
    cfile << params.shape(0) << " " << params.shape(1) << " " << params.shape(2) << std::endl;
    cfile << std::fixed; // fixed instead of scientific format
    for(unsigned int j = 0; j < coords.size(); j++) {
        for(it2=coords[j].begin(); it2 != coords[j].end(); it2++) {
            numSpots++;
            const Coord<float>& c = *it2;
            cfile << std::setprecision(4) << (float)c.x/params.getFactor() * params.getPixelSize() << " " << (float)c.y/params.getFactor() * params.getPixelSize() << " "
                << j << " " << std::setprecision(1) << c.val << " " << std::setprecision(3) << c.asymmetry << std::endl;
        }
    }
    return numSpots;
}

/**
 * finds the value, so that the given percentage of pixels is above / below that value.
 */
template<class T, class Iterator> //can't use Iterator::value_type for MultiArray
void findMinMaxPercentile(Iterator begin, Iterator end, double minPerc, double &minVal, double maxPerc, double &maxVal) {
    std::vector<T> v(begin, end);
    std::sort(v.begin(),v.end());
    minVal=v[(int)(v.size()*minPerc)];
    maxVal=v[(int)(v.size()*maxPerc)];
}
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
        const DataParams &params) {
    determineAsymmetry(s.first, s.second, s.third, coords, params);
}

template <class SrcIterator, class SrcAccessor, class T>
void determineAsymmetry(SrcIterator srcUpperLeft,
        SrcIterator srcLowerRight,
        SrcAccessor acc,
        std::set<Coord<T> >& coords,
        const DataParams &params) {
    int factor = params.getFactor();
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
        if (params.getDoAsymmetryCheck()) {
            if(ev1/ev2<params.getAsymmetryThreshold() and ev2/ev1 < params.getAsymmetryThreshold()) {
                Coord<float> cc (c.x, c.y, c.val, ev1/ev2);
                newcoords.insert(cc); // copy for now. Hack hack hack...
            }
        }
        else {
            Coord<float> cc (c.x, c.y, c.val, ev1/ev2);
            newcoords.insert(cc); // copy for now. Hack hack hack..
        }

    }
    coords=newcoords;
}

template <class SrcIterator, class SrcAccessor, class T>
inline void determineSNR(triple<SrcIterator, SrcIterator, SrcAccessor> s,
        std::set<Coord<T> >& coords,
        const int factor ) {
    determineSNR(s.first, s.second, s.third, coords, factor );
}

/*!
Determines signal-to-noise ratio based on the assumption that the variance of the background is equal to one.
*/
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

/*!
Functor to apply Anscombe transformation
*/
class transformationFunctor {
public:
	transformationFunctor(float a, float intercept, float minInt = 0): a(a), intercept(intercept), C(-2/a*std::sqrt(a*minInt+ intercept)){}
	float operator()(float poissInt) const {return 2/a*std::sqrt(a*poissInt + intercept)+ C;}

private:
	float a;
	float intercept;
	float C;
};

/*!
Fit a straight line to the selected points using R
*/
template <class T>
void fitSkellamPoints(DataParams &params,T meanValues[],T skellamParameters[],int numberPoints){
    int nbins = 10;
    std::cout<<std::endl;
    for (int i =0; i< numberPoints; ++i){
        std::cout<<meanValues[i]<<", ";
    }
    std::cout<<std::endl;
    for (int i =0; i< numberPoints; ++i){
        std::cout<<skellamParameters[i]<<", ";
    }
    std::cout<<std::endl;
    wienerStorm_R_mutex.lock();
    SEXP fun, t, tmp;
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
    SETCAR(t, Rf_install("fit.skellam.points"));
    t = CDR(t);
    SETCAR(t, tmp);
    t = CDR(t);
    SETCAR(t, bins);
    PROTECT(tmp = Rf_eval(fun, R_GlobalEnv));

    if (tmp == R_NilValue) {
        std::cerr << "Fit could not be performed, seting slope to one and intercept to zero..." << std::endl;
        params.setSlope(1);
        params.setIntercept(0);
    } else {
        double *coefs = REAL(tmp);
        params.setSlope(coefs[1]);
        params.setIntercept(-coefs[0] / coefs[1]);
    }
    UNPROTECT(4);
    R_gc();
    wienerStorm_R_mutex.unlock();
}

/*!
Calculates a mask that contains the information whether or not a pixel belongs to background or is signal for eachh pixel of the current frame,
 based on cumulative distribution. The asumption is that the background has a mean value of zeros and a variance of one.
 To get rid of noise just consider larger spots
 */
template <class T>
void getMask(const DataParams &params, const BasicImage<T>& array, int framenumber, MultiArray<2,T>& mask){
	//double cdf = params.getMaskThreshold();
    double cdf = qnorm(params.getAlpha(), 0, 1, 0, 0);
//     vigra::exportImage(srcImageRange(array),"/home/herrmannsdoerfer/tmpOutput/array.tif");
    vigra::transformImage(srcImageRange(array), destImage(mask), [&cdf](T p) {return p >= cdf ? 1 : 0;});
//     char name[1000];
//     sprintf(name, "/home/herrmannsdoerfer/tmpOutput/frameData/maskBeforeCC%d.tif", framenumber);
//     vigra::exportImage(srcImageRange(mask), name);
    vigra::IImage labels(array.width(), array.height());
    unsigned int nbrCC = vigra::labelImageWithBackground(srcImageRange(mask), destImage(labels), false, 0);
    std::valarray<int> bins(0, nbrCC + 1);
    auto fun = [&bins](int32_t p){++bins[p];};
    vigra::inspectImage(srcImageRange(labels), fun);
//     vigra::transformImage(srcImageRange(labels), destImage(mask), [&params, &bins] (T p) -> T {if(!p || bins[p] < std::max(3.,0.25*3.14 * std::pow(params.getSigma(), 2))) return 0; else return 1;});
    vigra::transformImage(srcImageRange(labels), destImage(mask), [&params, &bins] (T p) -> T {if(!p || bins[p] < 3) return 0; else return 1;});
//     char name2[1000];
//     sprintf(name2, "/home/herrmannsdoerfer/tmpOutput/frameData/maskAfterCC%d.tif", framenumber);
//     vigra::exportImage(srcImageRange(mask), name2);
}

//*! Estimates values for camera gain and offset using a mean-variance plot*/
template <class T>
void estimateCameraParameters(DataParams &params, ProgressFunctor &progressFunc) {
    bool needSkellam = !((params.getIgnoreSkellamFramesSaved() or params.getSkellamFramesSaved()) && params.getSlopeSaved() && params.getInterceptSaved());
    if (params.getIgnoreSkellamFramesSaved()) {
        std::cout<<"Values from GUI dialog:"<<std::endl;
        std::cout<<"Gain: "<<params.getSlope()<<" Offset: "<<params.getIntercept()<<std::endl;
        return;
    }
    if (!needSkellam) {
        std::cout<<"Values from settings-file:"<<std::endl;
        std::cout<<"Gain: "<<params.getSlope()<<" Offset: "<<params.getIntercept()<<std::endl;
        return;
    }
    unsigned int stacksize = params.getSkellamFrames();
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);

    progressFunc.setStage(CameraParameters);
    progressFunc.setStackSize(stacksize);
    int maxNumberPoints = 2000;

    T minVal = std::numeric_limits<T>::max(), maxVal = 0;
    MultiArray<3, T> meanArr(Shape3(w, h, 1)), *img = new MultiArray<3, T>(Shape3(w, h, 1)), *lastVal = new MultiArray<3, T>(Shape3(w, h, 1));
    MultiArray<2, vigra::acc::AccumulatorChain<T, vigra::acc::Select<vigra::acc::Sum, vigra::acc::Variance>>> skellamArr(Shape2(w, h));
    unsigned int passes = skellamArr(0, 0).passesRequired();
    int *listPixelCoordinates = (int*)std::malloc(maxNumberPoints * sizeof(int));
    T *meanValues = (T*)std::malloc(maxNumberPoints * sizeof(T));
    T *skellamParameters = (T*)std::malloc(maxNumberPoints * sizeof(T));

    for(int f = 0; f< stacksize; ++f) {
        if (progressFunc.getAbort())
            return;
        params.readBlock(Shape3(0,0,f),Shape3(w,h,1), *img);
        vigra::combineTwoMultiArrays(srcMultiArrayRange(*img), srcMultiArray(meanArr), destMultiArray(meanArr), std::plus<T>());
        if (f > 0) {
            FindMinMax<T> minmax;
            inspectMultiArray(srcMultiArrayRange(*img), minmax);
            if(minVal > minmax.min)
                minVal = minmax.min;
            if(maxVal < minmax.max)
                maxVal = minmax.max;
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
        progressFunc.frameFinished(f);
    }
    std::cout<<std::endl;
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
    fitSkellamPoints(params, meanValues, skellamParameters, intervalCounter);
    if(params.getIntercept() > 0)
        params.setIntercept(std::min(minVal, params.getIntercept()));
    else
        params.setIntercept(minVal);
    if(params.getSlope() <= 0)
        params.setSlope(1);
    std::cout<<"slope: "<< params.getSlope()<<" x0: "<<params.getIntercept() << std::endl;
    std::cout<<"min: "<<minVal<<" max: "<<maxVal<<std::endl;
//     if(params.getSlope()>(maxVal-minVal)/10.) //if the data has no beads the slope might get very high values, if so it is set to a quarter of the range.
//         params.setSlope((maxVal-minVal)/20.);
    delete lastVal;
    delete img;
    std::cout<<"Estimated values:"<<std::endl;
    std::cout<<"Gain: "<<params.getSlope()<<" Offset: "<<params.getIntercept()<<std::endl;
}

/*!
Take median of each chunk and store it for later interpolation to full image resolution
*/
template <class T>
void getPoissonMeansForChunk(const DataParams &params, int tChunkSize,const MultiArrayView<3, T> &img, MultiArrayView<2, T> &regionMeans) {
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);
    unsigned int xyChunkSize = params.getXYChunkSize();
    unsigned int xChunks = std::ceil(w / (float)xyChunkSize);
    unsigned int yChunks = std::ceil(h / (float)xyChunkSize);
    for (int x = 0, n = 0; x < xChunks; ++x) {
        for (int y = 0; y < yChunks; ++y, ++n) {
            vigra::Shape3 index(std::min((x + 1) * xyChunkSize, w), std::min((y + 1) * xyChunkSize, h), tChunkSize);
            auto nthroi = img.subarray(vigra::Shape3(x * xyChunkSize, y * xyChunkSize, 0), index);
            std::vector<T> vec(nthroi.begin(), nthroi.end());
            std::nth_element(vec.begin(), vec.begin() + vec.size()/2, vec.end());
            regionMeans(x, y) = vec[vec.size()/2];
        }
    }
}

template <class T, class Func>
void processChunk(const DataParams &params, MultiArray<3, T> &srcImage,
                  MultiArrayView<3, T> &poissonMeans, int &currframe, int middleChunk,
                  Func& functor, ProgressFunctor &progressFunc) {
    unsigned int middleChunkFrame = middleChunk * params.getTChunkSize();
    #pragma omp parallel for schedule(runtime) shared(srcImage, poissonMeans, functor, progressFunc)
    for (int f = 0; f < srcImage.shape()[2]; ++f) {
        auto currSrc = srcImage.bindOuter(f);
        auto currPoisson = poissonMeans.bindOuter(middleChunkFrame + f);
        vigra::combineTwoMultiArrays(srcMultiArrayRange(currSrc), srcMultiArray(currPoisson), destMultiArray(currSrc),
                                     [](T srcPixel, T poissonPixel) -> T {T val = srcPixel - poissonPixel; return (std::isnan(val)) ? 0 : val;});
//         char name[1000];
//         sprintf(name, "/home/herrmannsdoerfer/tmpOutput/frameData/frame%d.tif", middleChunkFrame+f);
//         vigra::exportImage(srcImageRange(currSrc), name);
        functor(params, currSrc, currframe + f);
        progressFunc.frameFinished(currframe + f);
    }
    currframe += srcImage.shape()[2];
}

template <class T, class F>
void readChunk(const DataParams &params, MultiArray<3, T>** srcImage,
               MultiArrayView<3, T> &poissonMeansRaw, MultiArrayView<3, T> &poissonMeans,
               int lastChunkSize, int chunk, F &tF) {
    unsigned int xChunks = std::ceil(params.shape(0) / (float)params.getXYChunkSize());
    unsigned int yChunks = std::ceil(params.shape(1) / (float)params.getXYChunkSize());
    unsigned int chunksInMemory = params.getChunksInMemory();
    unsigned int middleChunk = std::floor(0.5 * chunksInMemory);
    auto *tmp = srcImage[0];
    auto *tmp2 = srcImage[0];
    for (int i = 0; i < middleChunk; ++i) {
        srcImage[i] = srcImage[i + 1];
    }
    srcImage[middleChunk] = tmp;
    params.readBlock(Shape3(0, 0, chunk * params.getTChunkSize()), tmp->shape(), *tmp);
//     char name[1000];
//     sprintf(name, "/home/herrmannsdoerfer/tmpOutput/NotAtAllCorrectedFrame%d.tif", chunk*lastChunkSize);
//     vigra::exportImage(srcImageRange(tmp->bindOuter(0)),name);
    vigra::transformMultiArray(srcMultiArrayRange(*tmp), destMultiArrayRange(*tmp), [&params](T p){return (p - params.getIntercept()) / params.getSlope();});
    char name2[1000];
//     sprintf(name, "/home/herrmannsdoerfer/tmpOutput/PoissonCorrectedFrame%d.tif", chunk*lastChunkSize);
//     vigra::exportImage(srcImageRange(tmp->bindOuter(0)),name);
    vigra::transformMultiArray(srcMultiArrayRange(*tmp), destMultiArrayRange(*tmp), [&tF](T p){return tF(p);});
//     char name3[1000];
//     sprintf(name, "/home/herrmannsdoerfer/tmpOutput/AnscombeCorrectedFrame%d.tif", chunk*lastChunkSize);
//     vigra::exportImage(srcImageRange(tmp->bindOuter(0)),name);
    //vigra::exportImage(srcImageRange(tmp->bindOuter(0)),"/home/herrmannsdoerfer/tmpOutput/skellamCorrectedFrame.tif");
    for (int z = 0; z < chunksInMemory - 1; ++z) {
        for (int x = 0; x < xChunks; ++x) {
            for (int y = 0; y < yChunks; ++y) {
                poissonMeansRaw(x, y, z) = poissonMeansRaw(x, y, z + 1);
            }
        }
    }
    auto currRawMean = poissonMeansRaw.bindOuter(chunksInMemory - 1);
    getPoissonMeansForChunk(params, lastChunkSize, *tmp, currRawMean);
    vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(poissonMeansRaw), destMultiArrayRange(poissonMeans), vigra::BSpline<3>());
}

/*!
This function organizes everything, chunkwise loading the data, handling first and last chunk
*/
template <class T, class Func>
void processStack(const DataParams &params, Func& functor, ProgressFunctor &progressFunc, unsigned int stacksize = 0) {

    if (!stacksize)
        stacksize = params.shape(2);
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
     *     h elper::rangeSplit(params*.getFrameRange(), i_beg, i_end, i_stride);
     *        if(i_beg < 0) i_end = stacksize+i_beg; // allow counting backwards from the end
     *        if(i_end < 0) i_end = stacksize+i_end; // allow counting backwards from the end
     *        if(params.verbose) std::cout << "processing frames [" << i_beg << ":"
     *            << i_end << ":" << i_stride << "]" << std::endl;
}*/

    // TODO: Precondition: res must have size (params.getFactor()*(w-1)+1, params.getFactor()*(h-1)+1)
    // filter must have the size of input
    progressFunc.setStackSize(stacksize);

    #ifdef OPENMP_FOUND
    omp_set_schedule(omp_sched_dynamic, omp_get_num_threads() / params.getTChunkSize());
    #endif

    transformationFunctor tF(1, 3./8,0);

    MultiArray<3, T>** srcImage = (MultiArray<3, T>**)std::malloc(chunksInMemory * sizeof(MultiArray<3, T>*));
    for (int i = 0; i < chunksInMemory; ++i) {
        srcImage[i] = new MultiArray<3, T>(Shape3(w, h, params.getTChunkSize()));
    }
    MultiArray<3, T> poissonMeansRaw(Shape3(xChunks, yChunks, chunksInMemory));
    MultiArray<3, T> poissonMeans(Shape3(w, h, params.getTChunkSize() * chunksInMemory));

    int currframe = 0, chunk;

    for (chunk = 0; chunk < chunksInMemory; ++chunk) {
        params.readBlock(Shape3(0, 0, chunk * params.getTChunkSize()), srcImage[chunk]->shape(), *srcImage[chunk]);
        auto currRawMean = poissonMeansRaw.bindOuter(chunk);
        vigra::transformMultiArray(srcMultiArrayRange(*srcImage[chunk]), destMultiArrayRange(*srcImage[chunk]), [&params, &tF](T p){return tF((p - params.getIntercept()) / params.getSlope());});
        getPoissonMeansForChunk(params, params.getTChunkSize(), *srcImage[chunk], currRawMean);

    }
    vigra::resizeMultiArraySplineInterpolation(srcMultiArrayRange(poissonMeansRaw), destMultiArrayRange(poissonMeans), vigra::BSpline<3>());
    for (int c = 0; c <= middleChunk; ++c) {
        if (progressFunc.getAbort())
            return;
        processChunk(params, *srcImage[c], poissonMeans, currframe, c, functor, progressFunc);
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
        if (progressFunc.getAbort())
            return;
        readChunk(params, srcImage, poissonMeansRaw, poissonMeans, params.getTChunkSize(), chunk, tF);
        processChunk(params, *srcImage[0], poissonMeans, currframe, middleChunk, functor, progressFunc);
    }
    if (lastChunkSize) {
        if (progressFunc.getAbort())
            return;
        srcImage[0]->reshape(Shape3(w, h, lastChunkSize));
        readChunk(params, srcImage, poissonMeansRaw, poissonMeans, lastChunkSize, chunk, tF);
        processChunk(params, *srcImage[0], poissonMeans, currframe, middleChunk, functor, progressFunc);
    }
    delete srcImage[0];
    for (int c = middleChunk + 1; c < chunksInMemory; ++c) {
        if (progressFunc.getAbort())
            return;
        int cIndex = c - middleChunk;
        processChunk(params, *srcImage[cIndex], poissonMeans, currframe, c, functor, progressFunc);
        delete srcImage[cIndex];
    }
    std::free(srcImage);
}

/*!
adds new powerspecra for each frame to the previous ones
*/
template <class T, class S>
void accumulatePowerSpectrum(const DataParams &params, const FFTWPlan<2, S> &fplan, MultiArrayView<2, T>& in, MultiArrayView<2, double> &ps, int roiwidth, int nbrRoisPerFrame, int &rois) {

    int w = params.shape(0), h = params.shape(1);
    int roiwidth2 = roiwidth / 2;
    BasicImageView<T> input = makeBasicImageView(in);
    std::vector<Coord<T> > maxima;
    typename BasicImageView<T>::traverser it = input.upperLeft();
    VectorPushAccessor<Coord<T>, typename BasicImageView<T>::traverser> maxima_acc(maxima, it);
    vigra::localMaxima(srcImageRange(input), destImage(input, maxima_acc));
    if (maxima.empty()) {
        return;
    }

    T min = std::numeric_limits<T>::max(), max = std::numeric_limits<T>::min();
    int nbrMaxima = 0;
    for (auto i = maxima.begin(); i != maxima.end(); ++i) {
        if (i->x < roiwidth2 || i->x > w - roiwidth2  || i->y < roiwidth2 || i->y > h - roiwidth2 )
            continue;
        if (i->val < min)
            min = i->val;
        if (i->val > max)
            max = i->val;
        nbrMaxima += 1;
    }
    int roi = 0;
    T thresh = 0.7 * (max - min) + min;
    std::vector<size_t> maxima_indices(maxima.size());
    std::iota(maxima_indices.begin(), maxima_indices.end(), 0);
    std::random_shuffle(maxima_indices.begin(), maxima_indices.end());

    nbrRoisPerFrame = std::min(nbrRoisPerFrame, nbrMaxima);

    for (auto i = maxima_indices.begin(); roi < nbrRoisPerFrame && i != maxima_indices.end(); ++i) {
        /*const*/ Coord<T> &maximum = maxima[*i];
        if (maximum.x < roiwidth2+1 || maximum.x > w - roiwidth2-1 || maximum.y < roiwidth2+1 || maximum.y > h - roiwidth2-1 || maximum.val < thresh)
            continue;
        Shape2 roi_ul(maximum.x - roiwidth2, maximum.y - roiwidth2);
        Shape2 roi_lr(maximum.x - roiwidth2 + roiwidth, maximum.y - roiwidth2 + roiwidth);
        MultiArray<2, FFTWComplex<S>> fourier(ps.shape());

        MultiArray<2, FFTWComplex<S>> workImage(ps.shape()); //libfftw needs continuous memory
        vigra::copyMultiArray(srcMultiArrayRange(in.subarray(roi_ul,  roi_lr)), destMultiArray(workImage, FFTWWriteRealAccessor<S>()));
        fplan.execute(workImage, fourier);
        vigra::combineTwoMultiArrays(srcMultiArrayRange(ps),
                                     srcMultiArray(fourier, FFTWSquaredMagnitudeAccessor<double>()),
                                     destMultiArray(ps), std::plus<double>());
        ++roi;
        ++rois;
    }

}

void fitPSF(DataParams&, MultiArray<2, double>&);

/*!
Estimates background variance from a fit of the histogram of intensities for each frame. The assumption is that the backgroun pixels intensity values
form a gaussian in the histogram, with a variance that can be used to correct the gain.
*/
template <class T>
void getBGVariance(const DataParams &params, const MultiArrayView<2, T> &img, std::vector<T> &BGVar, int currframe) {
    vigra::acc::AccumulatorChain<T, vigra::acc::Select<vigra::acc::AutoRangeHistogram<0>>> accChain;
    auto iter = img.begin();
    auto iterEnd = iter.getEndIterator();
    vigra::FindMinMax<T> imgMinMax;
    inspectImage(srcImageRange(img), imgMinMax);
//     char namePic[1000];
//     sprintf(namePic, "/home/herrmannsdoerfer/tmpOutput/bildforGetBgVar%d.tif",currframe);
//     vigra::exportImage(srcImageRange(img), namePic);
    double varBG = 0;
    int numberBins = 100;
    if (int(imgMinMax.max - imgMinMax.min)>0) {
        wienerStorm_R_mutex.lock();
        vigra::HistogramOptions histogram_opt;
        histogram_opt = histogram_opt.setBinCount(numberBins);
        accChain.setHistogramOptions(histogram_opt);
        vigra::acc::extractFeatures(iter, iterEnd, accChain);
        vigra::MultiArray<1, double> hist2 = get<vigra::acc::AutoRangeHistogram<0>>(accChain);
//         std::cout<<std::endl;
//         for (auto it =hist2.begin(); it != hist2.end(); it++){
//             std::cout<<*it<<", ";
//         }
//         std::cout<<std::endl;
        SEXP vec, minimum, maximum, nbrbins, fun, t;
        PROTECT(vec = Rf_allocVector(REALSXP, hist2.size()));
        std::memcpy(REAL(vec), hist2.data(), hist2.size() * sizeof(double));
        PROTECT(minimum = Rf_ScalarReal(imgMinMax.min));
        PROTECT(maximum = Rf_ScalarReal(imgMinMax.max));
        PROTECT(nbrbins = Rf_allocVector(INTSXP, 1));
        std::memcpy(INTEGER(nbrbins), &numberBins, 1 * sizeof(int));
        PROTECT(fun = t = Rf_allocList(5));
        SET_TYPEOF(fun, LANGSXP);
        SETCAR(t, Rf_install("fit.BG2"));
        t = CDR(t);
        SETCAR(t, vec);
        t = CDR(t);
        SETCAR(t, minimum);
        t = CDR(t);
        SETCAR(t, maximum);
        t = CDR(t);
        SETCAR(t, nbrbins);

        PROTECT(t = Rf_eval(fun, R_GlobalEnv));
        varBG = *REAL(t);
        UNPROTECT(6);
        wienerStorm_R_mutex.unlock();
    }
    BGVar[currframe] = varBG;

}

/*!
The gain factor is adjusted iteratively until the backgrounds variance is equal to one
*/
template <class T>
void checkCameraParameters(DataParams &params, ProgressFunctor &progressFunc) {
    bool needSkellam = !((params.getIgnoreSkellamFramesSaved() or params.getSkellamFramesSaved()) && params.getSlopeSaved() && params.getInterceptSaved());
    if (!needSkellam)
        return;
    unsigned int stacksize = params.getSkellamFrames();
    std::vector<T> BGVars(stacksize);
    auto func = [&params, &BGVars](const DataParams &params, const MultiArrayView<2, T> &currSrc, int currframe) {getBGVariance(params, currSrc, BGVars, currframe);};
    progressFunc.setStage(ParameterCheck);
    T medBGVar;
    T medBGVar2;
    int maxIterLimit = 10;
    int counter = 0;
    std::vector<T> slopeResults(maxIterLimit+1);
    T initialSlope = params.getSlope();
    slopeResults[0] = initialSlope;
    while (true) {
        counter += 1;
        processStack<T>(params, func, progressFunc, stacksize);
        std::nth_element(BGVars.begin(), BGVars.begin() + BGVars.size()/2, BGVars.end());
        medBGVar = BGVars[BGVars.size()/2];
        medBGVar2 = medBGVar * medBGVar;
        if (std::abs(medBGVar - 1)< 0.1 ){
            std::cout<<std::endl<<"changing slope from: "<< params.getSlope()<<" to "<<params.getSlope()*medBGVar2<<" based on estimated background variance of: "<<medBGVar<<std::endl;
            params.setSlope(params.getSlope()*medBGVar2);
            break;
        }
        std::cout<<std::endl<<"changing slope from: "<< params.getSlope()<<" to "<<params.getSlope()*medBGVar2<<" based on estimated background variance of: "<<medBGVar<<std::endl;
        params.setSlope(params.getSlope()*medBGVar2);
        slopeResults[counter] = params.getSlope();
        BGVars.clear();
        if (counter == maxIterLimit){
            std::cout<<"maximal number of iterations for parameter check is reached without converging variance. The slope is set to the average of the last 2 guesses."<<std::endl;
            params.setSlope((slopeResults[counter-1]+slopeResults[counter-2])/2.);
            break;
        }
    }
}


/*! A multivariate gaussian is fitted to the mean power spectrum, after that corresponding sigma in spatial domain is calculated*/
template <class T>
void estimatePSFParameters(DataParams &params, ProgressFunctor &progressFunc) {
    bool needFilter = !(params.getSkellamFramesSaved() && params.getSigmaSaved());
    if (!needFilter) {
        std::cout<<"Values from settings-file:"<<std::endl;
        std::cout<<"Sigma: "<<params.getSigma();
        return;}
    progressFunc.setStage(PSFWidth);
    unsigned int stacksize = params.getSkellamFrames();
    int roiwidth = 3 * params.getRoilen();
    int nbrRoisPerFrame = 10;
    int rois = 0;

    MultiArray<2, double> ps(Shape2(roiwidth, roiwidth));
    ps.init(0.0);
    MultiArray<2, FFTWComplex<double>> filterInit(Shape2(roiwidth, roiwidth));
    FFTWPlan<2, double> fplan(filterInit, filterInit, FFTW_FORWARD, FFTW_MEASURE);
    auto func = [&fplan, &ps, &roiwidth, &nbrRoisPerFrame, &rois](const DataParams &params, MultiArrayView<2, T> &currSrc, int currframe){accumulatePowerSpectrum(params, fplan, currSrc, ps, roiwidth, nbrRoisPerFrame, rois);};
    processStack<T>(params, func, progressFunc, stacksize);
    if (progressFunc.getAbort())
        return;
    moveDCToCenter(ps);
    vigra::transformMultiArray(srcMultiArrayRange(ps), destMultiArray(ps),
                          [&stacksize, &roiwidth, &rois](double p){return p / (rois * roiwidth * roiwidth);});
    fitPSF(params, ps);
    std::cout<<"Estimated value:"<<std::endl;
    std::cout<<"Sigma: "<<params.getSigma()<<std::endl;
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
void wienerStorm(DataParams &params, std::vector<std::set<Coord<T> > >& maxima_coords, ProgressFunctor &progressFunc) {
    estimateCameraParameters<T>(params, progressFunc);
    if (progressFunc.getAbort()) {
        progressFunc.setFinished();
        return;
    }
    checkCameraParameters<T>(params, progressFunc);
    estimatePSFParameters<T>(params, progressFunc);
    if (progressFunc.getAbort()) {
        progressFunc.setFinished();
        return;
    }
    auto func = [&maxima_coords](const DataParams &params, const MultiArrayView<2, T> &currSrc, int currframe) {wienerStormSingleFrame(params, currSrc, maxima_coords[currframe], currframe);};
    progressFunc.setStage(Localization);
    maxima_coords.resize(params.shape(2));
    processStack<T>(params, func, progressFunc);
    progressFunc.setFinished();
}

template <class T>
void wienerStormAsync(DataParams &params, std::vector<std::set<Coord<T> > >& maxima_coords, ProgressFunctor &progressFunc) {
    std::thread(wienerStorm<T>, std::ref(params), std::ref(maxima_coords), std::ref(progressFunc)).detach();
}

template <class T>
void wienerStormSingleFrame(const DataParams &params, const MultiArrayView<2, T>& in, std::set<Coord<T> >& maxima_coords, int framenumber)
{
    int w = params.shape(0); // width
    int h = params.shape(1); // height

    BasicImage<T> filtered(w,h);

    int factor = params.getFactor();
    int mylen = params.getRoilen();

    const int mylen2 = mylen/2;
    unsigned int w_roi = factor*(mylen-1)+1;
    unsigned int h_roi = factor*(mylen-1)+1;
    BasicImage<T> im_xxl(w_roi, h_roi);
    BasicImage<T> unfiltered(w,h);
    BasicImageView<T> input = makeBasicImageView(in);  // access data as BasicImage
    BasicImageView<T> filteredView(filtered.data(), filtered.size());

    vigra::copyImage(srcImageRange(input), destImage(unfiltered));
    float kernelWidth = params.getSigma()*params.getPrefactorSigma();
    gaussianSmoothing(srcImageRange(input), destImage(filteredView), kernelWidth);
    //gaussianSharpening(srcImageRange(input), destImage(filteredView), 0.85, kernelWidth);

    MultiArray<2, T> mask(Shape2(w, h));
    getMask(params, unfiltered, framenumber, mask);
    std::set<Coord<T> > maxima_candidates_vect;  // we use a set for the coordinates to automatically squeeze duplicates
                                                 // (from overlapping ROIs)
    SetPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_candidates(maxima_candidates_vect, filtered.upperLeft(), 1, mask);
    vigra::localMaxima(srcImageRange(filtered), destImage(filtered, maxima_candidates), vigra::LocalMinmaxOptions().neighborhood(4));
//     char name2[1000];
//     sprintf(name2, "/home/herrmannsdoerfer/tmpOutput/frameData/imgBeforeMaximaDetection%d.tif", framenumber);
//     vigra::exportImage(srcImageRange(filtered), name2);

    SetPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_acc(maxima_coords, im_xxl.upperLeft(), factor, mask);
    //upscale filtered image regions with spline interpolation
    std::set<Coord<float> >::iterator it2;

    for(it2=maxima_candidates_vect.begin(); it2 != maxima_candidates_vect.end(); it2++) {
            Coord<float> c = *it2;
            if(unfiltered(c.x,c.y)<1 or mask(c.x,c.y) == 0.0) { // skip very low signals with SNR lower 3 or maxima not covered by the mask
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
            // find local maxima.
            // here we include only internal pixels, no border
            // to get every maximum only once, the maxima are pushed into a std::set
            maxima_acc.setOffset(Diff2D(factor*(c.x-mylen2), factor*(c.y-mylen2)));
            vigra::localMaxima(srcIterRange(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), im_xxl.lowerRight()+xxl_lr-Diff2D(factor,factor)),
                               destIter(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), maxima_acc), vigra::LocalMinmaxOptions().neighborhood(8));
    }
    determineAsymmetry(srcImageRange(unfiltered), maxima_coords, params);
    determineSNR(srcImageRange(unfiltered), maxima_coords, factor);
}

bool initR(int argc, char **argv, bool withRestart = true);

void endR();

#endif
