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

//#include <algorithm>
#ifdef OPENMP_FOUND
    #include <omp.h>
#endif //OPENMP_FOUND
#include "util.h"
#include "dataparams.h"

#define R_INTERFACE_PTRS
#define R_NO_REMAP
#include <Rembedded.h>
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

enum WienerStormStage{CameraParameters, PSFWidth, Localization};
class ProgressFunctor
{
public:
    ProgressFunctor() : m_abort(false) {};
    virtual ~ProgressFunctor(){};
    virtual void setStage(WienerStormStage) = 0;
    virtual void setStackSize(int) = 0;
    virtual void setFrame(int) = 0;
    void abort() {m_abort = true;};
    bool getAbort() {return m_abort;};

protected:
    std::atomic<bool> m_abort;
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
void fitSkellamPoints(DataParams &params,T meanValues[],T skellamParameters[],int numberPoints){
    int nbins = 10;

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

//get Mask calculates the probabilities for each pixel of the current frame to be foreground,
//there are two different methods available: logliklihood, based on cumulative distribution function
template <class T>
void getMask(const DataParams &params, const BasicImage<T>& array, int framenumber, MultiArray<2,T>& mask){
	double cdf = params.getMaskThreshold();
    vigra::transformImage(srcImageRange(array), destImage(mask), [&cdf](T p) {return p >= cdf ? 1 : 0;});

//     vigra::IImage labels(array.width(), array.height());
//     unsigned int nbrCC = vigra::labelImageWithBackground(srcImageRange(mask), destImage(labels), false, 0);
//     std::valarray<int> bins(0, nbrCC + 1);
//     auto fun = [&bins](int32_t p){++bins[p];};
//     vigra::inspectImage(srcImageRange(labels), fun);
//     vigra::transformImage(srcImageRange(labels), destImage(mask), [&params, &bins](T p) {if(!p || bins[p] < /*3.14 * std::pow(params.getSigma(), 2)*/3) return 0; else return 1;});
}

//To estimate the gain factor points with different mean intensities are needed. This functions searches for
//good candidates, it tries to find pixels with as much different mean values as possible.
template <class T>
void estimateCameraParameters(DataParams &params, ProgressFunctor &progressFunc) {
    bool needSkellam = !(params.getSkellamFramesSaved() && params.getSlopeSaved() && params.getInterceptSaved());
    if (!needSkellam)
        return;
    unsigned int stacksize = params.getSkellamFrames();
    unsigned int w = params.shape(0);
    unsigned int h = params.shape(1);

    progressFunc.setStage(CameraParameters);
    progressFunc.setStackSize(stacksize);
    int maxNumberPoints = 2000;

    T minVal = std::numeric_limits<T>::max();
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
        progressFunc.setFrame(f);
    }
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
    delete lastVal;
    delete img;
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
            regionMeans(x, y) = vigra::acc::get<vigra::acc::StandardQuantiles<vigra::acc::AutoRangeHistogram<0>>>(accChain, i)[3];
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
                                     [](T srcPixel, T poissonPixel){T val = srcPixel - poissonPixel; return (std::isnan(val)) ? 0 : val;});
        functor(params, currSrc, currframe + f);
        progressFunc.setFrame(currframe + f);
    }
    currframe += srcImage.shape()[2];
}

template <class T, class L, class F>
void readChunk(const DataParams &params, MultiArray<3, T>** srcImage,
               MultiArrayView<3, T> &poissonMeansRaw, MultiArrayView<3, T> &poissonMeans,
               const MultiArrayView<3, L> &poissonLabels, int chunk, F &tF) {
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
    vigra::transformMultiArray(srcMultiArrayRange(*tmp), destMultiArrayRange(*tmp), [&params, &tF](T p){return tF((p - params.getIntercept()) / params.getSlope());});
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
    MultiArray<3, int> poissonLabels;
    MultiArray<3, T> poissonMeansRaw(Shape3(xChunks, yChunks, chunksInMemory));
    MultiArray<3, T> poissonMeans(Shape3(w, h, params.getTChunkSize() * chunksInMemory));
    getPoissonLabelsArray(params, poissonLabels);

    int currframe = 0, chunk;

    for (chunk = 0; chunk < chunksInMemory; ++chunk) {
        params.readBlock(Shape3(0, 0, chunk * params.getTChunkSize()), srcImage[chunk]->shape(), *srcImage[chunk]);
        auto currRawMean = poissonMeansRaw.bindOuter(chunk);
        vigra::transformMultiArray(srcMultiArrayRange(*srcImage[chunk]), destMultiArrayRange(*srcImage[chunk]), [&params, &tF](T p){return tF((p - params.getIntercept()) / params.getSlope());});
        getPoissonMeansForChunk(params, poissonLabels, *srcImage[chunk], currRawMean);

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
        readChunk(params, srcImage, poissonMeansRaw, poissonMeans, poissonLabels, chunk, tF);
        processChunk(params, *srcImage[0], poissonMeans, currframe, middleChunk, functor, progressFunc);
    }
    if (lastChunkSize) {
        if (progressFunc.getAbort())
            return;
        srcImage[0]->reshape(Shape3(w, h, lastChunkSize));
        Shape3 labelsShape = poissonLabels.shape();
        labelsShape[2] = lastChunkSize;
        auto lastPoissonLabelsView = poissonLabels.subarray(Shape3(0, 0, 0), labelsShape);
        readChunk(params, srcImage, poissonMeansRaw, poissonMeans, lastPoissonLabelsView, chunk, tF);
        processChunk(params, *srcImage[0], poissonMeans, currframe, middleChunk, functor, progressFunc);
    }
    delete srcImage[0];
    for (int c = middleChunk + 1; c < chunksInMemory; ++c) {
        if (progressFunc.getAbort())
            return;
        int cIndex = c - middleChunk;
        processChunk(params, *srcImage[cIndex], poissonMeans, currframe, c, functor, progressFunc);
        helper::progress(currframe, i_end);
        delete srcImage[cIndex];
    }
    std::free(srcImage);
}

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
    for (auto i = maxima.begin(); i != maxima.end(); ++i) {
        if (i->x < roiwidth2 || i->x > w - roiwidth2 || i->y < roiwidth2 || i->y > h - roiwidth2)
            continue;
        if (i->val < min)
            min = i->val;
        if (i->val > max)
            max = i->val;
    }
    int roi = 0;
    T thresh = 0.7 * (max - min) + min;
    std::vector<size_t> maxima_indices(maxima.size());
    std::iota(maxima_indices.begin(), maxima_indices.end(), 0);
    std::random_shuffle(maxima_indices.begin(), maxima_indices.end());

    for (auto i = maxima_indices.begin(); roi < nbrRoisPerFrame && i != maxima_indices.end(); ++i) {
        const Coord<T> &maximum = maxima[*i];
        if (maximum.x < roiwidth2 || maximum.x > w - roiwidth2 || maximum.y < roiwidth2 || maximum.y > h - roiwidth2 || maximum.val < thresh)
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

void fitPSF(DataParams &params, MultiArray<2, double> &ps) {
    SEXP mat, fun, t;
    PROTECT(mat = Rf_allocMatrix(REALSXP, ps.shape(1), ps.shape(0)));
    std::memcpy(REAL(mat), ps.data(), ps.size() * sizeof(double));
    PROTECT(fun = t = Rf_allocList(2));
    SET_TYPEOF(fun, LANGSXP);
    SETCAR(t, Rf_install("fit.filter"));
    t = CDR(t);
    SETCAR(t, mat);

    PROTECT(t = Rf_eval(fun, R_GlobalEnv));
    double *sigmas = REAL(t);
    double sigmax = ps.shape(0) / (2 * std::sqrt(2) * M_PI * sigmas[0]);
    double sigmay = ps.shape(1) / (2 * std::sqrt(2) * M_PI * sigmas[1]);
    params.setSigma((sigmax + sigmay) / 2);

    UNPROTECT(3);
}

template <class T>
void estimatePSFParameters(DataParams &params, ProgressFunctor &progressFunc) {
    std::srand(42);
    bool needFilter = !(params.getSkellamFramesSaved() && params.getSigmaSaved());
    if (!needFilter)
        return;
    progressFunc.setStage(PSFWidth);
    unsigned int stacksize = params.getSkellamFrames();
    int roiwidth = 3 * params.getRoilen();
    int nbrRoisPerFrame = 20;
    int rois = 0;

    MultiArray<2, double> ps(Shape2(roiwidth, roiwidth));
    ps.init(0.0);
    MultiArray<2, FFTWComplex<double>> filterInit(Shape2(roiwidth, roiwidth));
    FFTWPlan<2, double> fplan(filterInit, filterInit, FFTW_FORWARD, FFTW_MEASURE);
    auto func = [&fplan, &ps, &roiwidth, &nbrRoisPerFrame, &rois](const DataParams &params, MultiArrayView<2, T> &currSrc, int currframe){accumulatePowerSpectrum(params, fplan, currSrc, ps, roiwidth, nbrRoisPerFrame, rois);};
    processStack<T>(params, func, progressFunc, stacksize);
    moveDCToCenter(ps);
    vigra::transformMultiArray(srcMultiArrayRange(ps), destMultiArray(ps),
                          [&stacksize, &roiwidth, &rois](double p){return p / (rois * roiwidth * roiwidth);});
    fitPSF(params, ps);
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
    estimatePSFParameters<T>(params, progressFunc);
    auto func = [&maxima_coords](const DataParams &params, const MultiArrayView<2, T> &currSrc, int currframe) {wienerStormSingleFrame(params, currSrc, maxima_coords[currframe], currframe);};
    progressFunc.setStage(Localization);
    processStack<T>(params, func, progressFunc);
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
    float kernelWidth = params.getSigma() < 0.85 ? 0.0 : std::sqrt(std::pow(params.getSigma(), 2)-std::pow(0.85, 2));
    gaussianSmoothing(srcImageRange(input), destImage(filteredView), kernelWidth);

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
    SetPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_candidates(maxima_candidates_vect, filtered.upperLeft(), 1, mask);
    vigra::localMaxima(srcImageRange(filtered), destImage(filtered, maxima_candidates), vigra::LocalMinmaxOptions().neighborhood(4));

    SetPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_acc(maxima_coords, im_xxl.upperLeft(), factor, mask);
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
                               destIter(im_xxl.upperLeft()+xxl_ul+Diff2D(factor,factor), maxima_acc), vigra::LocalMinmaxOptions().neighborhood(8));
    }
    determineAsymmetry(srcImageRange(unfiltered), maxima_coords, factor);
    determineSNR(srcImageRange(unfiltered), maxima_coords, factor);
}

void preventRConsoleWrite(const char* buf, int buflen)
{}

bool initR(const StormParams &params, int argc, char **argv, bool withRestart = true)
{
    if (std::getenv("R_HOME") == nullptr) {
        if (!withRestart)
            return false;
        char **args = (char**)std::malloc((argc + 3) * sizeof(char*));
        args[0] = (char*)"R";
        args[1] = (char*)"CMD";
        for (int i = 0, j = 2; i < argc; ++i, ++j)
            args[j] = argv[i];
        args[argc + 2] = nullptr;
        int ret = execvp(args[0], args);
        /*std::string reason;
        switch (errno) {
            case ENOENT:
                reason << "ENOENT";
                break;
            case ENOTDIR:
                reason << "ENOTDIR";
                break;
            case E2BIG:
                reason << "E2BIG";
                break;
            case EACCES:
                reason << "EACCES";
                break;
            case EINVAL:
                reason << "EINVAL";
                break;
            case ELOOP:
                reason << "ELOOP";
                break;
            case ENOMEM:
                reason << "ENOMEM";
                break;
            case ETXTBSY:
                reason << "ETXTBSY";
                break;
            default:
                reason << "unknown";
                break;
        }*/
        std::free(args);
        return false;
    }
    char *Rargv[] = {(char*)"REmbeddedStorm", (char*)"--silent", (char*)"--no-save"};
    R_SignalHandlers = FALSE;
    Rf_initEmbeddedR(sizeof(Rargv) / sizeof(Rargv[0]), Rargv);

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
    return true;
}

void endR()
{
    // prevent printing of R warnings
    void (*ptr_R_WriteConsole_old)(const char *, int) = ptr_R_WriteConsole;
    FILE *R_Consolefile_old = R_Consolefile;
    ptr_R_WriteConsole = preventRConsoleWrite;
    R_Consolefile = NULL;
    Rf_endEmbeddedR(0);
    ptr_R_WriteConsole = ptr_R_WriteConsole_old;
    R_Consolefile = R_Consolefile_old;
}
