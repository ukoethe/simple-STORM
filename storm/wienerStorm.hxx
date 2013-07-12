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
#include <vigra/linear_solve.hxx>

#include <iostream>


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

#ifndef Q_MOC_RUN
#include <boost/math/distributions/normal.hpp>
#endif // Q_MOC_RUN


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
void fitGaussian1D(double data[], int size, T &sigma = 2, T &scale = 0.5, T &offset = 0, T &x0 = 0){
    try
    {
        vigra::linalg::Matrix<double> points(Shape2(size, 2), data);

        //double sigma = 2.0, scale = 0.5, offset = 0.0;

        double t = 1.4, l = 0.1;

        for(int k=0; k< 10; ++k)
        {
            vigra::linalg::Matrix<double> jr(4,1), jj(4,4), j(4,1);
            double tr = 0.0;

            for(int i=0; i<size; ++i)
            {
                double xs = sq((points(i,0) - x0) / sigma);
                double e = std::exp(-0.5 * xs);
                double r = points(i, 1) - (scale * e + offset);

                j(0,0) = scale * e * xs / sigma;
                j(1,0) = e;
                j(2,0) = 1.0;
                j(3,0) = scale * e * (points(i,0) - x0)/sq(sigma);

                jr += r * j;
                jj += j * transpose(j);
                tr += sq(r);
            }

            vigra::linalg::Matrix<double> jj1(jj), jj2(jj), d1(4,1), d2(4,1);

            jj1.diagonal() *= 1.0 + l;
            jj2.diagonal() *= 1.0 + l / t;

            vigra::linearSolve(jj1, jr, d1);
            vigra::linearSolve(jj2, jr, d2);

            double si1 = sigma + d1(0,0), s1 = scale + d1(1,0), o1 = offset + d1(2,0), c1 = x0 + d1(3,0);
            double si2 = sigma + d2(0,0), s2 = scale + d2(1,0), o2 = offset + d2(2,0), c2 = x0 + d2(3,0);
            double tr1 = 0.0, tr2 = 0.0;

            for(int i=0; i<size; ++i)
            {
                double r1 = points(i, 1) - (s1 * std::exp(-0.5 * sq((points(i,0) - c1) / si1)) + o1);
                double r2 = points(i, 1) - (s2 * std::exp(-0.5 * sq((points(i,0) - c2) / si2)) + o2);
                tr1 += sq(r1);
                tr2 += sq(r2);
            }

            if(tr1 < tr2)
            {
                if(tr1 < tr)
                {
                    sigma = si1;
                    scale = s1;
                    offset = o1;
                    x0 = c1;
                }
                else
                {
                    l *= t;
                }
            }
            else
            {
                if(tr2 < tr)
                {
                    sigma = si2;
                    scale = s2;
                    offset = o2;
                    x0 = c2;
                    l /= t;
                }
                else
                {
                    l *= t;
                }
            }
            //std::cerr << "tr: " << (0.5*tr) << " sigma: " << sigma << " scale: " << scale << " offset: " << offset << " center: "<< x0 << "\n";
            if(std::abs((tr - std::min(tr1, tr2)) / tr) < 1e-15)
                break;
        }
        //std::cerr << "sigma: " << sigma << " scale: " << scale << " offset: " << offset << "\n";
    }
    catch (std::exception & e)
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
    }
}

template <class T>
void fitGaussian(double data[], int size, T &sigma = 2, T &scale = 0.5, T &offset = 0){
    try
    {
        vigra::linalg::Matrix<double> points(Shape2(size, 2), data);

        //double sigma = 2.0, scale = 0.5, offset = 0.0;

        double t = 1.4, l = 0.1;

        for(int k=0; k< 30; ++k)
        {
            vigra::linalg::Matrix<double> jr(3,1), jj(3,3), j(3,1);
            double tr = 0.0;

            for(int i=0; i<size; ++i)
            {
                double xs = sq(points(i,0) / sigma);
                double e = std::exp(-0.5 * xs);
                double r = points(i, 1) - (scale * e + offset);
                j(0,0) = scale * e * xs / sigma;
                j(1,0) = e;
                j(2,0) = 1.0;

                jr += r * j;
                jj += j * transpose(j);
                tr += sq(r);
            }

            vigra::linalg::Matrix<double> jj1(jj), jj2(jj), d1(3,1), d2(3,1);

            jj1.diagonal() *= 1.0 + l;
            jj2.diagonal() *= 1.0 + l / t;

            vigra::linearSolve(jj1, jr, d1);
            vigra::linearSolve(jj2, jr, d2);

            double si1 = sigma + d1(0,0), s1 = scale + d1(1,0), o1 = offset + d1(2,0);
            double si2 = sigma + d2(0,0), s2 = scale + d2(1,0), o2 = offset + d2(2,0);
            double tr1 = 0.0, tr2 = 0.0;

            for(int i=0; i<size; ++i)
            {
                double r1 = points(i, 1) - (s1 * std::exp(-0.5 * sq(points(i,0) / si1)) + o1);
                double r2 = points(i, 1) - (s2 * std::exp(-0.5 * sq(points(i,0) / si2)) + o2);
                tr1 += sq(r1);
                tr2 += sq(r2);
            }

            if(tr1 < tr2)
            {
                if(tr1 < tr)
                {
                    sigma = si1;
                    scale = s1;
                    offset = o1;
                }
                else
                {
                    l *= t;
                }
            }
            else
            {
                if(tr2 < tr)
                {
                    sigma = si2;
                    scale = s2;
                    offset = o2;
                    l /= t;
                }
                else
                {
                    l *= t;
                }
            }
            std::cerr << "tr: " << (0.5*tr) << " sigma: " << sigma << " scale: " << scale << " offset: " << offset << "\n";
            if(std::abs((tr - std::min(tr1, tr2)) / tr) < 1e-15)
                break;
        }
        std::cerr << "sigma: " << sigma << " scale: " << scale << " offset: " << offset << "\n";
    }
    catch (std::exception & e)
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
    }
}


template <class T>
void doRansac(DataParams &params,T meanValues[],T skellamParameters[],int numberPoints){
    int minNumberRequired = 2;
    int numberIterations = 10000;
    int numberClose = numberPoints/10;
    int threshold = (skellamParameters[numberPoints-2]-skellamParameters[0])/20;
    T bestm =0, bestc=0;
    T lasterror = 9999999999;
    for(int i=0;i<numberIterations;++i){
        std::vector<int> iter(numberPoints);
        linearSequence(iter.begin(), iter.end());
        std::random_shuffle(iter.begin(), iter.end());
        T m = (skellamParameters[iter[1]]-skellamParameters[iter[0]])/(meanValues[iter[1]]-meanValues[iter[0]]);
        T c = skellamParameters[iter[0]]-meanValues[iter[0]]*m;
        int numberPointsInRange = 0;
        T currentError = 0;
        std::vector<T> consSetX, consSetY;
        for (int j=0; j< numberPoints; ++j){
//             currentError += std::pow(skellamParameters[j]-meanValues[j]*m+c,2);
            if((std::abs(skellamParameters[j]-meanValues[j]*m+c))<threshold){
                consSetX.push_back(meanValues[j]);
                consSetY.push_back(skellamParameters[j]);
                numberPointsInRange+=1;
            }
        }
        std::vector<int> indices(numberPointsInRange);
        linearSequence(indices.begin(), indices.end());
        fitLine(consSetX, consSetY, indices, numberPointsInRange, m,c,currentError);

//         std::cout<<"error: "<<currentError<<" m:"<<m<<"c: "<<c<<" numberPoints: "<<numberPointsInRange<<std::endl;

        if (numberPointsInRange>numberClose and currentError < lasterror){
            lasterror = currentError;
            bestm = m;
            bestc = c;
        }

    }
    if (bestm!=0 or bestc!=0){
        params.setSlope(bestm);
        params.setIntercept(-bestc/bestm);
    }
    else{
        std::cout<<"no fit found"<<std::endl;
    }


}

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
    cfile << params.shape(0) << " " << params.shape(1) << " " << params.shape(2) << " "<< params.getPixelSize()<< " " <<params.getFactor()<< " " << params.getSigma()<< " " <<params.getPrefactorSigma()<< std::endl;
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
Calculates a mask that contains the information whether or not a pixel belongs to background or is signal for eachh pixel of the current frame,
 based on cumulative distribution. The asumption is that the background has a mean value of zeros and a variance of one.
 To get rid of noise just consider larger spots
 */
template <class T>
void getMask(const DataParams &params, const BasicImage<T>& array, int framenumber, MultiArray<2,T>& mask){
	//double cdf = params.getMaskThreshold();
	boost::math::normal dist(0.0, 1.0);
    double cdf = quantile(dist, 1-params.getAlpha());//qnorm(params.getAlpha(), 0, 1, 0, 0);
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
//     if (params.getIgnoreSkellamFramesSaved()) {
//         std::cout<<"Values from GUI dialog:"<<std::endl;
//         std::cout<<"Gain: "<<params.getSlope()<<" Offset: "<<params.getIntercept()<<std::endl;
//         return;
//     }
    if (params.getSlopeSaved() and params.getInterceptSaved()) {
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
    MultiArray<2, vigra::acc::AccumulatorChain<T, vigra::acc::Select<vigra::acc::Sum, vigra::acc::Variance>>> varianceArr(Shape2(w, h));
    unsigned int passes = skellamArr(0, 0).passesRequired();
    int *listPixelCoordinates = (int*)std::malloc(maxNumberPoints * sizeof(int));
    T *meanValues = (T*)std::malloc(maxNumberPoints * sizeof(T));
    T *skellamParameters = (T*)std::malloc(maxNumberPoints * sizeof(T));

    bool useSkellam = false;

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
            if (useSkellam) {
                vigra::combineTwoMultiArrays(srcMultiArrayRange(*lastVal), srcMultiArrayRange(*img), destMultiArrayRange(*lastVal), std::minus<T>());
                auto imgIter = lastVal->begin(), imgEnd = lastVal->end();
                auto skellamIter = skellamArr.begin(), skellamEnd = skellamArr.end();
                for (; imgIter != imgEnd && skellamIter != skellamEnd; ++imgIter, ++skellamIter) {
                    for (int n = 1; n <= passes; ++n) {
                        skellamIter->updatePassN(*imgIter, n);
                    }
                }
            }
        }
        if (useSkellam) {
            auto *tmp = lastVal;
            lastVal = img;
            img = tmp;
        }

        if (not useSkellam) {
            auto imgIter = img->begin(), imgEnd = img->end();
            auto varIter = varianceArr.begin(), varEnd = varianceArr.end();
            for (; imgIter != imgEnd && varIter != varEnd; ++imgIter, ++varIter) {
                for (int n = 1; n <= passes; ++n) {
                    varIter->updatePassN(*imgIter, n);
                }
            }
        }
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
    indexSort(meanArr.begin(), meanArr.end(), iter.begin()); //means are sorted iter stores the permutation

    int intervalCounter = 0;
    T intervalDistance = (minmax.max - minmax.min)/ maxNumberPoints;

    for(int i = 0; i< w * h && intervalCounter < maxNumberPoints; i++) {
        if(meanArr[iter[i]] > minmax.min + intervalDistance * intervalCounter) {
            listPixelCoordinates[intervalCounter] = iter[i];
            meanValues[intervalCounter] = meanArr[iter[i]];
            if (useSkellam) {
                skellamParameters[intervalCounter] = 0.5 * (vigra::acc::get<vigra::acc::Sum>(skellamArr[iter[i]]) / stacksize + vigra::acc::get<vigra::acc::Variance>(skellamArr[iter[i]]));
            }
            else {
                skellamParameters[intervalCounter] =  vigra::acc::get<vigra::acc::Variance>(varianceArr[iter[i]]);
            }
            ++intervalCounter;
        }
    }

    std::cout<<std::endl;
    for (int i =0; i< intervalCounter; ++i){
        std::cout<<meanValues[i]<<", ";
    }
    std::cout<<std::endl;
    for (int i =0; i< intervalCounter; ++i){
        std::cout<<skellamParameters[i]<<", ";
    }

    doRansac(params, meanValues, skellamParameters, intervalCounter);
    std::cout<<"y4 = "<< params.getSlope()<<" * x1 + "<<-params.getSlope()*params.getIntercept() <<"#(real RANSAC)"<< std::endl;

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
    vigra::transformMultiArray(srcMultiArrayRange(in), destMultiArray(in),
                               [](double p){return std::pow((p+3./8.)/2.,2)-3./8.;}); //inverse Anscombe transform to avoid spread of PSF
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

        MultiArray<2, T> expImage(ps.shape());
        vigra::copyMultiArray(srcMultiArrayRange(in.subarray(roi_ul,  roi_lr)), destMultiArray(expImage));

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
void getBGVariance(const DataParams &params, const MultiArrayView<2, T> &img, std::vector<T> &BGStd, int currframe) {
    vigra::acc::AccumulatorChain<T, vigra::acc::Select<vigra::acc::AutoRangeHistogram<0>>> accChain;
    auto iter = img.begin();
    auto iterEnd = iter.getEndIterator();
    vigra::FindMinMax<T> imgMinMax;
    inspectImage(srcImageRange(img), imgMinMax);
//     std::cout<<"minimage: "<<imgMinMax.min<<"max image: "<<imgMinMax.max<<std::endl;
//     char namePic[1000];
//     sprintf(namePic, "/home/herrmannsdoerfer/tmpOutput/bildforGetBgVar%d.tif",currframe);
//     vigra::exportImage(srcImageRange(img), namePic);
    double stdBG = 0;
    int numberBins = 100;
    T sigma = 1;
    if (int(imgMinMax.max - imgMinMax.min)>0) {
        vigra::HistogramOptions histogram_opt;
        histogram_opt = histogram_opt.setBinCount(numberBins);
        accChain.setHistogramOptions(histogram_opt);
        vigra::acc::extractFeatures(iter, iterEnd, accChain);
        vigra::MultiArray<1, double> hist2 = get<vigra::acc::AutoRangeHistogram<0>>(accChain);

        double data[2*numberBins];
        T minimum = imgMinMax.min, maximum = imgMinMax.max;

        T delta = (maximum - minimum)/ (1.*numberBins);
//         std::cout<<std::endl;
        T mindata = 9999999, maxdata = 0;
        for (auto counter = 0; counter < numberBins; ++counter){
            data[2*counter] = minimum + counter * delta;
            data[2*counter+1] = hist2(counter);
//             std::cout<<hist2(counter)<<",";
            if (hist2(counter)>maxdata){maxdata = hist2(counter);}
            if (hist2(counter)<mindata){mindata = hist2(counter);}
        }
        T scale = maxdata - mindata, offset = mindata, center = 0;
//         std::cout<<std::endl;
//         std::cout<<maximum<<" "<<minimum<<" "<<scale<<" "<<offset<<" "<<delta<<std::endl;
        fitGaussian1D(data, numberBins, sigma, scale, offset, center);
    }
    BGStd[currframe] = sigma;
//     std::cout<<"cur stdBG: "<<stdBG<<" BGStd[currframe]: "<<BGStd[currframe]<<" currframe: "<<currframe<<std::endl;
//     std::cout<<BGStd[currframe]<<" currframe: "<<currframe<<" ";
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
    std::vector<T> BGStds(stacksize);
    auto func = [&params, &BGStds](const DataParams &params, const MultiArrayView<2, T> &currSrc, int currframe) {getBGVariance(params, currSrc, BGStds, currframe);};
    progressFunc.setStage(ParameterCheck);
    T medBGStd;
    T medBGStd2;
    int maxIterLimit =10;
    int counter = 0;
    std::vector<T> slopeResults(maxIterLimit+1);
    T initialSlope = params.getSlope();
    slopeResults[0] = initialSlope;
    while (true) {
        counter += 1;
        processStack<T>(params, func, progressFunc, stacksize);
//         for(auto it= BGStds.begin(); it!=BGStds.end(); ++it){
//             std::cout<<*it<<" ";
//         }
//         std::cout<<std::endl;
        std::nth_element(BGStds.begin(), BGStds.begin() + BGStds.size()/2, BGStds.end());
        medBGStd = BGStds[BGStds.size()/2];
        medBGStd2 = medBGStd * medBGStd;
        if (std::abs(medBGStd- 1)< 0.1 ){
            std::cout<<std::endl<<"changing slope from: "<< params.getSlope()<<" to "<<params.getSlope()*medBGStd2<<" based on estimated background standard deviation of: "<<medBGStd<<std::endl;
            params.setSlope(params.getSlope()*medBGStd2);
            break;
        }
        std::cout<<std::endl<<"changing slope from: "<< params.getSlope()<<" to "<<params.getSlope()*medBGStd2<<" based on estimated background standard deviation of: "<<medBGStd<<std::endl;
        params.setSlope(params.getSlope()*medBGStd2);
        slopeResults[counter] = params.getSlope();
//         BGStds.clear();
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
    std::cout<<params.getSkellamFramesSaved()<<" "<<params.getSigmaSaved()<<" "<<params.getIgnoreSkellamFramesSaved()<<std::endl;
    if ( params.getSigmaSaved()) {
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
//     BasicImage<T> output(w,h);
//     vigra::copyImage(srcImageRange(input), destImage(output));
//     char name2[1000];
//     sprintf(name2, "/home/herrmannsdoerfer/tmpOutput/frameData/imgBeforeMaximaDetection%d.tif", framenumber);
//     vigra::exportImage(srcImageRange(filtered), name2);

    MultiArray<2, T> mask(Shape2(w, h));
    getMask(params, unfiltered, framenumber, mask);
    std::set<Coord<T> > maxima_candidates_vect;  // we use a set for the coordinates to automatically squeeze duplicates
                                                 // (from overlapping ROIs)
    SetPushAccessor<Coord<T>, T, typename BasicImage<T>::const_traverser> maxima_candidates(maxima_candidates_vect, filtered.upperLeft(), 1, mask);
    vigra::localMaxima(srcImageRange(filtered), destImage(filtered, maxima_candidates), vigra::LocalMinmaxOptions().neighborhood(4));


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



