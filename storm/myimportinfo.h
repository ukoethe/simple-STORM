/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*      Copyright 2010-2011 by Joachim Schleicher                       */
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

#include <string>
#include <vigra/impex.hxx>
#include <vigra/sifImport.hxx>
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#ifndef MYIMPORTINFO_H
#define MYIMPORTINFO_H

#define MYIMPORT_N 3 // could eventually be a template parameter later on

#include "program_options_getopt.h"
#include "configVersion.hxx"

inline double convertToDouble(const char* const s) {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     throw BadConversion("convertToDouble(\"\")");
   return x;
}

enum FileType { UNDEFINED, TIFF, HDF5, SIF };

using namespace vigra;

class MyImportInfo {
    typedef vigra::MultiArrayShape<MYIMPORT_N>::type Shape;
  public:
    MyImportInfo(int argc, char **argv);
    ~MyImportInfo();

    void printUsage(const char* prog);
    void setDefaults();
    int parseProgramOptions(int argc, char **argv);

    int getFactor() const;
    int getRoilen() const;
    float getThreshold() const;
    float getPixelsize() const;
    const std::string& getInfile() const;
    const std::string& getOutfile() const;
    const std::string& getCoordsfile() const;
    const std::string& getFilterfile() const;
    const std::string& getFrameRange() const;
    char getVerbose() const;

    const Shape & shape() const;
    vigra::MultiArrayIndex shape(const int dim) const;
    vigra::MultiArrayIndex shapeOfDimension(const int dim) const;
    FileType type() const;
    const std::string& executableDir() const;

    char verbose;

    void * ptr; // hack

  private:
    std::string m_filename;
    Shape m_shape;
    FileType m_type;
    std::string m_executableDir;
    int m_factor;
    int m_roilen;
    float m_threshold;
    float m_pixelsize;
    std::string m_infile;
    std::string m_outfile;
    std::string m_coordsfile;
    std::string m_filterfile;
    std::string m_frames;


};

template <class  T>
void readVolume(MyImportInfo & info, MultiArrayView<MYIMPORT_N, T> & array) {
    std::string filename = info.getInfile();
    switch(info.type()) {
        case TIFF:
        {
            ImageImportInfo* info2 = reinterpret_cast<ImageImportInfo*>(info.ptr);
            vigra_precondition(array.size(2)==info2->numImages(),"array shape and number of images in tiff file differ.");
            for(int i = 0; i < info2->numImages(); ++i) {
                MultiArrayView <2, T> img = array.bindOuter(i);
                BasicImageView <T> v = makeBasicImageView(img);
                info2->setImageIndex(i);
                importImage(*info2, destImage(v));
            }
            break;
        }
        case SIF:
        {
            SIFImportInfo* info2 = reinterpret_cast<SIFImportInfo*>(info.ptr);
            readSIF(*info2, array);
            break;
        }
        #ifdef HDF5_FOUND
        case HDF5:
        {
            HDF5File* info2 = reinterpret_cast<HDF5File*>(info.ptr);
            info2->read("/data", array);
            break;
        }
        #endif // HDF5_FOUND
        default:
            vigra_fail("decoder for type not implemented.");
    }
}

template <class  T>
void readBlock(const MyImportInfo & info,
            const MultiArrayShape<MYIMPORT_N>::type& blockOffset,
            const MultiArrayShape<MYIMPORT_N>::type& blockShape,
            MultiArrayView<MYIMPORT_N, T> & array)
{
    std::string filename = info.getInfile();
    switch(info.type()) {
        case TIFF:
        {
            ImageImportInfo* info2 = reinterpret_cast<ImageImportInfo*>(info.ptr);
            vigra_precondition(blockOffset[0]==0 && blockOffset[1]==0 &&
                    blockShape[0]==info.shapeOfDimension(0) && blockShape[1]==info.shapeOfDimension(1),
                    "for Tiff images only complete Frames are currently supported as ROIs");
            vigra_precondition(array.size(2)==blockShape[2],"array shape and number of images in ROI for tiff file differ.");
            vigra_precondition(blockShape[2] <= info2->numImages(), "block shape larger than number of frames in the image");
            for(int i = 0; i < blockShape[2]; ++i) {
                MultiArrayView <2, T> img = array.bindOuter(i);
                BasicImageView <T> v = makeBasicImageView(img);
                info2->setImageIndex(i+blockOffset[2]);
                importImage(*info2, destImage(v));
            }
            break;
        }
        case SIF:
        {
            SIFImportInfo* info2 = reinterpret_cast<SIFImportInfo*>(info.ptr);
            readSIFBlock(*info2, blockOffset, blockShape, array);
            break;
        }
        #ifdef HDF5_FOUND
        case HDF5:
        {
            HDF5File* info2 = reinterpret_cast<HDF5File*>(info.ptr);
            info2->readBlock("/data", blockOffset, blockShape, array);
            break;
        }
        #endif // HDF5_FOUND
        default:
            vigra_fail("decoder for type not implemented.");
    }
}

#endif // MYIMPORTINFO_H
