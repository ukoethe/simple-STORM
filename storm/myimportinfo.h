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

#include "configVersion.hxx"

enum FileType { UNDEFINED, TIFF, HDF5, SIF };

class BadConversion : public std::runtime_error {
public:
    BadConversion(std::string const& s)
    : std::runtime_error(s)
    { }
};

class MyImportInfo {
public:
    typedef vigra::MultiArrayShape<MYIMPORT_N>::type Shape;

    MyImportInfo(int argc, char **argv);
    ~MyImportInfo();

    void printUsage() const;
    void printVersion() const;
    int getFactor() const;
    int getRoilen() const;
    float getThreshold() const;
    float getPixelsize() const;
    unsigned int getSkellamFrames() const;
    unsigned int getXYChunkSize() const;
    unsigned int getTChunkSize() const;
    unsigned int getChunksInMemory() const;
    const std::string& getInfile() const;
    const std::string& getOutfile() const;
    const std::string& getCoordsfile() const;
    const std::string& getFilterfile() const;
    const std::string& getFrameRange() const;
    char getVerbose() const;

    const Shape & shape() const;
    vigra::MultiArrayIndex shape(const int dim) const;
    FileType type() const;
    const std::string& executableDir() const;

    template <typename  T>
    void readVolume(vigra::MultiArrayView<MYIMPORT_N, T> &) const;
    template <typename  T>
    void readBlock(const Shape&, const Shape&, vigra::MultiArrayView<MYIMPORT_N, T>&) const;

    char verbose;

    mutable void * ptr; // hack

private:
    int parseProgramOptions(int argc, char **argv);
    void setDefaults();

    std::string m_filename;
    Shape m_shape;
    FileType m_type;
    std::string m_executableDir;
    std::string m_executableName;
    int m_factor;
    int m_roilen;
    float m_threshold;
    float m_pixelsize;
    unsigned int m_skellamFrames;
    unsigned int m_xyChunkSize;
    unsigned int m_tChunkSize;
    unsigned int m_chunksInMemory;
    std::string m_infile;
    std::string m_outfile;
    std::string m_coordsfile;
    std::string m_filterfile;
    std::string m_frames;
};

extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, int8_t>&) const;
extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, int16_t>&) const;
extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, int32_t>&) const;
extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, unsigned int8_t>&) const;
extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, unsigned int16_t>&) const;
extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, unsigned int32_t>&) const;
template<>
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, float>&) const;
extern template
void MyImportInfo::readVolume(vigra::MultiArrayView<MYIMPORT_N, double>&) const;

extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, int8_t>&) const;
extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, int16_t>&) const;
extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, int32_t>&) const;
extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, unsigned int8_t>&) const;
extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, unsigned int16_t>&) const;
extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, unsigned int32_t>&) const;
template<>
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, float>&) const;
extern template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
                const MyImportInfo::Shape&,
                vigra::MultiArrayView<MYIMPORT_N, double>&) const;

#endif // MYIMPORTINFO_H
