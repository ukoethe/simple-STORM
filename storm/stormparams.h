/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*      Copyright 2010-2013 by Joachim Schleicher, Ilia Kats            */
/*				and Frank Herrmannsdoerfer
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

#ifndef STORMPARAMS_H
#define STORMPARAMS_H

#include <string>
#include <set>
#include <vigra/impex.hxx>
#include <vigra/sifImport.hxx>
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#define STORMPARAMS_N 3 // could eventually be a template parameter later on

namespace rude {
class Config;
}

enum FileType { UNDEFINED, TIFF, HDF5, SIF };

class BadConversion : public std::runtime_error {
public:
    BadConversion(std::string const& s)
    : std::runtime_error(s)
    { }
};

/*! 
Class that contains all needed parameters, settings and also provides the information about the input file
*/
class StormParams {
public:
    typedef vigra::MultiArrayShape<STORMPARAMS_N>::type Shape;

    StormParams();
    StormParams(const StormParams&);
    StormParams(int argc, char **argv);
    virtual ~StormParams();
    StormParams& operator=(const StormParams&);

    void printUsage() const;
    void printVersion() const;
    int getFactor() const;
    void setFactor(int);
    bool getFactorSaved() const;
    int getRoilen() const;
    void setRoilen(int);
    bool getRoilenSaved() const;
    float getPixelSize() const;
    void setPixelSize(float);
    bool getPixelSizeSaved() const;
    unsigned int getSkellamFrames() const;
    void setSkellamFrames(unsigned int);
    bool getSkellamFramesSaved() const;
    unsigned int getXYChunkSize() const;
    void setXYChunkSize(unsigned int);
    bool getXYChunkSizeSaved() const;
    unsigned int getTChunkSize() const;
    void setTChunkSize(unsigned int);
    bool getTChunkSizeSaved() const;
    unsigned int getChunksInMemory() const;
    void setChunksInMemory(unsigned int);
    bool getChunksInMemorySaved() const;
    const std::string& getInFile() const;
    void setInFile(const std::string&, bool forceDefaults = false);
    const std::string& getOutFile() const;
    void setOutFile(const std::string&);
    const std::string& getCoordsFile() const;
    void setCoordsFile(const std::string&);
    const std::string& getSettingsFile() const;
    void setSettingsFile(const std::string&);
    const std::string& getFrameRange() const;
    void setFrameRange(const std::string&);
    bool getFrameRangeSaved() const;
    float getAlpha() const;
    void setAlpha(float);
    double getMaskThreshold() const;
    bool getDoAsymmetryCheck() const;
    void setDoAsymmetryCheck(bool);
    float getAsymmetryThreshold() const;
    void setAsymmetryThreshold(float);
    bool getDoAsymmetryCheckSaved() const;
    int getMaxXyChunksize() const;
    int getMinXyChunksize() const;
    int getMaxTChunkSize() const;
    int getMinTChunkSize() const;
    float getMaxAsymmetryThreshold() const;
    float getMinAsymmetryThreshold() const;
    bool getVerbose() const;
    void setVerbose(bool);

    const Shape & shape() const;
    vigra::MultiArrayIndex shape(const int dim) const;
    FileType type() const;
    const std::set<std::string>& acceptedFileTypes();

    template <typename  T>
    void readVolume(vigra::MultiArrayView<STORMPARAMS_N, T> &) const;
    template <typename  T>
    void readBlock(const Shape&, const Shape&, vigra::MultiArrayView<STORMPARAMS_N, T>&) const;

    virtual void save() const;
    void load(bool propagate = true);
    void doSanityChecks();

protected:
    mutable rude::Config *m_config;
    virtual void loadSettings(bool);

private:
    int parseProgramOptions(int argc, char **argv);
    void setDefaultFileNames(bool);
    void setDefaults();

    mutable void * ptr; // hack
    Shape m_shape;
    FileType m_type;
    int m_factor;
    bool m_factorSaved;
    int m_roilen;
    bool m_roilenSaved;
    float m_pixelsize;
    bool m_pixelsizeSaved;
    unsigned int m_skellamFrames;
    bool m_skellamFramesSaved;
    unsigned int m_xyChunkSize;
    bool m_xyChunkSizeSaved;
    unsigned int m_tChunkSize;
    bool m_tChunkSizeSaved;
    unsigned int m_chunksInMemory;
    bool m_chunksInMemorySaved;
    bool m_framesSaved;
    float m_alpha;
    double m_thresholdMask;
    float m_asymmetryThreshold;
    bool m_doAsymmetryCheck;
    bool m_doAsymmetryCheckSaved;
    const int m_minXyChunksize = 3;
    const int m_maxXyChunksize = 100;
    const int m_minTChunksize = 3;
    const int m_maxTChunksize = 100;
    const float m_minAsymmetryThreshold = 0;
    const float m_maxAsymmetryThreshold = 5;
    bool m_verbose;
    std::string m_infile;
    std::string m_outfile;
    std::string m_coordsfile;
    std::string m_settingsfile;
    std::string m_frames;

    std::set<std::string> m_acceptedFileTypes;

    static const std::string s_section;
};

extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, int8_t>&) const;
extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, int16_t>&) const;
extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, int32_t>&) const;
extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, unsigned int8_t>&) const;
extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, unsigned int16_t>&) const;
extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, unsigned int32_t>&) const;
template<>
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, float>&) const;
extern template
void StormParams::readVolume(vigra::MultiArrayView<STORMPARAMS_N, double>&) const;

extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, int8_t>&) const;
extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, int16_t>&) const;
extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, int32_t>&) const;
extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, unsigned int8_t>&) const;
extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, unsigned int16_t>&) const;
extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, unsigned int32_t>&) const;
template<>
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, float>&) const;
extern template
void StormParams::readBlock(const StormParams::Shape&,
                const StormParams::Shape&,
                vigra::MultiArrayView<STORMPARAMS_N, double>&) const;

#endif // STORMPARAMS_H
