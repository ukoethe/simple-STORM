/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*      Copyright 2010-2013 by Joachim Schleicher, Ilia Kats            */
/*                          and Frank Herrmannsdoerfer					*/
/*															            */
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

#include "stormparams.h"
#include "version.h"

#include <cstdint>
#include <string>
#include <iostream>
#ifndef EMULATE_GETOPT
#include <getopt.h>
#else
#include "getoptMSVC.h"
#endif // EMULATE_GETOPT

#include <vigra/impex.hxx>
#include <vigra/sifImport.hxx>
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif
#include <rude/config.h>
#include <Rmath.h>

using namespace vigra;

inline double convertToDouble(const char* const s) {
    char *endptr;
    double x = std::strtod(s, &endptr);
    if (endptr == s)
        throw BadConversion("convertToDouble(\"\")");
    return x;
}

inline float convertToFloat(const char* const s) {
    char *endptr;
    double x = std::strtof(s, &endptr);
    if (endptr == s)
        throw BadConversion("convertToFloat(\"\")");
    return x;
}

inline long convertToLong(const char* const s) {
    char *endptr;
    long x = std::strtol(s, &endptr, 0);
    if (endptr == s)
        throw BadConversion("convertToLong(\"\")");
    return x;
}

inline unsigned long convertToULong(const char* const s) {
    char *endptr;
    unsigned long x = std::strtoul(s, &endptr, 0);
    if (endptr == s)
        throw BadConversion("convertToULong(\"\")");
    return x;
}

const std::string StormParams::s_section = "stormparams";

StormParams::StormParams()
: m_config(new rude::Config()) {
    setDefaults();
}

StormParams::StormParams(const StormParams &other)
: m_config(new rude::Config()), m_shape(other.m_shape), m_type(other.m_type), m_factor(other.m_factor), m_factorSaved(other.m_factorSaved), m_roilen(other.m_roilen), m_roilenSaved(other.m_roilenSaved), m_pixelsize(other.m_pixelsize), m_pixelsizeSaved(other.m_pixelsize), m_skellamFrames(other.m_skellamFrames), m_skellamFramesSaved(other.m_skellamFramesSaved), m_xyChunkSize(other.m_xyChunkSize), m_xyChunkSizeSaved(other.m_xyChunkSizeSaved), m_tChunkSize(other.m_tChunkSize), m_tChunkSizeSaved(other.m_tChunkSizeSaved), m_chunksInMemory(other.m_chunksInMemory), m_chunksInMemorySaved(other.m_chunksInMemorySaved), m_framesSaved(other.m_framesSaved), m_alpha(other.m_alpha), m_thresholdMask(other.m_thresholdMask), m_doAsymmetryCheck(other.m_doAsymmetryCheck), m_doAsymmetryCheckSaved(other.m_doAsymmetryCheckSaved), m_verbose(other.m_verbose), m_outfile(other.m_outfile), m_coordsfile(other.m_coordsfile), m_settingsfile(other.m_settingsfile), m_frames(other.m_frames), m_acceptedFileTypes(other.m_acceptedFileTypes) {
    setInFile(other.m_infile, false);
    m_config->setConfigFile(m_settingsfile.c_str());
    m_config->load();
}

StormParams& StormParams::operator=(const StormParams &other)
{
    m_shape = other.m_shape;
    m_type = other.m_type;
    m_factor = other.m_factor;
    m_factorSaved = other.m_factorSaved;
    m_roilen = other.m_roilen;
    m_roilenSaved = other.m_roilenSaved;
    m_pixelsize = other.m_pixelsize;
    m_pixelsizeSaved = other.m_pixelsizeSaved;
    m_skellamFrames = other.m_skellamFrames;
    m_skellamFramesSaved = other.m_skellamFramesSaved;
    m_xyChunkSize = other.m_xyChunkSize;
    m_xyChunkSizeSaved = other.m_xyChunkSizeSaved;
    m_tChunkSize = other.m_tChunkSize;
    m_tChunkSizeSaved = other.m_tChunkSizeSaved;
    m_chunksInMemory = other.m_chunksInMemory;
    m_chunksInMemorySaved = other.m_chunksInMemorySaved;
    m_framesSaved = other.m_framesSaved;
    m_alpha = other.m_alpha;
    m_thresholdMask = other.m_thresholdMask;
    m_verbose = other.m_verbose;
    m_outfile = other.m_outfile;
    m_coordsfile = other.m_coordsfile;
    m_settingsfile = other.m_settingsfile;
    m_frames = other.m_frames;
    m_acceptedFileTypes = other.m_acceptedFileTypes;

    setInFile(other.m_infile, false);
    m_config->clear();
    m_config->setConfigFile(m_settingsfile.c_str());
    m_config->load();
}

StormParams::StormParams(int argc, char **argv)
: m_config(new rude::Config()) {
    setDefaults();
    parseProgramOptions(argc, argv);
}

StormParams::~StormParams() {
    switch(m_type) {
        case TIFF:
            delete (ImageImportInfo*)ptr;
            break;
        case SIF:
            delete (SIFImportInfo*)ptr;
            break;
        #ifdef HDF5_FOUND
        case HDF5:
            delete (HDF5File*)ptr;
            break;
        #endif // HDF5_FOUND
        default:
            break;
    }
    delete m_config;
}

int StormParams::getFactor() const {
    return m_factor;
}
void StormParams::setFactor(int factor) {
    if (factor != m_factor) {
        m_factor = factor;
        m_factorSaved = false;
    }
}
bool StormParams::getFactorSaved() const {
    return m_factorSaved;
}

int StormParams::getRoilen() const {
    return m_roilen;
}
void StormParams::setRoilen(int roilen) {
    if (roilen != m_roilen) {
        m_roilen = roilen;
        m_roilenSaved = false;
    }
}
bool StormParams::getRoilenSaved() const {
    return m_roilenSaved;
}

float StormParams::getPixelSize() const {
    return m_pixelsize;
}
void StormParams::setPixelSize(float pixelsize) {
    if (pixelsize != m_pixelsize) {
        m_pixelsize = pixelsize;
        m_pixelsizeSaved = false;
    }
}
bool StormParams::getPixelSizeSaved() const {
    return m_pixelsizeSaved;
}

unsigned int StormParams::getSkellamFrames() const {
    return m_skellamFrames;
}
void StormParams::setSkellamFrames(unsigned int frames) {
    if (frames != m_skellamFrames) {
        m_skellamFrames = frames;
        m_skellamFramesSaved = false;
    }
}
bool StormParams::getSkellamFramesSaved() const {
    return m_skellamFramesSaved;
}

unsigned int StormParams::getXYChunkSize() const {
    return m_xyChunkSize;
}
void StormParams::setXYChunkSize(unsigned int chunksize) {
    if (chunksize != m_xyChunkSize) {
        m_xyChunkSize = chunksize;
        m_xyChunkSizeSaved = false;
    }
}
bool StormParams::getXYChunkSizeSaved() const {
    return m_xyChunkSizeSaved;
}

unsigned int StormParams::getTChunkSize() const {
    return m_tChunkSize;
}
void StormParams::setTChunkSize(unsigned int chunksize) {
    if (chunksize != m_tChunkSize) {
        m_tChunkSize = chunksize;
        m_tChunkSizeSaved = false;
    }
}
bool StormParams::getTChunkSizeSaved() const {
    return m_tChunkSizeSaved;
}

unsigned int StormParams::getChunksInMemory() const {
    return m_chunksInMemory;
}
void StormParams::setChunksInMemory(unsigned int chunks) {
    if (chunks != m_chunksInMemory) {
        m_chunksInMemory = chunks;
        m_chunksInMemorySaved = false;
    }
}
bool StormParams::getChunksInMemorySaved() const {
    return m_chunksInMemorySaved;
}
const std::string& StormParams::getInFile() const {
    return m_infile;
}
void StormParams::setInFile(const std::string& infile, bool forceDefaults) {
    m_infile = infile;

    std::string extension = m_infile.substr( m_infile.find_last_of('.'));
    if(extension==".tif" || extension==".tiff") {
        m_type = TIFF;
        vigra::ImageImportInfo* info = new vigra::ImageImportInfo(m_infile.c_str());
        ptr = (void*) info;
        m_shape = Shape(info->width(), info->height(), info->numImages());
    }
    else if(extension==".sif") {
        m_type = SIF;
        vigra::SIFImportInfo* info = new vigra::SIFImportInfo (m_infile.c_str());
        ptr = (void*) info;
        //m_shape = Shape(info->shape()[0],info->shape()[1],100);
        m_shape = Shape(info->shape()[0],info->shape()[1],info->shape()[2]);
    }
    #ifdef HDF5_FOUND
    else if (extension==".h5" || extension==".hdf" || extension==".hdf5") {
        m_type = HDF5;
        vigra::HDF5File* h5file = new vigra::HDF5File(m_infile.c_str(), HDF5File::Open);
        ArrayVector<hsize_t> shape = h5file->getDatasetShape("/data");
        m_shape = Shape(shape[0],shape[1],shape[2]);
        ptr = (void*) h5file;
    }
    #endif // HDF5_FOUND
    else {
        vigra_precondition(false, "Wrong infile-extension given. Currently supported: .sif .h5 .hdf .hdf5 .tif .tiff");
    }

    if (forceDefaults)
        setDefaults();
    setDefaultFileNames(forceDefaults);
}

const std::string& StormParams::getOutFile() const {
    return m_outfile;
}
void StormParams::setOutFile(const std::string &file) {
    m_outfile = file;
}

const std::string& StormParams::getCoordsFile() const {
    return m_coordsfile;
}
void StormParams::setCoordsFile(const std::string &file) {
    m_coordsfile = file;
}

const std::string& StormParams::getSettingsFile() const {
    return m_settingsfile;
}
void StormParams::setSettingsFile(const std::string &file) {
    if (file != m_settingsfile) {
        m_settingsfile = file;
        m_config->clear();
        m_config->setConfigFile(m_settingsfile.c_str());
        load();
    }
}

const std::string& StormParams::getFrameRange() const {
    return m_frames;
}
void StormParams::setFrameRange(const std::string &frames) {
    if (frames != m_frames) {
        m_frames = frames;
        m_framesSaved = false;
    }
}

bool StormParams::getFrameRangeSaved() const {
    return m_framesSaved;
}

float StormParams::getAlpha() const {
    return m_alpha;
}
void StormParams::setAlpha(float alpha) {
    if (alpha != m_alpha) {
        m_alpha = alpha;
        m_thresholdMask = qnorm(m_alpha, 0, 1, 0, 0);
    }
}

double StormParams::getMaskThreshold() const {
    return m_thresholdMask;
}

void StormParams::setAsymmetryThreshold(float asymmetryThreshold) {
    m_asymmetryThreshold = asymmetryThreshold;
}

float StormParams::getAsymmetryThreshold() const
{
    return m_asymmetryThreshold;
}

void StormParams::setDoAsymmetryCheck(bool doAsymmetryCheck) {
    m_doAsymmetryCheck = doAsymmetryCheck;
}

bool StormParams::getDoAsymmetryCheck() const {
    return m_doAsymmetryCheck;
}

int StormParams::getMaxTChunkSize() const
{
    return m_maxTChunksize;
}

int StormParams::getMinTChunkSize() const
{
    return m_minTChunksize;
}

int StormParams::getMaxXyChunksize() const
{
    return m_maxXyChunksize;
}

int StormParams::getMinXyChunksize() const
{
    return m_minXyChunksize;
}

float StormParams::getMaxAsymmetryThreshold() const
{
    return m_maxAsymmetryThreshold;
}

float StormParams::getMinAsymmetryThreshold() const
{
    return m_minAsymmetryThreshold;
}

bool StormParams::getVerbose() const {
    return m_verbose;
}
void StormParams::setVerbose(bool verbose) {
    m_verbose = verbose;
}

const StormParams::Shape & StormParams::shape() const {
    return m_shape;
}
vigra::MultiArrayIndex StormParams::shape(const int dim) const {
    return m_shape[dim];
}
FileType StormParams::type() const { return m_type; };

const std::set<std::string>& StormParams::acceptedFileTypes()
{
    if (m_acceptedFileTypes.empty()) {
        m_acceptedFileTypes.insert(".tif");
        m_acceptedFileTypes.insert(".tiff");
        m_acceptedFileTypes.insert(".sif");
#ifdef HDF5_FOUND
        m_acceptedFileTypes.insert(".h5");
        m_acceptedFileTypes.insert(".hdf");
        m_acceptedFileTypes.insert(".hdf5");
#endif
    }
    return m_acceptedFileTypes;
}

void StormParams::printUsage() const {
	std::cout << "Usage: storm [Options] infile.sif [outfile.png]" << std::endl
	<< "Allowed Options: " << std::endl
	<< "  --help                 Print this help message" << std::endl
	<< "  -v or --verbose        verbose message output" << std::endl
	<< "  --factor=Arg           Resize factor equivalent to the subpixel-precision" << std::endl
	<< "  --cam-param-frames=Arg Number of frames to use for estimation of gain and offset." << std::endl
    << "                         Set to 0 to use the whole stack." << std::endl
	<< "  --coordsfile=Arg       filename for output of the found Coordinates" << std::endl
	<< "  --pixelsize=Arg        Pixel size in nanometers. If set, the coordinates" << std::endl
	<< "                         will be in nanometers, otherwise in pixels" << std::endl
	<< "  --filter=Arg           Text file with filter width (in pixels) for filtering in the" << std::endl
    << "                         FFT domain. If the file does not exist, generate a new filter" << std::endl
	<< "                         from the data" << std::endl
	<< "  --roi-len=Arg          size of the roi around maxima candidates" << std::endl
	<< "  --frames=Arg           run only on a subset of the stack (frames=start:end)" << std::endl
	<< "  --a                    enables check for asymmetry of points to skip asymmetric ones" <<std::endl
	<< "  --version              print version information and exit" << std::endl
	;
}

void StormParams::printVersion() const {
    std::cout << "simpleSTORM " << STORM_VERSION_STRING << std::endl
    << STORM_AUTHORS << std::endl
    << "This is free software; see the source for copying conditions.  There is NO" << std::endl
    << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl
    ;
}

/**
 * Defaults for unset variables are important
 */
void StormParams::setDefaults() {
    m_factor = 8;
    m_factorSaved = true;
    m_roilen = 9;
    m_roilenSaved = true;
    m_pixelsize = 1;
    m_pixelsizeSaved = true;
    m_skellamFrames = 200;
    m_skellamFramesSaved = true;
    m_xyChunkSize = 10;
    m_xyChunkSizeSaved = true;
    m_tChunkSize = 10;
    m_tChunkSizeSaved = true;
    m_chunksInMemory = 5;
    m_chunksInMemorySaved = true;
    m_thresholdMask = 0;
    setAlpha(0.001);
    setDoAsymmetryCheck(false);
    m_verbose = false;
}

void StormParams::setDefaultFileNames(bool force) {
    // defaults: save out- and coordsfile into the same folder as input stack
	size_t pos = m_infile.find_last_of('.');
    if (m_outfile.empty() || force) {
    	m_outfile = m_infile;
    	m_outfile.replace(pos, 255, ".png"); // replace extension
	}
    if (m_coordsfile.empty() || force) {
    	m_coordsfile = m_infile;
    	m_coordsfile.replace(pos, 255, ".txt"); // replace extension
	}
    if( m_settingsfile.empty() || force) {
        std::string tmp = m_infile;
    	tmp.replace(pos, 255, "_settings.txt"); // replace extension
        setSettingsFile(tmp);
	}
}

void StormParams::doSanityChecks() {
    if (m_skellamFrames > m_shape[2] || m_skellamFrames <= 0) {
        std::cout<<"selected value for skellamFrames ("<<m_skellamFrames<<") is either larger than the maximal stacksize ("<<m_shape[2]<<") or smaller than zero! It is set to: "<<m_shape[2]<<std::endl;
        m_skellamFrames = m_shape[2];}
    if (m_chunksInMemory > m_shape[2]) {
        m_chunksInMemory = m_shape[2];
    }
    if (m_chunksInMemory < 3) {
        std::cout<<"selected value for chunksInMemory ("<<m_chunksInMemory<<") is smaller than 3, it is therefore set to 3"<<std::endl;
        m_chunksInMemory = 3;}
    if (m_xyChunkSize > m_maxXyChunksize) {
        std::cout<<"selected value for xyChunkSize ("<<m_xyChunkSize<<") is smaller than 3, it is therefore set to 3"<<std::endl;
        m_xyChunkSize = m_maxXyChunksize;}
    if (m_xyChunkSize < m_minXyChunksize)
        m_xyChunkSize = m_minXyChunksize;
    if (m_tChunkSize < m_minTChunksize)
        m_tChunkSize = m_minTChunksize;
    if (m_tChunkSize >m_skellamFrames/m_chunksInMemory)
        m_tChunkSize = m_skellamFrames/m_chunksInMemory;

}

/**
 * Parse Options from argv into parameter-maps
 */
int StormParams::parseProgramOptions(int argc, char **argv)
{
	int c;
	int digit_optind = 0;

	while (1) {
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		static struct option long_options[] = {
			{"help",             no_argument,       0,  '?'},
			{"verbose",          no_argument,       0,  'v'},
			{"version",          no_argument,       0,  'V'},
			{"factor",           required_argument, 0,  'g'},
            {"cam-param-frames", required_argument, 0,  'P'},
			{"coordsfile",       required_argument, 0,  'c'},
            {"pixelsize",        required_argument, 0,  'p'},
			{"settings",         required_argument, 0,  's'},
			{"roi-len",          required_argument, 0,  'm'},
			{"frames",           required_argument, 0,  'F'},
            {"doAsymmetryCheck", required_argument, 0,  'a'},
			{0, 0, 0, 0 }

		};

		// valid options: "vc:" => -v option without parameter, c flag requires parameter
		c = getopt_long(argc, argv, "?vVg:P:t:c:s:p:m:F:a",
				long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case 'g': // factor
            setFactor(convertToLong(optarg));
			break;
        case 'P': // cam-param-frames
            setSkellamFrames(convertToULong(optarg));
            break;
		case 'm': // roi-len
			setRoilen(convertToLong(optarg));
			break;
        case 'p': // pixelsize
        	setPixelSize(convertToFloat(optarg));
        	break;
		case 'c': // coordsfile
            setCoordsFile(optarg);
			break;
		case 's': // settingsfile
			setSettingsFile(optarg);
			break;
		case 'F': // frames
			setFrameRange(optarg);
			break;
        case 'a':
            setDoAsymmetryCheck(true);
            break;
		case 'v':
			setVerbose(true); // verbose mode
			break;
		// Option -? and in case of unknown option or missing argument
		case '?':
			printUsage();
			std::exit(0);
		case 'V':
			printVersion();
			std::exit(0);
		default:
			std::cout << "?? getopt returned character code 0%o ??\n" << std::endl;
		}
	}

	if (optind < argc)
        setInFile(argv[optind++]);
    if (optind < argc)
        m_outfile = argv[optind++];
    if (optind < argc)
        std::cout << "unrecognized non-option Argument: " << argv[optind++] << std::endl;

	// if no input file given
	if(m_infile.empty()) {
		std::cerr << "error: no input file given" << std::endl;
		printUsage();
		std::exit(-1);
	}
	return 0;
}

template <class T>
bool readTVolume(MultiArrayView<STORMPARAMS_N, T> & array, FileType type, void *ptr) {
    switch(type) {
        case TIFF:
        {
            ImageImportInfo* info = reinterpret_cast<ImageImportInfo*>(ptr);
            vigra_precondition(array.size(2)==info->numImages(),"array shape and number of images in tiff file differ.");
            for(int i = 0; i < info->numImages(); ++i) {
                MultiArrayView <2, T> img = array.bindOuter(i);
                BasicImageView <T> v = makeBasicImageView(img);
                info->setImageIndex(i);
                vigra::importImage(*info, destImage(v));
            }
            return true;
        }
        #ifdef HDF5_FOUND
        case HDF5:
        {
            HDF5File* info = reinterpret_cast<HDF5File*>(ptr);
            info->read("/data", array);
            return true;
        }
        #endif // HDF5_FOUND
        default:
            return false;
    }
}

template <class T>
bool readTBlock(const StormParams::Shape& blockOffset,
               const StormParams::Shape& blockShape,
                MultiArrayView<STORMPARAMS_N, T> & array, FileType type,
               const StormParams::Shape &shape, void *ptr) {
    switch(type) {
        case TIFF:
        {
            ImageImportInfo* info = reinterpret_cast<ImageImportInfo*>(ptr);
            vigra_precondition(blockOffset[0]==0 && blockOffset[1]==0 &&
            blockShape[0]==shape[0] && blockShape[1]==shape[1],
            "for Tiff images only complete Frames are currently supported as ROIs");
            vigra_precondition(array.size(2)==blockShape[2],"array shape and number of images in ROI for tiff file differ.");
            vigra_precondition(blockShape[2] <= info->numImages(), "block shape larger than number of frames in the image");
            for(int i = 0; i < blockShape[2]; ++i) {
                MultiArrayView <2, T> img = array.bindOuter(i);
                BasicImageView <T> v = makeBasicImageView(img);
                info->setImageIndex(i+blockOffset[2]);
                vigra::importImage(*info, destImage(v));
            }
            return true;
        }
        #ifdef HDF5_FOUND
        case HDF5:
        {
            HDF5File* info = reinterpret_cast<HDF5File*>(ptr);
            info->readBlock("/data", blockOffset, blockShape, array);
            return true;
        }
        #endif // HDF5_FOUND
        default:
            return false;
    }
}

template <class  T>
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, T> & array) const {
    if (!readTVolume(array, m_type, ptr))
        vigra_fail("decoder for type not implemented.");
}



template <class  T>
void StormParams::readBlock(const StormParams::Shape& blockOffset,
               const StormParams::Shape& blockShape,
               MultiArrayView<STORMPARAMS_N, T> & array) const
{
    if (!readTBlock(blockOffset, blockShape, array, m_type, m_shape, ptr))
        vigra_fail("decoder for type not implemented.");
}



template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, int8_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, int16_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, int32_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, uint8_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, uint16_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, uint32_t>&) const;
template<>
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, float>& array) const {
    if (!readTVolume(array, m_type, ptr)) {
        if (m_type == SIF) {
            SIFImportInfo* info = reinterpret_cast<SIFImportInfo*>(ptr);
            readSIF(*info, array);
        } else
            vigra_fail("decoder for type not implemented.");
    }
}
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, double>&) const;

template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, int8_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, int16_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, int32_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, uint8_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, uint16_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, uint32_t>&) const;
template<>
void StormParams::readBlock(const StormParams::Shape& blockOffset,
               const StormParams::Shape& blockShape,
               MultiArrayView<STORMPARAMS_N, float>& array) const {
    if (!readTBlock(blockOffset, blockShape, array, m_type, m_shape, ptr)) {
        if (m_type == SIF) {
            SIFImportInfo* info = reinterpret_cast<SIFImportInfo*>(ptr);
            readSIFBlock(*info, blockOffset, blockShape, array);
        } else
            vigra_fail("decoder for type not implemented.");
    }
}
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, double>&) const;

void StormParams::save() const
{
    m_config->setSection(s_section.c_str());
    m_config->setIntValue("factor", m_factor);
    m_config->setIntValue("roilen", m_roilen);
    m_config->setDoubleValue("pixelsize", m_pixelsize);
    m_config->setIntValue("skellamFrames", m_skellamFrames);
    m_config->setIntValue("xyChunkSize", m_xyChunkSize);
    m_config->setIntValue("tChunkSize", m_tChunkSize);
    m_config->setIntValue("chunksInMemory", m_chunksInMemory);
    m_config->setDoubleValue("alpha", m_alpha);
    m_config->setBoolValue("doAsymmetryCheck", m_doAsymmetryCheck);
    m_config->save();
}

void StormParams::load(bool propagate)
{
    m_config->load();
    loadSettings(propagate);
}

void StormParams::loadSettings(bool)
{
    m_config->setSection(s_section.c_str());
    if (m_factorSaved && m_config->exists("factor"))
        m_factor = m_config->getIntValue("factor");
    else
        m_factorSaved = false;
    if (m_roilenSaved && m_config->exists("roilen"))
        m_roilenSaved = m_config->getIntValue("roilen");
    else
        m_roilenSaved = false;
    if (m_pixelsizeSaved && m_config->exists("pixelsize"))
        m_pixelsize = m_config->getDoubleValue("pixelsize");
    else
        m_pixelsizeSaved = false;
    if (m_skellamFramesSaved && m_config->exists("skellamFrames"))
        m_skellamFrames = m_config->getIntValue("skellamFrames");
    else
        m_skellamFramesSaved = false;
    if (m_xyChunkSizeSaved && m_config->exists("xyChunkSize"))
        m_xyChunkSize = m_config->getIntValue("xyChunkSize");
    else
        m_xyChunkSizeSaved = false;
    if (m_tChunkSizeSaved && m_config->exists("tChunkSize"))
        m_tChunkSize = m_config->getIntValue("tChunkSize");
    else
        m_tChunkSizeSaved = false;
    if (m_chunksInMemorySaved && m_config->exists("chunksInMemory"))
        m_chunksInMemory = m_config->getIntValue("chunksInMemory");
    else
        m_chunksInMemorySaved = false;
    if (m_config->exists("alpha"))
        m_alpha = m_config->getDoubleValue("alpha");
    if (m_doAsymmetryCheckSaved && m_config->exists("doAsymmetryCheck"))
        m_doAsymmetryCheck = m_config->getBoolValue("doAsymmetryCheck");
    else
        m_doAsymmetryCheckSaved = false;
}
