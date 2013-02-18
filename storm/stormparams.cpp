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
#include <iostream>
#ifndef __WIN__
#include <libgen.h>
#else
#include <cstdlib>
#endif
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

#include "stormparams.h"

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

StormParams::StormParams(int argc, char **argv)
: m_factor(8), m_factorSaved(true), m_roilen(9), m_roilenSaved(true), m_threshold(250),
  m_thresholdSaved(true), m_pixelsize(1), m_pixelsizeSaved(true), m_skellamFrames(200),
  m_skellamFramesSaved(true), m_xyChunkSize(10), m_xyChunkSizeSaved(true),
  m_tChunkSize(10), m_tChunkSizeSaved(true), m_chunksInMemory(5), m_chunksInMemorySaved(true),
  m_verbose(false) {
#ifndef __WIN__
    m_executableDir.append(dirname(argv[0]));
    m_executableName.append(basename(argv[0]));
#else
    char drive[_MAX_DRIVE];
    char dir[_MAXDIR];
    char fname[ _MAX_FNAME];
    char ext[_MAX_EXT];
    _splitpath(argv0, drive, dir, fname, ext);
    m_executableDir.append(drive).append(dir);
    m_executableName.append(fname).append(ext);
#endif

    parseProgramOptions(argc, argv);

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

    setDefaults();
    m_config = new rude::Config();
    m_config->setConfigFile(m_settingsfile.c_str());
    load();
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
bool StormParams::getFactorSaved() const {
    return m_factorSaved;
}
int StormParams::getRoilen() const {
    return m_roilen;
}
bool StormParams::getRoilenSaved() const {
    return m_roilenSaved;
}
float StormParams::getThreshold() const {
    return m_threshold;
}
bool StormParams::getThresholdSaved() const {
    return m_thresholdSaved;
}
float StormParams::getPixelSize() const {
    return m_pixelsize;
}
bool StormParams::getPixelSizeSaved() const {
    return m_pixelsizeSaved;
}
unsigned int StormParams::getSkellamFrames() const {
    return m_skellamFrames;
}
bool StormParams::getSkellamFramesSaved() const {
    return m_skellamFramesSaved;
}
unsigned int StormParams::getXYChunkSize() const {
    return m_xyChunkSize;
}
bool StormParams::getXYChunkSizeSaved() const {
    return m_xyChunkSizeSaved;
}
unsigned int StormParams::getTChunkSize() const {
    return m_tChunkSize;
}
bool StormParams::getTChunkSizeSaved() const {
    return m_tChunkSizeSaved;
}
unsigned int StormParams::getChunksInMemory() const {
    return m_chunksInMemory;
}
bool StormParams::getChunksInMemorySaved() const {
    return m_chunksInMemorySaved;
}
const std::string& StormParams::getInFile() const {
    return m_infile;
}
const std::string& StormParams::getOutFile() const {
    return m_outfile;
}
const std::string& StormParams::getCoordsFile() const {
    return m_coordsfile;
}
const std::string& StormParams::getSettingsFile() const {
    return m_settingsfile;
}
const std::string& StormParams::getFrameRange() const {
    return m_frames;
}
bool StormParams::getFrameRangeSaved() const {
    return m_framesSaved;
}
bool StormParams::getVerbose() const {
    return m_verbose;
}
const StormParams::Shape & StormParams::shape() const {
    return m_shape;
}
vigra::MultiArrayIndex StormParams::shape(const int dim) const {
    return m_shape[dim];
}
FileType StormParams::type() const { return m_type; };

const std::string& StormParams::executableDir() const {
    return m_executableDir;
}

void StormParams::printUsage() const {
	std::cout << "Usage: " << m_executableName << " [Options] infile.sif [outfile.png]" << std::endl
	<< "Allowed Options: " << std::endl
	<< "  --help                 Print this help message" << std::endl
	<< "  -v or --verbose        verbose message output" << std::endl
	<< "  --factor=Arg           Resize factor equivalent to the subpixel-precision" << std::endl
	<< "  --cam-param-frames=Arg Number of frames to use for estimation of gain and offset." << std::endl
    << "                         Set to 0 to use the whole stack." << std::endl
	<< "  --threshold=Arg        Threshold for background suppression" << std::endl
	<< "  --coordsfile=Arg       filename for output of the found Coordinates" << std::endl
	<< "  --pixelsize=Arg        Pixel size in nanometers. If set, the coordinates" << std::endl
	<< "                         will be in nanometers, otherwise in pixels" << std::endl
	<< "  --filter=Arg           Text file with filter width (in pixels) for filtering in the" << std::endl
    << "                         FFT domain. If the file does not exist, generate a new filter" << std::endl
	<< "                         from the data" << std::endl
	<< "  --roi-len=Arg          size of the roi around maxima candidates" << std::endl
	<< "  --frames=Arg           run only on a subset of the stack (frames=start:end)" << std::endl
	<< "  --version              print version information and exit" << std::endl
	;
}

void StormParams::printVersion() const {
    std::cout << "STORM analysis software version " << versionString() << std::endl
    << "Copyright (C) 2011 Joachim Schleicher and Ullrich Koethe" << std::endl
    << "This is free software; see the source for copying conditions.  There is NO" << std::endl
    << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl
    ;
}

/**
 * Defaults for unset variables are important
 */
void StormParams::setDefaults() {
    // defaults: save out- and coordsfile into the same folder as input stack
	size_t pos = m_infile.find_last_of('.');
    if (m_outfile.empty()) {
    	m_outfile = m_infile;
    	m_outfile.replace(pos, 255, ".png"); // replace extension
	}
    if (m_coordsfile.empty()) {
    	m_coordsfile = m_infile;
    	m_coordsfile.replace(pos, 255, ".txt"); // replace extension
	}
    if( m_settingsfile.empty()) {
    	m_settingsfile = m_infile;
    	m_settingsfile.replace(pos, 255, "_settings.txt"); // replace extension
	}
    if (m_skellamFrames > m_shape[2] || m_skellamFrames <= 0)
        m_skellamFrames = m_shape[2];
    if (m_chunksInMemory < 4)
        m_chunksInMemory = 3;
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
			{"threshold",        required_argument, 0,  't'},
			{"coordsfile",       required_argument, 0,  'c'},
            {"pixelsize",        required_argument, 0,  'p'},
			{"filter",           required_argument, 0,  'f'},
			{"roi-len",          required_argument, 0,  'm'},
			{"frames",           required_argument, 0,  'F'},
			{0, 0, 0, 0 }

		};

		// valid options: "vc:" => -v option without parameter, c flag requires parameter
		c = getopt_long(argc, argv, "?vVg:P:t:c:f:p:m:F:",
				long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case 't': // threshold
			m_threshold = convertToFloat(optarg);
            m_thresholdSaved = false;
			break;
		case 'g': // factor
			m_factor = convertToLong(optarg);
            m_factorSaved = false;
			break;
        case 'P': // cam-param-frames
            m_skellamFrames = convertToULong(optarg);
            m_skellamFramesSaved = false;
            break;
		case 'm': // roi-len
			m_roilen = convertToLong(optarg);
            m_roilenSaved = false;
			break;
        case 'p': // pixelsize
        	m_pixelsize = convertToFloat(optarg);
            m_pixelsizeSaved = false;
        	break;
		case 'c': // coordsfile
			m_coordsfile = optarg;
			break;
		case 'f': // filter
			m_settingsfile = optarg;
            m_pixelsizeSaved = false;
			break;
		case 'F': // frames
			m_frames = optarg;
            m_framesSaved = false;
			break;
		case 'v':
			m_verbose = true; // verbose mode
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

	while (optind < argc) {
		if(m_infile == "") m_infile = argv[optind++];
		else if (m_outfile == "") m_outfile = argv[optind++];
		else std::cout << "unrecognized non-option Argument: " << argv[optind++] << std::endl;
	}

	// if no input file given
	if(m_infile == "" ) {
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
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, unsigned int8_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, unsigned int16_t>&) const;
template
void StormParams::readVolume(MultiArrayView<STORMPARAMS_N, unsigned int32_t>&) const;
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
               MultiArrayView<STORMPARAMS_N, unsigned int8_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, unsigned int16_t>&) const;
template
void StormParams::readBlock(const StormParams::Shape&,
               const StormParams::Shape&,
               MultiArrayView<STORMPARAMS_N, unsigned int32_t>&) const;
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
    m_config->setDoubleValue("threshold", m_threshold);
    m_config->setDoubleValue("pixelsize", m_pixelsize);
    m_config->setIntValue("skellamFrames", m_skellamFrames);
    m_config->setIntValue("xyChunkSize", m_xyChunkSize);
    m_config->setIntValue("tChunkSize", m_tChunkSize);
    m_config->setIntValue("chunksInMemory", m_chunksInMemory);
    m_config->save();
}

void StormParams::load()
{
    m_config->load();
    m_config->setSection(s_section.c_str());
    if (m_factorSaved && m_config->exists("factor"))
        m_factor = m_config->getIntValue("factor");
    else
        m_factorSaved = false;
    if (m_roilenSaved && m_config->exists("roilen"))
        m_roilenSaved = m_config->getIntValue("roilen");
    else
        m_roilenSaved = false;
    if (m_thresholdSaved && m_config->exists("threshold"))
        m_threshold = m_config->getDoubleValue("threshold");
    else
        m_thresholdSaved = false;
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
}
