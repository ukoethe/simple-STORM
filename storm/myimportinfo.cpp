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

#include "myimportinfo.h"

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

MyImportInfo::MyImportInfo(int argc, char **argv)
: m_factor(8), m_roilen(9), m_threshold(250), m_pixelsize(1), m_skellamFrames(200) {
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
}

MyImportInfo::~MyImportInfo() {
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
}

int MyImportInfo::getFactor() const {return m_factor;}
int MyImportInfo::getRoilen() const {return m_roilen;}
float MyImportInfo::getThreshold() const {return m_threshold;}
float MyImportInfo::getPixelsize() const {return m_pixelsize;}
unsigned int MyImportInfo::getSkellamFrames() const {return m_skellamFrames;}
const std::string& MyImportInfo::getInfile() const {return m_infile;}
const std::string& MyImportInfo::getOutfile() const {return m_outfile;}
const std::string& MyImportInfo::getCoordsfile() const {return m_coordsfile;}
const std::string& MyImportInfo::getFilterfile() const {return m_filterfile;}
const std::string& MyImportInfo::getFrameRange() const {return m_frames;}
char MyImportInfo::getVerbose() const {return verbose;}

const MyImportInfo::Shape & MyImportInfo::shape() const { return m_shape; }
vigra::MultiArrayIndex MyImportInfo::shape(const int dim) const { return m_shape[dim]; }
FileType MyImportInfo::type() const { return m_type; };

const std::string& MyImportInfo::executableDir() const
{
    return m_executableDir;
}

void MyImportInfo::printUsage() const {
	std::cout << "Usage: " << m_executableName << " [Options] infile.sif [outfile.png]" << std::endl
	<< "Allowed Options: " << std::endl
	<< "  --help                 Print this help message" << std::endl
	<<	"  -v or --verbose        verbose message output" << std::endl
	<< "  --factor=Arg           Resize factor equivalent to the subpixel-precision" << std::endl
	<< "  --cam-param-frames=Arg Number of frames to use for estimation of gain and offset" << std::endl
	<< "  --threshold=Arg        Threshold for background suppression" << std::endl
	<< "  --coordsfile=Arg       filename for output of the found Coordinates" << std::endl
	<< "  --pixelsize=Arg        Pixel size in nanometers. If set, the coordinates" << std::endl
	<< "                         will be in nanometers, otherwise in pixels" << std::endl
	<< "  --filter=Arg           Text file with filter width (in pixels) for filtering in the" << std::endl
    << "                        FFT domain. If the file does not exist, generate a new filter" << std::endl
	<< "                        from the data" << std::endl
	<< "  --roi-len=Arg         size of the roi around maxima candidates" << std::endl
	<< "  --frames=Arg          run only on a subset of the stack (frames=start:end)" << std::endl
	<< "  --version             print version information and exit" << std::endl
	;
}

void MyImportInfo::printVersion() const {
    std::cout << "STORM analysis software version " << versionString() << std::endl
    << "Copyright (C) 2011 Joachim Schleicher and Ullrich Koethe" << std::endl
    << "This is free software; see the source for copying conditions.  There is NO" << std::endl
    << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl
    ;
}

/**
 * Defaults for unset variables are important
 */
void MyImportInfo::setDefaults() {
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
    if( m_filterfile.empty()) {
    	m_filterfile = m_infile;
    	m_filterfile.replace(pos, 255, "_filter.txt"); // replace extension
	}
    if (m_skellamFrames > m_shape[2])
        m_skellamFrames = m_shape[2];
}

/**
 * Parse Options from argv into parameter-maps
 */
int MyImportInfo::parseProgramOptions(int argc, char **argv)
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
			break;
		case 'g': // factor
			m_factor = convertToLong(optarg);
			break;
        case 'P': // cam-param-frames
            m_skellamFrames =convertToULong(optarg);
            break;
		case 'm': // roi-len
			m_roilen = convertToLong(optarg);
			break;
        case 'p': // pixelsize
        	m_pixelsize = convertToFloat(optarg);
        	break;
		case 'c': // coordsfile
			m_coordsfile = optarg;
			break;
		case 'f': // filter
			m_filterfile = optarg;
			break;
		case 'F': // frames
			m_frames = optarg;
			break;
		case 'v':
			verbose = 1; // verbose mode
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
bool readTVolume(MultiArrayView<MYIMPORT_N, T> & array, FileType type, void *ptr) {
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
bool readTBlock(const MyImportInfo::Shape& blockOffset,
               const MyImportInfo::Shape& blockShape,
               MultiArrayView<MYIMPORT_N, T> & array, FileType type,
               const MyImportInfo::Shape &shape, void *ptr) {
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
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, T> & array) const {
    if (!readTVolume(array, m_type, ptr))
        vigra_fail("decoder for type not implemented.");
}



template <class  T>
void MyImportInfo::readBlock(const MyImportInfo::Shape& blockOffset,
               const MyImportInfo::Shape& blockShape,
               MultiArrayView<MYIMPORT_N, T> & array) const
{
    if (!readTBlock(blockOffset, blockShape, array, m_type, m_shape, ptr))
        vigra_fail("decoder for type not implemented.");
}



template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, int8_t>&) const;
template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, int16_t>&) const;
template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, int32_t>&) const;
template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, unsigned int8_t>&) const;
template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, unsigned int16_t>&) const;
template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, unsigned int32_t>&) const;
template<>
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, float>& array) const {
    if (!readTVolume(array, m_type, ptr)) {
        if (m_type == SIF) {
            SIFImportInfo* info = reinterpret_cast<SIFImportInfo*>(ptr);
            readSIF(*info, array);
        } else
            vigra_fail("decoder for type not implemented.");
    }
}
template
void MyImportInfo::readVolume(MultiArrayView<MYIMPORT_N, double>&) const;

template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, int8_t>&) const;
template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, int16_t>&) const;
template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, int32_t>&) const;
template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, unsigned int8_t>&) const;
template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, unsigned int16_t>&) const;
template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, unsigned int32_t>&) const;
template<>
void MyImportInfo::readBlock(const MyImportInfo::Shape& blockOffset,
               const MyImportInfo::Shape& blockShape,
               MultiArrayView<MYIMPORT_N, float>& array) const {
    if (!readTBlock(blockOffset, blockShape, array, m_type, m_shape, ptr)) {
        if (m_type == SIF) {
            SIFImportInfo* info = reinterpret_cast<SIFImportInfo*>(ptr);
            readSIFBlock(*info, blockOffset, blockShape, array);
        } else
            vigra_fail("decoder for type not implemented.");
    }
}
template
void MyImportInfo::readBlock(const MyImportInfo::Shape&,
               const MyImportInfo::Shape&,
               MultiArrayView<MYIMPORT_N, double>&) const;
