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

#ifndef __WIN__
#include <libgen.h>
#else
#include <stdio.h>
#include <cstdlib>
#endif

#include <vigra/impex.hxx>
#include <vigra/sifImport.hxx>
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#include "myimportinfo.h"

using namespace vigra;


MyImportInfo::MyImportInfo(int argc, char **argv)
{
	parseProgramOptions(argc, argv);
	setDefaults();

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
#ifndef __WIN__
    m_executableDir.append(dirname(argv[0]));
#else
    char drive[_MAX_DRIVE];
    char dir[_MAXDIR];
    _splitpath(argv0, drive, dir);
    m_executableDir.append(drive).append(dir);
#endif

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

typedef vigra::MultiArrayShape<MYIMPORT_N>::type Shape;
int MyImportInfo::getFactor() const {return m_factor;}
int MyImportInfo::getRoilen() const {return m_roilen;}
float MyImportInfo::getThreshold() const {return m_threshold;}
float MyImportInfo::getPixelsize() const {return m_pixelsize;}
const std::string& MyImportInfo::getInfile() const {return m_infile;}
const std::string& MyImportInfo::getOutfile() const {return m_outfile;}
const std::string& MyImportInfo::getCoordsfile() const {return m_coordsfile;}
const std::string& MyImportInfo::getFilterfile() const {return m_filterfile;}
const std::string& MyImportInfo::getFrameRange() const {return m_frames;}
char MyImportInfo::getVerbose() const {return verbose;}

const Shape & MyImportInfo::shape() const { return m_shape; }
vigra::MultiArrayIndex MyImportInfo::shape(const int dim) const { return m_shape[dim]; }
vigra::MultiArrayIndex MyImportInfo::shapeOfDimension(const int dim) const { return m_shape[dim]; }
FileType MyImportInfo::type() const { return m_type; };

const std::string& MyImportInfo::executableDir() const
{
    return m_executableDir;
}
void MyImportInfo::printUsage(const char* prog) {
	std::cout << "Usage: " << prog << " [Options] infile.sif [outfile.png]" << std::endl
	 << "Allowed Options: " << std::endl
	 << "  --help           Print this help message" << std::endl
	 <<	"  -v or --verbose  verbose message output" << std::endl
	 << "  --factor=Arg     Resize factor equivalent to the subpixel-precision" << std::endl
	 << "  --threshold=Arg  Threshold for background suppression" << std::endl
	 << "  --coordsfile=Arg filename for output of the found Coordinates" << std::endl
	 << "  --pixelsize=Arg  Pixel size in nanometers. If set, the coordinates" << std::endl
	 << "                   will be in nanometers, otherwise in pixels" << std::endl
	 << "  --filter=Arg     tif input for filtering in fft domain. If the file" << std::endl
	 << "                   does not exist, generate a new filter from the data" << std::endl
	 << "  --roi-len=Arg    size of the roi around maxima candidates" << std::endl
	 << "  --frames=Arg     run only on a subset of the stack (frames=start:end)" << std::endl
	 << "  --version        print version information and exit" << std::endl
	 ;
}

/**
 * Defaults for unset variables are important
 */
void MyImportInfo::setDefaults() {
    // defaults:
	m_factor = 8;
	m_roilen = 9;
	m_threshold = 250;
	m_pixelsize = 100;

    // defaults: save out- and coordsfile into the same folder as input stack
	size_t pos = m_infile.find_last_of('.');
    if(m_outfile=="") {
    	m_outfile = m_infile;
    	m_outfile.replace(pos, 255, ".png"); // replace extension
	}
    if(m_coordsfile=="") {
    	m_coordsfile = m_infile;
    	m_coordsfile.replace(pos, 255, ".txt"); // replace extension
	}
    if(m_filterfile=="") {
    	m_filterfile = m_infile;
    	m_filterfile.replace(pos, 255, "_filter.txt"); // replace extension
	}

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
			{"help",     no_argument, 0,  '?' },
			{"verbose",     no_argument, 0,  'v' },
			{"version",     no_argument, 0,  'V' },
			{"factor",  required_argument,       0,  'g' },
			{"threshold",  required_argument, 0,  't' },
			{"coordsfile",    required_argument, 0,  'c'},
            {"pixelsize",    required_argument, 0,  'p'},
			{"filter",    required_argument, 0,  'f' },
			{"roi-len",    required_argument, 0,  'm' },
			{"frames",    required_argument, 0,  'F' },
			{0,         0,                 0,  0 }

		};

		// valid options: "vc:" => -v option without parameter, c flag requires parameter
		c = getopt_long(argc, argv, "?vVt:c:f:p:m:F:",
				long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case 't': // threshold
			m_threshold = convertToDouble(optarg);
			break;
		case 'g': // factor
			m_factor = convertToDouble(optarg);
			break;
		case 'm': // roi-len
			m_roilen = convertToDouble(optarg);
			break;
        case 'p': // pixelsize
        	m_pixelsize = convertToDouble(optarg);
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
			printUsage(argv[0]);
			return -1;
			break;

		case 'V':
			std::cout << "STORM analysis software version " << versionString() << std::endl
			 << "Copyright (C) 2011 Joachim Schleicher and Ullrich Koethe" << std::endl
			 << "This is free software; see the source for copying conditions.  There is NO" << std::endl
			 << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl
			 ;
			return -1;
			break;

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
		printUsage(argv[0]);
		return -1;
	}

	setDefaults();

	return 0;
}


