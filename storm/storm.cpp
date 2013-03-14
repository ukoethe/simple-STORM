/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*         Copyright 2010 by Joachim Schleicher and Ullrich Koethe      */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*    joachim.schleicher@iwr.uni-heidelberg.de                          */
/************************************************************************/

#define CHUNKSIZE 10

#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>

#include <cstdlib>
#ifdef __WIN__
#include <process.h>
#else
#include <unistd.h>
#endif

#define R_INTERFACE_PTRS
#define R_NO_REMAP
#include <Rembedded.h>
#include <Rinterface.h>

#include "wienerStorm.hxx"
#include "configVersion.hxx"

#include <vigra/impex.hxx>
#include "dataparams.h"
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#include <vigra/timing.hxx>

// MAIN
int main(int argc, char* argv[]) {
    try
    {
        DataParams params(argc, argv);
        if (!initR(params, argc, argv)) {
            std::cerr << "Could not initialize R" << std::endl
            << "You probably do not have R installed or do not have it in your PATH." << std::endl;
            return 1;
        }

        //~ in.reshape(info.shape());
        //~ readVolume(info, in);
        int stacksize = params.shape(2);
        Size2D size2 (params.shape(0), params.shape(1)); // isnt' there a slicing operator?


        if(params.getVerbose()) {
            std::cout << "Images with Shape: " << params.shape() << std::endl;
            std::cout << "Processing a stack of " << stacksize << " images..." << std::endl;
        }


        // found spots. One Vector over all images in stack
        // the inner set contains all spots in the image
        std::vector<std::set<Coord<float> > > res_coords(stacksize);

        DImage res((size2-Diff2D(1,1))*params.getFactor()+Diff2D(1,1));
        // check if outfile is writable, otherwise throw error -> exit
        exportImage(srcImageRange(res), ImageExportInfo(params.getOutFile().c_str()));
        std::ofstream cf;
        if(!params.getCoordsFile().empty()) {
            cf.open(params.getCoordsFile());
            vigra_precondition(cf.is_open(), "Could not open coordinate-file for writing.");
        }

        USE_NESTED_TICTOC;
        //USETICTOC;
        TICPUSH;  // measure the time

        // STORM Algorithmut

        wienerStorm(params, res_coords);
        std::cout<<"a: "<<params.getSlope()<<" b: "<<params.getIntercept() << " sigma: " << params.getSigma()<<std::endl;

        // resulting image
        drawCoordsToImage<Coord<float> >(res_coords, res);

        int numSpots = 0;
        if(cf.is_open()) {
            numSpots = saveCoordsFile(params, cf, res_coords);
            cf.close();
        }

        // end: done.
        std::cout << std::endl << TOCS << std::endl;
        std::cout << "detected " << numSpots << " spots." << std::endl;

        // some maxima are very strong so we scale the image as appropriate :
        double maxlim = 0., minlim = 0;
        findMinMaxPercentile(res, 0., minlim, 0.996, maxlim);
        std::cout << "cropping output values to range [" << minlim << ", " << maxlim << "]" << std::endl;
        if(maxlim > minlim) {
            transformImage(srcImageRange(res), destImage(res), ifThenElse(Arg1()>Param(maxlim), Param(maxlim), Arg1()));
        }
        exportImage(srcImageRange(res), ImageExportInfo(params.getOutFile().c_str()));
        if (!params.getSettingsFile().empty())
            params.save();


    }
    catch (vigra::StdException & e)
    {
        std::cout<<"There was an error:"<<std::endl;
        std::cout << e.what() << std::endl;
        return 1;
    }

    endR();

    return 0;
}
