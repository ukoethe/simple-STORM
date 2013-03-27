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

#include "wienerStorm.hxx"

#include <vigra/impex.hxx>
#include "dataparams.h"
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#include <vigra/timing.hxx>

class CliProgressFunctor : public ProgressFunctor
{
public:
    CliProgressFunctor() : ProgressFunctor(), m_stacksize(0), m_frame(0) {};
    ~CliProgressFunctor() {};
    virtual void setStage(WienerStormStage stage)
    {
        m_stage = stage;
        switch (m_stage) {
            case CameraParameters:
                std::cout << std::endl << "Estimating camera gain and offset..." << std::endl;
                break;
            case PSFWidth:
                std::cout << std::endl << "Estimating PSF width..." << std::endl;
                break;
            case ParameterCheck:
                std::cout << std::endl << "Checking parameter..." << std::endl;
                break;
            case Localization:
                std::cout << std::endl << "Localizing molecules..." << std::endl;
                break;
        }
    }

    virtual void setStackSize(int stacksize)
    {
        m_stacksize = stacksize;
        m_frame = 0;
        helper::progress(-1, -1);
    }

    virtual void frameFinished(int frame)
    {
        ++m_frame;
#ifdef OPENMP_FOUND
        if(omp_get_thread_num()==0) { // master thread
            helper::progress(m_frame, m_stacksize);
        }
#else
        helper::progress(m_frame, m_stacksize);
#endif
    }
private:
    WienerStormStage m_stage;
    int m_stacksize;
    std::atomic<int> m_frame;
};

// MAIN
int main(int argc, char* argv[]) {
    try
    {
        DataParams params(argc, argv);
        if (!initR(argc, argv)) {
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

        CliProgressFunctor func;
        wienerStorm(params, res_coords, func);

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
