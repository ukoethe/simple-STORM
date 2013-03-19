#include "analysisworker.h"

QtProgressFunctor::QtProgressFunctor(QObject *parent)
: QObject(parent), ProgressFunctor(), m_stacksize(0), m_frame(0) {}

QtProgressFunctor::~QtProgressFunctor(){}

void QtProgressFunctor::setStage(WienerStormStage stage)
{
    m_stage = stage;
    m_frame = 0;
    switch (m_stage) {
        case CameraParameters:
            emit stageChanged(stage, "Estimating camera gain and offset: %v / %m (%p%)");
            break;
        case PSFWidth:
            emit stageChanged(stage, "Estimating PSF width: %v / %m (%p%)");
            break;
        case Localization:
            emit stageChanged(stage, "Localizing molecules: %v / %m (%p%)");
            break;
    }
}

WienerStormStage QtProgressFunctor::getStage() const
{
    return m_stage;
}

void QtProgressFunctor::setStackSize(int stacksize)
{
    m_stacksize = stacksize;
    emit stackSizeChanged(m_stacksize);
}

void QtProgressFunctor::frameFinished(int frame)
{
    ++m_frame;
    emit progress(m_frame);
    emit frameCompleted(frame);
}

void QtProgressFunctor::abort()
{
    ProgressFunctor::abort();
}

AnalysisWorker::AnalysisWorker(QtProgressFunctor &functor, DataParams &params, std::vector<std::set<Coord<float>>> &res, QObject *parent)
: QThread(parent), m_functor(functor), m_params(params), m_res(res)
{}

AnalysisWorker::~AnalysisWorker() {}

void AnalysisWorker::run()
{
    wienerStorm(m_params, m_res, m_functor);
}
