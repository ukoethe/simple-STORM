#ifndef ANALYSISWORKER_H
#define ANALYSISWORKER_H

#include "dataparams.h"
#include "wienerStorm.hxx"

#include <QThread>

class QtProgressFunctor : public QObject, public ProgressFunctor
{
Q_OBJECT
public:
    QtProgressFunctor(QObject *parent = 0);
    ~QtProgressFunctor();
    virtual void setStage(WienerStormStage stage);
    WienerStormStage getStage() const;
    virtual void setStackSize(int stacksize);
    virtual void frameFinished(int frame);

public Q_SLOTS:
    virtual void abort();

Q_SIGNALS:
    void stageChanged(int, const QString&);
    void stackSizeChanged(int);
    void progress(int);
    void frameCompleted(int);

private:
    WienerStormStage m_stage;
    int m_stacksize;
    std::atomic<int> m_frame;
};

class AnalysisWorker : public QThread
{
Q_OBJECT
public:
    AnalysisWorker(QtProgressFunctor&, DataParams&, std::vector<std::set<Coord<float>>>&, QObject *parent);
    ~AnalysisWorker();

protected:
    void run();
private:
    QtProgressFunctor &m_functor;
    DataParams &m_params;
    std::vector<std::set<Coord<float>>> &m_res;
};

#endif
