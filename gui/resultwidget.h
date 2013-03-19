#ifndef RESULTWIDGET_H
#define RESULTWIDGET_H

#include "wienerStorm.hxx"
#include "guiparams.h"

#include <atomic>

#include <QWidget>

namespace Ui
{
class ResultWidget;
}
class QtProgressFunctor;
class AnalysisWorker;

class ResultWidget : public QWidget
{
Q_OBJECT
public:
    ResultWidget(QWidget *parent = 0);
    ~ResultWidget();
    void start(const GuiParams&);
    bool canBeDeleted();
    bool canBeClosed();
    void waitFor();

Q_SIGNALS:
    void userAttentionRequired(QWidget*);
    void finished(ResultWidget*);

private:
    void saveCoordinates(const QString&);

    Ui::ResultWidget *m_ui;
    QtProgressFunctor *m_functor;
    AnalysisWorker *m_worker;
    std::vector<std::set<Coord<float>>> m_result;
    GuiParams m_params;
    bool m_coordinatesSaved;
    bool m_reconstructionSaved;

private Q_SLOTS:
    void stageChanged(int, const QString&);
    void workerFinished();
    void saveClicked();
};

#endif
