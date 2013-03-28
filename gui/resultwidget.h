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
    QString makeToolTip();

Q_SIGNALS:
    void userAttentionRequired(QWidget*);
    void finished(ResultWidget*);
    void toolTip(const QString&, ResultWidget*);

protected:
    virtual void resizeEvent(QResizeEvent*);

private:
    void saveCoordinates(const QString&);

    Ui::ResultWidget *m_ui;
    QtProgressFunctor *m_functor;
    AnalysisWorker *m_worker;
    std::vector<std::set<Coord<float>>> m_result;
    GuiParams m_params;
    bool m_coordinatesSaved;
    bool m_reconstructionSaved;
    bool m_autoZoom;
    bool m_previewEnabled;
    int m_minZoom;
    int m_maxZoom;
    int m_zoomStep;
    int m_zoomPageStep;

private Q_SLOTS:
    void stageChanged(int, const QString&);
    void workerFinished();
    void saveClicked();
    void scaleChanged(float, bool);
    void scaleChanged(int, bool);
    void zoomSliderValueChanged(int);
    void zoomSliderActionTriggered(int);
    void zoomValueChanged();
    void setPreviewEnabled(bool);
};

#endif
