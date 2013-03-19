#ifndef PREVIEWWIDGET_H
#define PREVIEWWIDGET_H

#include "wienerStorm.hxx"
#include "guiparams.h"

#include <chrono>
#include <set>

#include <vigra/multi_array.hxx>

#include <QWidget>
#include <QPainter>

namespace Ui
{
class PreviewWidget;
}
class QResizeEvent;

class PreviewWidget : public QWidget
{
Q_OBJECT
public:
    PreviewWidget(QWidget *parent = 0);
    ~PreviewWidget();
    void setResults(const std::vector<std::set<Coord<float>>>*);
    void setUpdateInterval(const std::chrono::milliseconds&);
    void setParams(const GuiParams *params);
    void saveImage(const QString &file);

Q_SIGNALS:
    void detections(const QString&);

public Q_SLOTS:
    void frameCompleted(int);
    void update();

private Q_SLOTS:
    void autoZoom();
    void zoomInClicked();
    void zoomOutClicked();
    void zoomSliderValueChanged(int);
    void zoomValueChanged();

protected:
    virtual void resizeEvent(QResizeEvent*);

private:
    Ui::PreviewWidget *m_ui;
    std::chrono::steady_clock::time_point m_lastProcessed;
    std::chrono::milliseconds m_updateInterval;
    unsigned long int m_detections;
    bool m_paramsSet;
    float m_zoom;
    const std::vector<std::set<Coord<float>>> *m_results; // have to be pointers in order to have
                                                          // constructor without arguments and use
                                                          // this in Qt Designer
    std::vector<int> m_unprocessed;
    vigra::MultiArray<2, double> m_result;
    std::multiset<double> m_resultsForScaling;
    size_t m_quantileIndex;
    const GuiParams *m_params;
    float m_factor;
    QPixmap m_pixmap;
    QPainter m_painter;
};

#endif
