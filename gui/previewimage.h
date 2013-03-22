#ifndef PREVIEWIMAGE_H
#define PREVIEWIMAGE_H

#include "wienerStorm.hxx"

#include <chrono>
#include <vector>
#include <set>

#include <vigra/multi_array.hxx>

#include <QSet>
#include <QPair>
#include <QWidget>
#include <QPainter>

class GuiParams;

class PreviewImage : public QWidget
{
Q_OBJECT
public:
    PreviewImage(QWidget *parent = 0);
    ~PreviewImage();
    void setResults(const std::vector<std::set<Coord<float>>>*);
    void setUpdateInterval(const std::chrono::milliseconds&);
    void setParams(const GuiParams *params);
    void setIntensityScaleFactor(float);
    void saveImage(const QString &file);
    virtual QSize sizeHint() const;
    float scale() const;

public Q_SLOTS:
    void updateImage();

Q_SIGNALS:
    void detections(const QString&);
    void initialized();

protected:
    virtual void paintEvent(QPaintEvent *);
    virtual void resizeEvent(QResizeEvent *);
    virtual void moveEvent(QMoveEvent *);

private Q_SLOTS:
    void frameCompleted(int);
    void initializationFinished();

private:
    void updatePixmap(const QRect&);
    void init(const QRect&);

    std::chrono::steady_clock::time_point m_lastProcessed;
    std::chrono::milliseconds m_updateInterval;
    unsigned long int m_detections;
    float m_intensityScaleFactor;
    float m_scale;
    QRect m_geometry;
    QImage m_pixmap;
    QRect m_pixmapGeometry;
    QSize m_size;
    bool m_needRepaint;
    QPainter m_painter;
    const std::vector<std::set<Coord<float>>> *m_results; // have to be pointers in order to have
                                                          // constructor without arguments and use
                                                          // this in Qt Designer
    const GuiParams *m_params;
    bool m_initialized;
    bool m_needCompleteRepaint;
    std::vector<int> m_unprocessed;
    QSet<QPair<int, int>> m_toPaint;
    vigra::MultiArray<2, float> m_result;
    QSize m_resultSize;
    std::multiset<float> m_resultsForScaling;
    float m_sizeFactor;
    float m_intensityFactor;
    std::pair<float, float> m_limits;
};

#endif
