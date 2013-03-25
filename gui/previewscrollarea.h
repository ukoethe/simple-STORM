#ifndef PREVIEWSCROLLAREA_H
#define PREVIEWSCROLLAREA_H

#include <QScrollArea>

class PreviewScrollArea : public QScrollArea
{
Q_OBJECT
public:
    PreviewScrollArea(QWidget *parent = 0);
    ~PreviewScrollArea();
    float scale() const;

public Q_SLOTS:
    void scrollBy(int, int);
    void setMinScale(float);
    void setMaxScale(float);
    void setScale(float);
    void setScaleStep(float);
    void autoZoom();
    void zoomIn();
    void zoomOut();

Q_SIGNALS:
    void scaleChanged(float, bool);

protected:
    virtual void mousePressEvent(QMouseEvent*);
    virtual void mouseMoveEvent(QMouseEvent*);
    virtual void mouseReleaseEvent(QMouseEvent*);
    virtual void wheelEvent(QWheelEvent*);
    virtual void resizeEvent(QResizeEvent*);

private Q_SLOTS:
    void horizontalScrollBarRangeChanged(int, int);
    void verticalScrollBarRangeChanged(int, int);

private:
    bool hasHorizontalScrollBar();
    bool hasVerticalScrollBar();
    void zoomIn(const QPoint&);
    void zoomOut(const QPoint&);
    void setScale(float, const QPoint&);
    int m_dx;
    int m_dy;
    bool m_dragging;
    QPoint m_lastDragPoint;
    float m_scale;
    float m_minScale;
    float m_maxScale;
    float m_scaleStep;
};

#endif
