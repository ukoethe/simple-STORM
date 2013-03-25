#include "previewscrollarea.h"

#include <QScrollBar>
#include <QMouseEvent>
#include <QWheelEvent>
#include <iostream>
PreviewScrollArea::PreviewScrollArea(QWidget *parent)
: QScrollArea(parent), m_dragging(false), m_scale(1), m_minScale(0.1), m_maxScale(5.), m_scaleStep(0.1)
{
    connect(horizontalScrollBar(), SIGNAL(rangeChanged(int, int)), this, SLOT(horizontalScrollBarRangeChanged(int, int)));
    connect(verticalScrollBar(), SIGNAL(rangeChanged(int, int)), this, SLOT(verticalScrollBarRangeChanged(int, int)));
}

PreviewScrollArea::~PreviewScrollArea()
{}

void PreviewScrollArea::scrollBy(int dx, int dy)
{
    horizontalScrollBar()->setValue(horizontalScrollBar()->value() + dx);
    verticalScrollBar()->setValue(verticalScrollBar()->value() + dy);
}

void PreviewScrollArea::horizontalScrollBarRangeChanged(int min, int max)
{
    horizontalScrollBar()->setValue(horizontalScrollBar()->value() + m_dx);
    m_dx = 0;
    if (max - min > 0)
        setCursor(Qt::OpenHandCursor);
    else
        setCursor(Qt::ArrowCursor);
}

void PreviewScrollArea::verticalScrollBarRangeChanged(int min, int max)
{
    verticalScrollBar()->setValue(verticalScrollBar()->value() + m_dy);
    m_dy = 0;
    if (max - min > 0)
        setCursor(Qt::OpenHandCursor);
    else
        setCursor(Qt::ArrowCursor);
}

void PreviewScrollArea::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton && (hasHorizontalScrollBar() || hasVerticalScrollBar())) {
        m_lastDragPoint = event->pos();
        m_dragging = true;
        setCursor(Qt::ClosedHandCursor);
        event->accept();
    } else
        QScrollArea::mousePressEvent(event);
}

void PreviewScrollArea::mouseMoveEvent(QMouseEvent *event)
{
    QScrollArea::mouseMoveEvent(event);
    if (m_dragging && (event->buttons() & Qt::LeftButton)) {
        QPoint p = (m_lastDragPoint - event->pos());
        scrollBy(p.x(), p.y());
        m_lastDragPoint = event->pos();
        event->accept();
    }
}

void PreviewScrollArea::mouseReleaseEvent(QMouseEvent *event)
{
    if (m_dragging && event->button() == Qt::LeftButton) {
        QPoint p = event->pos() - m_lastDragPoint;
        scrollBy(p.x(), p.y());
        m_dragging = false;
        setCursor(Qt::OpenHandCursor); // don't need to check for scroll bars here, as drag events are only
                                       // accepted when scroll bars are present
        event->accept();
    } else
        QScrollArea::mouseReleaseEvent(event);
}

void PreviewScrollArea::wheelEvent(QWheelEvent *event)
{
    if (event->modifiers().testFlag(Qt::ControlModifier)) {
        event->accept();
        if (event->delta() > 0)
            zoomIn(event->pos() - widget()->pos());
        else
            zoomOut(event->pos() - widget()->pos());
    } else
        QScrollArea::wheelEvent(event);
}

void PreviewScrollArea::resizeEvent(QResizeEvent *event)
{
    QScrollArea::resizeEvent(event);
    QSize s = widget()->size();
    m_scale = s.width() /(float)widget()->baseSize().width();
    if (!hasHorizontalScrollBar() && !hasVerticalScrollBar() && (s.width() == maximumViewportSize().width() || s.height() == maximumViewportSize().height()))
        emit scaleChanged(m_scale, true);
}

bool PreviewScrollArea::hasHorizontalScrollBar()
{
    return horizontalScrollBar()->maximum() - horizontalScrollBar()->minimum() > 0;
}

bool PreviewScrollArea::hasVerticalScrollBar()
{
    return verticalScrollBar()->maximum() - verticalScrollBar()->minimum() > 0;
}


void PreviewScrollArea::setMinScale(float s)
{
    m_minScale = s;
}

void PreviewScrollArea::setMaxScale(float s)
{
    m_maxScale = s;
}

void PreviewScrollArea::setScale(float sc)
{
    QRect r = contentsRect();
    r.moveTo(-widget()->pos());
    setScale(sc, r.center());
}

float PreviewScrollArea::scale() const
{
    return m_scale;
}

void PreviewScrollArea::setScaleStep(float s)
{
    m_scaleStep = s;
}

void PreviewScrollArea::zoomIn()
{
    QRect r = contentsRect();
    r.moveTo(-widget()->pos());
    zoomIn(r.center());
}

void PreviewScrollArea::zoomOut()
{
    QRect r = contentsRect();
    r.moveTo(-widget()->pos());
    zoomOut(r.center());
}

void PreviewScrollArea::zoomIn(const QPoint &p)
{
    if (m_scale < m_maxScale) {
        float scale = m_scale + m_scaleStep;
        if (scale > m_maxScale)
            scale = m_maxScale;
            setScale(scale, p);
    }
}

void PreviewScrollArea::zoomOut(const QPoint &p)
{
    if (m_scale > m_minScale) {
        float scale = m_scale - m_scaleStep;
        if (scale < m_minScale)
            scale = m_minScale;
            setScale(scale, p);
    }
}

void PreviewScrollArea::setScale(float sc, const QPoint &p)
{
    QPoint c = p / m_scale * sc - p;
    m_scale = sc;
    QSize s = widget()->baseSize() * sc;
    widget()->resize(s);
    emit scaleChanged(m_scale, !hasHorizontalScrollBar() && !hasVerticalScrollBar() && (s.width() == maximumViewportSize().width() || s.height() == maximumViewportSize().height()));
    if (hasHorizontalScrollBar())
        m_dx += c.x();
    else
        m_dx = c.x();
    if (hasVerticalScrollBar())
        m_dy += c.y();
    else
        m_dy = c.y();
}

void PreviewScrollArea::autoZoom()
{
    setScale(std::min(maximumViewportSize().width() / (float)widget()->baseSize().width(), maximumViewportSize().height() / (float)widget()->baseSize().height()));
}