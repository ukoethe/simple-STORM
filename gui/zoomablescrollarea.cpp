#include "zoomablescrollarea.h"

#include <QWheelEvent>

ZoomableScrollArea::ZoomableScrollArea(QWidget *parent)
: QScrollArea(parent)
{}

ZoomableScrollArea::~ZoomableScrollArea()
{}

void ZoomableScrollArea::wheelEvent(QWheelEvent *event)
{
    if (event->modifiers().testFlag(Qt::ControlModifier)) {
        event->accept();
        if (event->delta() > 0)
            emit zoomIn();
        else
            emit zoomOut();
    } else
        QScrollArea::wheelEvent(event);
}
