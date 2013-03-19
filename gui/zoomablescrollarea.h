#ifndef ZOOMABLESCROLLAREA_H
#define ZOOMABLESCROLLAREA_H

#include <QScrollArea>

class ZoomableScrollArea : public QScrollArea
{
Q_OBJECT
public:
    ZoomableScrollArea(QWidget *parent = 0);
    ~ZoomableScrollArea();

Q_SIGNALS:
    void zoomIn();
    void zoomOut();

protected:
    virtual void wheelEvent(QWheelEvent *);
};

#endif
