#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "guiparams.h"

#include <QMainWindow>
#include <QSet>

class QTabWidget;
class QCloseEvent;

class ResultWidget;

class MainWindow : public QMainWindow
{
Q_OBJECT
public:
    MainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
    ~MainWindow();

protected:
    virtual void closeEvent(QCloseEvent *event);

private:
    QTabWidget *m_widget;
    QSet<ResultWidget*> m_removedJobs;

private Q_SLOTS:
    void aboutClicked();
    void run(const GuiParams&);
    bool tabClosed(int);
    void deletedResultFinished(ResultWidget*);
    void onApplicationClose();
};

#endif
