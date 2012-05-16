#ifndef MAINWINDOWSTRUCTUREANALYSIS_H
#define MAINWINDOWSTRUCTUREANALYSIS_H

#include <QMainWindow>

namespace Ui {
    class MainWindowStructureAnalysis;
}

class MainWindowStructureAnalysis : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindowStructureAnalysis(QWidget *parent = 0);
    ~MainWindowStructureAnalysis();

private:
    Ui::MainWindowStructureAnalysis *ui;
};

#endif // MAINWINDOWSTRUCTUREANALYSIS_H
