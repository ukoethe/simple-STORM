#include "mainwindowstructureanalysis.h"
#include "ui_mainwindowstructureanalysis.h"

MainWindowStructureAnalysis::MainWindowStructureAnalysis(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindowStructureAnalysis)
{
    ui->setupUi(this);
}

MainWindowStructureAnalysis::~MainWindowStructureAnalysis()
{
    delete ui;
}
