#include "mainwindow.h"

#include "wienerStorm.hxx"

#include <clocale>

#include <QApplication>
#include <QTextCodec>
#include <QTranslator>
#include <QLocale>
#include <QLibraryInfo>

int main(int argc, char* argv[])
{
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
    QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8"));
    QApplication app(argc, argv);
    initR(app.applicationDirPath().toStdString(), argc, argv);
    std::setlocale(LC_NUMERIC,"C");

    MainWindow mainWindow(0, Qt::Window);
    mainWindow.setWindowTitle("simpleSTORM");
    //mainWindow.setWindowState(Qt::WindowMaximized);
    mainWindow.show();
    int ret = app.exec();
    endR();
    return ret;
}
