#include "globals.h"
#include "mainwindow.h"

#include "wienerStorm.hxx"

#include <clocale>

#include <QApplication>
#include <QTextCodec>
#include <QTranslator>
#include <QLocale>
#include <QLibraryInfo>
#include <QDir>

QString fileOpenDialogDirectory;
QString fileSaveDialogDirectory;

int main(int argc, char* argv[])
{
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
    QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8"));

    fileOpenDialogDirectory = QDir::homePath();
    fileSaveDialogDirectory = QDir::homePath();

    QApplication app(argc, argv);
    initR(argc, argv);
    std::setlocale(LC_NUMERIC,"C");

    MainWindow mainWindow(0, Qt::Window);
    mainWindow.setWindowTitle("simpleSTORM");
    mainWindow.setWindowIcon(QIcon(":/StormIcon.png"));
    //mainWindow.setWindowState(Qt::WindowMaximized);
    mainWindow.show();
    int ret = app.exec();
    endR();
    return ret;
}
