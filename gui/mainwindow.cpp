#include "mainwindow.h"
#include "inputwidget.h"
#include "resultwidget.h"
#include "version.h"

#include "util.h"

#include <QDesktopWidget>
#include <QApplication>
#include <QMessageBox>
#include <QCloseEvent>
#include <QPushButton>
#include <QTabWidget>
#include <QFileInfo>
#include <QTabBar>
#include <QStyle>

MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags)
: QMainWindow(parent, flags), m_widget(new QTabWidget(this))
{
    setCentralWidget(m_widget);
    m_widget->setTabsClosable(true);
    m_widget->setDocumentMode(true);
    connect(m_widget, SIGNAL(tabCloseRequested(int)), this, SLOT(tabClosed(int)));
    InputWidget *in = new InputWidget(this);
    m_widget->addTab(in, "I&nput");
    connect(in, SIGNAL(run(const GuiParams&)), this, SLOT(run(const GuiParams&)));

    QPushButton *aboutButton = new QPushButton(this);
    aboutButton->setFlat(true);
    aboutButton->setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed));
    aboutButton->setToolTip("About simpleSTORM");
#ifdef Q_WS_X11
    aboutButton->setIcon(QIcon::fromTheme("help-about"));
#else
    aboutButton->setIcon(QIcon(":/help-about.svgz"));
#endif
    m_widget->setCornerWidget(aboutButton, Qt::TopRightCorner);
    connect(aboutButton, SIGNAL(clicked()), this, SLOT(aboutClicked()));

    QTabBar *tabBar = m_widget->findChild<QTabBar *>();
    tabBar->adjustSize();
    tabBar->setTabButton(0, (QTabBar::ButtonPosition)tabBar->style()->styleHint(QStyle::SH_TabBar_CloseButtonPosition), 0);
    connect(QCoreApplication::instance(), SIGNAL(aboutToQuit()), this, SLOT(onApplicationClose()));
    adjustSize();
    QSize s = size();
    s.setWidth(s.height());
    setGeometry(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, s, QApplication::desktop()->availableGeometry()));

}

MainWindow::~MainWindow()
{
}

void MainWindow::aboutClicked()
{
    QMessageBox::about(this, "About simpleSTORM",
                       QString("simpleSTORM %1\n"
                       "simpleSTORM GUI %2\n"
                       "GUI Version %3, STORM Version %4").arg(QString::fromStdString(wienerStormAuthors())).arg(STORMGUI_AUTHORS).arg(STORMGUI_VERSION_STRING).arg(QString::fromStdString(wienerStormVersion())));
}

void MainWindow::run(const GuiParams &params)
{
    ResultWidget *result = new ResultWidget(this);
    m_widget->setCurrentIndex(m_widget->addTab(result, QFileInfo(QString::fromStdString(params.getInFile())).fileName()));
    connect(result, SIGNAL(userAttentionRequired(QWidget*)), m_widget, SLOT(setCurrentWidget(QWidget*)));
    result->start(params);
}

bool MainWindow::tabClosed(int index)
{
    ResultWidget *result = static_cast<ResultWidget*>(m_widget->widget(index));
    if (result->canBeClosed()) {
        m_widget->removeTab(index);
        if (result->canBeDeleted())
            delete result;
        else {
            connect(result, SIGNAL(finished(ResultWidget*)), this, SLOT(deletedResultFinished(ResultWidget*)));
            m_removedJobs.insert(result);
        }
        return true;
    }
    return false;
}

void MainWindow::deletedResultFinished(ResultWidget *result)
{
    m_removedJobs.remove(result);
    delete result;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    for (int i = 1; i < m_widget->count(); ++i) {
        if (!tabClosed(i)){
            event->ignore();
            return;
        }
    }
    event->accept();
}

void MainWindow::onApplicationClose()
{
    for (ResultWidget *job : m_removedJobs) {
        job->waitFor();
        delete job;
    }
}
