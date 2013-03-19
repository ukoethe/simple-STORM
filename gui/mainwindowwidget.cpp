#include "mainwindowwidget.h"
#include "ui_mainwindowwidget.h"
#include "inputwidget.h"
#include "previewwidget.h"

MainWindowWidget::MainWindowWidget(QWidget *parent)
: QWidget(parent), m_ui(new Ui::MainWindowWidget()), m_input(new InputWidget(this))
{
    m_ui->setupUi(this);
#ifdef Q_WS_X11
    m_ui->btn_run->setIcon(QIcon::fromTheme("system-run"));
    m_ui->btn_abort->setIcon(QIcon::fromTheme("dialog-cancel"));
    m_ui->btn_about->setIcon(QIcon::FromTheme("dialog-about"));
#else
    m_ui->btn_run->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogOkButton));
    m_ui->btn_abort->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogCancelButton));
    m_ui->btn_about->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogHelpButton));
#endif
    m_ui->btn_abort->setVisible(false);
    m_ui->prgrs_bar->setVisible(false);

    connect(m_input, SIGNAL(valid(bool)), this, SLOT(validInput(bool)));
    connect(m_ui->btn_run, SIGNAL(clicked()), this, SLOT(runClicked()));
}

MainWindowWidget::~MainWindowWidget()
{
    delete m_ui;
}

void MainWindowWidget::validInput(bool valid)
{
    m_ui->btn_run->setEnabled(valid);
}

void MainWindowWidget::runClicked()
{

}


