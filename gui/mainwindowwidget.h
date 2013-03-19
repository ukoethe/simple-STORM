#ifndef MAINWINDOWWIDGET_H
#define MAINWINDOWWIDGET_H

#include <QWidget>

namespace Ui
{
class MainWindowWidget;
}
class InputWidget;
class PreviewWidget;

class MainWindowWidget : public QWidget
{
Q_OBJECT
public:
    MainWindowWidget(QWidget *parent = 0);
    ~MainWindowWidget();

private:
    Ui::MainWindowWidget *m_ui;
    InputWidget *m_input;

private Q_SLOTS:
    void runClicked();
    void abortClicked();
    void aboutClicked();
};

#endif
