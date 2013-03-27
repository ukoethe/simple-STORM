#ifndef INPUTWIDGET_H
#define INPUTWIDGET_H

#include "guiparams.h"

#include <QWidget>

namespace Ui
{
class InputWidget;
class AdvancedSettingsGroupBox;
class BackgroundLevelGroupBox;
}
class QGroupBox;

class InputWidget : public QWidget
{
Q_OBJECT
public:
    InputWidget(QWidget *parent = 0);
    ~InputWidget();

Q_SIGNALS:
    void run(const GuiParams&);

private:
    void setFieldsFromDefaults();
    void enableInput(bool);

    Ui::InputWidget *m_ui;
    Ui::AdvancedSettingsGroupBox *m_uiAdvancedSettings;
    Ui::BackgroundLevelGroupBox *m_uiBackgroundLevel;
    GuiParams m_params;

private Q_SLOTS:
    void inputFileButtonClicked();
    void inputFileEdited(const QString&);
    void settingsFileButtonClicked();
    void settingsFileEdited(const QString&);
    void advancedSettingsToggled(bool);
    void runClicked();
    void factorEdited(int);
    void pixelSizeEdited(double);
    void reconstructionResolutionEdited(int);
};

#endif
