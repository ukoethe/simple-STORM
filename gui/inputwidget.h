#ifndef INPUTWIDGET_H
#define INPUTWIDGET_H

#include "guiparams.h"

#include <QWidget>

class OptionWidget;

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
	void enableInput(bool);
	GuiParams m_params;
	OptionWidget *opt;
	Ui::InputWidget *m_ui;
	
Q_SIGNALS:
    void run(const GuiParams&);

private:
    void setFieldsFromDefaults();
    
	
    
    Ui::AdvancedSettingsGroupBox *m_uiAdvancedSettings;
    Ui::BackgroundLevelGroupBox *m_uiBackgroundLevel;

private Q_SLOTS:
    void inputFileButtonClicked();
    void inputFileEdited(const QString&);
    void settingsFileButtonClicked();
    void settingsFileEdited(const QString&);
    void advancedSettingsToggled(bool);
    void doAsymmetryCheckToggled(bool);
	void optionsClicked();
    void runClicked();
    void factorEdited();
    void pixelSizeEdited();
    void reconstructionResolutionEdited();
};

#endif
