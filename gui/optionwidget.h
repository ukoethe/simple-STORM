#ifndef OPTIONWIDGET_H
#define OPTIONWIDGET_H

#include <QDialog>
#include "inputwidget.h"

typedef std::basic_string<TCHAR> tstring;

namespace rude{
	class Config;
}

namespace Ui {
class OptionWidget;
}

class OptionWidget : public QDialog
{
    Q_OBJECT
    
public:
    explicit OptionWidget(InputWidget *parent = 0);
    ~OptionWidget();
	InputWidget *parenthandle;

private:
	
    Ui::OptionWidget *ui;
	rude::Config *m_config;
	tstring defaultDirectory();
private Q_SLOTS:
	void SaveButtonClicked();
	void CloseWindow();
	void loadSavedDefaults();
	void AdjustReconstructionResolution();
};

#endif // OPTIONWIDGET_H