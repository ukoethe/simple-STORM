#ifndef OPTIONWIDGET_H
#define OPTIONWIDGET_H

#include <QDialog>
#include "inputwidget.h"

#ifdef WIN32
	typedef std::basic_string<TCHAR> tstring;
#endif
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
#ifdef WIN32
	tstring defaultDirectory();
#else
	std::string defaultDirectory();
#endif
private Q_SLOTS:
	void SaveButtonClicked();
	void CloseWindow();
	void loadSavedDefaults();
	void AdjustReconstructionResolution();
};

#endif // OPTIONWIDGET_H
