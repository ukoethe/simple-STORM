#include "optionwidget.h"
#include "ui_optionwidget.h"
#include "inputwidget.h"
#include <iostream>
#include <rude/config.h>

#ifdef WIN32
	#include <shlobj.h>
	#include <shobjidl.h>
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif


OptionWidget::OptionWidget(InputWidget *parent) :
    QDialog(parent),
    ui(new Ui::OptionWidget), m_config(new rude::Config())
{
    ui->setupUi(this);
	parenthandle = parent;
	connect(ui->btn_saveSettings, SIGNAL(clicked()), this, SLOT(SaveButtonClicked()));
	connect(ui->btn_cancel,SIGNAL(clicked()), this, SLOT(CloseWindow()));
	connect(ui->spn_default_factor, SIGNAL(editingFinished()), this, SLOT(AdjustReconstructionResolution()));
	connect(ui->spn_default_input_resolution, SIGNAL(editingFinished()), this, SLOT(AdjustReconstructionResolution()));
	std::string defaultdir = defaultDirectory();
	std::string fname = (defaultdir+"/SimpleSTORMsettings.txt");
	m_config->setConfigFile((fname).c_str());
	m_config->load();
	loadSavedDefaults();
	//int a = parent->m_ui->spn_factor()->value();
	//ui->txtlbl_current_factor->setText(parent->m_ui->spn_factor()->value());
	//ui->txtlbl_current_input_resolution->setText(QString::number(parent->m_ui->));
}

OptionWidget::~OptionWidget()
{
    delete ui;
	delete m_config;
}

void OptionWidget::SaveButtonClicked()
{
	std::string dir;
	dir = defaultDirectory();
	m_config->setIntValue("factor", ui->spn_default_factor->value());
	m_config->setDoubleValue("pixelsize", ui->spn_default_input_resolution->value());
	m_config->setDoubleValue("reconstructionresolution", ui->lbl_default_reconstruction_resolution->text().toDouble());
	m_config->setIntValue("skellamFrames", ui->spn_default_skellam->value());
	m_config->setDoubleValue("alpha", ui->spn_default_sensitivity->value()/100);
	m_config->save();
}

void OptionWidget::CloseWindow()
{
	this->close();
}

void OptionWidget::loadSavedDefaults()
{
	try
	{
		ui->spn_default_factor->setValue(m_config->getIntValue("factor"));
		ui->spn_default_input_resolution->setValue(m_config->getDoubleValue("pixelsize"));
		ui->lbl_default_reconstruction_resolution->setText(QString::number(m_config->getDoubleValue("reconstructionresolution")));
		ui->spn_default_skellam->setValue(m_config->getIntValue("skellamFrames"));
		ui->spn_default_sensitivity->setValue(m_config->getDoubleValue("alpha")*100);
	}
	catch (int e)
	{std::cout<<"Not standard values found!"<<e<<std::endl;}
}

void OptionWidget::AdjustReconstructionResolution()
{
	try
	{
		ui->lbl_default_reconstruction_resolution->setText(QString::number(ui->spn_default_input_resolution->value()/ui->spn_default_factor->value()));
	}
	catch (int e)
	{
		ui->lbl_default_reconstruction_resolution->setText("0");
		std::cout<<e<<std::endl;
	}
}

tstring OptionWidget::defaultDirectory()
{
	#ifdef WIN32
	{
		char szPath[MAX_PATH];
		SHGetFolderPath(NULL,CSIDL_APPDATA,0,NULL,szPath);
		tstring a(szPath);
		return a;
	}
	#endif
	char c[FILENAME_MAX];
	GetCurrentDir(c, sizeof(c));
	return c;
}

std::string TCharToString(LPCTSTR t)
{
	// Handy for converting TCHAR to std::string (char)
	// If the conversion fails, an empty string is returned.
	std::string str;
	#ifdef UNICODE
		// calculate the size of the char string required
		// Note: If wcstombs encounters a wide character it cannot convert
		//       to a multibyte character, it returns –1.
		int len = 1 + wcstombs(0, t, 0);
		if (0 == len) return str;

		char* c = new char[len];
		if (NULL == c) throw std::bad_alloc();
		c[0] = '\0';

		wcstombs(c, t, len);
		str = c;
		delete []c;
		#else
		str = t;
	#endif
	return str;
}