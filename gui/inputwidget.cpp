#include "inputwidget.h"
#include "ui_inputwidget.h"
#include "ui_advancedsettingsgroupbox.h"
#include "ui_backgroundlevelgroupbox.h"
#include "globals.h"

#include <QGroupBox>
#include <QFileDialog>

InputWidget::InputWidget(QWidget *parent)
: QWidget(parent), m_ui(new Ui::InputWidget()), m_uiAdvancedSettings(new Ui::AdvancedSettingsGroupBox()), m_uiBackgroundLevel(new Ui::BackgroundLevelGroupBox())
{
    m_ui->setupUi(this);
    QGroupBox *tmp = new QGroupBox(this);
    m_uiBackgroundLevel->setupUi(tmp);
    m_ui->stck_advancedSettings->addWidget(tmp);
    tmp = new QGroupBox(this);
    m_uiAdvancedSettings->setupUi(tmp);
    m_ui->stck_advancedSettings->addWidget(tmp);
#ifdef Q_WS_X11
    m_ui->btn_run->setIcon(QIcon::fromTheme("system-run"));
#else
    m_ui->btn_run->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogOkButton));
#endif
    connect(m_ui->btn_inputFile, SIGNAL(clicked()), this, SLOT(inputFileButtonClicked()));
    connect(m_ui->btn_settingsFile, SIGNAL(clicked()), this, SLOT(settingsFileButtonClicked()));
    connect(m_ui->lne_inputFile, SIGNAL(textEdited(const QString&)), this, SLOT(inputFileEdited(const QString&)));
    connect(m_ui->lne_settingsFile, SIGNAL(textEdited(const QString&)), this, SLOT(settingsFileEdited(const QString&)));
    connect(m_ui->btn_run, SIGNAL(clicked()), this, SLOT(runClicked()));
    connect(m_ui->chk_advancedSettings, SIGNAL(toggled(bool)), this, SLOT(advancedSettingsToggled(bool)));
    advancedSettingsToggled(m_ui->chk_advancedSettings->isChecked());
    enableInput(false);
}

InputWidget::~InputWidget()
{
    delete m_ui;
    delete m_uiAdvancedSettings;
    delete m_uiBackgroundLevel;
}

void InputWidget::inputFileButtonClicked()
{
    QString file = QFileDialog::getOpenFileName(this, "Select input file", fileOpenDialogDirectory, "Accepted Files (*.sif *.h5 *.hdf *.hdf5 *.tif *.tiff)");
    if (!file.isEmpty()) {
        fileOpenDialogDirectory = QFileInfo(file).absoluteDir().absolutePath();
        inputFileEdited(file);
    }
}

void InputWidget::settingsFileButtonClicked()
{
    QString file = QFileDialog::getOpenFileName(this, "Select settings file", fileOpenDialogDirectory, "INI-style text files (*.txt *.ini)");
    if (!file.isEmpty()) {
        fileOpenDialogDirectory = QFileInfo(file).absoluteDir().absolutePath();
        m_ui->lne_inputFile->setText(file);
        m_params.setSettingsFile(file.toStdString());
        m_ui->lne_settingsFile->setText(file);
        setFieldsFromDefaults();
    }
}

void InputWidget::inputFileEdited(const QString &file)
{
    if (QFileInfo(file).exists()) {
        m_ui->lne_inputFile->setText(file);
        m_params.setInFile(file.toStdString(), true);
        m_ui->lne_settingsFile->setText(QString::fromStdString(m_params.getSettingsFile()));
        setFieldsFromDefaults();
        enableInput(true);
        m_ui->btn_run->setEnabled(true);
    } else {
        enableInput(false);
        m_ui->btn_run->setEnabled(false);
    }
}

void InputWidget::settingsFileEdited(const QString &file) {
    m_params.setSettingsFile(file.toStdString());
    setFieldsFromDefaults();
}

void InputWidget::advancedSettingsToggled(bool toggled)
{
    m_ui->stck_advancedSettings->setCurrentIndex(toggled ? 1 : 0);
    m_params.setAdvancedSettingsEnabled(toggled);
}

void InputWidget::runClicked()
{
    m_params.setFactor(m_ui->spn_factor->value());
    m_params.setPixelSize(m_ui->spn_pixelSize->value());
    m_params.setReconstructionResolution(m_ui->spn_reconstructionRes->value());
    m_params.setSkellamFrames(m_ui->spn_skellamFrames->value());
    m_params.setAlpha(m_ui->spn_alpha->value());

    if (m_ui->chk_advancedSettings->isChecked()) {
        m_params.setRoilen(m_uiAdvancedSettings->spn_roilen->value());
        m_params.setXYChunkSize(m_uiAdvancedSettings->spn_xyChunkSize->value());
        m_params.setTChunkSize(m_uiAdvancedSettings->spn_tChunkSize->value());
        m_params.setChunksInMemory(m_uiAdvancedSettings->spn_chunksInMemory->value());
        m_params.setDoAsymmetryCheck(m_uiAdvancedSettings->cB_doAsymmetryCheck->isChecked());
    } else {
        int resx = m_params.shape()[0], resy = m_params.shape()[1], resz = m_params.shape()[2];
        int currValueXYSlider = m_uiBackgroundLevel->sldr_xyBackgroundLevel->value();
        int maxXYSize = std::min(std::min(resx/2-1,resy/2-1), m_params.getMaxXyChunksize()), minXYSize = m_params.getMinXyChunksize();
        int maxValueXYSlider = m_uiBackgroundLevel->sldr_xyBackgroundLevel->maximum();
        m_params.setXYChunkSize((int)((maxValueXYSlider-currValueXYSlider)/(float)maxValueXYSlider*(maxXYSize - minXYSize) + minXYSize));
        int currValueTSlider = m_uiBackgroundLevel->sldr_tBackgroundLevel->value();
        int maxTSize = std::min((int)(resz/m_params.getChunksInMemory()), m_params.getMaxTChunkSize()), minTSize = m_params.getMinTChunkSize();
        int maxValueTSlider = m_uiBackgroundLevel->sldr_tBackgroundLevel->maximum();
        m_params.setTChunkSize((int)((maxValueTSlider-currValueTSlider)/(float)maxValueTSlider*(maxTSize - minTSize) + minTSize));
        int currValueSpotSlider = m_uiBackgroundLevel->sldr_spotDensity->value();
        int maxAsymmVal = m_params.getMaxAsymmetryThreshold(), minAsymmVal = m_params.getMinAsymmetryThreshold();
        int maxValueSpotSlider = m_uiBackgroundLevel->sldr_tBackgroundLevel->maximum();
        m_params.setDoAsymmetryCheck(true);
        m_params.setAsymmetryThreshold(currValueSpotSlider/(float)maxValueSpotSlider*(maxAsymmVal-minAsymmVal)+minAsymmVal);
        std::cout<<"xyChunksize: "<<(int)((maxValueXYSlider-currValueXYSlider)/(float)maxValueXYSlider*(maxXYSize - minXYSize) + minXYSize)<<std::endl
        <<"tChunksize: "<< (int)((maxValueTSlider-currValueTSlider)/(float)maxValueTSlider*(maxTSize - minTSize) + minTSize)<<std::endl
        <<"Asymmetry Threshold: "<<currValueSpotSlider/(float)maxValueSpotSlider*(maxAsymmVal-minAsymmVal)+minAsymmVal<<std::endl;
        m_params.setChunksInMemory(m_uiAdvancedSettings->spn_chunksInMemory->value());
        m_params.setRoilen(m_uiAdvancedSettings->spn_roilen->value());
    }
    if (m_ui->spn_gain->value() > 0) {
        m_params.setSlope(m_ui->spn_gain->value());
        m_params.setUseSavedSlope(true);
    } else
        m_params.setUseSavedSlope(false);
    if (m_ui->spn_offset->value() > 0) {
        m_params.setIntercept(m_ui->spn_offset->value());
        m_params.setUseSavedIntercept(true);
    } else
        m_params.setUseSavedIntercept(false);
    if (m_ui->spn_sigma->value() > 0) {
        m_params.setSigma(m_ui->spn_sigma->value());
        m_params.setUseSavedSigma(true);
    } else
        m_params.setUseSavedSigma(false);
    m_params.doSanityChecks();
    emit run(m_params);
}

void InputWidget::setFieldsFromDefaults()
{
    m_ui->spn_factor->setValue(m_params.getFactor());
    m_ui->spn_pixelSize->setValue(m_params.getPixelSize());
    m_ui->spn_reconstructionRes->setValue(m_params.getReconstructionResolution());
    m_ui->spn_skellamFrames->setValue(m_params.getSkellamFrames());
    m_ui->spn_alpha->setValue(m_params.getAlpha());

    m_ui->chk_advancedSettings->setChecked(m_params.getAdvancedSettingsEnabled());

    m_uiAdvancedSettings->spn_roilen->setValue(m_params.getRoilen());
    m_uiAdvancedSettings->spn_xyChunkSize->setValue(m_params.getXYChunkSize());
    m_uiAdvancedSettings->spn_tChunkSize->setValue(m_params.getTChunkSize());
    m_uiAdvancedSettings->spn_chunksInMemory->setValue(m_params.getChunksInMemory());

    if (m_params.getSlopeSaved())
        m_ui->spn_gain->setValue(m_params.getSlope());
    if (m_params.getInterceptSaved())
        m_ui->spn_offset->setValue(m_params.getIntercept());
    if (m_params.getSigmaSaved())
        m_ui->spn_sigma->setValue(m_params.getSigma());
}

void InputWidget::enableInput(bool enable)
{
    m_ui->btn_run->setEnabled(enable);
    m_ui->lne_settingsFile->setEnabled(enable);
    m_ui->btn_settingsFile->setEnabled(enable);
    m_ui->chk_advancedSettings->setEnabled(enable);
    m_ui->grp_general->setEnabled(enable);
    m_ui->stck_advancedSettings->setEnabled(enable);
    m_ui->grp_data->setEnabled(enable);
}
