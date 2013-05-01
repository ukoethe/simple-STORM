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
    QStringList ftypes;
    QString validator("^.*(");
    for (const std::string &ftype : m_params.acceptedFileTypes()) {
        ftypes.append("\\" + QString::fromStdString(ftype));
    }
    validator.append(ftypes.join("|")).append(")$");
    m_ui->lne_inputFile->setValidator(new QRegExpValidator(QRegExp(validator, Qt::CaseInsensitive), this));
    QGroupBox *tmp = new QGroupBox(this);
    m_uiBackgroundLevel->setupUi(tmp);
    m_ui->stck_advancedSettings->addWidget(tmp);
    tmp = new QGroupBox(this);
    m_uiAdvancedSettings->setupUi(tmp);
    m_ui->stck_advancedSettings->addWidget(tmp);
#ifdef Q_WS_X11
    m_ui->btn_run->setIcon(QIcon::fromTheme("system-run"));
#else
    m_ui->btn_run->setIcon(QIcon(":/system-run.svgz"));
#endif
    connect(m_ui->btn_inputFile, SIGNAL(clicked()), this, SLOT(inputFileButtonClicked()));
    connect(m_ui->btn_settingsFile, SIGNAL(clicked()), this, SLOT(settingsFileButtonClicked()));
    connect(m_ui->lne_inputFile, SIGNAL(textEdited(const QString&)), this, SLOT(inputFileEdited(const QString&)));
    connect(m_ui->lne_settingsFile, SIGNAL(textEdited(const QString&)), this, SLOT(settingsFileEdited(const QString&)));
    connect(m_ui->btn_run, SIGNAL(clicked()), this, SLOT(runClicked()));
    connect(m_ui->chk_advancedSettings, SIGNAL(toggled(bool)), this, SLOT(advancedSettingsToggled(bool)));
    advancedSettingsToggled(m_ui->chk_advancedSettings->isChecked());
    connect(m_ui->spn_factor, SIGNAL(editingFinished()), this, SLOT(factorEdited()));
    connect(m_ui->spn_pixelSize, SIGNAL(editingFinished()), this, SLOT(pixelSizeEdited()));
    connect(m_ui->spn_reconstructionRes, SIGNAL(editingFinished()), this, SLOT(reconstructionResolutionEdited()));
    connect(m_uiAdvancedSettings->chk_doAsymmetryCheck, SIGNAL(toggled(bool)), this, SLOT(doAsymmetryCheckToggled(bool)));
    doAsymmetryCheckToggled(m_uiAdvancedSettings->chk_doAsymmetryCheck->isChecked());
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
    QStringList accepted;
    for (const std::string &ftype : m_params.acceptedFileTypes()) {
        accepted.append("*" + QString::fromStdString(ftype));
    }
    QString file = QFileDialog::getOpenFileName(this, "Select input file", fileOpenDialogDirectory, "Accepted Files (" + accepted.join(" ") + ")");
    if (!file.isEmpty()) {
        fileOpenDialogDirectory = QFileInfo(file).absoluteDir().absolutePath();
        m_ui->lne_inputFile->setText(file);
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
    QFileInfo finfo(file);
    if (m_ui->lne_inputFile->hasAcceptableInput() && finfo.isFile() && finfo.exists()) {
        fileOpenDialogDirectory = finfo.absoluteDir().absolutePath();
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

void InputWidget::doAsymmetryCheckToggled(bool toggled)
{
    m_uiAdvancedSettings->spn_AsymmetryThreshold->setEnabled(toggled);
    m_uiAdvancedSettings->label_AsymThreshold->setEnabled(toggled);
}

void InputWidget::runClicked()
{
    m_params.setFactor(m_ui->spn_factor->value());
    m_params.setPixelSize(m_ui->spn_pixelSize->value());
    m_params.setReconstructionResolution(m_ui->spn_reconstructionRes->value());
    m_params.setAlpha((m_ui->spn_alpha->value())/100.);

    if (m_ui->chk_advancedSettings->isChecked()) {
        m_params.setRoilen(m_uiAdvancedSettings->spn_roilen->value());
        m_params.setXYChunkSize(m_uiAdvancedSettings->spn_xyChunkSize->value());
        m_params.setTChunkSize(m_uiAdvancedSettings->spn_tChunkSize->value());
        m_params.setChunksInMemory(m_uiAdvancedSettings->spn_chunksInMemory->value());
        m_params.setDoAsymmetryCheck(m_uiAdvancedSettings->chk_doAsymmetryCheck->isChecked());
        m_params.setPrefactorSigma(m_uiAdvancedSettings->spn_factorSigma->value());
        if(m_uiAdvancedSettings->chk_doAsymmetryCheck->isChecked()) {
            m_params.setAsymmetryThreshold(m_uiAdvancedSettings->spn_AsymmetryThreshold->value());
            m_params.setDoAsymmetryCheck(true);
        }
    }
    else {
        int resx = m_params.shape()[0], resy = m_params.shape()[1], resz = m_params.shape()[2];
        int currValueXYSlider = m_uiBackgroundLevel->sldr_xyBackgroundLevel->value();
        int maxXYSize = std::min(std::min(resx/2-1,resy/2-1), m_params.getMaxXyChunksize()), minXYSize = m_params.getMinXyChunksize();
        int maxValueXYSlider = m_uiBackgroundLevel->sldr_xyBackgroundLevel->maximum();
        m_params.setXYChunkSize((int)((maxValueXYSlider-currValueXYSlider)/(float)maxValueXYSlider*(maxXYSize - minXYSize) + minXYSize));
        int currValueTSlider = m_uiBackgroundLevel->sldr_tBackgroundLevel->value();
        int maxTSize = std::min((int)(m_ui->spn_skellamFrames->value()/m_params.getChunksInMemory()), m_params.getMaxTChunkSize()), minTSize = m_params.getMinTChunkSize();
        int maxValueTSlider = m_uiBackgroundLevel->sldr_tBackgroundLevel->maximum();
        m_params.setTChunkSize((int)((maxValueTSlider-currValueTSlider)/(float)maxValueTSlider*(maxTSize - minTSize) + minTSize));
        int currValueSpotSlider = m_uiBackgroundLevel->sldr_spotDensity->value();
        int maxAsymmVal = m_params.getMaxAsymmetryThreshold(), minAsymmVal = m_params.getMinAsymmetryThreshold();
        int maxValueSpotSlider = m_uiBackgroundLevel->sldr_tBackgroundLevel->maximum();
        m_params.setDoAsymmetryCheck(true);
        m_params.setAsymmetryThreshold(currValueSpotSlider/(float)maxValueSpotSlider*(maxAsymmVal-minAsymmVal)+minAsymmVal);
        float maxPSFFactorVal = 1, minPSFFactorVal = 0;
        m_params.setPrefactorSigma((maxValueSpotSlider-currValueSpotSlider)/(float)maxValueSpotSlider*(maxPSFFactorVal-minPSFFactorVal)+minPSFFactorVal);
        m_params.setChunksInMemory(m_uiAdvancedSettings->spn_chunksInMemory->value());
        m_params.setRoilen(m_uiAdvancedSettings->spn_roilen->value());
    }
    m_params.setSlope(m_ui->spn_gain->value());
    m_params.setUseSavedSlope(m_ui->spn_gain->value() > 0);
    m_params.setIntercept(m_ui->spn_offset->value());
    m_params.setUseSavedIntercept(m_ui->spn_offset->value() > 0);
    m_params.setSigma(m_ui->spn_sigma->value());
    m_params.setUseSavedSigma(m_ui->spn_sigma->value() > 0);
    std::cout<<"getSlopeSaved:"<<m_params.getSlopeSaved()<<std::endl;
    if (!(m_params.getSlopeSaved() and m_params.getInterceptSaved())) {
        m_params.setSkellamFrames(m_ui->spn_skellamFrames->value());
    }
    else
        m_params.setIgnoreSkellamFramesSaved(true);
    std::cout<<"IgnoreSkellamFramesSaved: "<<m_params.getIgnoreSkellamFramesSaved()<<std::endl;
    m_params.doSanityChecks();

    emit run(m_params);
}

void InputWidget::setFieldsFromDefaults()
{
    m_ui->spn_factor->setValue(m_params.getFactor());
    m_ui->spn_pixelSize->setValue(m_params.getPixelSize());
    m_ui->spn_reconstructionRes->setValue(m_params.getReconstructionResolution());
    m_ui->spn_skellamFrames->setValue(m_params.getSkellamFrames());
    m_ui->spn_alpha->setValue(m_params.getAlpha()*100);

    m_ui->chk_advancedSettings->setChecked(m_params.getAdvancedSettingsEnabled());

    m_uiAdvancedSettings->spn_roilen->setValue(m_params.getRoilen());
    m_uiAdvancedSettings->spn_xyChunkSize->setValue(m_params.getXYChunkSize());
    m_uiAdvancedSettings->spn_tChunkSize->setValue(m_params.getTChunkSize());
    m_uiAdvancedSettings->spn_chunksInMemory->setValue(m_params.getChunksInMemory());


    if (m_params.getSlopeSaved())
        m_ui->spn_gain->setValue(m_params.getSlope());
    else
        m_ui->spn_gain->setValue(0);
    if (m_params.getInterceptSaved())
        m_ui->spn_offset->setValue(m_params.getIntercept());
    else
        m_ui->spn_offset->setValue(0);
    if (m_params.getSigmaSaved())
        m_ui->spn_sigma->setValue(m_params.getSigma());
    else
        m_ui->spn_offset->setValue(0);
    if (m_params.getDoAsymmetryCheck())
        m_uiAdvancedSettings->spn_AsymmetryThreshold->setValue(m_params.getAsymmetryThreshold());
    m_uiAdvancedSettings->spn_factorSigma->setValue(m_params.getPrefactorSigma());
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
    if (m_ui->chk_advancedSettings->isChecked() && m_params.getDoAsymmetryCheckSaved())
        m_uiAdvancedSettings->chk_doAsymmetryCheck->setCheckState(Qt::Checked);
}

void InputWidget::factorEdited()
{
    if (m_ui->spn_reconstructionRes->value() * m_ui->spn_factor->value() < m_ui->spn_pixelSize->value())
        m_ui->spn_reconstructionRes->setValue(m_ui->spn_pixelSize->value() / m_ui->spn_factor->value());
}

void InputWidget::pixelSizeEdited()
{
    if (m_ui->spn_reconstructionRes->value() * m_ui->spn_factor->value() < m_ui->spn_pixelSize->value())
        m_ui->spn_reconstructionRes->setValue(m_ui->spn_pixelSize->value() / m_ui->spn_factor->value());
}

void InputWidget::reconstructionResolutionEdited()
{
    if (m_ui->spn_pixelSize->value() / m_ui->spn_reconstructionRes->value() > m_ui->spn_factor->value())
        m_ui->spn_factor->setValue(std::ceil(m_ui->spn_pixelSize->value() / m_ui->spn_reconstructionRes->value()));
}
