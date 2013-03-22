#include "resultwidget.h"
#include "ui_resultwidget.h"
#include "analysisworker.h"
#include "globals.h"

#include "wienerStorm.hxx"

#include <iostream>

#include <QMessageBox>
#include <QScrollBar>
#include <QFileDialog>
#include <QFileInfo>
#include <QDir>

ResultWidget::ResultWidget(QWidget *parent)
: QWidget(parent), m_ui(new Ui::ResultWidget()), m_functor(0), m_worker(0), m_coordinatesSaved(false), m_reconstructionSaved(false), m_autoZoom(true), m_minZoom(10), m_maxZoom(500), m_zoomStep(10), m_zoomPageStep(5 * m_zoomStep)
{
    m_ui->setupUi(this);
#ifdef Q_WS_X11
    m_ui->btn_abort->setIcon(QIcon::fromTheme("dialog-cancel"));
    m_ui->btn_save->setIcon(QIcon::fromTheme("document-save"));
    m_ui->btn_scaleToFit->setIcon(QIcon::fromTheme("zoom-fit-best"));
    m_ui->btn_scaleToFit->setToolTip("Scale to fit");
    m_ui->btn_zoomIn->setIcon(QIcon::fromTheme("zoom-in"));
    m_ui->btn_zoomOut->setIcon(QIcon::fromTheme("zoom-out"));
#else
    m_ui->btn_abort->setIcon(QApplication::style()->standardIcon(QStyle::SP_DialogCancelButton));
    m_ui->btn_save->setIcon(QApplication::style()->standardIcon(QStyle::SP_DriveFDIcon));
    m_ui->btn_scaleToFit->setText("Scale to fit");
    m_ui->btn_zoomIn->setText("+");
    m_ui->btn_zoomOut->setText("-");
#endif
    m_ui->btn_abort->setVisible(false);
    m_ui->btn_save->setVisible(false);
    m_ui->progressBar->setVisible(false);
    setPreviewEnabled(false);
    m_ui->sldr_zoom->setMinimum(m_minZoom);
    m_ui->sldr_zoom->setMaximum(m_maxZoom);
    m_ui->sldr_zoom->setSingleStep(m_zoomStep);
    m_ui->sldr_zoom->setPageStep(m_zoomPageStep);
    m_ui->spn_zoom->setMinimum(m_minZoom);
    m_ui->spn_zoom->setMaximum(m_maxZoom);
    m_ui->scrl_preview->setMinScale(m_minZoom / 100.);
    m_ui->scrl_preview->setMaxScale(m_maxZoom / 100.);
    m_ui->scrl_preview->setScaleStep(m_zoomStep / 100.);
    m_ui->preview->setIntensityScaleFactor(0.1);
    connect(m_ui->btn_zoomIn, SIGNAL(clicked()), m_ui->scrl_preview, SLOT(zoomIn()));
    connect(m_ui->btn_zoomOut, SIGNAL(clicked()), m_ui->scrl_preview, SLOT(zoomOut()));
    connect(m_ui->btn_scaleToFit, SIGNAL(clicked()), m_ui->scrl_preview, SLOT(autoZoom()));
    connect(m_ui->spn_zoom, SIGNAL(editingFinished()), this, SLOT(zoomValueChanged()));
    connect(m_ui->sldr_zoom, SIGNAL(sliderMoved(int)), this, SLOT(zoomSliderValueChanged(int)));
    connect(m_ui->sldr_zoom, SIGNAL(actionTriggered(int)), this, SLOT(zoomSliderActionTriggered(int)));
    connect(m_ui->scrl_preview, SIGNAL(scaleChanged(float, bool)), this, SLOT(scaleChanged(float, bool)));
}

ResultWidget::~ResultWidget()
{
    delete m_ui;
    delete m_worker;
    delete m_functor;
}

void ResultWidget::start(const GuiParams &params)
{
    if (m_functor || m_worker)
        return;
    m_params = params;
    m_ui->preview->setResults(&m_result);
    m_ui->btn_abort->setVisible(true);
    m_ui->progressBar->setVisible(true);
    m_functor = new QtProgressFunctor(this);
    m_worker = new AnalysisWorker(*m_functor, m_params, m_result, this);
    connect(m_functor, SIGNAL(stageChanged(int, const QString&)), this, SLOT(stageChanged(int, const QString&)));
    connect(m_functor, SIGNAL(stackSizeChanged(int)), m_ui->progressBar, SLOT(setMaximum(int)));
    connect(m_functor, SIGNAL(progress(int)), m_ui->progressBar, SLOT(setValue(int)));
    connect(m_ui->btn_abort, SIGNAL(clicked()), m_functor, SLOT(abort()));
    connect(m_worker, SIGNAL(finished()), this, SLOT(workerFinished()));
    connect(m_ui->preview, SIGNAL(detections(const QString&)), m_ui->lbl_detections, SLOT(setText(const QString&)));
    m_ui->bottomLayout->removeItem(m_ui->spcr);
    m_worker->setPriority(QThread::LowPriority);
    connect(m_ui->preview, SIGNAL(initialized()), m_worker, SLOT(start()));
    m_ui->progressBar->setRange(0, 0);
    m_ui->preview->setParams(&params);
    m_ui->scrl_preview->autoZoom();
}

bool ResultWidget::canBeClosed()
{
    if (m_worker->isRunning()) {
        emit userAttentionRequired(this);
        if(QMessageBox::warning(this, "Abort job", "This job is still running. Do you want to abort it?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
            return true;
        else
            return false;
    } else if (!m_coordinatesSaved || !m_reconstructionSaved) {
        emit userAttentionRequired(this);
        QStringList notSaved;
        if (!m_coordinatesSaved)
            notSaved.append("the coordinates file");
        if (!m_reconstructionSaved)
            notSaved.append("the reconstructed image");
        if (QMessageBox::warning(this, "Abort job", QString("You have not saved %1. Do you want to close this tab and discard the results?").arg(notSaved.join(" or ")), QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
            return true;
        else
            return false;
    } else
        return true;
}

bool ResultWidget::canBeDeleted()
{
    if (m_worker->isFinished()) {
       return true;
    } else {
        m_functor->abort();
        return false;
    }
}

void ResultWidget::waitFor()
{
    m_worker->wait();
}

void ResultWidget::stageChanged(int stage, const QString &format)
{
    m_ui->progressBar->setFormat(format);
    m_ui->progressBar->setValue(0);
    if ((WienerStormStage)stage == Localization) {
        setPreviewEnabled(true);
        connect(m_functor, SIGNAL(frameCompleted(int)), m_ui->preview, SLOT(frameCompleted(int)));
    } else
        setPreviewEnabled(false);
}

void ResultWidget::workerFinished()
{
    m_ui->preview->update();
    m_ui->bottomLayout->insertSpacerItem(0, m_ui->spcr);
    m_ui->btn_abort->setVisible(false);
    m_ui->progressBar->setVisible(false);
    m_ui->btn_save->setVisible(true);
    connect(m_ui->btn_save, SIGNAL(clicked()), this, SLOT(saveClicked()));
    if (m_functor->getStage() == Localization)
        m_params.save();
    emit finished(this);
}

void ResultWidget::saveClicked()
{
    QStringList filters;
    filters.append("Coordinate file and reconstructed image (*.txt *.tif)");
    filters.append("Reconstructed image (*. tif)");
    filters.append("Coordinate file (*.txt)");
    QString selectedFilter;
    QString file = QFileDialog::getSaveFileName(this, "Save as", QFileInfo(fileSaveDialogDirectory, QFileInfo(QString::fromStdString(m_params.getCoordsFile())).fileName()).absoluteFilePath(), filters.join(";;"), &selectedFilter);
    QFileInfo info(file);
    fileSaveDialogDirectory = info.absolutePath();
    if (selectedFilter == filters[0]) {
        m_ui->preview->saveImage(info.absolutePath() + '/' + info.baseName() + ".tif");
        saveCoordinates(info.absolutePath() + '/' + info.baseName() + ".txt");
        m_coordinatesSaved = m_reconstructionSaved = true;
    } else if (selectedFilter == filters[1]){
        m_ui->preview->saveImage(file);
        m_reconstructionSaved = true;
    }
    else{
        saveCoordinates(file);
        m_coordinatesSaved = true;
    }
}

void ResultWidget::saveCoordinates(const QString &file)
{
    std::ofstream out(QDir::toNativeSeparators(file).toStdString());
    saveCoordsFile(m_params, out, m_result);
    out.close();
}

void ResultWidget::zoomSliderValueChanged(int val)
{
    m_ui->scrl_preview->setScale(val / 100.);
}

void ResultWidget::zoomSliderActionTriggered(int action)
{
    switch (action) {
        case QAbstractSlider::SliderSingleStepAdd:
            m_ui->scrl_preview->zoomIn();
            break;
        case QAbstractSlider::SliderSingleStepSub:
            m_ui->scrl_preview->zoomOut();
            break;
        case QAbstractSlider::SliderPageStepAdd:
        case QAbstractSlider::SliderPageStepSub:
            m_ui->scrl_preview->setScale(m_ui->sldr_zoom->sliderPosition() / 100.);
            break;
    }
}

void ResultWidget::zoomValueChanged()
{
    m_ui->scrl_preview->setScale(m_ui->spn_zoom->value() / 100.);
}

void ResultWidget::scaleChanged(float scale, bool autoZoom)
{
    scaleChanged((int)(scale * 100), autoZoom);
}

void ResultWidget::scaleChanged(int scale, bool autoZoom)
{
    m_autoZoom = autoZoom;
    m_ui->sldr_zoom->setValue(scale);
    m_ui->spn_zoom->setValue(scale);
    m_ui->btn_zoomOut->setEnabled(m_previewEnabled && scale > m_minZoom);
    m_ui->btn_zoomIn->setEnabled(m_previewEnabled && scale < m_maxZoom);
    m_ui->btn_scaleToFit->setEnabled(m_previewEnabled && !autoZoom);
}

void ResultWidget::resizeEvent(QResizeEvent *event)
{
    QSize ps = m_ui->preview->size(), vs = m_ui->scrl_preview->maximumViewportSize();
    if (vs.width() > ps.width() && vs.height() > ps.height())
        m_autoZoom = true;
    if (m_autoZoom)
        m_ui->scrl_preview->autoZoom();
    QWidget::resizeEvent(event);
}

void ResultWidget::setPreviewEnabled(bool enabled)
{
    m_previewEnabled = enabled;
    m_ui->btn_scaleToFit->setEnabled(enabled && !m_autoZoom);
    m_ui->btn_zoomIn->setEnabled(enabled && m_ui->scrl_preview->scale() * 100 < m_maxZoom);
    m_ui->btn_zoomOut->setEnabled(enabled && m_ui->scrl_preview->scale() * 100 > m_minZoom);
    m_ui->sldr_zoom->setEnabled(enabled);
    m_ui->spn_zoom->setEnabled(enabled);
    m_ui->scrl_preview->setEnabled(enabled);
}
