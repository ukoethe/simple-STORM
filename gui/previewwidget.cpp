#include "previewwidget.h"
#include "ui_previewwidget.h"

#include <QDir>
#include <QPainter>

PreviewWidget::PreviewWidget(QWidget *parent)
: QWidget(parent), m_ui(new Ui::PreviewWidget()), m_lastProcessed(std::chrono::steady_clock::now()), m_updateInterval(std::chrono::milliseconds(1000)), m_detections(0), m_paramsSet(false), m_zoom(-1)
{
    m_ui->setupUi(this);
    m_ui->btn_scaleToFit->setEnabled(false);
    m_ui->sldr_zoom->setMinimum(10);
    m_ui->sldr_zoom->setMaximum(500);
    m_ui->sldr_zoom->setSingleStep(10);
    m_ui->spn_zoom->setMinimum(10);
    m_ui->spn_zoom->setMaximum(500);
    connect(m_ui->btn_zoomIn, SIGNAL(clicked()), this, SLOT(zoomInClicked()));
    connect(m_ui->btn_zoomOut, SIGNAL(clicked()), this, SLOT(zoomOutClicked()));
    connect(m_ui->btn_scaleToFit, SIGNAL(clicked()), this, SLOT(autoZoom()));
    connect(m_ui->spn_zoom, SIGNAL(editingFinished()), this, SLOT(zoomValueChanged()));
    connect(m_ui->scrl_preview, SIGNAL(zoomIn()), this, SLOT(zoomInClicked()));
    connect(m_ui->scrl_preview, SIGNAL(zoomOut()), this, SLOT(zoomOutClicked()));
#ifdef Q_WS_X11
    m_ui->btn_scaleToFit->setIcon(QIcon::fromTheme("zoom-fit-best"));
    m_ui->btn_scaleToFit->setToolTip("Scale to fit");
    m_ui->btn_zoomIn->setIcon(QIcon::fromTheme("zoom-in"));
    m_ui->btn_zoomOut->setIcon(QIcon::fromTheme("zoom-out"));
#else
    m_ui->btn_scaleToFit->setText("Scale to fit");
    m_ui->btn_zoomIn->setText("+");
    m_ui->btn_zoomOut->setText("-");
#endif
}

PreviewWidget::~PreviewWidget()
{
    delete m_ui;
}

void PreviewWidget::setResults(const std::vector<std::set<Coord<float>>> *results)
{
    m_results = results;
}

void PreviewWidget::setUpdateInterval(const std::chrono::milliseconds &interval)
{
    m_updateInterval = interval;
}

void PreviewWidget::setParams(const GuiParams *params)
{
    m_params = params;
    m_paramsSet = true;

    m_result.reshape(vigra::Shape2(m_params->shape(0) * m_params->getPixelSize() / m_params->getReconstructionResolution(), m_params->shape(1) * m_params->getPixelSize() / m_params->getReconstructionResolution()), 0.0);
    m_factor = m_params->getPixelSize() / (m_params->getReconstructionResolution() * m_params->getFactor());
    m_quantileIndex = m_result.size() * 0.996;
    for (const double &d : m_result)
        m_resultsForScaling.insert(d);
    m_pixmap = QPixmap(m_result.shape(0), m_result.shape(1));
    m_painter.begin(&m_pixmap);
    autoZoom();
}

void PreviewWidget::autoZoom()
{
    m_zoom = -1;
    disconnect(m_ui->sldr_zoom, SIGNAL(valueChanged(int)), this, SLOT(zoomSliderValueChanged(int)));
    m_ui->btn_scaleToFit->setEnabled(false);
    QSize size(m_result.shape(0), m_result.shape(1));
    size.scale(m_ui->scrl_preview->maximumViewportSize(), Qt::KeepAspectRatio);
    float zoom = size.height() / (float)m_result.shape(1);
    m_ui->sldr_zoom->setValue(zoom * 100);
    m_ui->spn_zoom->setValue(zoom * 100);
    m_ui->lbl_preview->setPixmap(m_pixmap.scaled(size, Qt::KeepAspectRatio, Qt::FastTransformation));
    connect(m_ui->sldr_zoom, SIGNAL(valueChanged(int)), this, SLOT(zoomSliderValueChanged(int)));
}

void PreviewWidget::frameCompleted(int frame)
{
    m_unprocessed.push_back(frame);
    std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    if (now - m_lastProcessed >= m_updateInterval) {
        update();
    }
}

void PreviewWidget::update()
{
    if (!m_paramsSet)
        return;
    m_lastProcessed = std::chrono::steady_clock::now();
    for (int frame : m_unprocessed) {
        for (const Coord<float> &c : m_results->at(frame)) {
            double val = m_result(std::round(c.x * m_factor), std::round(c.y * m_factor));
            m_resultsForScaling.erase(m_resultsForScaling.find(val));
            val += c.val;
            m_resultsForScaling.insert(val);
            m_result(std::round(c.x * m_factor), std::round(c.y * m_factor)) = val;
            ++m_detections;
        }
    }
    m_unprocessed.clear();
    emit detections(QString("%1 detections").arg(m_detections));

    // some maxima are very strong so we scale the image as appropriate :
    double minlim = *(m_resultsForScaling.begin()), maxlim = 0;
    auto it = m_resultsForScaling.end();
    for (size_t i = m_resultsForScaling.size(); i >=m_quantileIndex; --i, --it)
        maxlim = *it;
    double factor = 255 / (maxlim - minlim);
    for (int y = 0; y < m_result.shape(1); ++y) {
        for (int x = 0; x < m_result.shape(0); ++x) {
            double p = m_result(x, y);
            uchar v = (p > maxlim) ? 255 : p * factor;
            m_painter.fillRect(x, y, 1, 1, QColor(v, v, v));
        }
    }
    QSize scaled = (m_zoom == -1) ? m_ui->scrl_preview->maximumViewportSize() : (m_pixmap.size() * m_zoom);
    m_ui->lbl_preview->setPixmap(m_pixmap.scaled(scaled, Qt::KeepAspectRatio, Qt::FastTransformation));
}

void PreviewWidget::saveImage(const QString &file)
{
    double maxval = std::numeric_limits<double>::min();
    for (const double &val : m_result) {
        if (val > maxval)
            maxval = val;
    }
    vigra::MultiArray<2, float> scaled(m_result.shape());
    if (maxval > std::numeric_limits<float>::max()) {
        float factor = std::numeric_limits<float>::max() / maxval;
        vigra::transformMultiArray(srcMultiArrayRange(m_result), destMultiArray(scaled), [&factor](double p){return p * factor;});
    } else
        vigra::transformMultiArray(srcMultiArrayRange(m_result), destMultiArray(scaled), [](double p){return (float)p;});
    vigra::exportImage(srcImageRange(scaled), QDir::toNativeSeparators(file).toStdString().c_str());
}

void PreviewWidget::resizeEvent(QResizeEvent *event)
{
    Q_UNUSED(event);
    if (m_zoom == -1)
        autoZoom();
}

void PreviewWidget::zoomInClicked()
{
    m_ui->sldr_zoom->triggerAction(QAbstractSlider::SliderSingleStepAdd);
}

void PreviewWidget::zoomOutClicked()
{
    m_ui->sldr_zoom->triggerAction(QAbstractSlider::SliderSingleStepSub);
}

void PreviewWidget::zoomSliderValueChanged(int val)
{
    m_ui->spn_zoom->setValue(val);
    if (val == m_ui->sldr_zoom->minimum())
        m_ui->btn_zoomOut->setEnabled(false);
    else
        m_ui->btn_zoomOut->setEnabled(true);
    if (val == m_ui->sldr_zoom->maximum())
        m_ui->btn_zoomIn->setEnabled(false);
    else
        m_ui->btn_zoomOut->setEnabled(true);
    m_zoom = val / 100.;
    m_ui->lbl_preview->setPixmap(m_pixmap.scaled(m_pixmap.size() * m_zoom, Qt::KeepAspectRatio, Qt::FastTransformation));
    m_ui->btn_scaleToFit->setEnabled(true);
}

void PreviewWidget::zoomValueChanged()
{
    m_ui->sldr_zoom->setValue(m_ui->spn_zoom->value());
}
