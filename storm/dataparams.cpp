#include "dataparams.h"
#include <iostream>
#include <rude/config.h>

const std::string DataParams::s_section = "dataparams";

DataParams::DataParams()
: StormParams(), m_slopeSaved(false), m_interceptSaved(false), m_sigmaSaved(false)
{}

DataParams::DataParams(int argc, char **argv)
: StormParams(argc, argv), m_slopeSaved(false), m_interceptSaved(false), m_sigmaSaved(false)
{
std::cout<<"load triggered"<<std::endl;
    load(false);
	std::cout<<"load triggered"<<std::endl;
}

DataParams::~DataParams() {}

float DataParams::getSlope() const {
    return m_slope;
}
bool DataParams::getSlopeSaved() const {
    return m_slopeSaved;
}
void DataParams::setSlope(float slope) {
    m_slope = slope;
}
void DataParams::setUseSavedSlope(bool saved) {
    m_slopeSaved = saved;
}
float DataParams::getIntercept() const {
    return m_intercept;
}
bool DataParams::getInterceptSaved() const {
    return m_interceptSaved;
}
void DataParams::setIntercept(float intercept) {
    m_intercept = intercept;
}
void DataParams::setUseSavedIntercept(bool saved) {
    m_interceptSaved = saved;
}
float DataParams::getSigma() const {
    return m_sigma;
}
bool DataParams::getSigmaSaved() const {
    return m_sigmaSaved;
}
void DataParams::setSigma(float sigma) {
    m_sigma = sigma;
}
void DataParams::setUseSavedSigma(bool saved) {
    m_sigmaSaved = saved;
}

void DataParams::save() const {
    m_config->setSection(s_section.c_str());
    m_config->setDoubleValue("slope", m_slope);
    m_config->setDoubleValue("intercept", m_intercept);
    m_config->setDoubleValue("sigma", m_sigma);
    StormParams::save();
}

void DataParams::loadSettings(bool propagate) {
    m_config->setSection(s_section.c_str());
    if (m_config->exists("slope")) {
        m_slope = m_config->getDoubleValue("slope");
        m_slopeSaved = true;
    }
    if (m_config->exists("intercept")) {
        m_intercept = m_config->getDoubleValue("intercept");
        m_interceptSaved = true;
    }
    if (m_config->exists("sigma")) {
        m_sigma = m_config->getDoubleValue("sigma");
        m_sigmaSaved = true;
    }
    if (propagate)
        StormParams::loadSettings(propagate);
}
