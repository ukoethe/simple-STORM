#include "guiparams.h"
#include <iostream>
#include <rude/config.h>

const std::string GuiParams::s_section = "guiparams";

GuiParams::GuiParams()
: DataParams(), m_advancedSettingsEnabled(false), m_reconstructionResolution((float)13.4)
{}

GuiParams::~GuiParams() {}

bool GuiParams::getAdvancedSettingsEnabled() const
{
    return m_advancedSettingsEnabled;
}

void GuiParams::setAdvancedSettingsEnabled(bool enabled)
{
    m_advancedSettingsEnabled = enabled;
}

float GuiParams::getReconstructionResolution() const
{
    return m_reconstructionResolution;
}

void GuiParams::setReconstructionResolution(float res)
{
    m_reconstructionResolution = res;
}

float GuiParams::getReconstructionResolutionDefaults() const
{
	return m_reconstructionResolutionDefaults;
}

bool GuiParams::getReconstructionResolutionDefaultsSet() const
{
	return m_reconstructionResolutionDefaultsSet;
}

bool GuiParams::getReconstructionResolutionFromSettingsFile() const
{
	return m_reconstructionResolutionFromSettingsFile;
}

void GuiParams::save() const {
    m_config->setSection(s_section.c_str());
    m_config->setBoolValue("advancedsettings", m_advancedSettingsEnabled);
    m_config->setDoubleValue("reconstructionresolution", m_reconstructionResolution);
    DataParams::save();
}

void GuiParams::loadSettings(bool propagate) {
    m_config->setSection(s_section.c_str());
    if (m_config->exists("advancedsettings")) {
        m_advancedSettingsEnabled = m_config->getBoolValue("advancedsettings");
    }
    if (m_config->exists("reconstructionresolution")) {
        m_reconstructionResolution = m_config->getDoubleValue("reconstructionresolution");
		m_reconstructionResolutionFromSettingsFile = true;
    }
	else
		m_reconstructionResolutionFromSettingsFile = false;
    if (propagate)
        DataParams::loadSettings(propagate);
	
	m_configDefaults->setConfigFile((getDefaultsFileFilename()).c_str());
	m_configDefaults->load();
	if (m_configDefaults->exists("reconstructionresolution"))
	{	
		m_reconstructionResolutionDefaults = m_configDefaults->getDoubleValue("reconstructionresolution");
		m_reconstructionResolutionDefaultsSet = true;
	}
	else
		m_reconstructionResolutionDefaultsSet = false;
	
}
