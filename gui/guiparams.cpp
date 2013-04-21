#include "guiparams.h"

#include <rude/config.h>

const std::string GuiParams::s_section = "guiparams";

GuiParams::GuiParams()
: DataParams(), m_advancedSettingsEnabled(false), m_reconstructionResolution(1)
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

void GuiParams::save() const {
    m_config->setSection(s_section.c_str());
    m_config->setBoolValue("advancedsettings", m_advancedSettingsEnabled);
    m_config->setIntValue("reconstructionresolution", m_reconstructionResolution);
    DataParams::save();
}

void GuiParams::loadSettings(bool propagate) {
    m_config->setSection(s_section.c_str());
    if (m_config->exists("advancedsettings")) {
        m_advancedSettingsEnabled = m_config->getBoolValue("advancedsettings");
    }
    if (m_config->exists("reconstructionresolution")) {
        m_reconstructionResolution = m_config->getIntValue("reconstructionresolution");
    }
    if (propagate)
        DataParams::loadSettings(propagate);
}
