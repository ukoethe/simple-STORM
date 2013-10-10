#ifndef GUIPARAMS_H
#define GUIPARAMS_H

#include "dataparams.h"

class GuiParams : public DataParams {
public:
    GuiParams();
    ~GuiParams();

    bool getAdvancedSettingsEnabled() const;
    void setAdvancedSettingsEnabled(bool);
    float getReconstructionResolution() const;
    void setReconstructionResolution(float);
	float getReconstructionResolutionDefaults() const;
    virtual void save() const;
	bool getReconstructionResolutionDefaultsSet() const;
	bool getReconstructionResolutionFromSettingsFile() const;

protected:
    virtual void loadSettings(bool);

private:
    bool m_advancedSettingsEnabled;
    float m_reconstructionResolution;
	float m_reconstructionResolutionDefaults;
	bool m_reconstructionResolutionDefaultsSet;
	bool m_reconstructionResolutionFromSettingsFile;
    static const std::string s_section;
};

#endif
