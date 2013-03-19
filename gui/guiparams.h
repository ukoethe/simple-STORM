#ifndef GUIPARAMS_H
#define GUIPARAMS_H

#include "dataparams.h"

class GuiParams : public DataParams {
public:
    GuiParams();
    ~GuiParams();

    bool getAdvancedSettingsEnabled() const;
    void setAdvancedSettingsEnabled(bool);
    int getReconstructionResolution() const;
    void setReconstructionResolution(int);

    virtual void save() const;

protected:
    virtual void loadSettings(bool);

private:
    bool m_advancedSettingsEnabled;
    int m_reconstructionResolution;

    static const std::string s_section;
};

#endif
