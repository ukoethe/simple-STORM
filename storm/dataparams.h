#ifndef DATAPARAMS_H
#define DATAPARAMS_H

#include "stormparams.h"

class DataParams : public StormParams {
public:
    DataParams();
    DataParams(int argc, char **argv);
    ~DataParams();

    float getSlope() const;
    bool getSlopeSaved() const;
    void setSlope(float);
    float getIntercept() const;
    bool getInterceptSaved() const;
    void setIntercept(float);
    float getSigma() const;
    bool getSigmaSaved() const;
    void setSigma(float);

    virtual void save() const;
    void load();

private:
    float m_slope;
    bool m_slopeSaved;
    float m_intercept;
    bool m_interceptSaved;
    float m_sigma;
    bool m_sigmaSaved;

    static const std::string s_section;
};

#endif
