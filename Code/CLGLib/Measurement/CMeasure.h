//=============================================================================
// FILENAME : CMeasure.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASURE_H_
#define _CMEASURE_H_

__BEGIN_NAMESPACE

class CLGAPI CMeasure : public CBase
{
public:
    CMeasure() : m_pOwner(NULL), m_fLastRealResult(F(0.0)) {}

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
    {
        m_pOwner = pOwner;
        m_pLatticeData = pLatticeData;
        m_byId = byId;
    }

    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) = 0;
    virtual void Average(UINT uiConfigurationCount) = 0;
    virtual void Report() = 0;
    virtual void Reset() = 0;

protected:

    class CMeasurementManager* m_pOwner;
    class CLatticeData* m_pLatticeData;
    BYTE m_byId;

public:
    //============================================================
    //some simple measurement only produce real or complex results
    Real m_fLastRealResult;
    CLGComplex m_cLastComplexResult;
    
};

__END_NAMESPACE

#endif //#ifndef _CMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================