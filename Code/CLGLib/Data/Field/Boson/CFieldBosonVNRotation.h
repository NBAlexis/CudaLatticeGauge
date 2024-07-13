//=============================================================================
// FILENAME : CFieldBosonVNRotation.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [07/13/2024 nbale]
//=============================================================================
#include "CFieldBosonVN.h"

#ifndef _CFIELDBOSONVNROTATION_H_
#define _CFIELDBOSONVNROTATION_H_

__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge>
class __DLL_EXPORT CFieldBosonVNRotation : public CFieldBosonVN<deviceDataBoson, deviceDataGauge>
{
public:
    CFieldBosonVNRotation();
    ~CFieldBosonVNRotation();

    void CopyTo(CField* U) const override;

    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const override;

    void InitialOtherParameters(CParameters& param) override;
    CCString GetInfos(const CCString& tab) const override;

    Real m_fOmega;
    INT* m_pDevicePath;
    _deviceCoeffFunctionPointer m_pfCx;

protected:

    void DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
};

__CLG_REGISTER_HELPER_HEADER(CFieldBosonVNRotationU1)
class CLGAPI CFieldBosonVNRotationU1 : public CFieldBosonVNRotation<CLGComplex, CLGComplex>
{ 
__CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldBosonVNRotationU1)
public: 
    EFieldType GetFieldType() const override { return EFT_BosonComplex; }
    UINT VectorN() const override { return 1; } 
    UINT FloatN() const override { return 2; } 
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVNROTATION_H_

//=============================================================================
// END OF FILE
//=============================================================================