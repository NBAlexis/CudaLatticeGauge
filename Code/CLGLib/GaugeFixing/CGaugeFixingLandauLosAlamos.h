//=============================================================================
// FILENAME : CGaugeFixingLandauLosAlamos.h
// 
// DESCRIPTION:
//
// Note!! this does not pass the test,
// I must have some misunderstanding of the Los Alamos gauge fixing!!
//
// REVISION:
//  [09/21/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGLANDAUCLOSALAMOS_H_
#define _CGAUGEFIXINGLANDAUCLOSALAMOS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingLandauLosAlamos)

class CLGAPI CGaugeFixingLandauLosAlamos : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingLandauLosAlamos)
public:

    CGaugeFixingLandauLosAlamos()
    : CGaugeFixing()
    , m_fOmega(F(1.0))
    , m_iCheckErrorStep(1000)
    , m_pG(NULL)
    , m_pA11(NULL)
    , m_pA12(NULL)
    , m_pA13(NULL)
    , m_pA22(NULL)
    , m_pA23(NULL)
    {
    }

    ~CGaugeFixingLandauLosAlamos()
    {
        cudaSafeFree(m_pG);
        cudaSafeFree(m_pA11);
        cudaSafeFree(m_pA12);
        cudaSafeFree(m_pA13);
        cudaSafeFree(m_pA22);
        cudaSafeFree(m_pA23);
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeFixing(CFieldGauge* pResGauge) override;
    Real CheckRes(const CFieldGauge* pGauge) override;
    CCString GetInfos(const CCString& sTab) const override;

    Real m_fOmega;
    UINT m_iCheckErrorStep;
    deviceSU3* m_pG;
    Real* m_pA11;
    CLGComplex* m_pA12;
    CLGComplex* m_pA13;
    Real* m_pA22;
    CLGComplex* m_pA23;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGLANDAUCLOSALAMOS_H_

//=============================================================================
// END OF FILE
//=============================================================================