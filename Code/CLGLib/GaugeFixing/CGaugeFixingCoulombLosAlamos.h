//=============================================================================
// FILENAME : CGaugeFixingCoulombLosAlamos.h
// 
// DESCRIPTION:
//
// 
// 
//
// REVISION:
//  [09/23/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGCOULOMBLOSALAMOS_H_
#define _CGAUGEFIXINGCOULOMBLOSALAMOS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingCoulombLosAlamos)

class CLGAPI CGaugeFixingCoulombLosAlamos : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingCoulombLosAlamos)
public:

    CGaugeFixingCoulombLosAlamos()
    : CGaugeFixing()
    , m_pDDecomp(NULL)
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

    ~CGaugeFixingCoulombLosAlamos()
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
    void GaugeFixingForT(deviceSU3* pResGauge, SBYTE uiT);
    Real CheckRes(const CFieldGauge* pGauge) override;
    Real CheckResDeviceBuffer(const deviceSU3* __restrict__ pGauge);
    Real CheckResDeviceBufferOnlyT(const deviceSU3* __restrict__ pGauge, SBYTE uiT);
    CCString GetInfos(const CCString& sTab) const override;

    UINT m_pHDecomp[6];
    UINT* m_pDDecomp;

    Real m_fOmega;
    UINT m_iCheckErrorStep;
    deviceSU3* m_pG;
    Real* m_pA11;
    CLGComplex* m_pA12;
    CLGComplex* m_pA13;
    Real* m_pA22;
    CLGComplex* m_pA23;

    TArray<INT> m_lstDims;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGCOULOMBLOSALAMOS_H_

//=============================================================================
// END OF FILE
//=============================================================================