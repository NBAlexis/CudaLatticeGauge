//=============================================================================
// FILENAME : CGaugeFixingLandauCornell.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/18/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGLANDAUCORNELL_H_
#define _CGAUGEFIXINGLANDAUCORNELL_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingLandauCornell)

class CLGAPI CGaugeFixingLandauCornell : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingLandauCornell)
public:

    CGaugeFixingLandauCornell()
    : CGaugeFixing()
    , m_fAlpha(F(0.08))
    , m_pA11(NULL)
    , m_pA12(NULL)
    , m_pA13(NULL)
    , m_pA22(NULL)
    , m_pA23(NULL)
    , m_pGamma11(NULL)
    , m_pGamma12(NULL)
    , m_pGamma13(NULL)
    , m_pGamma22(NULL)
    , m_pGamma23(NULL)
    , m_pG(NULL)
    , m_pMomentumTable(NULL)
    , m_pTempFFTBuffer(NULL)
    , m_bFA(TRUE)
    {
    }

    ~CGaugeFixingLandauCornell()
    {
        cudaSafeFree(m_pA11);
        cudaSafeFree(m_pA12);
        cudaSafeFree(m_pA13);
        cudaSafeFree(m_pA22);
        cudaSafeFree(m_pA23);
        cudaSafeFree(m_pGamma11);
        cudaSafeFree(m_pGamma12);
        cudaSafeFree(m_pGamma13);
        cudaSafeFree(m_pGamma22);
        cudaSafeFree(m_pGamma23);
        cudaSafeFree(m_pG);
        cudaSafeFree(m_pMomentumTable);
        cudaSafeFree(m_pTempFFTBuffer);
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeFixing(CFieldGauge* pResGauge) override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CheckRes(const CFieldGauge* pGauge) override;
#else
    Real CheckRes(const CFieldGauge* pGauge) override;
#endif
    CCString GetInfos(const CCString& sTab) const override;

#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fAlpha;

    //device SU3 is not alligned, therefor use CLGComplex*
    DOUBLE* m_pA11;
    cuDoubleComplex* m_pA12;
    cuDoubleComplex* m_pA13;
    DOUBLE* m_pA22;
    cuDoubleComplex* m_pA23;
    DOUBLE* m_pGamma11;
    cuDoubleComplex* m_pGamma12;
    cuDoubleComplex* m_pGamma13;
    DOUBLE* m_pGamma22;
    cuDoubleComplex* m_pGamma23;
    deviceSU3* m_pG;
    DOUBLE* m_pMomentumTable;
    cuDoubleComplex* m_pTempFFTBuffer;
#else
    Real m_fAlpha;

    //device SU3 is not alligned, therefor use CLGComplex*
    Real* m_pA11;
    CLGComplex* m_pA12;
    CLGComplex* m_pA13;
    Real* m_pA22;
    CLGComplex* m_pA23;
    Real* m_pGamma11;
    CLGComplex* m_pGamma12;
    CLGComplex* m_pGamma13;
    Real* m_pGamma22;
    CLGComplex* m_pGamma23;
    deviceSU3* m_pG;
    Real* m_pMomentumTable;
    CLGComplex* m_pTempFFTBuffer;
#endif

    //FFT accelaration
    UBOOL m_bFA;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGLANDAUCORNELL_H_

//=============================================================================
// END OF FILE
//=============================================================================