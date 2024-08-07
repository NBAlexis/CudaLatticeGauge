//=============================================================================
// FILENAME : CGaugeFixingCoulombCornell.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/20/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGCOULOMBCORNELL_H_
#define _CGAUGEFIXINGCOULOMBCORNELL_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingCoulombCornell)

class CLGAPI CGaugeFixingCoulombCornell : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingCoulombCornell)
public:

    CGaugeFixingCoulombCornell()
        : CGaugeFixing()
        , m_fAlpha(F(0.08))
        , m_pDDecomp(NULL)
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

    ~CGaugeFixingCoulombCornell()
    {
        cudaSafeFree(m_pDDecomp);
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
    DOUBLE CheckRes(const CFieldGauge* pGauge) override;
    void GaugeFixingOneTimeSlice(deviceSU3* pResGauge, SBYTE uiT, BYTE byFieldId);
    
    CCString GetInfos(const CCString& sTab) const override;

    DOUBLE m_fAlpha;

    //device SU3 is not alligned, therefor use CLGComplex*
    UINT m_pHDecomp[6];
    UINT* m_pDDecomp;

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

    //FFT accelaration not support now
    UBOOL m_bFA;
    TArray<INT> m_lstDims;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGCOULOMBCORNELL_H_

//=============================================================================
// END OF FILE
//=============================================================================