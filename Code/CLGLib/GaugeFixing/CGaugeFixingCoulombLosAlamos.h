//=============================================================================
// FILENAME : CGaugeFixingCoulombLosAlamos.h
// 
// DESCRIPTION:
//
// It seems the Los Alamos cannot fix gauge with logarithm definition
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
#if _CLG_MULTI_GPU
    , m_pSavedG(NULL)
    , m_pSavedA11(NULL)
    , m_pSavedA12(NULL)
    , m_pSavedA13(NULL)
    , m_pSavedA22(NULL)
    , m_pSavedA23(NULL)
#endif
    , m_bMixed(FALSE)
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
    void GaugeFixingForT(deviceSU3* pResGauge, SCHAR uiT, BYTE byFieldId);

    DOUBLE CheckRes(const CFieldGauge* pGauge) override;
    DOUBLE CheckResDeviceBuffer(const deviceSU3* pGauge, BYTE byFieldId);
    DOUBLE CheckResDeviceBufferOnlyT(const deviceSU3* pGauge, SCHAR uiT, BYTE byFieldId);

    CCString GetInfos(const CCString& sTab) const override;

#if _CLG_MULTI_GPU
    //P4-2.2: temporary global-lattice fixing buffers (rank 0 only). The fixing
    //buffers below are sized to the LOCAL lattice in Initial(); under the
    //temporary GLOBAL context (see CGaugeFixing::MGEnterGlobalFixerContext) the
    //kernels sweep the whole lattice, so the buffers must be re-allocated to the
    //global 3D volume and restored afterwards.
    void ResizeBuffersToGlobal();
    void RestoreLocalBuffers();
    deviceSU3* m_pSavedG;
    Real* m_pSavedA11;
    CLGComplex* m_pSavedA12;
    CLGComplex* m_pSavedA13;
    Real* m_pSavedA22;
    CLGComplex* m_pSavedA23;
#endif

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

    UBOOL m_bMixed;

    TArray<INT> m_lstDims;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGCOULOMBLOSALAMOS_H_

//=============================================================================
// END OF FILE
//=============================================================================