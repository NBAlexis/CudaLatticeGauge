//=============================================================================
// FILENAME : CGaugeFixingLandauLosAlamos.h
// 
// DESCRIPTION:
// It seems the Los Alamos cannot fix gauge with logarithm definition
// 
// 
//
// REVISION:
//  [09/21/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGLANDAULOSALAMOS_H_
#define _CGAUGEFIXINGLANDAULOSALAMOS_H_

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
#if _CLG_MULTI_GPU
    , m_pSavedG(NULL)
    , m_pSavedA11(NULL)
    , m_pSavedA12(NULL)
    , m_pSavedA13(NULL)
    , m_pSavedA22(NULL)
    , m_pSavedA23(NULL)
#endif
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
    DOUBLE CheckRes(const CFieldGauge* pGauge) override;

    CCString GetInfos(const CCString& sTab) const override;

    //P4-2.4: core deviation loop over one gauge buffer (local or gathered-global).
    DOUBLE CheckResLocal(const deviceSU3* pGaugeData, BYTE byFieldId);

    //P4-2.4: the 4D iteration loop, run over a local or a gathered-global buffer.
    void GaugeFixingLoop(deviceSU3* pDeviceBufferPointer, BYTE byFieldId);

    Real m_fOmega;
    UINT m_iCheckErrorStep;
    deviceSU3* m_pG;
    Real* m_pA11;
    CLGComplex* m_pA12;
    CLGComplex* m_pA13;
    Real* m_pA22;
    CLGComplex* m_pA23;

#if _CLG_MULTI_GPU
    //P4-2.4: temporary global-lattice fixing buffers (rank 0 only), mirroring
    //P4-2.2/P4-2.3: the 4D-volume buffers above are sized to the LOCAL lattice in
    //Initial(); under the temporary GLOBAL context they are re-allocated to the
    //global volume and restored on exit.
    void ResizeBuffersToGlobal();
    void RestoreLocalBuffers();
    deviceSU3* m_pSavedG;
    Real* m_pSavedA11;
    CLGComplex* m_pSavedA12;
    CLGComplex* m_pSavedA13;
    Real* m_pSavedA22;
    CLGComplex* m_pSavedA23;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGLANDAULOSALAMOS_H_

//=============================================================================
// END OF FILE
//=============================================================================