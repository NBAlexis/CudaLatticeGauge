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
#if _CLG_MULTI_GPU
    , m_pSavedA11(NULL)
    , m_pSavedA12(NULL)
    , m_pSavedA13(NULL)
    , m_pSavedA22(NULL)
    , m_pSavedA23(NULL)
    , m_pSavedGamma11(NULL)
    , m_pSavedGamma12(NULL)
    , m_pSavedGamma13(NULL)
    , m_pSavedGamma22(NULL)
    , m_pSavedGamma23(NULL)
    , m_pSavedG(NULL)
    , m_pSavedMomentumTable(NULL)
    , m_pSavedTempFFTBuffer(NULL)
#endif
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
    DOUBLE CheckRes(const CFieldGauge* pGauge) override;

    CCString GetInfos(const CCString& sTab) const override;

    //P4-2.4: core deviation loop over one gauge buffer (local or gathered-global).
    DOUBLE CheckResLocal(const deviceSU3* pGaugeData, BYTE byFieldId);

    //P4-2.4: the 4D iteration loop, run over a local or a gathered-global buffer.
    void GaugeFixingLoop(deviceSU3* pDeviceBufferPointer, BYTE byFieldId);

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


    //FFT accelaration
    UBOOL m_bFA;

#if _CLG_MULTI_GPU
    //P4-2.4: temporary global-lattice fixing buffers (rank 0 only), mirroring
    //P4-2.3/P4-2.4: the 4D-volume buffers above are sized to the LOCAL lattice in
    //Initial(); under the temporary GLOBAL context they are re-allocated to the
    //global volume, the momentum table is re-baked, then restored on exit.
    void ResizeBuffersToGlobal();
    void RestoreLocalBuffers();
    DOUBLE* m_pSavedA11;
    cuDoubleComplex* m_pSavedA12;
    cuDoubleComplex* m_pSavedA13;
    DOUBLE* m_pSavedA22;
    cuDoubleComplex* m_pSavedA23;
    DOUBLE* m_pSavedGamma11;
    cuDoubleComplex* m_pSavedGamma12;
    cuDoubleComplex* m_pSavedGamma13;
    DOUBLE* m_pSavedGamma22;
    cuDoubleComplex* m_pSavedGamma23;
    deviceSU3* m_pSavedG;
    DOUBLE* m_pSavedMomentumTable;
    cuDoubleComplex* m_pSavedTempFFTBuffer;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGLANDAUCORNELL_H_

//=============================================================================
// END OF FILE
//=============================================================================