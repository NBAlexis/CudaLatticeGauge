//=============================================================================
// FILENAME : CGaugeFixingRandom.h
//
// DESCRIPTION:
//
// This is a gauge transform to test gauge invarience of observables.
//
//
// REVISION:
//  [09/25/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXINGRANDOM_H_
#define _CGAUGEFIXINGRANDOM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CGaugeFixingRandom)

class CLGAPI CGaugeFixingRandom : public CGaugeFixing
{
    __CLGDECLARE_CLASS(CGaugeFixingRandom)
public:

    CGaugeFixingRandom()
    : CGaugeFixing()
    , m_pGSU2(NULL)
    , m_pG(NULL)
#if _CLG_MULTI_GPU
    , m_pSavedGSU2(NULL)
    , m_pSavedG(NULL)
#endif
    {
    }

    ~CGaugeFixingRandom()
    {
        cudaSafeFree(m_pGSU2);
        cudaSafeFree(m_pG);
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeFixing(CFieldGauge* pResGauge) override;
    DOUBLE CheckRes(const CFieldGauge* ) override { return 0.0; }

    CCString GetInfos(const CCString& sTab) const override;

    /**
     * Call fix fermion just after gaugefixing, will use the same gauge transform as before
     * When never gaugefixing is called the gauge transform is randomized
     */
    void AlsoFixingFermion(CFieldFermion* pFermion) const;
    void AlsoFixingFermionWilsonSU3(CFieldFermionWilsonSquareSU3 * pFermion) const;
    void AlsoFixingFermionKSSU3(CFieldFermionKSSU3* pFermion) const;
    void AlsoFixingAphys(CFieldGauge* pGauge) const;
    deviceSU2* m_pGSU2;
    deviceSU3* m_pG;

#if _CLG_MULTI_GPU
    //P4-2.5 fix: temporary global-lattice buffers (rank 0 only). m_pGSU2/m_pG are
    //sized to the LOCAL lattice in Initial(); under the temporary GLOBAL context
    //they are re-allocated to the global volume and restored on exit. The RNG's
    //curand states are switched in parallel via CRandom::EnterGlobalContext()
    //(see Tools/Math/Random.cu), because _deviceRandomF indexes them by GLOBAL
    //site index there.
    void ResizeBuffersToGlobal();
    void RestoreLocalBuffers();
    deviceSU2* m_pSavedGSU2;
    deviceSU3* m_pSavedG;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGRANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================
