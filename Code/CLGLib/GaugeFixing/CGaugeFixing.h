//=============================================================================
// FILENAME : CGaugeFixing.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/18/2019 nbale]
//=============================================================================

#ifndef _CGAUGEFIXING_H_
#define _CGAUGEFIXING_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EGaugeFixingType,
    EGFT_LandauCornell,
    EGFT_LandauCornellFFT,
    EGFT_CoulombCornellFFT,

    EGFT_ForceDWORD = 0x7fffffff,
    )


class CLGAPI CGaugeFixing : public CBase
{
public:

    CGaugeFixing()
    : m_pOwner(NULL)
#if !_CLG_DOUBLEFLOAT
    , m_fAccuracy(F(0.00001))
#else
    , m_fAccuracy(F(0.00000000001))
#endif
    , m_iIterate(0)
    , m_iMaxIterate(1000000)
    , m_iLinearIterate(20)
    , m_iShowErrorStep(1000)
    {
        
    }

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params) = 0;

    /**
     * pResGauge should be a copy of the gauge, and will be changed.
     */
    virtual void GaugeFixing(CFieldGauge* pResGauge) = 0;
    virtual CCString GetInfos(const CCString& sTab) const = 0;
    virtual DOUBLE CheckRes(const CFieldGauge* pGauge) = 0;

#if _CLG_MULTI_GPU
    /**
     * P4-2.1: rank0 global-lattice fixer context helpers.
     *
     * Gauge fixing (esp. the cuFFT variants) needs the WHOLE lattice, so on a
     * decomposed lattice the field is gathered to rank 0 and the existing
     * single-GPU fixing loop runs there under a TEMPORARY GLOBAL lattice
     * context (lattice constants + index cache sized to the global lattice).
     *
     * MGEnterGlobalFixerContext: on rank 0, saves the local (decomposed)
     * lattice constants and index cache, switches the constants to the GLOBAL
     * lattice (from CLGComm::GlobalLattice) and bakes a global index cache, so
     * the subsequent fixing loop and CheckRes run over the whole lattice. On
     * non-root ranks it returns FALSE immediately (they block on the scatter
     * call in the caller). Returns TRUE only on rank 0 (the rank that must run
     * the fixing loop); single-GPU builds (no comm) return TRUE as a no-op.
     *
     * MGExitGlobalFixerContext: rank 0 restores the local constants and index
     * cache and frees the temporary global ones.
     */
    UBOOL MGEnterGlobalFixerContext();
    void MGExitGlobalFixerContext();
    UINT m_uiSavedConstIntegers[128];
    class CIndexData* m_pSavedGlobalIndexCache;
    class CIndex* m_pSavedGlobalIndex;
#endif

    class CLatticeData* m_pOwner;
    Real m_fAccuracy;
    UINT m_iIterate;
    UINT m_iMaxIterate;
    UINT m_iLinearIterate;
    UINT m_iShowErrorStep;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXING_H_

//=============================================================================
// END OF FILE
//=============================================================================