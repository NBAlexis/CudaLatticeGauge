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
    , m_pG(NULL)
    {
    }

    ~CGaugeFixingRandom()
    {
        cudaSafeFree(m_pG);
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    void GaugeFixing(CFieldGauge* pResGauge) override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CheckRes(const CFieldGauge* ) override { return 0.0; }
#else
    Real CheckRes(const CFieldGauge* pGauge) override { return F(0.0); }
#endif
    CCString GetInfos(const CCString& sTab) const override;

    /**
     * Call fix fermion just after gaugefixing, will use the same gauge transform as before
     * When never gaugefixing is called the gauge transform is randomized
     */
    void AlsoFixingFermion(CFieldFermion* pFermion) const;
    void AlsoFixingFermionWilsonSU3(CFieldFermionWilsonSquareSU3 * pFermion) const;
    void AlsoFixingFermionKSSU3(CFieldFermionKSSU3* pFermion) const;
    void AlsoFixingAphys(CFieldGauge* pGauge) const;
    deviceSU3* m_pG;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGRANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================