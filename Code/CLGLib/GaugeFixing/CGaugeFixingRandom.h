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
    Real CheckRes(const CFieldGauge* pGauge) override { return F(0.0); }
    CCString GetInfos(const CCString& sTab) const override;
    void AlsoFixingFermion(CFieldFermionWilsonSquareSU3 * pFermion) const;
    void AlsoFixingAphys(CFieldGauge* pGauge) const;
    deviceSU3* m_pG;
};

__END_NAMESPACE

#endif //#ifndef _CGAUGEFIXINGRANDOM_H_

//=============================================================================
// END OF FILE
//=============================================================================