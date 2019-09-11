//=============================================================================
// FILENAME : CHMC.h
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================

#ifndef _CHMC_H_
#define _CHMC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CHMC)

class CLGAPI CHMC : public CUpdator
{
    __CLGDECLARE_CLASS(CHMC)

public:

    CHMC() : CUpdator(), m_pIntegrator(NULL), m_bMetropolis(FALSE) {  }
    ~CHMC();
    UINT Update(UINT iSteps, UBOOL bMeasure) override;
    Real CalculateEnergy() override { return 0.0f; }
    EUpdatorType GetUpdatorType() const override { return EUT_HMC; }

    class CIntegrator* m_pIntegrator;

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    CCString GetInfos(const CCString &tab) const override;

protected:

    UBOOL m_bMetropolis;
};

__END_NAMESPACE

#endif //#ifndef _CHMC_H_

//=============================================================================
// END OF FILE
//=============================================================================