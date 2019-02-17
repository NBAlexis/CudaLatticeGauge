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

class CLGAPI CHMC : public CUpdator
{
    __CLGDECLARE_CLASS(CHMC)

public:

    CHMC() : CUpdator(), m_pIntegrator(NULL), m_bMetropolis(FALSE) { ; }
    ~CHMC();
    virtual UINT Update(UINT iSteps, UBOOL bMeasure);
    virtual Real CalculateEnergy() { return 0.0f; }
    virtual EUpdatorType GetUpdatorType() const { return EUT_HMC; }

    class CIntegrator* m_pIntegrator;

    virtual void Initial(class CLatticeData* pOwner, const CParameters& params);
    virtual CCString GetInfos(const CCString &tab) const;

protected:

    UBOOL m_bMetropolis;
};

__END_NAMESPACE

#endif //#ifndef _CHMC_H_

//=============================================================================
// END OF FILE
//=============================================================================