//=============================================================================
// FILENAME : CUpdator.h
// 
// DESCRIPTION:
// This is the class for update the field
// It can be continous for U(1), SU(2), SU(3) gauge, using hybrid Monte Carlo
// It can also be descrete for Z2, Zn, tetrahydraul, octahydraul, icosahydraul, using Heatbath
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================

#ifndef _CUPDATOR_H_
#define _CUPDATOR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EUpdatorType,
    EUT_HMC,
    EUT_Max,

    EUT_ForceDWORD = 0x7fffffff,
)

class CLGAPI CUpdator : public CBase
{
public:

    CUpdator() : m_pOwner(NULL), m_uiUpdateCall(0), m_bSaveConfigurations(FALSE) { ; }
    virtual UINT Update(UINT iSteps, UBOOL bMeasure) = 0;
    virtual Real CalculateEnergy() = 0;

    virtual EUpdatorType GetUpdatorType() const { return EUT_Max; }
    virtual void Initial(class CLatticeData* pOwner, const CParameters& params) = 0;
    virtual CCString GetInfos(const CCString& sTab) const = 0;

    class CLatticeData* m_pOwner;

    void SetConfigurationCount(UINT uiCount) { m_iAcceptedConfigurationCount = uiCount; }
    UINT GetConfigurationCount() const { return m_iAcceptedConfigurationCount; }
    void OnConfigurationAccepted();
    void SaveConfiguration(UINT uiUpdateStep);

protected:

    UINT m_uiUpdateCall;
    UINT m_iAcceptedConfigurationCount;
    UBOOL m_bSaveConfigurations;
    CCString m_sConfigurationPrefix;
};

__END_NAMESPACE

#endif //#ifndef _CUPDATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================