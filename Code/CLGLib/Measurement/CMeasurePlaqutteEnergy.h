//=============================================================================
// FILENAME : CMeasurePlaqutteEnergy.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASUREPLAQUTTEENERGY_H_
#define _CMEASUREPLAQUTTEENERGY_H_

__BEGIN_NAMESPACE

class CLGAPI CMeasurePlaqutteEnergy : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePlaqutteEnergy)
public:
    CMeasurePlaqutteEnergy() : CMeasure() {}
    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) 
    { 
        m_pOwner = pOwner; 
        m_pLatticeData = pLatticeData;
        m_byId = byId; 
    }

    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

protected:

    TArray<Real> m_lstData;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPLAQUTTEENERGY_H_

//=============================================================================
// END OF FILE
//=============================================================================