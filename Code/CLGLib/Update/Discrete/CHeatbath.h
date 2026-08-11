//=============================================================================
// FILENAME : CHeatbath.h
//
// DESCRIPTION:
// Heat bath updater for discrete gauge groups (Z_N and D_N)
// Supports multiple gauge fields simultaneously -- each field is paired
// with its CActionDiscreteGauge action (matched by field ID via GaugeFields config).
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================

#ifndef _CHEATBATH_H_
#define _CHEATBATH_H_

#include "Data/Action/CActionDiscreteGauge.h"
#include "Update/CUpdator.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CHeatbath)

class CLGAPI CHeatbath : public CUpdator
{
    __CLGDECLARE_CLASS(CHeatbath)

public:
    CHeatbath();
    ~CHeatbath();

    EUpdatorType GetUpdatorType() const override { return EUT_Heatbath; }

    UINT Update(UINT iSteps, UBOOL bMeasure) override;
    Real CalculateEnergy() override
    {
        appCrucial(_T("CHeatbath does not support energy calculation.\n"));
        return F(0.0);
    }

    void Initial(class CLatticeData* pOwner, const CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;
    void SetAutoCorrection(UBOOL ) override
    {
        appGeneral(_T("CHeatbath: SetAutoCorrection has no effect for discrete gauge heatbath\n"));
    }

protected:

    void UpdateOneParity(UBOOL bEven);

    template<typename deviceGauge>
    void UpdateOneField(deviceGauge* buffer, deviceGauge* staple, DeviceDiscreteLinkActionFunc* pSweapFuncion, DOUBLE fBeta, UBOOL bEven, BYTE byFieldId);

};

__END_NAMESPACE

#endif //#ifndef _CHEATBATH_H_

//=============================================================================
// END OF FILE
//=============================================================================
