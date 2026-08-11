//=============================================================================
// FILENAME : CActionDiscreteGauge.h
//
// DESCRIPTION:
// Base class for discrete gauge actions (Z_N, D8, etc.)
// Stores a device function pointer for link-level action weight computation
// so the heatbath kernel can call different actions through the same interface.
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONDISCRETEGAUGE_H_
#define _CACTIONDISCRETEGAUGE_H_

#include "Data/Action/CAction.h"

__BEGIN_NAMESPACE

using DeviceDiscreteLinkActionFunc = DOUBLE(*)(const void* __restrict__ pGauge, const void* __restrict__ pStaple, UINT linkIndex, UINT k, DOUBLE beta);

class CLGAPI CActionDiscreteGauge : public CAction
{

public:
    CActionDiscreteGauge();
    virtual ~CActionDiscreteGauge();

    CActionDiscreteGauge(const CActionDiscreteGauge&) = delete;
    CActionDiscreteGauge& operator=(const CActionDiscreteGauge&) = delete;

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString& tab) const override;

    DeviceDiscreteLinkActionFunc* GetDeviceFuncPtr() const { return m_pDeviceFuncPtr; }
    DOUBLE GetBeta() const { return m_fBeta; }
    BYTE GetPrimaryGaugeFieldId() const { return m_byGaugeFieldIds.Num() > 0 ? m_byGaugeFieldIds[0] : 1; }

    UBOOL IsDiscreteGauge() const override { return TRUE; }

protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;

    virtual DeviceDiscreteLinkActionFunc GetDeviceFunc() const = 0;

    DOUBLE m_fBeta;
    DeviceDiscreteLinkActionFunc* m_pDeviceFuncPtr;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONDISCRETEGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================
