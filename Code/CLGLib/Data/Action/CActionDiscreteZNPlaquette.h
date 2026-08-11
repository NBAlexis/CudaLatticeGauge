//=============================================================================
// FILENAME : CActionDiscreteZNPlaquette.h
//
// DESCRIPTION:
// Z_N Wilson plaquette action for discrete gauge groups
// Sets the device function pointer to _ZNPlaquetteLinkActionWeight<N>
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONDISCRETEZNPLAQUETTE_H_
#define _CACTIONDISCRETEZNPLAQUETTE_H_

#include "Data/Action/CActionDiscreteGauge.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionDiscreteZNPlaquette)

class CLGAPI CActionDiscreteZNPlaquette : public CActionDiscreteGauge
{
    __CLGDECLARE_CLASS(CActionDiscreteZNPlaquette)

public:
    CActionDiscreteZNPlaquette();

protected:
    DeviceDiscreteLinkActionFunc GetDeviceFunc() const override;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONDISCRETEZNPLAQUETTE_H_

//=============================================================================
// END OF FILE
//=============================================================================
