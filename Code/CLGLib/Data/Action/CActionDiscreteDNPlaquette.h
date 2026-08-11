//=============================================================================
// FILENAME : CActionDiscreteDNPlaquette.h
//
// DESCRIPTION:
// D_N Wilson plaquette action for dihedral gauge groups (D3, D4, D8).
// Extends CActionDiscreteGauge for infrastructure.
//
// DeviceDiscreteLinkActionFunc passes staples as const void* (typed in the
// kernel via static_cast), so D_N receives full 2x2 deviceDN<N> matrices.
// The same function-pointer interface works for both Z_N and D_N.
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONDISCRETEDNPLAQUETTE_H_
#define _CACTIONDISCRETEDNPLAQUETTE_H_

#include "Data/Action/CActionDiscreteGauge.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionDiscreteDNPlaquette)

class CLGAPI CActionDiscreteDNPlaquette : public CActionDiscreteGauge
{
    __CLGDECLARE_CLASS(CActionDiscreteDNPlaquette)

public:
    CActionDiscreteDNPlaquette();

protected:
    DeviceDiscreteLinkActionFunc GetDeviceFunc() const override;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONDISCRETEDNPLAQUETTE_H_

//=============================================================================
// END OF FILE
//=============================================================================
