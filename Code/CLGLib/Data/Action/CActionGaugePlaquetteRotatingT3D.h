//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingT3D.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [mm/dd/yy]
//  [07/14/2024 nbale]
//=============================================================================
#pragma once

#include "CActionGaugePlaquetteRotatingT.h"

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGT3D_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATINGT3D_H_

__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CActionGaugePlaquetteRotatingT3D : public CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>
{
protected:

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple = NULL) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const override;
};

#if !_CLG_WIN
extern template class CActionGaugePlaquetteRotatingT3D<CLGComplex, 1>;
extern template class CActionGaugePlaquetteRotatingT3D<deviceSU3, 3>;
#endif

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotatingU1_3D)
class CLGAPI CActionGaugePlaquetteRotatingU1_3D : public CActionGaugePlaquetteRotatingT3D<CLGComplex, 1>
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotatingU1_3D)
};

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotating3D)

class CLGAPI CActionGaugePlaquetteRotating3D : public CActionGaugePlaquetteRotatingT3D<deviceSU3, 3>
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotating3D)
};

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATINGT3D_H_

//=============================================================================
// END OF FILE
//=============================================================================