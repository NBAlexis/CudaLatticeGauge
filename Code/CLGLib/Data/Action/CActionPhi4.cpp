//=============================================================================
// FILENAME : CActionPhi4.cpp
// 
// DESCRIPTION:
// (nabla phi)*(nabla phi) + m phi* phi + lambda (phi*phi)^2
//
// REVISION:
//  [06/13/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionPhi4)

void CActionPhi4::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) 
{
    CAction::Initial(pOwner, param, byId);

    param.FetchValueReal(_T("Mass"), m_fM);
    param.FetchValueReal(_T("Lambda"), m_fLambda);
}

/**
* calculate Dphi
* calculate |Dphi|^2
* calculate |phi|^2
*/
DOUBLE CActionPhi4::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields)
{
    if (1 != m_byBosonFieldIds.Num())
    {
        appCrucial(_T("there must be boson field!\n"));
        return 0.0;
    }
    INT ibosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, bosonFields, m_byBosonFieldIds[0]);
    const CFieldBoson* bosonfield = bosonFields[ibosonidx];
    CFieldBoson* dboson = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledCopy(bosonfield));
    dboson->D(gaugeNum, bosonNum, gaugeFields, bosonFields);
    checkCudaErrors(cudaDeviceSynchronize());
    cuDoubleComplex partiald = bosonfield->Dot(dboson);
    checkCudaErrors(cudaDeviceSynchronize());
    cuDoubleComplex phi2 = bosonfield->Dot(bosonfield);
    checkCudaErrors(cudaDeviceSynchronize());
    bosonfield->CopyTo(dboson);
    dboson->Mul(bosonfield);
    checkCudaErrors(cudaDeviceSynchronize());
    cuDoubleComplex phi4 = dboson->Dot(dboson);
    checkCudaErrors(cudaDeviceSynchronize());
    m_fLastEnergy = (8.0 + m_fM) * phi2.x + m_fLambda * phi4.x - partiald.x;
    dboson->Return();
    return m_fLastEnergy;
}

/**
* calculate force on boson
* calculate force on gauge
*/
UBOOL CActionPhi4::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    if (1 != m_byBosonFieldIds.Num())
    {
        appCrucial(_T("there must be boson field!\n"));
        return FALSE;
    }

    //Boson force
    INT ibosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, bosonFields, m_byBosonFieldIds[0]);
    const CFieldBoson* bosonfield = bosonFields[ibosonidx];
    //bosonfield->DebugPrintMe();
    CFieldBoson* dboson = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledCopy(bosonfield));
    dboson->D(gaugeNum, bosonNum, gaugeFields, bosonFields);
    checkCudaErrors(cudaDeviceSynchronize());
    
    bosonForces[ibosonidx]->Axpy(_make_cuComplex(F(1.0), F(0.0)), dboson);
    checkCudaErrors(cudaDeviceSynchronize());
    bosonForces[ibosonidx]->Axpy(_make_cuComplex(-F(1.0) * (F(8.0) + m_fM), F(0.0)), bosonfield);
    checkCudaErrors(cudaDeviceSynchronize());
    bosonfield->CopyTo(dboson);
    dboson->Mul(bosonfield);
    checkCudaErrors(cudaDeviceSynchronize());
    dboson->Mul(bosonfield, FALSE);
    checkCudaErrors(cudaDeviceSynchronize());
    bosonForces[ibosonidx]->Axpy(_make_cuComplex(-F(2.0) * m_fLambda, F(0.0)), dboson);
    checkCudaErrors(cudaDeviceSynchronize());
    dboson->Return();
    //Gauge force
    //bosonfield->DebugPrintMe();
    return TRUE;
}

void CActionPhi4::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    if (0 == iUpdateIterate)
    {
        m_fLastEnergy = Energy(TRUE, gaugeNum, bosonNum, gaugeFields, bosonFields, NULL);
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================