//=============================================================================
// FILENAME : CHMC.cpp
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CHMC)


CHMC::~CHMC()
{
    appSafeDelete(m_pIntegrator);
}

void CHMC::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
}

void CHMC::Update(UINT iSteps)
{
    UBOOL bAccepted = FALSE;
    Real fEnergy = 0.0f;
    for (UINT i = 0; i < iSteps; ++i)
    {
        m_pIntegrator->Prepare(bAccepted);
        if (0 == i)
        {
            fEnergy = m_pIntegrator->GetEnergy();
        }
        m_pIntegrator->Evaluate();
        Real fEnergyNew = m_pIntegrator->GetEnergy();

        //Metropolis
        appGeneral(_T(" HMC Update: step = %d, energy before update = %f, energy after update = %f"),
            i + 1, fEnergy, fEnergyNew);
        Real diff_H = _exp(fEnergy - fEnergyNew);  // Delta H (SA)
        Real rand = GetRandomReal();
        if (rand <= diff_H)
        {
            appGeneral(_T("  Accepted\n"));
            m_pIntegrator->Accept();
            bAccepted = TRUE;
            fEnergy = fEnergyNew;
        }
        else
        {
            appGeneral(_T("  Rejected\n"));
            bAccepted = FALSE;
        }
    }

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================