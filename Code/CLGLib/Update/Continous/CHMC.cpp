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
    m_iAcceptedConfigurationCount = 0;
}

UINT CHMC::Update(UINT iSteps, UBOOL bMeasure)
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
        Real fDiff = fEnergy - fEnergyNew;
        //Metropolis
        appGeneral(_T(" HMC: step = %d, H (before, after, diff) = (%f, %f, %f)\n"),
            i + 1, fEnergy, fEnergyNew, fDiff);
        Real diff_H = _exp(fDiff);  // Delta H (SA)
        Real rand = GetRandomReal();
        if (TRUE)//rand <= diff_H)
        {
            ++m_iAcceptedConfigurationCount;
            appGeneral(_T("  Accepted (accepted:%d)\n"), m_iAcceptedConfigurationCount);
            m_pIntegrator->Accept(); //Here we copy the gauge field back
            bAccepted = TRUE;
            fEnergy = fEnergyNew;
        }
        else
        {
            appGeneral(_T("  Rejected (accepted:%d)\n"), m_iAcceptedConfigurationCount);
            bAccepted = FALSE;
        }

        //If rejected, just accept the old configuration and trigger the measure
        if (bMeasure)
        {
            m_pOwner->OnUpdatorConfigurationAccepted();
        }
    }

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    m_pOwner->OnUpdatorFinished();
    return m_iAcceptedConfigurationCount;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================