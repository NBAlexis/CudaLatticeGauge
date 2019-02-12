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

    INT iMetro = 0;
    params.FetchValueINT(_T("Metropolis"), iMetro);
    m_bMetropolis = (0 != iMetro);
}

UINT CHMC::Update(UINT iSteps, UBOOL bMeasure)
{
    UBOOL bAccepted = FALSE;

    Real fEnergy = F(0.0);
    Real fEnergyNew = F(0.0);

    for (UINT i = 0; i < iSteps; ++i)
    {
        m_pIntegrator->Prepare(bAccepted, i);
        if (m_bMetropolis)
        {
            fEnergy = m_pIntegrator->GetEnergy(TRUE);
        }
        m_pIntegrator->Evaluate();
        if (m_bMetropolis)
        {
            fEnergyNew = m_pIntegrator->GetEnergy(FALSE);
        }

        Real diff_H = F(1.0);
        Real rand = F(0.0);

        if (m_bMetropolis)
        {
            Real fDiff = fEnergy - fEnergyNew;
            //Metropolis
            appDetailed(_T(" HMC: step = %d, H_dff = %f\n"),
                i + 1, 
                fDiff);
            diff_H = _hostexp(fDiff);  // Delta H (SA)
            rand = GetRandomReal();
        }

        if (rand <= diff_H)
        {
            ++m_iAcceptedConfigurationCount;
            appGeneral(_T("  Accepted (accepted:%d)\n"), m_iAcceptedConfigurationCount);
            bAccepted = TRUE;
        }
        else
        {
            appGeneral(_T("  Rejected (accepted:%d)\n"), m_iAcceptedConfigurationCount);
            bAccepted = FALSE;
        }
        m_pIntegrator->OnFinishTrajectory(bAccepted); //Here we copy the gauge field back

        //If rejected, just accept the old configuration and trigger the measure
        if (bMeasure)
        {
            m_pOwner->OnUpdatorConfigurationAccepted();
        }
    }

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    m_pOwner->OnUpdatorFinished(bMeasure);
    return m_iAcceptedConfigurationCount;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================