//=============================================================================
// FILENAME : CIntegratorNested11Stage.cpp
// 
// DESCRIPTION:
// This is the Approximate force gradient integrator for HMC
//
// REVISION:
//  [08/19/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIntegratorNested11Stage)

void CIntegratorNested11Stage::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CNestedIntegrator::Initial(pOwner, pLattice, params);
    INT iPosition = 0;
    params.FetchValueINT(_T("Position"), iPosition);
    m_bPosition = (0 != iPosition);

    if (!params.FetchValueReal(_T("Rho"), m_fRho))
    {
        m_fRho = m_bPosition ? _11Stage_RhoP : _11Stage_Rho;
    }

    if (!params.FetchValueReal(_T("Theta"), m_fTheta))
    {
        m_fTheta = m_bPosition ? _11Stage_ThetaP : _11Stage_Theta;
    }

    if (!params.FetchValueReal(_T("VarTheta"), m_fVarTheta))
    {
        m_fVarTheta = m_bPosition ? _11Stage_VarThetaP : _11Stage_VarTheta;
    }

    if (!params.FetchValueReal(_T("Lambda"), m_fLambda))
    {
        m_fLambda = m_bPosition ? _11Stage_LambdaP : _11Stage_Lambda;
    }
}

void CIntegratorNested11Stage::Evaluate()
{
    if (m_lstActions.Num() < 2)
    {
        appCrucial(_T("Nested Updator only work with actions more than 2!"));
        _FAIL_EXIT;
    }

    const Real fVarThetaT = m_fVarTheta * m_fEStep;
    const Real f2VarThetaT = fVarThetaT * F(2.0);
    const Real fThetaT = m_fTheta * m_fEStep;
    const Real fRhoT = m_fRho * m_fEStep;
    const Real fLambdaT = m_fLambda * m_fEStep;
    const Real fOneMinus2LambdaTheta = F(0.5) * (F(1.0) - F(2.0) * (m_fLambda + m_fVarTheta)) * m_fEStep;
    const Real fOneMinus2ThetaRho = (F(1.0) - F(2.0) * (m_fTheta + m_fRho)) * m_fEStep;

    appDetailed("  11Stage sub step 0\n");

    UpdatePF(fVarThetaT, ESP_StartTrajectory);

    for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
    {
        Nested(fRhoT, FALSE);
        UpdatePF(fLambdaT, ESP_InTrajectory);
        Nested(fThetaT, FALSE);
        UpdatePF(fOneMinus2LambdaTheta, ESP_InTrajectory);
        Nested(fOneMinus2ThetaRho, FALSE);
        UpdatePF(fOneMinus2LambdaTheta, ESP_InTrajectory);
        Nested(fThetaT, FALSE);
        UpdatePF(fLambdaT, ESP_InTrajectory);
        

        if (uiStep < m_uiStepCount)
        {
            Nested(fRhoT, FALSE);
            appDetailed("  11Stage sub step %d\n", uiStep);
            UpdatePF(f2VarThetaT, ESP_InTrajectory);
        }
        else
        {
            Nested(fRhoT, TRUE);
            appDetailed("  11Stage last step %d\n", uiStep);
            UpdatePF(fVarThetaT, ESP_EndTrajectory);
        }
    }

    FinishEvaluate();
}

void CIntegratorNested11Stage::NestedEvaluate(UBOOL bLast)
{
    const Real fVarThetaT = m_fVarTheta * m_fNestedStepLength;
    const Real f2VarThetaT = fVarThetaT * F(2.0);
    const Real fThetaT = m_fTheta * m_fNestedStepLength;
    const Real fRhoT = m_fRho * m_fNestedStepLength;
    const Real fLambdaT = m_fLambda * m_fNestedStepLength;
    const Real fOneMinus2LambdaTheta = F(0.5) * (F(1.0) - F(2.0) * (m_fLambda + m_fVarTheta)) * m_fNestedStepLength;
    const Real fOneMinus2ThetaRho = (F(1.0) - F(2.0) * (m_fTheta + m_fRho)) * m_fNestedStepLength;

    appDetailed("  11Stage nested sub step 0\n");
    UpdatePG(fVarThetaT, FALSE);

    for (UINT uiStep = 1; uiStep < m_uiNestedStep + 1; ++uiStep)
    {
        UpdateU(fRhoT);
        UpdatePG(fLambdaT, FALSE);
        UpdateU(fThetaT);
        UpdatePG(fOneMinus2LambdaTheta, FALSE);
        UpdateU(fOneMinus2ThetaRho);
        UpdatePG(fOneMinus2LambdaTheta, FALSE);
        UpdateU(fThetaT);
        UpdatePG(fLambdaT, FALSE);
        UpdateU(fRhoT);

        if (uiStep < m_uiNestedStep)
        {
            appDetailed("  11Stage nested sub step %d\n", uiStep);
            UpdatePG(f2VarThetaT, FALSE);
        }
        else
        {
            appDetailed("  11Stage nested last step %d\n", uiStep);
            UpdatePG(fVarThetaT, bLast);
        }
    }
}

CCString CIntegratorNested11Stage::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    sRet = sTab + _T("Name : Nested 11 Stage\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appFloatToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appIntToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appFloatToString(m_fEStep * m_uiStepCount) + _T("\n");
    sRet = sRet + GetNestedInfo(sTab);
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================