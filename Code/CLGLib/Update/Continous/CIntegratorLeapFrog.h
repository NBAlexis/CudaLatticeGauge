//=============================================================================
// FILENAME : CIntegratorLeapFrog.h
// 
// DESCRIPTION:
// This is the Leap frog integrator for HMC
//
// REVISION:
//  [12/11/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATORLEAPFROG_H_
#define _CINTEGRATORLEAPFROG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorLeapFrog)

class CLGAPI CIntegratorLeapFrog : public CIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorLeapFrog)
public:

    /**
    * The gauge field is copied
    * The random momentum field is generated
    */
    void Evaluate() override
    {
        const Real fHalfPstep = F(0.5) * m_fEStep;
        UpdateP(fHalfPstep, FALSE, ESP_StartTrajectory);
        appDetailed("  leap frog sub step 0\n");

        for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
        {
            UpdateU(m_fEStep);

            if (uiStep < m_uiStepCount)
            {
                UpdateP(m_fEStep, FALSE, ESP_InTrajectory);
                appDetailed("  leap frog sub step %d\n", uiStep);
            }
            else 
            {
                UpdateP(fHalfPstep, TRUE, ESP_EndTrajectory);
                appDetailed("  leap frog last step %d\n", uiStep);
            }
        }

        FinishEvaluate();
    }

    CCString GetInfos(const CCString& sTab) const override
    {
        CCString sRet;
        sRet = sTab + _T("Name : LeapFrog\n");
        sRet = sRet + sTab + _T("Epsilon : ") + appToString(m_fEStep) + _T("\n");
        sRet = sRet + sTab + _T("Step : ") + appToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
        sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
        sRet = sRet + sTab + _T("Tau : ") + appToString(m_fEStep * m_uiStepCount) + _T("\n");
        return sRet;
    }
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORLEAPFROG_H_

//=============================================================================
// END OF FILE
//=============================================================================