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

class CLGAPI CIntegratorLeapFrog : public CIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorLeapFrog)
public:

    /**
    * The gauge field is copied
    * The random momentum field is generated
    */
    virtual void Evaluate()
    {
        Real fHalfPstep = 0.5f * m_fEStep;

        if (m_uiStepCount > 0)
        {
            appGeneral("  leap frog sub step 0\n");
            UpdateP(fHalfPstep);
        }

        for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
        {
            appGeneral("  leap frog sub step %d\n", uiStep);

            UpdateU(m_fEStep);

            if (uiStep < m_uiStepCount)
            {
                UpdateP(m_fEStep);
            }
            else 
            {
                UpdateP(fHalfPstep);
            }
        }
    }
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORLEAPFROG_H_

//=============================================================================
// END OF FILE
//=============================================================================