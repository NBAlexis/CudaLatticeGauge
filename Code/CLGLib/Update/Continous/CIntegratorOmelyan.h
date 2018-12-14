//=============================================================================
// FILENAME : CIntegratorOmelyan.h
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [12/12/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATOROMELYAN_H_
#define _CINTEGRATOROMELYAN_H_

#define OmelyanLambda (0.19318332750378364f)

__BEGIN_NAMESPACE

class CLGAPI CIntegratorOmelyan : public CIntegrator
{
public:

    /**
    * The gauge field is copied
    * The random momentum field is generated
    */
    virtual void Evaluate()
    {
        FLOAT fHalfEstep = 0.5f * m_fEStep;
        appGeneral("  Omelyan sub step 0\n");
        UpdateP(OmelyanLambda * m_fEStep);

        for (UINT uiStep = 1; uiStep < m_uiStepCount; ++uiStep)
        {
            appGeneral("  Omelyan sub step %d\n", uiStep);

            UpdateU(fHalfEstep);
            UpdateP(m_fEStep * (1.0f - 2.0f * OmelyanLambda));
            UpdateU(fHalfEstep);

            if (uiStep < m_uiStepCount)
            {
                UpdateP(m_fEStep * 2.0f * OmelyanLambda);
            }
            else
            {
                UpdateP(OmelyanLambda * m_fEStep);
            }
        }
    }
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOROMELYAN_H_

//=============================================================================
// END OF FILE
//=============================================================================