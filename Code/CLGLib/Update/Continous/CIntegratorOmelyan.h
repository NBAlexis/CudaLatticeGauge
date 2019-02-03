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

#define OmelyanLambda (F(0.19318332750378364))

__BEGIN_NAMESPACE

class CLGAPI CIntegratorOmelyan : public CIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorOmelyan)
public:

    /**
    * The gauge field is copied
    * The random momentum field is generated
    */
    virtual void Evaluate()
    {
        Real fHalfEstep = F(0.5) * m_fEStep;
        appDetailed("  Omelyan sub step 0\n");
        UpdateP(OmelyanLambda * m_fEStep, FALSE);

        for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
        {
            UpdateU(fHalfEstep);
            UpdateP(m_fEStep * (F(1.0) - F(2.0) * OmelyanLambda), FALSE);
            UpdateU(fHalfEstep);

            if (uiStep < m_uiStepCount)
            {
                appDetailed("  Omelyan sub step %d\n", uiStep);
                UpdateP(m_fEStep * F(2.0) * OmelyanLambda, FALSE);
            }
            else
            {
                appDetailed("  Omelyan last step %d\n", uiStep);
                UpdateP(OmelyanLambda * m_fEStep, TRUE);
            }
        }

        m_pGaugeField->ElementNormalize();
    }
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOROMELYAN_H_

//=============================================================================
// END OF FILE
//=============================================================================