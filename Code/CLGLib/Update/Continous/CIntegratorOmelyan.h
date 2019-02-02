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
        UpdateP(OmelyanLambda * m_fEStep);

        for (UINT uiStep = 1; uiStep < m_uiStepCount; ++uiStep)
        {
            appDetailed("  Omelyan sub step %d\n", uiStep);

            UpdateU(fHalfEstep);
            UpdateP(m_fEStep * (F(1.0) - F(2.0) * OmelyanLambda));
            UpdateU(fHalfEstep);

            if (uiStep < m_uiStepCount)
            {
                UpdateP(m_fEStep * F(2.0) * OmelyanLambda);
            }
            else
            {
                UpdateP(OmelyanLambda * m_fEStep);
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