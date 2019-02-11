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

//#define OmelyanLambda (F(0.19318332750378364))

//typically, 0.3-0.5 - arXiv:002.4232
#define OmelyanLambda2 (F(0.38636665500756728))

__BEGIN_NAMESPACE

class CLGAPI CIntegratorOmelyan : public CIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorOmelyan)

public:

    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);

    /**
    * The gauge field is copied
    * The random momentum field is generated
    */
    virtual void Evaluate()
    {
        Real fHalfEstep = F(0.5) * m_fEStep;
        appDetailed("  Omelyan sub step 0\n");
        UpdateP(m_f2Lambda * fHalfEstep, FALSE);

        for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
        {
            UpdateU(fHalfEstep);
            UpdateP(m_fEStep * (F(1.0) - m_f2Lambda), FALSE);
            UpdateU(fHalfEstep);

            if (uiStep < m_uiStepCount)
            {
                appDetailed("  Omelyan sub step %d\n", uiStep);
                UpdateP(m_fEStep * m_f2Lambda, FALSE);
            }
            else
            {
                appDetailed("  Omelyan last step %d\n", uiStep);
                UpdateP(m_f2Lambda * fHalfEstep, TRUE);
            }
        }

        m_pGaugeField->ElementNormalize();
    }

protected:

    Real m_f2Lambda;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOROMELYAN_H_

//=============================================================================
// END OF FILE
//=============================================================================