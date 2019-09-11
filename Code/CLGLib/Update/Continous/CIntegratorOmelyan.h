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

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorOmelyan)

class CLGAPI CIntegratorOmelyan : public CIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorOmelyan)

public:

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;

    /**
    * The gauge field is copied
    * The random momentum field is generated
    */
    void Evaluate() override
    {
        const Real fHalfEstep = F(0.5) * m_fEStep;
        appDetailed("  Omelyan sub step 0\n");
        UpdateP(m_f2Lambda * fHalfEstep, FALSE, ESP_StartTrajectory);

        for (UINT uiStep = 1; uiStep < m_uiStepCount + 1; ++uiStep)
        {
            UpdateU(fHalfEstep);
            UpdateP(m_fEStep * (F(1.0) - m_f2Lambda), FALSE, ESP_InTrajectory);
            UpdateU(fHalfEstep);

            if (uiStep < m_uiStepCount)
            {
                appDetailed("  Omelyan sub step %d\n", uiStep);
                UpdateP(m_fEStep * m_f2Lambda, FALSE, ESP_InTrajectory);
            }
            else
            {
                appDetailed("  Omelyan last step %d\n", uiStep);
                UpdateP(m_f2Lambda * fHalfEstep, TRUE, ESP_EndTrajectory);
            }
        }

        FinishEvaluate();
    }

    CCString GetInfos(const CCString& sTab) const override;

protected:

    Real m_f2Lambda;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOROMELYAN_H_

//=============================================================================
// END OF FILE
//=============================================================================