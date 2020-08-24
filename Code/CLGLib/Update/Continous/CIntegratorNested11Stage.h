//=============================================================================
// FILENAME : CIntegratorNested11Stage.h
// 
// DESCRIPTION:
//  11 Stage, from 10.1016/s0010-4655(02)00754-3
//  Velocity version, Variant 8
//
// REVISION:
//  [08/19/2020 nbale]
//=============================================================================

#ifndef _CINTEGRATORNESTED11STAGE_H_
#define _CINTEGRATORNESTED11STAGE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorNested11Stage)

class CLGAPI CIntegratorNested11Stage : public CNestedIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorNested11Stage)

public:

    CIntegratorNested11Stage()
        : CNestedIntegrator()
        , m_fRho(_11Stage_Rho)
        , m_fTheta(_11Stage_Theta)
        , m_fVarTheta(_11Stage_VarTheta)
        , m_fLambda(_11Stage_Lambda)
        , m_bPosition(FALSE)
    {
        
    }

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;

    void Evaluate() override;

    void NestedEvaluate(UBOOL bLast);

    void Nested(Real fStep, UBOOL bLast)
    {
        m_fNestedStepLength = fStep / m_uiNestedStep;
        if (m_bInnerLeapFrog)
        {
            NestedEvaluateLeapfrog(bLast);
        }
        else
        {
            NestedEvaluate(bLast);
        }
    }

    CCString GetInfos(const CCString& sTab) const override;

protected:

    Real m_fRho;
    Real m_fTheta;
    Real m_fVarTheta;
    Real m_fLambda;

    UBOOL m_bPosition;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORNESTED11STAGE_H_

//=============================================================================
// END OF FILE
//=============================================================================