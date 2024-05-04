//=============================================================================
// FILENAME : CSolverTFQMR.cpp
// 
// DESCRIPTION:
// This is the class for TFQMR solver
//
// REVISION:
//  [20/06/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSolverTFQMR)

CSolverTFQMR::CSolverTFQMR()
    : CSLASolver()
    , m_uiReTry(1)
    , m_uiDevationCheck(10)
    , m_uiStepCount(20)
    , m_fAccuracy(F(0.000001))
{

}

CSolverTFQMR::~CSolverTFQMR()
{
    CSolverTFQMR::ReleaseBuffers();
}

void CSolverTFQMR::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;
    if (param.FetchValueINT(_T("DiviationStep"), iValue))
    {
        m_uiDevationCheck = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("MaxStep"), iValue))
    {
        m_uiStepCount = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("Restart"), iValue))
    {
        m_uiReTry = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("AbsoluteAccuracy"), iValue))
    {
        m_bAbsoluteAccuracy = (0 != iValue);
    }
    if (param.FetchValueReal(_T("Accuracy"), fValue))
    {
        m_fAccuracy = fValue;
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small (%2.18f), set to be %2.18f\n"), fValue, m_fAccuracy);
        }
    }
}

void CSolverTFQMR::AllocateBuffers(const CField* )
{

}

void CSolverTFQMR::ReleaseBuffers()
{

}

UBOOL CSolverTFQMR::Solve(CField* pFieldX, const CField* pFieldB, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pD = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pV = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pAU = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pU = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pRh = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
#if !_CLG_DOUBLEFLOAT
        fBLength = _cuCabsf(_cToFloat(pFieldB->Dot(pFieldB)));
#else
        fBLength = _cuCabsf(pFieldB->Dot(pFieldB));
#endif
    }

    appParanoiac(_T("-- CSolverTFQMR::Solve start operator: %s--\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str());

    //pFieldB->CopyTo(pB);

    //Using b as the guess, (Assuming A is near identity?)
    //r_0 = b - A x_0
    if (NULL == pStart)
    {
        pFieldB->CopyTo(pV);
        pFieldB->CopyTo(pX);
    }
    else 
    {
        pStart->CopyTo(pV);
        pStart->CopyTo(pX);
    }

    UBOOL bDone = FALSE;
    for (UINT i = 0; i < m_uiReTry; ++i)
    {
        pV->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields, EOCT_Minus); //-A x_0
        pV->AxpyPlus(pFieldB); //pR->AxpyPlus(pB); //b - A x_0
        pV->CopyTo(pRh);
        pRh->Dagger();

        pV->CopyTo(pU);
        pV->CopyTo(pW);
        pV->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields); //v0 = A u0
        pV->CopyTo(pAU);

        Real thetaSq = F(0.0);
#if !_CLG_DOUBLEFLOAT
        Real tau = static_cast<Real>(_sqrtd(pU->Dot(pU).x));
#else
        Real tau = _sqrt(pU->Dot(pU).x);
#endif

#if !_CLG_DOUBLEFLOAT
        CLGComplex rho = _cToFloat(pRh->Dot(pU));
#else
        CLGComplex rho = pRh->Dot(pU);
#endif
        CLGComplex alpha = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex eta = _make_cuComplex(F(0.0), F(0.0));

        for (UINT j = 0; j < m_uiStepCount * m_uiDevationCheck; ++j)
        {
            //==========
            //One step
            const Real fErrorRho = __cuCabsSqf(rho);
            if (fErrorRho < _CLG_FLT_MIN) //if rho = 0, failed (assume will not)    
            {
                appParanoiac(_T("CSolverTFQMR::rho too small:%0.28f\n"), fErrorRho);
                if (i < m_uiReTry - 1)
                {
                    i = m_uiReTry - 1;
                }
                break;
            }

            if (0 == (j & 1)) //even
            {
#if !_CLG_DOUBLEFLOAT
                alpha = _cuCdivf(rho, _cToFloat(pRh->Dot(pV)));
#else
                alpha = _cuCdivf(rho, pRh->Dot(pV));
#endif
            }

            pW->Axpy(_make_cuComplex(-alpha.x, -alpha.y), pAU);
            if (0 == j)
            {
                pU->CopyTo(pD);
            }
            else
            {
                pD->ScalarMultply(_cuCmulf(_cuCdivf(_make_cuComplex(thetaSq, F(0.0)), alpha), eta));
                pD->AxpyPlus(pU);
            }
#if !_CLG_DOUBLEFLOAT
            thetaSq = static_cast<Real>(pW->Dot(pW).x / (tau * tau));
#else
            thetaSq = pW->Dot(pW).x / (tau * tau);
#endif
            const Real cSq = F(1.0) / (1 + thetaSq);
            tau = tau * _sqrt(thetaSq * cSq);
            eta = cuCmulf_cr(alpha, cSq);

            pX->Axpy(eta, pD);

            if (1 == (j & 1)) //odd
            {
#if !_CLG_DOUBLEFLOAT
                const CLGComplex newRho = _cToFloat(pRh->Dot(pW));
#else
                const CLGComplex newRho = pRh->Dot(pW);
#endif
                CLGComplex beta = _cuCdivf(newRho, rho);
                rho = newRho;
                pU->ScalarMultply(beta);
                pU->AxpyPlus(pW);

                pV->ScalarMultply(beta);
                pV->AxpyPlus(pAU);
                pV->ScalarMultply(beta); //v = beta (Au(m) + beta v(m))

                pU->CopyTo(pAU);
                pAU->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
                pV->AxpyPlus(pAU); //v = Au(m+1) + beta (Au(m) + beta v(m))

            }
            else
            {
                pU->Axpy(_make_cuComplex(-alpha.x, -alpha.y), pV);
                pU->CopyTo(pAU);
                pAU->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
            }

            if (0 == (j + 1) % m_uiDevationCheck)
            {
#if !_CLG_DOUBLEFLOAT
                const Real fErrorV = _cuCabsf(_cToFloat(pRh->Dot(pV)));
#else
                const Real fErrorV = _cuCabsf(pRh->Dot(pV));
#endif
                appParanoiac(_T("CSolverTFQMR::r* x v = %8.18f\n"), fErrorV);
                if (fErrorV < m_fAccuracy * fBLength)
                {
                    if (i < m_uiReTry - 1)
                    {
                        i = m_uiReTry - 1;
                    }
                    bDone = TRUE;
                    break;
                }
            }
        }

        //we are here, means we do not converge.
        //we need to restart with a new guess, we use last X
        pX->CopyTo(pV);
    }

    //The solver failed.
    appGeneral(_T("CSolverTFQMR Stopped (If it is rho too small, maybe already done)"));

    pX->CopyTo(pFieldX);

    pX->Return();
    pD->Return();
    pV->Return();
    pAU->Return();
    pU->Return();
    pW->Return();
    pRh->Return();
    return bDone;
}



__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================