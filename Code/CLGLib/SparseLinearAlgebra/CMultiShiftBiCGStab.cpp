//=============================================================================
// FILENAME : CMultiShiftBiCGStab.cpp
// 
// DESCRIPTION:
// This is the class for Multi-Shift BiCGStab Solver
//
// REVISION:
//  [20/06/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMultiShiftBiCGStab)

CMultiShiftBiCGStab::CMultiShiftBiCGStab()
    : CMultiShiftSolver()
    , m_uiDevationCheck(10)
    , m_uiStepCount(20)
    , m_fAccuracy(F(0.000001))
{

}

CMultiShiftBiCGStab::~CMultiShiftBiCGStab()
{
    CMultiShiftBiCGStab::ReleaseBuffers();
}

void CMultiShiftBiCGStab::Configurate(const CParameters& param)
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
    if (param.FetchValueINT(_T("AbsoluteAccuracy"), iValue))
    {
        m_bAbsoluteAccuracy = (0 != iValue);
    }
    if (param.FetchValueReal(_T("Accuracy"), fValue))
    {
        m_fAccuracy = fValue;
    }
}

void CMultiShiftBiCGStab::AllocateBuffers(const CField*)
{

}

void CMultiShiftBiCGStab::ReleaseBuffers()
{

}

UBOOL CMultiShiftBiCGStab::Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    appSetLogDate(FALSE);
    TArray<CLGComplex> betas;
    TArray<CLGComplex> zeta;
    TArray<CLGComplex> zetaold;
    TArray<CLGComplex> rho;
    TArray<CLGComplex> chis;
    TArray<CLGComplex> alphas;
    TArray<Real> sl;
    TArray<CField*> pSsigma;
    for (INT i = 0; i < cn.Num(); ++i)
    {
        pFieldX[i]->InitialField(EFIT_Zero);
        betas.AddItem(_onec);
        zeta.AddItem(_onec);
        zetaold.AddItem(_onec);
        rho.AddItem(_onec);
        chis.AddItem(_zeroc);
        alphas.AddItem(_zeroc);
        sl.AddItem(F(1.0));

        CField* s = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        pFieldB->CopyTo(s);
        pSsigma.AddItem(s);
    }

    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = _sqrt(pFieldB->Dot(pFieldB).x);
    }

    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pW0 = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pS = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pSA = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pWA = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    pFieldB->CopyTo(pR);
    pFieldB->CopyTo(pW);
    pFieldB->CopyTo(pS);
    pFieldB->CopyTo(pSA);
    pW->Dagger();
    pW->CopyTo(pW0);
    CLGComplex delta = pW0->Dot(pR);
    pSA->ApplyOperator(uiM, pGaugeFeild);
    CLGComplex phi = _cuCdivf(pW0->Dot(pSA), delta);
    CLGComplex beta = _zeroc;
    CLGComplex alpha = _zeroc;

    UBOOL bDone = FALSE;
    for (UINT i = 0; i < m_uiStepCount * m_uiDevationCheck; ++i)
    {
        const CLGComplex newbeta = _cuCdivf(_make_cuComplex(-F(1.0), F(0.0)), phi);
        for (INT n = 0; n < cn.Num(); ++n)
        {
            if (sl[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            if (0 == i)
            {
                //because alpha=0, alpha_n^sigma is zero
                zeta[n] = _cuCdivf(_onec, _cuCsubf(_onec, _cuCmulf(cn[n], newbeta)));
                betas[n] = _cuCmulf(newbeta, zeta[n]);
            }
            else
            {
                const CLGComplex zeta_nm1 = zetaold[n];
                zetaold[n] = zeta[n];

                //beta_n = newbeta, beta_{n-1} = beta
                //zeta_n = zetaold[n], zeta_{n-1} = zeta_nm1, zeta_{n+1}=zeta[n]
                const CLGComplex zetabeta_nm1 = _cuCmulf(zeta_nm1, beta);
                zeta[n] = _cuCdivf(_cuCmulf(zetaold[n], zetabeta_nm1),
                    _cuCaddf(_cuCmulf(alpha, _cuCmulf(newbeta,_cuCsubf(zeta_nm1, zetaold[n]))),
                        _cuCmulf(zetabeta_nm1, _cuCsubf(_onec, _cuCmulf(cn[n], newbeta))))
                );

                //now, betas[n]=beta_n^{sigma}
                betas[n] = _cuCdivf(_cuCmulf(newbeta, zeta[n]), zetaold[n]);
            }
        }
        beta = newbeta;

        pR->CopyTo(pW);
        pW->Axpy(beta, pSA);
        pW->CopyTo(pWA);
        pWA->ApplyOperator(uiM, pGaugeFeild);
        const CLGComplex chi = cuCdivf_cr(pWA->Dot(pW), pWA->Dot(pWA).x);

        for (INT n = 0; n < cn.Num(); ++n)
        {
            if (sl[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            const CLGComplex dinom = _cuCaddf(_onec, _cuCmulf(cn[n], chi));
            chis[n] = _cuCdivf(chi, dinom);
            const CLGComplex chirho = _cuCmulf(chis[n], rho[n]);
            const CLGComplex chirhozetanp1 = _cuCmulf(chirho, zeta[n]);
            pFieldX[n]->Axpy(_make_cuComplex(-betas[n].x, -betas[n].y), pSsigma[n]);
            pFieldX[n]->Axpy(chirhozetanp1, pW);

            const CLGComplex chirhozetanp1OverBeta = _cuCdivf(chirhozetanp1, betas[n]);
            pSsigma[n]->Axpy(_make_cuComplex(-chirhozetanp1OverBeta.x, -chirhozetanp1OverBeta.y), pW);
            pSsigma[n]->Axpy(_cuCdivf(_cuCmulf(chirho, zetaold[n]), betas[n]), pR);

            rho[n] = _cuCdivf(rho[n], dinom);
        }
        pW->CopyTo(pR);
        pR->Axpy(_make_cuComplex(-chi.x, -chi.y), pWA);

        const CLGComplex newdelta = pW0->Dot(pR);
        alpha = _cuCdivf(_cuCmulf(beta, _make_cuComplex(-newdelta.x, -newdelta.y)), _cuCmulf(chi, delta));
        delta = newdelta;

        for (INT n = 0; n < cn.Num(); ++n)
        {
            if (sl[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            alphas[n] = _cuCdivf(_cuCmulf(_cuCmulf(alpha, zeta[n]), betas[n]),
                _cuCmulf(zetaold[n], beta));

            pSsigma[n]->ScalarMultply(alphas[n]);
            pSsigma[n]->Axpy(_cuCmulf(zeta[n], rho[n]), pR);
        }

        pS->Axpy(_make_cuComplex(-chi.x, -chi.y), pSA);
        pS->ScalarMultply(alpha);
        pS->AxpyPlus(pR);

        pS->CopyTo(pSA);
        pSA->ApplyOperator(uiM, pGaugeFeild);
        phi = _cuCdivf(pW0->Dot(pSA), delta);

        if (0 == (i + 1) % m_uiDevationCheck)
        {
            Real fMaxErro = _sqrt(pS->Dot(pS).x);
            for (INT n = 0; n < cn.Num(); ++n)
            {
                if (sl[n] < m_fAccuracy * fBLength)
                {
                    continue;
                }
                sl[n] = _sqrt(pSsigma[n]->Dot(pSsigma[n]).x);
                if (sl[n] > fMaxErro)
                {
                    fMaxErro = sl[n];
                }
            }

            appParanoiac(_T("CMultiShiftBiCGStab: diviation is %2.20f\n"), fMaxErro);
            if (fMaxErro < m_fAccuracy * fBLength)
            {
                bDone = TRUE;
                break;
            }
        }
    }

    for (INT n = 0; n < cn.Num(); ++n)
    {
        pSsigma[n]->Return();
    }
    pR->Return();
    pW->Return();
    pW0->Return();
    pS->Return();
    pSA->Return();
    pWA->Return();
    return bDone;
}




__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================