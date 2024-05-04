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
#if !_CLG_DOUBLEFLOAT
    DOUBLE dValue;
    if (param.FetchValueDOUBLE(_T("Accuracy"), dValue))
    {
        m_fAccuracy = dValue;
    }
#else
    Real fValue;
    if (param.FetchValueReal(_T("Accuracy"), fValue))
    {
        m_fAccuracy = fValue;
    }
#endif

}

void CMultiShiftBiCGStab::AllocateBuffers(const CField*)
{

}

void CMultiShiftBiCGStab::ReleaseBuffers()
{

}

UBOOL CMultiShiftBiCGStab::Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    appPushLogDate(FALSE);
    TArray<cuDoubleComplex> betas;
    TArray<cuDoubleComplex> zeta;
    TArray<cuDoubleComplex> zetaold;
    TArray<cuDoubleComplex> rho;
    TArray<cuDoubleComplex> chis;
    TArray<cuDoubleComplex> alphas;
    TArray<DOUBLE> sl;
    TArray<CField*> pSsigma;
    for (INT i = 0; i < cn.Num(); ++i)
    {
        pFieldX[i]->InitialField(EFIT_Zero);
        betas.AddItem(make_cuDoubleComplex(1.0, 0.0));
        zeta.AddItem(make_cuDoubleComplex(1.0, 0.0));
        zetaold.AddItem(make_cuDoubleComplex(1.0, 0.0));
        rho.AddItem(make_cuDoubleComplex(1.0, 0.0));
        chis.AddItem(make_cuDoubleComplex(0.0, 0.0));
        alphas.AddItem(make_cuDoubleComplex(0.0, 0.0));
        sl.AddItem(1.0);

        CField* s = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        pFieldB->CopyTo(s);
        pSsigma.AddItem(s);
    }

    DOUBLE fBLength = 1.0;
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = sqrt(pFieldB->Dot(pFieldB).x);
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
    cuDoubleComplex delta = pW0->Dot(pR);
    pSA->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
    cuDoubleComplex phi = cuCdiv(pW0->Dot(pSA), delta);
    cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex alpha = make_cuDoubleComplex(0.0, 0.0);

    UBOOL bDone = FALSE;
    for (UINT i = 0; i < m_uiStepCount * m_uiDevationCheck; ++i)
    {
        const cuDoubleComplex newbeta = cuCdiv(make_cuDoubleComplex(-1.0, 0.0), phi);
        for (INT n = 0; n < cn.Num(); ++n)
        {
            if (sl[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            if (0 == i)
            {
                //because alpha=0, alpha_n^sigma is zero
                zeta[n] = cuCdiv(make_cuDoubleComplex(1.0, 0.0), cuCsub(make_cuDoubleComplex(1.0, 0.0), cuCmul(_cToDouble(cn[n]), newbeta)));
                betas[n] = cuCmul(newbeta, zeta[n]);
            }
            else
            {
                const cuDoubleComplex zeta_nm1 = zetaold[n];
                zetaold[n] = zeta[n];

                //beta_n = newbeta, beta_{n-1} = beta
                //zeta_n = zetaold[n], zeta_{n-1} = zeta_nm1, zeta_{n+1}=zeta[n]
                const cuDoubleComplex zetabeta_nm1 = cuCmul(zeta_nm1, beta);
                zeta[n] = cuCdiv(cuCmul(zetaold[n], zetabeta_nm1),
                    cuCadd(cuCmul(alpha, cuCmul(newbeta, cuCsub(zeta_nm1, zetaold[n]))),
                        cuCmul(zetabeta_nm1, cuCsub(make_cuDoubleComplex(1.0, 0.0), cuCmul(_cToDouble(cn[n]), newbeta))))
                );

                //now, betas[n]=beta_n^{sigma}
                betas[n] = cuCdiv(cuCmul(newbeta, zeta[n]), zetaold[n]);
            }
        }
        beta = newbeta;

        pR->CopyTo(pW);
        pW->Axpy(_cToFloat(beta), pSA);
        pW->CopyTo(pWA);
        pWA->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
        const cuDoubleComplex chi = cuCdivf_cd_host(pWA->Dot(pW), pWA->Dot(pWA).x);
        for (INT n = 0; n < cn.Num(); ++n)
        {
            if (sl[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            const cuDoubleComplex dinom = cuCadd(make_cuDoubleComplex(1.0, 0.0), cuCmul(_cToDouble(cn[n]), chi));
            chis[n] = cuCdiv(chi, dinom);
            const cuDoubleComplex chirho = cuCmul(chis[n], rho[n]);
            const cuDoubleComplex chirhozetanp1 = cuCmul(chirho, zeta[n]);
            pFieldX[n]->Axpy(_make_cuComplex(-static_cast<Real>(betas[n].x), -static_cast<Real>(betas[n].y)), pSsigma[n]);
            pFieldX[n]->Axpy(_cToFloat(chirhozetanp1), pW);

            const CLGComplex chirhozetanp1OverBeta = _cToFloat(cuCdiv(chirhozetanp1, betas[n]));
            pSsigma[n]->Axpy(_make_cuComplex(-chirhozetanp1OverBeta.x, -chirhozetanp1OverBeta.y), pW);
            pSsigma[n]->Axpy(_cToFloat(cuCdiv(cuCmul(chirho, zetaold[n]), betas[n])), pR);

            rho[n] = cuCdiv(rho[n], dinom);
        }
        pW->CopyTo(pR);
        pR->Axpy(_make_cuComplex(-static_cast<Real>(chi.x), -static_cast<Real>(chi.y)), pWA);

        const cuDoubleComplex newdelta = pW0->Dot(pR);
        alpha = cuCdiv(cuCmul(beta, make_cuDoubleComplex(-newdelta.x, -newdelta.y)), cuCmul(chi, delta));
        delta = newdelta;

        for (INT n = 0; n < cn.Num(); ++n)
        {
            if (sl[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            alphas[n] = cuCdiv(cuCmul(cuCmul(alpha, zeta[n]), betas[n]),
                cuCmul(zetaold[n], beta));

            pSsigma[n]->ScalarMultply(_cToFloat(alphas[n]));
            pSsigma[n]->Axpy(_cToFloat(cuCmul(zeta[n], rho[n])), pR);
        }

        pS->Axpy(_make_cuComplex(-static_cast<Real>(chi.x), -static_cast<Real>(chi.y)), pSA);
        pS->ScalarMultply(_cToFloat(alpha));
        pS->AxpyPlus(pR);

        pS->CopyTo(pSA);
        pSA->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
        phi = cuCdiv(pW0->Dot(pSA), delta);

        if (0 == (i + 1) % m_uiDevationCheck)
        {
            DOUBLE fMaxErro = _sqrtd(pS->Dot(pS).x);
            for (INT n = 0; n < cn.Num(); ++n)
            {
                if (sl[n] < m_fAccuracy * fBLength)
                {
                    continue;
                }
                sl[n] = sqrt(pSsigma[n]->Dot(pSsigma[n]).x);
                if (sl[n] > fMaxErro)
                {
                    fMaxErro = sl[n];
                }
            }

            appParanoiac(_T("CMultiShiftBiCGStab: diviation is %2.20f\n"), fMaxErro);
            if (fMaxErro < m_fAccuracy * fBLength)
            {
                appParanoiac(_T("CMultiShiftBiCGStab: Done\n"));
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
    appPopLogDate();
    return bDone;
}




__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================