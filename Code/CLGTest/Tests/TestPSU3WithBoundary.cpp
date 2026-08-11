//=============================================================================
// FILENAME : TestPSU3WithBoundary.cpp
//
// DESCRIPTION:
// Tests for CActionGaugePlaquettePSU3WithBoundary and CFieldTensor2Z3.
//
// Coverage (see PSU3WithBoundaryImplementationPlan.md section 11):
//   1. field : factory creation, EFIT_Zero/Identity == root 0, EFIT_Random
//              only produces Z3 roots, copy consistency, file round-trip
//   2. energy: B=1 agrees with the fundamental Wilson dynamic energy
//   3. force : finite-difference validation with a fixed B
//   4. heatbath (AllowMonopole=1): conditional distribution on a cold gauge
//              background, P(k) proportional to exp((Beta/3) Re[zeta^k Tr(U_p)])
//   5. flat sweep (AllowMonopole=0): link-star + global sheet sweeps keep
//              dB=0 (checked by CheckMonopole=1) and move B away from identity
//   6. callback/cache: post-sweep energy cache matches a full recomputation
//
// REVISION:
//  [08/09/26]
//=============================================================================
#include "CLGTest.h"
#include <cmath>

//=============================================================================
// helpers
//=============================================================================

static inline double GetClgRandom11()
{
    return (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * 2.0 - 1.0;
}

/**
* index of the nearest Z3 root of a host complex value
*/
static UINT HostRootIndex(double dRe, double dIm)
{
    UINT kBest = 0;
    double dBest = -1.0;
    for (UINT k = 0; k < 3; ++k)
    {
        const double dTheta = 2.0 * 3.14159265358979323846 * static_cast<double>(k) / 3.0;
        const double dDr = dRe - cos(dTheta);
        const double dDi = dIm - sin(dTheta);
        const double dDist2 = dDr * dDr + dDi * dDi;
        if (dBest < 0.0 || dDist2 < dBest)
        {
            dBest = dDist2;
            kBest = k;
        }
    }
    return kBest;
}

/**
* read all Z3 root indices of a tensor2 field into pOut (size 6*volume)
* returns the element count read
*/
static UINT ReadBFieldRoots(const CFieldTensor2Z3* pB, UINT* pOut)
{
    UINT uiSize = 0;
    BYTE* pData = pB->CopyDataOut(uiSize);
    const Real* pReal = (const Real*)pData;
    const UINT uiCount = uiSize / (2 * sizeof(Real));
    for (UINT e = 0; e < uiCount; ++e)
    {
        pOut[e] = HostRootIndex(static_cast<double>(pReal[2 * e + 0]), static_cast<double>(pReal[2 * e + 1]));
    }
    free(pData);
    return uiCount;
}

//=============================================================================
// 1. field tests
//=============================================================================

UINT TestZ3FieldCreate(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    if (NULL == pB)
    {
        appCrucial(_T("TestZ3FieldCreate: no CFieldTensor2Z3 field at id 3\n"));
        return 1;
    }

    const UINT uiElementCount = pB->GetElementCount();
    if (uiElementCount != 6 * _HC_Volume)
    {
        ++uiErrors;
        appCrucial(_T("TestZ3FieldCreate: unexpected element count %u\n"), uiElementCount);
        return uiErrors;
    }

    TArray<UINT> roots;

    // EFIT_Identity -> all root 0 (never complex 0)
    pB->InitialField(EFIT_Identity);
    roots.SetSize(uiElementCount);
    ReadBFieldRoots(pB, roots.GetData());
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (0 != roots[e])
        {
            ++uiErrors;
            appGeneral(_T("TestZ3FieldCreate: EFIT_Identity element %u is root %u\n"), e, roots[e]);
            break;
        }
    }

    // EFIT_Zero -> all root 0 as well
    pB->InitialField(EFIT_Zero);
    ReadBFieldRoots(pB, roots.GetData());
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (0 != roots[e])
        {
            ++uiErrors;
            appGeneral(_T("TestZ3FieldCreate: EFIT_Zero element %u is root %u\n"), e, roots[e]);
            break;
        }
    }

    // EFIT_Random -> all elements are one of the three roots
    pB->InitialField(EFIT_Random);
    ReadBFieldRoots(pB, roots.GetData());
    UINT uiNonZero = 0;
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (roots[e] > 2)
        {
            ++uiErrors;
            appGeneral(_T("TestZ3FieldCreate: EFIT_Random element %u is root %u\n"), e, roots[e]);
            break;
        }
        if (0 != roots[e])
        {
            ++uiNonZero;
        }
    }
    if (0 == uiNonZero)
    {
        ++uiErrors;
        appCrucial(_T("TestZ3FieldCreate: EFIT_Random produced only identity, RNG broken?\n"));
    }

    // GetCopy / CopyTo elementwise consistency
    CFieldTensor2Z3* pCopy = dynamic_cast<CFieldTensor2Z3*>(pB->GetCopy());
    if (NULL == pCopy)
    {
        ++uiErrors;
        appCrucial(_T("TestZ3FieldCreate: GetCopy failed\n"));
        return uiErrors;
    }
    TArray<UINT> rootsCopy;
    rootsCopy.SetSize(uiElementCount);
    ReadBFieldRoots(pCopy, rootsCopy.GetData());
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (rootsCopy[e] != roots[e])
        {
            ++uiErrors;
            appGeneral(_T("TestZ3FieldCreate: copy mismatch at element %u (%u vs %u)\n"), e, rootsCopy[e], roots[e]);
            break;
        }
    }
    appSafeDelete(pCopy);

    // file round-trip preserves root indices
    pB->SaveToFile(_T("testZ3field.con"), EFFT_CLGBin);
    pCopy = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetPooledFieldById(3, _T(__FILE__), __LINE__));
    if (NULL == pCopy)
    {
        ++uiErrors;
        appCrucial(_T("TestZ3FieldCreate: pooled field allocation failed\n"));
        return uiErrors;
    }
    pCopy->InitialFieldWithFile(_T("testZ3field.con"), EFFT_CLGBin);
    rootsCopy.SetSize(uiElementCount);
    ReadBFieldRoots(pCopy, rootsCopy.GetData());
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (rootsCopy[e] != roots[e])
        {
            ++uiErrors;
            appGeneral(_T("TestZ3FieldCreate: round-trip mismatch at element %u (%u vs %u)\n"), e, rootsCopy[e], roots[e]);
            break;
        }
    }
    pCopy->Return();
    appSynchronize();

    appGeneral(_T("TestZ3FieldCreate: %u errors\n"), uiErrors);
    return uiErrors;
}

//=============================================================================
// 2. energy test: B=1 equals the fundamental Wilson dynamic energy
//=============================================================================

UINT TestPSU3Energy(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CActionGaugePlaquettePSU3WithBoundary* pAction =
        dynamic_cast<CActionGaugePlaquettePSU3WithBoundary*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pB || NULL == pAction)
    {
        appCrucial(_T("TestPSU3Energy: missing gauge (1) / boundary (3) / action (1)\n"));
        return 1;
    }

    DOUBLE fBeta = 6.0;
    if (params.Exist(_T("Action1")))
    {
        CParameters actionParam = params.GetParameter(_T("Action1"));
        actionParam.FetchValueDOUBLE(_T("Beta"), fBeta);
    }
    const DOUBLE fBetaOverN = fBeta / 3.0;

    // B = 1
    pB->InitialField(EFIT_Identity);

    const CFieldGauge* g[1] = { pGauge };
    const CFieldTensor2* t[1] = { pB };

    // fundamental Wilson dynamic energy: E_fund = sum_P (N - ReTr(U_P)) * Beta/N
    // our action with B=1: E = -(Beta/3) sum_P ReTr(U_P) = E_fund - 6 * vol * Beta
    const DOUBLE eFund = pGauge->CalculatePlaqutteEnergy(fBetaOverN);
    const DOUBLE eAction = pAction->Energy(FALSE, 1, 0, 1, g, NULL, t, NULL);
    const DOUBLE eExpect = eFund - 6.0 * fBeta * static_cast<DOUBLE>(_HC_Volume);

    const DOUBLE dDiff = fabs(eAction - eExpect);
    const DOUBLE dTol = 1.0e-3 * fmax(1.0, fabs(eExpect));
    appGeneral(_T("TestPSU3Energy: B=1 action=%f, fundamental-derived=%f, diff=%e (tol=%e)\n"),
        eAction, eExpect, dDiff, dTol);

    if (!std::isfinite(eAction) || dDiff > dTol)
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3Energy: B=1 energy does not match fundamental Wilson\n"));
    }

    // random valid B (AllowMonopole=1 mode) must also give a finite energy
    pB->InitialField(EFIT_Random);
    const DOUBLE eRandom = pAction->Energy(FALSE, 1, 0, 1, g, NULL, t, NULL);
    if (!std::isfinite(eRandom))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3Energy: non-finite energy for random B\n"));
    }

    appGeneral(_T("TestPSU3Energy: %u errors\n"), uiErrors);
    return uiErrors;
}

//=============================================================================
// 3. force finite difference with fixed B
//=============================================================================

UINT TestPSU3Force(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CActionGaugePlaquettePSU3WithBoundary* pAction =
        dynamic_cast<CActionGaugePlaquettePSU3WithBoundary*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pB || NULL == pAction)
    {
        appCrucial(_T("TestPSU3Force: missing gauge (1) / boundary (3) / action (1)\n"));
        return 1;
    }

    INT iLinkTrials = 4;
    INT iDirTrials = 2;
    Real fDelta = F(1.0e-1);
    INT iDeltaLevels = 6;
    Real fTolRel = F(0.30);
    Real fTolAbs = F(5.0e-3);
    params.FetchValueINT(_T("LinkTrials"), iLinkTrials);
    params.FetchValueINT(_T("DirTrials"), iDirTrials);
    params.FetchValueReal(_T("Delta"), fDelta);
    params.FetchValueINT(_T("DeltaLevels"), iDeltaLevels);
    params.FetchValueReal(_T("TolAbs"), fTolAbs);
    params.FetchValueReal(_T("TolRel"), fTolRel);

    // fixed B: identity for the B=1 comparison; a random B afterwards
    pB->InitialField(EFIT_Identity);

    // save U0
    CFieldGaugeSU3* pU0 = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
    CFieldGaugeSU3* pWork = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
    CFieldGaugeSU3* pRawForce = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
    CFieldGaugeSU3* pProjForce = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
    if (NULL == pU0 || NULL == pWork || NULL == pRawForce || NULL == pProjForce)
    {
        appCrucial(_T("TestPSU3Force: allocation failed\n"));
        appSafeDelete(pProjForce);
        appSafeDelete(pRawForce);
        appSafeDelete(pWork);
        appSafeDelete(pU0);
        return 1;
    }
    pGauge->CopyTo(pU0);
    appSynchronize();

    const UINT uiLinkCount = _HC_LinkCount;
    const UINT uiVolume = _HC_Volume;

    UINT uiPassed = 0;
    UINT uiFailed = 0;

    for (INT iTrial = 0; iTrial < iLinkTrials * iDirTrials; ++iTrial)
    {
        // random link, stratified by direction
        const BYTE byDir = static_cast<BYTE>(iTrial % 4);
        const UINT uiSiteRand = static_cast<UINT>(fabs(GetClgRandom11()) * static_cast<double>(uiVolume)) % uiVolume;
        const UINT uiLinkIndex = uiSiteRand * 4 + byDir;
        if (uiLinkIndex >= uiLinkCount)
        {
            continue;
        }

        // random anti-hermitian traceless direction X (only this link)
        CLGComplex xMatrix[9];
        // reuse the same construction as TestHmcDiagnostics
        {
            const double a1 = GetClgRandom11();
            const double a2 = GetClgRandom11();
            const double a3 = -(a1 + a2);
            const double x01_re = GetClgRandom11();
            const double x01_im = GetClgRandom11();
            const double x02_re = GetClgRandom11();
            const double x02_im = GetClgRandom11();
            const double x12_re = GetClgRandom11();
            const double x12_im = GetClgRandom11();
            xMatrix[0] = _make_cuComplex(F(0.0), static_cast<Real>(a1));
            xMatrix[1] = _make_cuComplex(static_cast<Real>(x01_re), static_cast<Real>(x01_im));
            xMatrix[2] = _make_cuComplex(static_cast<Real>(x02_re), static_cast<Real>(x02_im));
            xMatrix[3] = _make_cuComplex(static_cast<Real>(-x01_re), static_cast<Real>(x01_im));
            xMatrix[4] = _make_cuComplex(F(0.0), static_cast<Real>(a2));
            xMatrix[5] = _make_cuComplex(static_cast<Real>(x12_re), static_cast<Real>(x12_im));
            xMatrix[6] = _make_cuComplex(static_cast<Real>(-x02_re), static_cast<Real>(x02_im));
            xMatrix[7] = _make_cuComplex(static_cast<Real>(-x12_re), static_cast<Real>(x12_im));
            xMatrix[8] = _make_cuComplex(F(0.0), static_cast<Real>(a3));
            double fNormSq = 2.0 * (x01_re * x01_re + x01_im * x01_im
                + x02_re * x02_re + x02_im * x02_im
                + x12_re * x12_re + x12_im * x12_im)
                + (a1 * a1 + a2 * a2 + a3 * a3);
            const double fScale = 1.0 / sqrt(fNormSq);
            for (UINT i = 0; i < 9; ++i)
            {
                xMatrix[i].x = static_cast<Real>(xMatrix[i].x * fScale);
                xMatrix[i].y = static_cast<Real>(xMatrix[i].y * fScale);
            }
        }

        CFieldGaugeSU3* pX = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
        if (NULL == pX)
        {
            ++uiFailed;
            continue;
        }
        pX->Zero();
        appSetGaugeLink(pX, uiLinkIndex, xMatrix);
        appSynchronize();

        // force at U0 with fixed B
        pRawForce->Zero();
        const CFieldGauge* g0[1] = { pU0 };
        CFieldGauge* f0[1] = { pRawForce };
        if (!pAction->CalculateForce(1, 0, g0, NULL, f0, NULL, NULL, ESP_Once))
        {
            appCrucial(_T("TestPSU3Force: CalculateForce failed\n"));
            ++uiFailed;
            appSafeDelete(pX);
            continue;
        }
        appSynchronize();

        pRawForce->CopyTo(pProjForce);
        pProjForce->LeftMul(pU0, FALSE, TRUE);
        pProjForce->TA();
        appSynchronize();

        const DOUBLE fDforce = -2.0 * static_cast<double>(pProjForce->Dot(pX).x);

        // finite difference: record the error at every step size and verify
        // h-convergence (central-difference truncation error ~ C*h^2, so the
        // error must decrease as h is halved in the double-precision regime)
        TArray<DOUBLE> fdAbsErrors;
        TArray<DOUBLE> fdDfd;
        UBOOL bAnyNonFinite = FALSE;
        for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
        {
            const Real fH = fDelta / static_cast<Real>(1U << iLevel);

            const CFieldGauge* wg[1] = { pWork };
            const CFieldTensor2* t[1] = { pB };

            pU0->CopyTo(pWork);
            appSynchronize();
            pX->ExpMult(fH, pWork);
            appSynchronize();
            const DOUBLE fSPlus = pAction->Energy(FALSE, 1, 0, 1, wg, NULL, t, NULL);

            pU0->CopyTo(pWork);
            appSynchronize();
            pX->ExpMult(-fH, pWork);
            appSynchronize();
            const DOUBLE fSMinus = pAction->Energy(FALSE, 1, 0, 1, wg, NULL, t, NULL);

            const DOUBLE fDfd = (fSPlus - fSMinus) / (2.0 * static_cast<double>(fH));
            const DOUBLE fAbsErr = fabs(fDfd - fDforce);
            fdAbsErrors.AddItem(fAbsErr);
            fdDfd.AddItem(fDfd);
            if (!std::isfinite(fDfd) || !std::isfinite(fAbsErr))
            {
                bAnyNonFinite = TRUE;
            }
        }

        // h-convergence: the central-difference truncation error ~ C*h^2 must
        // decrease monotonically in the truncation regime (the first three
        // levels h, h/2, h/4). At smaller h the single-precision energy noise
        // floor (~1e-5 here) dominates, so monotonicity is only required in
        // the truncation regime, not down to the noise floor.
        UBOOL bConverged = FALSE;
        if (!bAnyNonFinite && iDeltaLevels >= 3)
        {
            bConverged = fdAbsErrors[0] > fdAbsErrors[1] && fdAbsErrors[1] > fdAbsErrors[2];
        }

        // any level within tolerance
        UBOOL bWithinTol = FALSE;
        for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
        {
            const DOUBLE fScale = fmax(fabs(fdDfd[iLevel]), fabs(fDforce));
            const DOUBLE fBound = static_cast<double>(fTolAbs) + static_cast<double>(fTolRel) * fScale;
            if (fdAbsErrors[iLevel] <= fBound)
            {
                bWithinTol = TRUE;
                break;
            }
        }

        CCString sErrSeq;
        for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
        {
            CCString sItem;
            sItem.Format(_T("h=%.1e err=%.3e%s"),
                static_cast<double>(fDelta / static_cast<Real>(1U << iLevel)),
                fdAbsErrors[iLevel],
                (iLevel + 1 < iDeltaLevels) ? _T(", ") : _T(""));
            sErrSeq = sErrSeq + sItem;
        }

        const UBOOL bTrialPassed = bConverged && bWithinTol;

        if (bTrialPassed)
        {
            ++uiPassed;
            appGeneral(_T("TestPSU3Force: trial %d (link %u dir %u) PASS, %s\n"),
                iTrial, uiLinkIndex, byDir, sErrSeq.c_str());
        }
        else
        {
            ++uiFailed;
            appGeneral(_T("TestPSU3Force: trial %d (link %u dir %u) FAIL%s%s, %s\n"),
                iTrial, uiLinkIndex, byDir,
                bConverged ? _T("") : _T(" (not h-converged)"),
                bWithinTol ? _T("") : _T(" (out of tolerance)"),
                sErrSeq.c_str());
        }

        appSafeDelete(pX);
    }

    // also validate a random B (contains non-trivial plaquette phases)
    pB->InitialField(EFIT_Random);
    UINT uiFailedRandom = 0;
    for (INT iTrial = 0; iTrial < iDirTrials; ++iTrial)
    {
        const BYTE byDir = static_cast<BYTE>(iTrial % 4);
        const UINT uiLinkIndex = (static_cast<UINT>(fabs(GetClgRandom11()) * static_cast<double>(uiVolume)) % uiVolume) * 4 + byDir;

        CLGComplex xMatrix[9];
        const double a1 = GetClgRandom11();
        const double a2 = GetClgRandom11();
        const double a3 = -(a1 + a2);
        xMatrix[0] = _make_cuComplex(F(0.0), static_cast<Real>(a1));
        xMatrix[4] = _make_cuComplex(F(0.0), static_cast<Real>(a2));
        xMatrix[8] = _make_cuComplex(F(0.0), static_cast<Real>(a3));
        xMatrix[1] = _make_cuComplex(static_cast<Real>(GetClgRandom11()), static_cast<Real>(GetClgRandom11()));
        xMatrix[2] = _make_cuComplex(static_cast<Real>(GetClgRandom11()), static_cast<Real>(GetClgRandom11()));
        xMatrix[3] = _make_cuComplex(-xMatrix[1].x, xMatrix[1].y);
        xMatrix[5] = _make_cuComplex(static_cast<Real>(GetClgRandom11()), static_cast<Real>(GetClgRandom11()));
        xMatrix[6] = _make_cuComplex(-xMatrix[2].x, xMatrix[2].y);
        xMatrix[7] = _make_cuComplex(-xMatrix[5].x, xMatrix[5].y);

        CFieldGaugeSU3* pX = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
        if (NULL == pX)
        {
            ++uiFailedRandom;
            continue;
        }
        pX->Zero();
        appSetGaugeLink(pX, uiLinkIndex, xMatrix);
        appSynchronize();

        pRawForce->Zero();
        const CFieldGauge* g0[1] = { pU0 };
        CFieldGauge* f0[1] = { pRawForce };
        pAction->CalculateForce(1, 0, g0, NULL, f0, NULL, NULL, ESP_Once);
        appSynchronize();
        pRawForce->CopyTo(pProjForce);
        pProjForce->LeftMul(pU0, FALSE, TRUE);
        pProjForce->TA();
        appSynchronize();
        const DOUBLE fDforce = -2.0 * static_cast<double>(pProjForce->Dot(pX).x);

        TArray<DOUBLE> fdAbsErrorsR;
        TArray<DOUBLE> fdDfdR;
        UBOOL bAnyNonFiniteR = FALSE;
        for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
        {
            const Real fH = fDelta / static_cast<Real>(1U << iLevel);
            const CFieldGauge* wg[1] = { pWork };
            const CFieldTensor2* t[1] = { pB };

            pU0->CopyTo(pWork);
            appSynchronize();
            pX->ExpMult(fH, pWork);
            appSynchronize();
            const DOUBLE fSPlus = pAction->Energy(FALSE, 1, 0, 1, wg, NULL, t, NULL);

            pU0->CopyTo(pWork);
            appSynchronize();
            pX->ExpMult(-fH, pWork);
            appSynchronize();
            const DOUBLE fSMinus = pAction->Energy(FALSE, 1, 0, 1, wg, NULL, t, NULL);

            const DOUBLE fDfd = (fSPlus - fSMinus) / (2.0 * static_cast<double>(fH));
            const DOUBLE fAbsErr = fabs(fDfd - fDforce);
            fdAbsErrorsR.AddItem(fAbsErr);
            fdDfdR.AddItem(fDfd);
            if (!std::isfinite(fDfd) || !std::isfinite(fAbsErr))
            {
                bAnyNonFiniteR = TRUE;
            }
        }

        UBOOL bConvergedR = FALSE;
        if (!bAnyNonFiniteR && iDeltaLevels >= 3)
        {
            bConvergedR = fdAbsErrorsR[0] > fdAbsErrorsR[1] && fdAbsErrorsR[1] > fdAbsErrorsR[2];
        }

        UBOOL bWithinTolR = FALSE;
        for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
        {
            const DOUBLE fScale = fmax(fabs(fdDfdR[iLevel]), fabs(fDforce));
            const DOUBLE fBound = static_cast<double>(fTolAbs) + static_cast<double>(fTolRel) * fScale;
            if (fdAbsErrorsR[iLevel] <= fBound)
            {
                bWithinTolR = TRUE;
                break;
            }
        }

        CCString sErrSeqR;
        for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
        {
            CCString sItem;
            sItem.Format(_T("h=%.1e err=%.3e%s"),
                static_cast<double>(fDelta / static_cast<Real>(1U << iLevel)),
                fdAbsErrorsR[iLevel],
                (iLevel + 1 < iDeltaLevels) ? _T(", ") : _T(""));
            sErrSeqR = sErrSeqR + sItem;
        }

        const UBOOL bPass = bConvergedR && bWithinTolR;
        if (bPass)
        {
            appGeneral(_T("TestPSU3Force: random-B trial %d (link %u dir %u) PASS, %s\n"),
                iTrial, uiLinkIndex, byDir, sErrSeqR.c_str());
        }
        else
        {
            ++uiFailedRandom;
            appGeneral(_T("TestPSU3Force: random-B trial %d (link %u dir %u) FAIL%s%s, %s\n"),
                iTrial, uiLinkIndex, byDir,
                bConvergedR ? _T("") : _T(" (not h-converged)"),
                bWithinTolR ? _T("") : _T(" (out of tolerance)"),
                sErrSeqR.c_str());
        }
        appSafeDelete(pX);
    }

    // restore the lattice gauge
    pU0->CopyTo(pGauge);
    appSynchronize();

    appSafeDelete(pProjForce);
    appSafeDelete(pRawForce);
    appSafeDelete(pWork);
    appSafeDelete(pU0);

    appGeneral(_T("TestPSU3Force: passed=%u failed=%u (identity B), failed=%u (random B)\n"),
        uiPassed, uiFailed, uiFailedRandom);

    if (0 == uiPassed || 0 != uiFailed || 0 != uiFailedRandom)
    {
        ++uiErrors;
    }
    return uiErrors;
}

//=============================================================================
// 4. unconstrained heatbath conditional distribution on a cold background
//=============================================================================

UINT TestPSU3Heatbath(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CActionGaugePlaquettePSU3WithBoundary* pAction =
        dynamic_cast<CActionGaugePlaquettePSU3WithBoundary*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pB || NULL == pAction)
    {
        appCrucial(_T("TestPSU3Heatbath: missing gauge (1) / boundary (3) / action (1)\n"));
        return 1;
    }

    INT iSweeps = 2000;
    params.FetchValueINT(_T("HeatbathSweeps"), iSweeps);
    Real fFreqTol = F(0.02);
    params.FetchValueReal(_T("FreqTol"), fFreqTol);
    DOUBLE fBeta = 2.0;
    if (params.Exist(_T("Action1")))
    {
        CParameters actionParam = params.GetParameter(_T("Action1"));
        actionParam.FetchValueDOUBLE(_T("Beta"), fBeta);
    }

    // cold gauge background: U_p = 1 for every plaquette
    pGauge->InitialField(EFIT_Identity);
    pB->InitialField(EFIT_Identity);
    appSynchronize();

    CFieldGauge* g[1] = { pGauge };
    CFieldTensor2* t[1] = { pB };

    // warm-up
    for (INT i = 0; i < 50; ++i)
    {
        pAction->OnFinishTrajectory(TRUE, 1, 0, 1, g, NULL, t);
        pAction->OnFinishTrajectory(TRUE);
    }
    appSynchronize();

    UINT uiCount[3] = { 0, 0, 0 };
    const UINT uiElementCount = pB->GetElementCount();
    TArray<UINT> roots;
    roots.SetSize(uiElementCount);

    for (INT i = 0; i < iSweeps; ++i)
    {
        pAction->OnFinishTrajectory(TRUE, 1, 0, 1, g, NULL, t);
        pAction->OnFinishTrajectory(TRUE);
        ReadBFieldRoots(pB, roots.GetData());
        for (UINT e = 0; e < uiElementCount; ++e)
        {
            ++uiCount[roots[e]];
        }
    }
    appSynchronize();

    const double dTotal = static_cast<double>(uiCount[0] + uiCount[1] + uiCount[2]);
    const double dP0 = static_cast<double>(uiCount[0]) / dTotal;
    const double dP1 = static_cast<double>(uiCount[1]) / dTotal;
    const double dP2 = static_cast<double>(uiCount[2]) / dTotal;

    // theory on a cold background: Tr(U_p) = 3, P(k) ~ exp(Beta * Re[zeta^k])
    const double dExp0 = exp(fBeta);
    const double dExp1 = exp(-fBeta / 2.0);
    const double dZ = dExp0 + 2.0 * dExp1;
    const double dT0 = dExp0 / dZ;
    const double dT1 = dExp1 / dZ;

    appGeneral(_T("TestPSU3Heatbath: P = (%f, %f, %f), theory = (%f, %f, %f), samples=%f\n"),
        dP0, dP1, dP2, dT0, dT1, dT1, dTotal);

    if (fabs(dP0 - dT0) > static_cast<double>(fFreqTol)
        || fabs(dP1 - dT1) > static_cast<double>(fFreqTol)
        || fabs(dP2 - dT1) > static_cast<double>(fFreqTol))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3Heatbath: conditional distribution mismatch\n"));
    }

    // AllowMonopole=1 mode must be able to produce non-zero cubes; a random
    // background produces plenty of monopole defects (not checked numerically
    // here, only that the sweep runs and the field stays on the Z3 roots)
    pGauge->InitialField(EFIT_Random);
    pB->InitialField(EFIT_Identity);
    appSynchronize();
    for (INT i = 0; i < 20; ++i)
    {
        pAction->OnFinishTrajectory(TRUE, 1, 0, 1, g, NULL, t);
        pAction->OnFinishTrajectory(TRUE);
    }
    appSynchronize();
    ReadBFieldRoots(pB, roots.GetData());
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (roots[e] > 2)
        {
            ++uiErrors;
            appCrucial(_T("TestPSU3Heatbath: invalid root %u after random-background sweep\n"), roots[e]);
            break;
        }
    }

    appGeneral(_T("TestPSU3Heatbath: %u errors\n"), uiErrors);
    return uiErrors;
}

//=============================================================================
// 5. flat sweep (AllowMonopole=0): link-star + global sheet keep dB=0
//=============================================================================

UINT TestPSU3FlatSweep(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CActionGaugePlaquettePSU3WithBoundary* pAction =
        dynamic_cast<CActionGaugePlaquettePSU3WithBoundary*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pB || NULL == pAction)
    {
        appCrucial(_T("TestPSU3FlatSweep: missing gauge (1) / boundary (3) / action (1)\n"));
        return 1;
    }

    INT iSweepCount = 100;
    params.FetchValueINT(_T("FlatSweepCount"), iSweepCount);
    INT iMinNonZero = 16;
    params.FetchValueINT(_T("MinNonZero"), iMinNonZero);

    // start from flat B = 1 on a cold background; CheckMonopole=1 makes the
    // action itself abort if any sweep creates a monopole
    pGauge->InitialField(EFIT_Identity);
    pB->InitialField(EFIT_Identity);
    appSynchronize();

    CFieldGauge* g[1] = { pGauge };
    CFieldTensor2* t[1] = { pB };

    for (INT i = 0; i < iSweepCount; ++i)
    {
        pAction->OnFinishTrajectory(TRUE, 1, 0, 1, g, NULL, t);
        pAction->OnFinishTrajectory(TRUE);
    }
    appSynchronize();

    // the sweeps (link-star + global sheet) must have moved B away from
    // identity: at Beta=0.5 a fraction of plaquettes carries a non-trivial root
    const UINT uiElementCount = pB->GetElementCount();
    TArray<UINT> roots;
    roots.SetSize(uiElementCount);
    ReadBFieldRoots(pB, roots.GetData());
    UINT uiNonZero = 0;
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (roots[e] > 2)
        {
            ++uiErrors;
            appCrucial(_T("TestPSU3FlatSweep: invalid root %u after flat sweep\n"), roots[e]);
            break;
        }
        if (0 != roots[e])
        {
            ++uiNonZero;
        }
    }
    appGeneral(_T("TestPSU3FlatSweep: non-identity plaquette count = %u / %u\n"), uiNonZero, uiElementCount);
    if (uiNonZero < static_cast<UINT>(iMinNonZero))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3FlatSweep: sweeps did not move B (all identity?)\n"));
    }

    // flatness is guaranteed by CheckMonopole=1 (any violation aborts), and
    // the energy must stay finite
    const CFieldGauge* g2[1] = { pGauge };
    const CFieldTensor2* t2[1] = { pB };
    const DOUBLE eFinal = pAction->Energy(FALSE, 1, 0, 1, g2, NULL, t2, NULL);
    if (!std::isfinite(eFinal))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3FlatSweep: non-finite final energy\n"));
    }

    appGeneral(_T("TestPSU3FlatSweep: %u errors\n"), uiErrors);
    return uiErrors;
}

//=============================================================================
// 6. trajectory callback / post-sweep energy cache
//=============================================================================

UINT TestPSU3CallbackCache(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CActionGaugePlaquettePSU3WithBoundary* pAction =
        dynamic_cast<CActionGaugePlaquettePSU3WithBoundary*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pB || NULL == pAction)
    {
        appCrucial(_T("TestPSU3CallbackCache: missing gauge (1) / boundary (3) / action (1)\n"));
        return 1;
    }

    DOUBLE fTol = 1.0e-4;
    params.FetchValueDOUBLE(_T("CacheTol"), fTol);

    CFieldGauge* g[1] = { pGauge };
    CFieldTensor2* t[1] = { pB };
    const CFieldGauge* gc[1] = { pGauge };
    const CFieldTensor2* tc[1] = { pB };

    // PrepareForHMC initializes m_fLastEnergy to the current state
    pAction->PrepareForHMC(1, 0, g, NULL, 0);
    appSynchronize();

    const DOUBLE e0 = pAction->Energy(FALSE, 1, 0, 1, gc, NULL, tc, NULL);
    const DOUBLE eCached = pAction->Energy(TRUE, 1, 0, 1, gc, NULL, tc, NULL);
    if (fabs(eCached - e0) > fTol * fmax(1.0, fabs(e0)))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3CallbackCache: PrepareForHMC cache mismatch (%e vs %e)\n"), eCached, e0);
    }

    // an accepted trajectory: sweep B, record post-sweep energy, publish it
    pAction->OnFinishTrajectory(TRUE, 1, 0, 1, g, NULL, t);
    pAction->OnFinishTrajectory(TRUE);
    appSynchronize();

    const DOUBLE eRecompute = pAction->Energy(FALSE, 1, 0, 1, gc, NULL, tc, NULL);
    const DOUBLE eCached2 = pAction->Energy(TRUE, 1, 0, 1, gc, NULL, tc, NULL);
    if (fabs(eCached2 - eRecompute) > fTol * fmax(1.0, fabs(eRecompute)))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3CallbackCache: post-sweep cache mismatch after accept (%e vs %e)\n"), eCached2, eRecompute);
    }

    // a rejected trajectory: the lattice (U,B) is kept, B is swept on it, and
    // the cache must follow the swept state
    pAction->OnFinishTrajectory(FALSE, 1, 0, 1, g, NULL, t);
    pAction->OnFinishTrajectory(FALSE);
    appSynchronize();

    const DOUBLE eRecompute2 = pAction->Energy(FALSE, 1, 0, 1, gc, NULL, tc, NULL);
    const DOUBLE eCached3 = pAction->Energy(TRUE, 1, 0, 1, gc, NULL, tc, NULL);
    if (fabs(eCached3 - eRecompute2) > fTol * fmax(1.0, fabs(eRecompute2)))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3CallbackCache: post-sweep cache mismatch after reject (%e vs %e)\n"), eCached3, eRecompute2);
    }

    appGeneral(_T("TestPSU3CallbackCache: %u errors\n"), uiErrors);
    return uiErrors;
}

//=============================================================================
// 7. tensor2 companion save/load round-trip through CUpdator
//=============================================================================

UINT TestPSU3Tensor2SaveLoad(CParameters& params)
{
    UINT uiErrors = 0;

    CFieldTensor2Z3* pB = dynamic_cast<CFieldTensor2Z3*>(appGetLattice()->GetFieldById(3));
    CUpdator* pUpdator = appGetLattice()->m_pUpdator;
    if (NULL == pB || NULL == pUpdator)
    {
        appCrucial(_T("TestPSU3Tensor2SaveLoad: missing boundary field (3) or updator\n"));
        return 1;
    }

    // random B, record the saved state
    pB->InitialField(EFIT_Random);
    appSynchronize();
    const UINT uiElementCount = pB->GetElementCount();
    TArray<UINT> savedRoots;
    savedRoots.SetSize(uiElementCount);
    ReadBFieldRoots(pB, savedRoots.GetData());
    appSynchronize();

    // save the configuration: gauge .con + tensor2 companion .con
    pUpdator->SetSaveConfiguration(TRUE, _T("testPSU3"), 0, EFFT_CLGBin);
    pUpdator->SetConfigurationCount(0);
    pUpdator->SaveConfiguration(0);
    appSynchronize();

    // the companion file must exist
    const CCString sT2File = _T("testPSU3_0_t3.con");
    if (!CFileSystem::IsFileExist(sT2File))
    {
        ++uiErrors;
        appCrucial(_T("TestPSU3Tensor2SaveLoad: companion file %s not written\n"), sT2File.c_str());
        return uiErrors;
    }

    // destroy B, then reload from the companion file
    pB->InitialField(EFIT_Identity);
    pUpdator->LoadTensor2Configuration(_T("testPSU3"), 0, EFFT_CLGBin);
    appSynchronize();

    TArray<UINT> loadedRoots;
    loadedRoots.SetSize(uiElementCount);
    ReadBFieldRoots(pB, loadedRoots.GetData());
    for (UINT e = 0; e < uiElementCount; ++e)
    {
        if (loadedRoots[e] != savedRoots[e])
        {
            ++uiErrors;
            appGeneral(_T("TestPSU3Tensor2SaveLoad: mismatch at element %u (%u vs %u)\n"), e, loadedRoots[e], savedRoots[e]);
            break;
        }
    }

    // cleanup the temporary configuration files
    remove("testPSU3_0_t3.con");
    remove("testPSU3_0.con");
    remove("testPSU3_0.txt");

    appGeneral(_T("TestPSU3Tensor2SaveLoad: %u errors\n"), uiErrors);
    return uiErrors;
}

__REGIST_TEST(TestPSU3Tensor2SaveLoad, PSU3, TestPSU3Tensor2SaveLoad, "PSU3WithBoundary tensor2 companion save/load")
__REGIST_TEST(TestZ3FieldCreate, PSU3, TestZ3FieldCreate, "Z3 tensor2 field creation/init/copy/IO")
__REGIST_TEST(TestPSU3Energy, PSU3, TestPSU3Energy, "PSU3WithBoundary energy (B=1 vs fundamental)")
__REGIST_TEST(TestPSU3Force, PSU3, TestPSU3Force, "PSU3WithBoundary force finite difference")
__REGIST_TEST(TestPSU3Heatbath, PSU3, TestPSU3Heatbath, "PSU3WithBoundary unconstrained heatbath")
__REGIST_TEST(TestPSU3FlatSweep, PSU3, TestPSU3FlatSweep, "PSU3WithBoundary flat sweep dB=0")
__REGIST_TEST(TestPSU3CallbackCache, PSU3, TestPSU3CallbackCache, "PSU3WithBoundary callback/cache")
