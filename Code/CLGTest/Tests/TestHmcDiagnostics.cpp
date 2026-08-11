//=============================================================================
// FILENAME : TestHmcDiagnostics.cpp
//
// DESCRIPTION:
// HMC diagnostic suite (revision v3).  The suite contains two diagnostics:
//   1. Molecular-dynamics reversibility for a fixed Hamiltonian state.
//   3. Finite-difference validation of the production force for each action.
//
// All diagnostics are driven by YAML parameters.
//=============================================================================

#include "CLGTest.h"
#include <cmath>

#pragma region helpers

/**
 * Returns a random double in [-1, 1) using the lattice's own random generator.
 */
static inline double GetClgRandom11()
{
    CRandom* pRand = appGetLattice()->m_pRandom;
    if (NULL == pRand)
    {
        return 0.0;
    }
    return static_cast<double>(pRand->GetRandomF()) * 2.0 - 1.0;
}

/**
 * Compute max and average Frobenius norm of a gauge-field matrix.
 * Returns FALSE if the norm cannot be computed; output values are left unset.
 */
static UBOOL SU3MaxAvgNorm(const CFieldGauge* pField, DOUBLE& fMaxNorm, DOUBLE& fAvgNorm)
{
    if (NULL == pField)
    {
        LastProbem(_T("SU3MaxAvgNorm: null field\n"));
        return FALSE;
    }

    const UINT uiLinkCount = pField->GetLinkCount();
    if (0 == uiLinkCount)
    {
        LastProbem(_T("SU3MaxAvgNorm: empty field\n"));
        return FALSE;
    }

    UINT uiSize = 0;
    BYTE* pData = pField->CopyDataOut(uiSize);
    if (NULL == pData || 0 == uiSize)
    {
        LastProbem(_T("SU3MaxAvgNorm: failed to copy field data\n"));
        return FALSE;
    }

    const CLGComplex* pM = reinterpret_cast<const CLGComplex*>(pData);

    DOUBLE fTotalSq = 0.0;
    DOUBLE fLocalMax = 0.0;
    for (UINT uiLink = 0; uiLink < uiLinkCount; ++uiLink)
    {
        const CLGComplex* pMe = pM + uiLink * 9;
        DOUBLE fSq = 0.0;
        for (UINT i = 0; i < 9; ++i)
        {
            fSq += static_cast<DOUBLE>(pMe[i].x) * pMe[i].x
                 + static_cast<DOUBLE>(pMe[i].y) * pMe[i].y;
        }
        fTotalSq += fSq;
        const DOUBLE fNorm = sqrt(fSq);
        if (fNorm > fLocalMax)
        {
            fLocalMax = fNorm;
        }
    }
    free(pData);

    fAvgNorm = sqrt(fTotalSq / static_cast<DOUBLE>(uiLinkCount));
    fMaxNorm = fLocalMax;

    if (!std::isfinite(fAvgNorm) || !std::isfinite(fMaxNorm)
     || fAvgNorm < 0.0 || fMaxNorm < 0.0)
    {
        LastProbem(_T("SU3MaxAvgNorm: non-finite or negative norm\n"));
        return FALSE;
    }

    return TRUE;
}

/**
 * Difference between two SU3 gauge fields: max and average Frobenius norm.
 * Returns FALSE if the difference norm cannot be computed.
 */
static UBOOL SU3DiffNorms(const CFieldGauge* pA, const CFieldGauge* pB, DOUBLE& fMaxNorm, DOUBLE& fAvgNorm)
{
    if (NULL == pA || NULL == pB)
    {
        LastProbem(_T("SU3DiffNorms: null input field\n"));
        return FALSE;
    }

    CFieldGauge* pDiff = dynamic_cast<CFieldGauge*>(pA->GetCopy());
    if (NULL == pDiff)
    {
        LastProbem(_T("SU3DiffNorms: failed to allocate diff field\n"));
        return FALSE;
    }
    pA->CopyTo(pDiff);
    pDiff->AxpyMinus(pB);

    const UBOOL bOk = SU3MaxAvgNorm(pDiff, fMaxNorm, fAvgNorm);
    appSafeDelete(pDiff);
    return bOk;
}

/**
 * Set a single link of a gauge field to a host-computed 3x3 complex matrix.
 * The matrix is supplied as 9 CLGComplex numbers in row-major order.
 */
static void SetSingleLinkMatrix(CFieldGauge* pField, UINT uiLinkIndex, const CLGComplex* pMatrix)
{
    appSetGaugeLink(pField, uiLinkIndex, pMatrix);
    appSynchronize();
}

/**
 * Generate a random anti-Hermitian traceless SU(3) matrix with unit Frobenius norm.
 * Returns 9 CLGComplex numbers in row-major order.
 */
static void RandomSu3Direction(CLGComplex* pOut)
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

    // Anti-Hermitian: X^dag = -X, traceless.
    pOut[0] = _make_cuComplex(static_cast<Real>(0.0), static_cast<Real>(a1));
    pOut[1] = _make_cuComplex(static_cast<Real>(x01_re), static_cast<Real>(x01_im));
    pOut[2] = _make_cuComplex(static_cast<Real>(x02_re), static_cast<Real>(x02_im));
    pOut[3] = _make_cuComplex(static_cast<Real>(-x01_re), static_cast<Real>(x01_im));
    pOut[4] = _make_cuComplex(static_cast<Real>(0.0), static_cast<Real>(a2));
    pOut[5] = _make_cuComplex(static_cast<Real>(x12_re), static_cast<Real>(x12_im));
    pOut[6] = _make_cuComplex(static_cast<Real>(-x02_re), static_cast<Real>(x02_im));
    pOut[7] = _make_cuComplex(static_cast<Real>(-x12_re), static_cast<Real>(x12_im));
    pOut[8] = _make_cuComplex(static_cast<Real>(0.0), static_cast<Real>(a3));

    double fNormSq = 2.0 * (x01_re * x01_re + x01_im * x01_im
                          + x02_re * x02_re + x02_im * x02_im
                          + x12_re * x12_re + x12_im * x12_im)
                   + (a1 * a1 + a2 * a2 + a3 * a3);
    double fScale = 1.0 / sqrt(fNormSq);
    for (UINT i = 0; i < 9; ++i)
    {
        pOut[i].x = static_cast<Real>(pOut[i].x * fScale);
        pOut[i].y = static_cast<Real>(pOut[i].y * fScale);
    }
}

/**
 * Reproduce the production ordering in CIntegrator::OnCacheAndSmearing(3).
 * This refreshes staple caches, gauge smearing, and post-smearing caches for
 * the supplied thin gauge field.  It must be called before every force or
 * energy evaluation and after restoring the real lattice gauge field.
 */
static void RefreshGaugeDependentState(CFieldGauge* pGauge)
{
    if (NULL == pGauge)
    {
        return;
    }

    CStapleCache* pCache = appGetStapleCache(pGauge->m_byFieldId);
    if (NULL != pCache)
    {
        pCache->Cache(pGauge, ECC_BeforeGaugeUpdate);
        pCache->Cache(pGauge, ECC_BeforeAllUpdateBeforeSmearing);
    }

    CGaugeSmearing* pSmear = appGetGaugeSmearing(pGauge->m_byFieldId);
    if (NULL != pSmear && pSmear->CalledWhenUpdate())
    {
        pSmear->GaugeSmearingC(pGauge);
    }

    if (NULL != pCache)
    {
        pCache->Cache(pGauge, ECC_BeforeAllUpdateAfterSmearing);
    }

    appSynchronize();
}

/**
 * Build a list of dynamical link indices for the supplied gauge field ID.
 * Excludes links whose SIndex tag contains _kDirichlet and directions masked
 * by the integrator BindDir.
 */
static UBOOL BuildDynamicalLinkList(BYTE byFieldId, BYTE byBindDir,
                                    TArray<UINT>& linkIndices, TArray<BYTE>& dirs)
{
    linkIndices.RemoveAll();
    dirs.RemoveAll();

    CIndexData* pIndexData = appGetLattice()->m_pIndexCache;
    if (NULL == pIndexData || NULL == pIndexData->m_pIndexLinkToSIndex[byFieldId])
    {
        LastProbem(_T("BuildDynamicalLinkList: index cache missing for field\n"));
        return FALSE;
    }

    const UINT uiBigX = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    const UINT uiBigY = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    const UINT uiBigZ = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    const UINT uiBigT = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
    const UINT uiBigYZT = uiBigY * uiBigZ * uiBigT;
    const UINT uiBigZT  = uiBigZ * uiBigT;
    const UINT uiBigTsingle = uiBigT;

    const UINT uiTableSize = uiBigX * uiBigYZT * _HC_Dir;
    SIndex* pHostTable = reinterpret_cast<SIndex*>(malloc(sizeof(SIndex) * uiTableSize));
    if (NULL == pHostTable)
    {
        LastProbem(_T("BuildDynamicalLinkList: host allocation failed\n"));
        return FALSE;
    }

    appSimpleCopyDH(pHostTable, pIndexData->m_pIndexLinkToSIndex[byFieldId],
                    sizeof(SIndex) * uiTableSize);
    appSynchronize();

    for (UINT uiSite = 0; uiSite < static_cast<UINT>(_HC_Volume); ++uiSite)
    {
        const SSmallInt4 coord = __hostSiteIndexToInt4(uiSite);
        const UINT uiBigIdx =
              (coord.x + CIndexData::kCacheIndexEdge) * uiBigYZT
            + (coord.y + CIndexData::kCacheIndexEdge) * uiBigZT
            + (coord.z + CIndexData::kCacheIndexEdge) * uiBigTsingle
            + (coord.w + CIndexData::kCacheIndexEdge);

        for (BYTE byDir = 0; byDir < static_cast<BYTE>(_HC_Dir); ++byDir)
        {
            if (0 != (byBindDir & (1 << byDir)))
            {
                continue;
            }

            const UINT uiTableIdx = uiBigIdx * _HC_Dir + byDir;
            const SIndex& idx = pHostTable[uiTableIdx];
            if (0 != (idx.m_byTag & _kDirichlet))
            {
                continue;
            }

            const UINT uiLinkIndex = uiSite * _HC_Dir + byDir;
            linkIndices.AddItem(uiLinkIndex);
            dirs.AddItem(byDir);
        }
    }

    free(pHostTable);

    if (0 == linkIndices.Num())
    {
        LastProbem(_T("BuildDynamicalLinkList: no dynamical links found\n"));
        return FALSE;
    }
    return TRUE;
}

#pragma endregion

#pragma region diagnostic 1: MD reversibility

UINT TestHmcReversibility(CParameters& sParam)
{
    UINT uiError = 0;

    CUpdator* pUpdator = appGetLattice()->m_pUpdator;
    CHMC* pHMC = dynamic_cast<CHMC*>(pUpdator);
    if (NULL == pHMC)
    {
        LastProbem(_T("TestHmcReversibility: updator is not CHMC\n"));
        return 1;
    }
    CIntegrator* pIntegrator = pHMC->m_pIntegrator;
    if (NULL == pIntegrator)
    {
        LastProbem(_T("TestHmcReversibility: integrator is NULL\n"));
        return 1;
    }

    pUpdator->SetAutoCorrection(FALSE);
    pUpdator->SetTestHdiff(FALSE);
    pIntegrator->FixStep(FALSE);

    Real fTolAvgU = F(1.0e-10);
    Real fTolMaxU = F(1.0e-9);
    Real fTolAvgP = F(1.0e-9);
    Real fTolMaxP = F(1.0e-8);
    sParam.FetchValueReal(_T("TolAvgU"), fTolAvgU);
    sParam.FetchValueReal(_T("TolMaxU"), fTolMaxU);
    sParam.FetchValueReal(_T("TolAvgP"), fTolAvgP);
    sParam.FetchValueReal(_T("TolMaxP"), fTolMaxP);

    pIntegrator->Prepare(FALSE, 0);
    appSynchronize();

    const INT iFieldCount = pIntegrator->m_pGaugeField.Num();
    if (iFieldCount < 1)
    {
        LastProbem(_T("TestHmcReversibility: no gauge fields\n"));
        return 1;
    }
    if (pIntegrator->m_pMomentumField.Num() != iFieldCount)
    {
        LastProbem(_T("TestHmcReversibility: gauge and momentum arrays have different lengths\n"));
        return 1;
    }

    TArray<CFieldGauge*> initialU;
    TArray<CFieldGauge*> initialP;

    auto CleanupSavedState = [&initialU, &initialP]()
    {
        for (INT j = 0; j < initialU.Num(); ++j)
        {
            appSafeDelete(initialU[j]);
        }
        for (INT j = 0; j < initialP.Num(); ++j)
        {
            appSafeDelete(initialP[j]);
        }
        initialU.RemoveAll();
        initialP.RemoveAll();
    };

    for (INT i = 0; i < iFieldCount; ++i)
    {
        CFieldGauge* pCurrentU = pIntegrator->m_pGaugeField[i];
        CFieldGauge* pCurrentP = pIntegrator->m_pMomentumField[i];
        if (NULL == pCurrentU || NULL == pCurrentP)
        {
            LastProbem(_T("TestHmcReversibility: null current gauge or momentum field\n"));
            CleanupSavedState();
            return 1;
        }

        CFieldGauge* pU0 = dynamic_cast<CFieldGauge*>(pCurrentU->GetCopy());
        if (NULL == pU0)
        {
            LastProbem(_T("TestHmcReversibility: failed to allocate saved gauge\n"));
            CleanupSavedState();
            return 1;
        }
        pCurrentU->CopyTo(pU0);
        initialU.AddItem(pU0);

        CFieldGauge* pP0 = dynamic_cast<CFieldGauge*>(pCurrentP->GetCopy());
        if (NULL == pP0)
        {
            LastProbem(_T("TestHmcReversibility: failed to allocate saved momentum\n"));
            CleanupSavedState();
            return 1;
        }
        pCurrentP->CopyTo(pP0);
        initialP.AddItem(pP0);
    }

    pIntegrator->Evaluate();
    appSynchronize();

    for (INT i = 0; i < pIntegrator->m_pMomentumField.Num(); ++i)
    {
        if (NULL != pIntegrator->m_pMomentumField[i])
        {
            pIntegrator->m_pMomentumField[i]->ScalarMultply(F(-1.0));
        }
    }
    appSynchronize();

    pIntegrator->Evaluate();
    appSynchronize();

    for (INT i = 0; i < iFieldCount; ++i)
    {
        CFieldGauge* pCurrentU = pIntegrator->m_pGaugeField[i];
        CFieldGauge* pCurrentP = pIntegrator->m_pMomentumField[i];
        CFieldGauge* pSavedU = initialU[i];
        CFieldGauge* pSavedP = initialP[i];

        if (NULL == pCurrentU || NULL == pCurrentP || NULL == pSavedU || NULL == pSavedP)
        {
            LastProbem(_T("TestHmcReversibility: missing field for comparison\n"));
            ++uiError;
            continue;
        }

        DOUBLE fMaxU = 0.0, fAvgU = 0.0;
        DOUBLE fMaxP = 0.0, fAvgP = 0.0;

        if (!SU3DiffNorms(pSavedU, pCurrentU, fMaxU, fAvgU))
        {
            CCString sProblem;
            sProblem.Format(_T("TestHmcReversibility: failed to compute gauge difference norm for field %d\n"), i);
            LastProbem(sProblem);
            ++uiError;
            continue;
        }

        CFieldGauge* pMDiff = dynamic_cast<CFieldGauge*>(pSavedP->GetCopy());
        if (NULL == pMDiff)
        {
            LastProbem(_T("TestHmcReversibility: failed to allocate momentum diff\n"));
            ++uiError;
            continue;
        }
        pSavedP->CopyTo(pMDiff);
        pMDiff->AxpyPlus(pCurrentP);
        const UBOOL bMomentumNormOk = SU3MaxAvgNorm(pMDiff, fMaxP, fAvgP);
        appSafeDelete(pMDiff);

        if (!bMomentumNormOk)
        {
            CCString sProblem;
            sProblem.Format(_T("TestHmcReversibility: failed to compute momentum norm for field %d\n"), i);
            LastProbem(sProblem);
            ++uiError;
            continue;
        }

        appGeneral(_T("Reversibility field %d: avgU=%.6e maxU=%.6e (tol avg/max %.1e/%.1e), avgP=%.6e maxP=%.6e (tol avg/max %.1e/%.1e)\n"),
                   i, fAvgU, fMaxU, fTolAvgU, fTolMaxU, fAvgP, fMaxP, fTolAvgP, fTolMaxP);

        if (fAvgU > static_cast<DOUBLE>(fTolAvgU) || fMaxU > static_cast<DOUBLE>(fTolMaxU)
         || fAvgP > static_cast<DOUBLE>(fTolAvgP) || fMaxP > static_cast<DOUBLE>(fTolMaxP))
        {
            CCString sProblem;
            sProblem.Format(_T("Reversibility failed for field %d: avgU=%.6e maxU=%.6e avgP=%.6e maxP=%.6e\n"),
                            i, fAvgU, fMaxU, fAvgP, fMaxP);
            LastProbem(sProblem);
            ++uiError;
        }
    }

    CleanupSavedState();

    return uiError;
}

#pragma endregion

#pragma region diagnostic 3: production-force finite difference

/**
 * RAII guard for the finite-difference diagnostic.
 * On destruction it restores the real lattice gauge from U0 and refreshes all
 * derived gauge-dependent state.
 */
struct FiniteDiffCleanupGuard
{
    CFieldGauge* pGauge;
    CFieldGauge* pU0;
    CFieldGauge* pWorkGauge;
    CFieldGauge* pRawForce;
    CFieldGauge* pProjForce;
    UBOOL bRestore;

    FiniteDiffCleanupGuard(CFieldGauge* pGaugeIn, CFieldGauge* pU0In)
        : pGauge(pGaugeIn)
        , pU0(pU0In)
        , pWorkGauge(NULL)
        , pRawForce(NULL)
        , pProjForce(NULL)
        , bRestore(TRUE)
    {
    }

    ~FiniteDiffCleanupGuard()
    {
        if (bRestore && NULL != pU0 && NULL != pGauge)
        {
            pU0->CopyTo(pGauge);
            RefreshGaugeDependentState(pGauge);
            appSynchronize();
        }
        appSafeDelete(pProjForce);
        appSafeDelete(pRawForce);
        appSafeDelete(pWorkGauge);
        appSafeDelete(pU0);
    }
};

static UINT TestActionForceFiniteDiff(CParameters& sParam, CAction* pAction, const TCHAR* szName)
{
    UINT uiError = 0;

    if (NULL == pAction || pAction->IsDiscreteGauge())
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): no continuous action selected\n"), szName);
        LastProbem(sProblem);
        return 1;
    }

    INT iLinkTrials = 3;
    INT iDirTrials = 2;
    Real fDelta = F(1.0e-2);
    INT iDeltaLevels = 3;
    Real fTolRel = F(0.30);
    Real fTolAbs = F(5.0e-3);
    sParam.FetchValueINT(_T("LinkTrials"), iLinkTrials);
    sParam.FetchValueINT(_T("DirTrials"), iDirTrials);
    sParam.FetchValueReal(_T("Delta"), fDelta);
    sParam.FetchValueINT(_T("DeltaLevels"), iDeltaLevels);
    sParam.FetchValueReal(_T("TolRel"), fTolRel);
    sParam.FetchValueReal(_T("TolAbs"), fTolAbs);

    if (iLinkTrials <= 0 || iDirTrials <= 0)
    {
        LastProbem(_T("TestActionForceFiniteDiff: LinkTrials and DirTrials must be positive\n"));
        return 1;
    }
    if (fDelta <= F(0.0) || iDeltaLevels < 2)
    {
        LastProbem(_T("TestActionForceFiniteDiff: Delta > 0 and DeltaLevels >= 2 required\n"));
        return 1;
    }
    if (iDeltaLevels > 32)
    {
        LastProbem(_T("TestActionForceFiniteDiff: DeltaLevels too large for step-halving\n"));
        return 1;
    }
    if (!std::isfinite(static_cast<double>(fTolAbs)) || !std::isfinite(static_cast<double>(fTolRel))
        || fTolAbs < F(0.0) || fTolRel < F(0.0))
    {
        LastProbem(_T("TestActionForceFiniteDiff: tolerances must be finite and nonnegative\n"));
        return 1;
    }

    CUpdator* pUpdator = appGetLattice()->m_pUpdator;
    CHMC* pHMC = dynamic_cast<CHMC*>(pUpdator);
    BYTE byBindDir = 0;
    if (NULL != pHMC && NULL != pHMC->m_pIntegrator)
    {
        byBindDir = pHMC->m_pIntegrator->GetBindDir();
    }

    CLatticeData* pLattice = appGetLattice();
    const TArray<BYTE>& gaugeFieldIds = pAction->GetGaugeFieldIds();

    if (gaugeFieldIds.Num() < 1)
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): action has no gauge field ids\n"), szName);
        LastProbem(sProblem);
        return 1;
    }

    CFieldGauge* pGauge = NULL;
    for (INT i = 0; i < gaugeFieldIds.Num(); ++i)
    {
        CField* pField = pLattice->GetFieldById(gaugeFieldIds[i]);
        CFieldGauge* pG = dynamic_cast<CFieldGauge*>(pField);
        if (NULL == pG)
        {
            CCString sProblem;
            sProblem.Format(_T("TestActionForceFiniteDiff (%s): required gauge field id %u missing or wrong type\n"), szName, gaugeFieldIds[i]);
            LastProbem(sProblem);
            return 1;
        }
        if (NULL == pGauge)
        {
            pGauge = pG;
        }
    }

    if (NULL == pGauge || pGauge->GetFieldType() != EFT_GaugeSU3)
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): primary gauge field is not SU3\n"), szName);
        LastProbem(sProblem);
        return 1;
    }

    TArray<UINT> dynamicalLinks;
    TArray<BYTE> dynamicalDirs;
    if (!BuildDynamicalLinkList(pGauge->m_byFieldId, byBindDir, dynamicalLinks, dynamicalDirs))
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): failed to build dynamical link list\n"), szName);
        LastProbem(sProblem);
        return 1;
    }
    const UINT uiDynamicalCount = dynamicalLinks.Num();

    CFieldGauge* pU0 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
    if (NULL == pU0)
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): failed to allocate U0 copy\n"), szName);
        LastProbem(sProblem);
        return 1;
    }
    pGauge->CopyTo(pU0);
    appSynchronize();

    FiniteDiffCleanupGuard guard(pGauge, pU0);

    CFieldGauge* pWorkGauge = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
    if (NULL == pWorkGauge)
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): failed to allocate work gauge\n"), szName);
        LastProbem(sProblem);
        return uiError + 1;
    }
    guard.pWorkGauge = pWorkGauge;

    CFieldGauge* pRawForce = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
    if (NULL == pRawForce)
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): failed to allocate raw force\n"), szName);
        LastProbem(sProblem);
        return uiError + 1;
    }
    pRawForce->Zero();
    guard.pRawForce = pRawForce;

    CFieldGauge* pProjForce = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
    if (NULL == pProjForce)
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): failed to allocate projected force\n"), szName);
        LastProbem(sProblem);
        return uiError + 1;
    }
    guard.pProjForce = pProjForce;

    TArray<const CFieldGauge*> gaugeArray;
    gaugeArray.AddItem(pU0);

    RefreshGaugeDependentState(pU0);
    if (pAction->IsFermion())
    {
        pAction->PrepareForHMC(1, 0, gaugeArray.GetData(), NULL, 0);
        appSynchronize();
    }

    RefreshGaugeDependentState(pU0);
    pRawForce->Zero();
    CFieldGauge* pRawForces[1] = { pRawForce };
    if (!pAction->CalculateForce(1, 0, gaugeArray.GetData(), NULL, pRawForces, NULL, NULL, ESP_Once))
    {
        CCString sProblem;
        sProblem.Format(_T("TestActionForceFiniteDiff (%s): CalculateForce failed\n"), szName);
        LastProbem(sProblem);
        return uiError + 1;
    }
    appSynchronize();

    pRawForce->CopyTo(pProjForce);
    pProjForce->LeftMul(pU0, FALSE, TRUE);
    pProjForce->TA();
    pProjForce->SetOneDirectionZero(byBindDir);
    appSynchronize();

    const UINT uiTotalTrials = static_cast<UINT>(iLinkTrials) * static_cast<UINT>(iDirTrials);
    UINT uiPassedTrials = 0;
    UINT uiFailedTrials = 0;

    //stratify the candidates by link direction, so that with LinkTrials >= 4
    //every direction (including the temporal one) is covered at least once
    TArray<UINT> dirCandidates[4];
    for (INT i = 0; i < dynamicalLinks.Num(); ++i)
    {
        dirCandidates[dynamicalDirs[i]].AddItem(dynamicalLinks[i]);
    }

    for (INT iLink = 0; iLink < iLinkTrials; ++iLink)
    {
        const BYTE byDirWanted = static_cast<BYTE>(iLink & 3);
        UINT uiLinkIndex;
        BYTE byDir;
        if (dirCandidates[byDirWanted].Num() > 0)
        {
            const UINT uiRandLink = static_cast<UINT>(fabs(GetClgRandom11()) * static_cast<double>(dirCandidates[byDirWanted].Num())) % static_cast<UINT>(dirCandidates[byDirWanted].Num());
            uiLinkIndex = dirCandidates[byDirWanted][uiRandLink];
            byDir = byDirWanted;
        }
        else
        {
            const UINT uiRandLink = static_cast<UINT>(fabs(GetClgRandom11()) * static_cast<double>(uiDynamicalCount)) % uiDynamicalCount;
            uiLinkIndex = dynamicalLinks[uiRandLink];
            byDir = dynamicalDirs[uiRandLink];
        }

        for (INT iDir = 0; iDir < iDirTrials; ++iDir)
        {
            CLGComplex xMatrix[9];
            RandomSu3Direction(xMatrix);

            CFieldGauge* pX = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            if (NULL == pX)
            {
                CCString sProblem;
                sProblem.Format(_T("TestActionForceFiniteDiff (%s): failed to allocate direction field\n"), szName);
                LastProbem(sProblem);
                ++uiFailedTrials;
                pU0->CopyTo(pGauge);
                RefreshGaugeDependentState(pGauge);
                appSynchronize();
                continue;
            }
            pX->Zero();
            SetSingleLinkMatrix(pX, uiLinkIndex, xMatrix);
            appSynchronize();

            const DOUBLE fDforce = -2.0 * static_cast<double>(pProjForce->Dot(pX).x);
            if (!std::isfinite(fDforce))
            {
                CCString sProblem;
                sProblem.Format(_T("TestActionForceFiniteDiff (%s): non-finite force prediction\n"), szName);
                LastProbem(sProblem);
                ++uiFailedTrials;
                pU0->CopyTo(pGauge);
                RefreshGaugeDependentState(pGauge);
                appSynchronize();
                appSafeDelete(pX);
                continue;
            }

            TArray<Real> fdSteps;
            TArray<DOUBLE> fdDfd;
            TArray<DOUBLE> fdAbsErrors;
            TArray<DOUBLE> fdBounds;

            UBOOL bAnyNonFinite = FALSE;
            for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
            {
                const Real fH = fDelta / static_cast<Real>(1U << iLevel);
                const CFieldGauge* pWorkGauges[1] = { pWorkGauge };

                pU0->CopyTo(pWorkGauge);
                appSynchronize();
                pX->ExpMult(fH, pWorkGauge);
                pWorkGauge->SetOneDirectionUnity(byBindDir);
                RefreshGaugeDependentState(pWorkGauge);
                const DOUBLE fSPlus = pAction->Energy(FALSE, 1, 0, 0, pWorkGauges, NULL, NULL, NULL);

                pU0->CopyTo(pWorkGauge);
                appSynchronize();
                pX->ExpMult(-fH, pWorkGauge);
                pWorkGauge->SetOneDirectionUnity(byBindDir);
                RefreshGaugeDependentState(pWorkGauge);
                const DOUBLE fSMinus = pAction->Energy(FALSE, 1, 0, 0, pWorkGauges, NULL, NULL, NULL);

                const DOUBLE fDfd = (fSPlus - fSMinus) / (2.0 * static_cast<double>(fH));
                const DOUBLE fAbsErr = fabs(fDfd - fDforce);
                const DOUBLE fScale = fmax(fabs(fDfd), fabs(fDforce));
                const DOUBLE fBound = static_cast<double>(fTolAbs) + static_cast<double>(fTolRel) * fScale;

                fdSteps.AddItem(fH);
                fdDfd.AddItem(fDfd);
                fdAbsErrors.AddItem(fAbsErr);
                fdBounds.AddItem(fBound);

                if (!std::isfinite(fSPlus) || !std::isfinite(fSMinus)
                 || !std::isfinite(fDfd) || !std::isfinite(fAbsErr) || !std::isfinite(fBound))
                {
                    bAnyNonFinite = TRUE;
                }
            }

            INT iBestLevel = -1;
            DOUBLE fBestAbsErr = 1.0e300;
            for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
            {
                if (std::isfinite(fdAbsErrors[iLevel]) && fdAbsErrors[iLevel] < fBestAbsErr)
                {
                    fBestAbsErr = fdAbsErrors[iLevel];
                    iBestLevel = iLevel;
                }
            }

            UBOOL bTrialPassed = FALSE;
            if (!bAnyNonFinite && std::isfinite(fDforce))
            {
                for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
                {
                    if (fdAbsErrors[iLevel] <= fdBounds[iLevel])
                    {
                        bTrialPassed = TRUE;
                        break;
                    }
                }
            }

            CCString sLevels;
            for (INT iLevel = 0; iLevel < iDeltaLevels; ++iLevel)
            {
                CCString sThis;
                sThis.Format(_T(" h=%.4e fd=%.6e err=%.6e bound=%.6e%s"),
                             fdSteps[iLevel], fdDfd[iLevel], fdAbsErrors[iLevel], fdBounds[iLevel],
                             (fdAbsErrors[iLevel] <= fdBounds[iLevel]) ? _T("*") : _T(""));
                sLevels += sThis;
            }

            const DOUBLE fBestBound = (iBestLevel >= 0) ? fdBounds[iBestLevel] : 0.0;
            appGeneral(_T("FD trial (%s) link=%u dir=%u pred=%.6e%s best=%d minerr=%.6e bound=%.6e %s\n"),
                       szName, uiLinkIndex, byDir, fDforce, sLevels.c_str(),
                       iBestLevel, (iBestLevel >= 0) ? fdAbsErrors[iBestLevel] : -1.0,
                       fBestBound,
                       bTrialPassed ? _T("PASS") : _T("FAIL"));

            if (bTrialPassed)
            {
                ++uiPassedTrials;
            }
            else
            {
                ++uiFailedTrials;
                CCString sProblem;
                sProblem.Format(_T("Finite-difference (%s) link=%u dir=%u failed: no in-tolerance level\n"),
                                szName, uiLinkIndex, byDir);
                LastProbem(sProblem);
            }

            pU0->CopyTo(pGauge);
            RefreshGaugeDependentState(pGauge);
            appSynchronize();

            appSafeDelete(pX);
        }
    }

    appGeneral(_T("Finite-difference force validation (%s): trials=%u passed=%u failed=%u\n"),
               szName, uiTotalTrials, uiPassedTrials, uiFailedTrials);

    return uiError + uiFailedTrials;
}

#pragma endregion

#pragma region top-level driver and registration

static void BuildSelectedActionList(CParameters& sParam, TArray<CAction*>& actions, TArray<INT>& oneBasedIndices)
{
    actions.RemoveAll();
    oneBasedIndices.RemoveAll();

    TArray<UINT> actionIndices;
    const UBOOL bHasActionIndices = sParam.FetchValueArrayUINT(_T("ActionIndices"), actionIndices);

    const INT iActionCount = appGetLattice()->m_pActionList.Num();
    for (INT i = 0; i < iActionCount; ++i)
    {
        CAction* pAction = appGetLattice()->m_pActionList[i];
        if (NULL == pAction || pAction->IsDiscreteGauge())
        {
            continue;
        }

        const INT iOneBased = i + 1;
        if (bHasActionIndices)
        {
            UBOOL bSelected = FALSE;
            for (INT j = 0; j < actionIndices.Num(); ++j)
            {
                if (static_cast<INT>(actionIndices[j]) == iOneBased)
                {
                    bSelected = TRUE;
                    break;
                }
            }
            if (!bSelected)
            {
                continue;
            }
        }

        actions.AddItem(pAction);
        oneBasedIndices.AddItem(iOneBased);
    }
}

UINT TestHmcDiagnostics(CParameters& sParam)
{
    UINT uiError = 0;

    INT iRunReversibility = 1;
    INT iRunFiniteDiff = 1;
    sParam.FetchValueINT(_T("RunReversibility"), iRunReversibility);
    sParam.FetchValueINT(_T("RunFiniteDiff"), iRunFiniteDiff);

    TArray<CAction*> selectedActions;
    TArray<INT> selectedOneBased;
    BuildSelectedActionList(sParam, selectedActions, selectedOneBased);

    if (0 == selectedActions.Num())
    {
        LastProbem(_T("TestHmcDiagnostics: no continuous actions selected\n"));
        return 1;
    }

    if (0 != iRunReversibility)
    {
        uiError += TestHmcReversibility(sParam);
    }

    if (0 != iRunFiniteDiff)
    {
        for (INT i = 0; i < selectedActions.Num(); ++i)
        {
            CCString sName;
            sName.Format(_T("Action%d"), selectedOneBased[i]);
            uiError += TestActionForceFiniteDiff(sParam, selectedActions[i], sName);
        }
    }

    return uiError;
}

#pragma region anisotropic acceptance checks

/**
 * Acceptance checks of the anisotropic tree-level Symanzik gauge action and
 * the aHISQ fermion field (task-anisotropy.md section 4.2):
 * 1. gauge isotropic limit (Xi kernel weight == 1 agrees with the isotropic path)
 * 2. rectangle-off limit (RectOverPlaq = 0 agrees with the plaquette-only anisotropy)
 * 4. fermion isotropic limit (XiF = 1 agrees with CFieldFermionHISQSU3)
 * 5. temporal-only weighting on the identity gauge background
 * 7. two flavors with different XiF, the combined force is the sum of the separate forces
 * 8. solver residual and anti-Hermiticity of the massless hopping operator
 */
UINT TestAnisotropicChecks(CParameters& sParam)
{
    UINT uiError = 0;
    DOUBLE fTol = F(0.001);
    sParam.FetchValueDOUBLE(_T("CheckTol"), fTol);

    CFieldGauge* pGauge = appGetLattice()->m_pGaugeField[0];
    CFieldGaugeSU3TreeImproved* pTreeImproved = dynamic_cast<CFieldGaugeSU3TreeImproved*>(pGauge);
    CFieldFermionHISQSU3Anisotropic* pAniFermion = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(appGetLattice()->GetFieldById(2));
    if (NULL == pTreeImproved || NULL == pAniFermion)
    {
        appCrucial(_T("TestAnisotropicChecks: need CFieldGaugeSU3TreeImproved and CFieldFermionHISQSU3Anisotropic (id 2)\n"));
        return 1;
    }
    const DOUBLE fBetaOverN = F(5.0) / F(3.0);

    RefreshGaugeDependentState(pGauge);
    CGaugeSmearing* pSmearing = appGetGaugeSmearing(pGauge->m_byFieldId);
    const CFieldGauge* pGaugeEff = pSmearing->GetEffectiveGauge();

    //================ check 1: gauge isotropic limit ================
    {
        const DOUBLE fEIso = pGauge->CalculatePlaqutteEnergy(fBetaOverN);
        const DOUBLE fEAni = pGauge->CalculatePlaqutteEnergyAnisotropy(fBetaOverN, 1.0);
        const DOUBLE fEIsoC = pGauge->CalculatePlaqutteEnergyUseClover(fBetaOverN);
        const DOUBLE fEAniC = pGauge->CalculatePlaqutteEnergyUseCloverAnisotropy(fBetaOverN, 1.0);

        CFieldGauge* pF1 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pF2 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pStaple = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        pF1->Zero();
        pF2->Zero();
        pGauge->CalculateForceAnisotropy(pF1, fBetaOverN, 1.0);
        pGauge->CalculateForceAndStaple(pF2, pStaple, static_cast<Real>(fBetaOverN));
        pF1->AxpyMinus(pF2);
        //gauge field GetLength subtracts the generator count, use Dot for the force difference
        const DOUBLE fForceDiff = pF1->Dot(pF1).x;
        const DOUBLE fForceNorm = pF2->Dot(pF2).x;

        appGeneral(_T("check1 gauge iso limit: dE=%2.6e dEClover=%2.6e dForce=%2.6e (|F|=%2.6e)\n"),
            abs(fEIso - fEAni), abs(fEIsoC - fEAniC), fForceDiff, fForceNorm);
        if (abs(fEIso - fEAni) > fTol * abs(fEIso) || abs(fEIsoC - fEAniC) > fTol * abs(fEIsoC)
         || fForceDiff > fTol * (fForceNorm + F(1.0)))
        {
            LastProbem(_T("check1 gauge isotropic limit failed"));
            ++uiError;
        }
        appSafeDelete(pF1);
        appSafeDelete(pF2);
        appSafeDelete(pStaple);
    }

    //================ check 1b: one-loop isotropic limit ================
    {
        CFieldGaugeSU3OneLoopImproved* pOneLoop = dynamic_cast<CFieldGaugeSU3OneLoopImproved*>(appCreate(_T("CFieldGaugeSU3OneLoopImproved")));
        if (NULL == pOneLoop)
        {
            appCrucial(_T("TestAnisotropicChecks: can not create CFieldGaugeSU3OneLoopImproved\n"));
            ++uiError;
        }
        else
        {
            pTreeImproved->CopyParamTo(pOneLoop);
            pOneLoop->m_byFieldId = pGauge->m_byFieldId;
            pOneLoop->m_pOwner = appGetLattice();
            pGauge->CopyBufferTo(pOneLoop);

            const DOUBLE fEIso = pOneLoop->CalculatePlaqutteEnergy(fBetaOverN);
            const DOUBLE fEAni = pOneLoop->CalculatePlaqutteEnergyAnisotropy(fBetaOverN, 1.0);

            CFieldGauge* pF1 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            CFieldGauge* pF2 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            CFieldGauge* pStaple = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            pF1->Zero();
            pF2->Zero();
            pOneLoop->CalculateForceAnisotropy(pF1, fBetaOverN, 1.0);
            pOneLoop->CalculateForceAndStaple(pF2, pStaple, static_cast<Real>(fBetaOverN));
            pF1->AxpyMinus(pF2);

            appGeneral(_T("check1b one-loop iso limit: dE=%2.6e dForce=%2.6e\n"), abs(fEIso - fEAni), pF1->Dot(pF1).x);
            if (abs(fEIso - fEAni) > fTol * abs(fEIso) || pF1->Dot(pF1).x > fTol * (pF2->Dot(pF2).x + F(1.0)))
            {
                LastProbem(_T("check1b one-loop isotropic limit failed"));
                ++uiError;
            }
            appSafeDelete(pF1);
            appSafeDelete(pF2);
            appSafeDelete(pStaple);
        }
        appSafeDelete(pOneLoop);
    }

    //================ check 2: rectangle-off limit ================
    {
        const DOUBLE fSaveRect = pTreeImproved->m_fRectOverPlaq;
        pTreeImproved->m_fRectOverPlaq = 0.0;
        const DOUBLE fKernelXi = 1.0 / 1.7;
        const DOUBLE fEFull = pGauge->CalculatePlaqutteEnergyAnisotropy(fBetaOverN, fKernelXi);
        //qualified call bypasses the virtual dispatch, this is the plaquette-only anisotropy
        const DOUBLE fEPlaqOnly = pTreeImproved->CFieldGaugeSU3::CalculatePlaqutteEnergyAnisotropy(fBetaOverN, fKernelXi);
        const DOUBLE fEFullC = pGauge->CalculatePlaqutteEnergyUseCloverAnisotropy(fBetaOverN, fKernelXi);
        const DOUBLE fEPlaqOnlyC = pTreeImproved->CFieldGaugeSU3::CalculatePlaqutteEnergyUseCloverAnisotropy(fBetaOverN, fKernelXi);
        pTreeImproved->m_fRectOverPlaq = fSaveRect;

        appGeneral(_T("check2 rectangle-off: dE=%2.6e dEClover=%2.6e\n"), abs(fEFull - fEPlaqOnly), abs(fEFullC - fEPlaqOnlyC));
        if (abs(fEFull - fEPlaqOnly) > fTol * abs(fEFull) || abs(fEFullC - fEPlaqOnlyC) > fTol * abs(fEFullC))
        {
            LastProbem(_T("check2 rectangle-off limit failed"));
            ++uiError;
        }
    }

    //================ check 4: fermion isotropic limit ================
    {
        const Real fSaveXiF = pAniFermion->m_fXiF;
        pAniFermion->m_fXiF = F(1.0);

        CFieldFermionHISQSU3* pRef = dynamic_cast<CFieldFermionHISQSU3*>(appCreate(_T("CFieldFermionHISQSU3")));
        if (NULL == pRef)
        {
            appCrucial(_T("TestAnisotropicChecks: can not create CFieldFermionHISQSU3\n"));
            return uiError + 1;
        }
        pAniFermion->CopyParamTo(pRef);
        pRef->m_byFieldId = pAniFermion->m_byFieldId;
        pRef->m_pOwner = appGetLattice();
        pAniFermion->CopyBufferTo(pRef);

        const CFieldGauge* gauges[1] = { pGaugeEff };
        CFieldFermion* pDAni = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        CFieldFermion* pDRef = dynamic_cast<CFieldFermion*>(pRef->GetCopy());
        pDAni->D(1, 0, 0, gauges, NULL, NULL);
        pDRef->D(1, 0, 0, gauges, NULL, NULL);
        pDAni->AxpyMinus(pDRef);
        const DOUBLE fDDiff = pDAni->GetLength();
        const DOUBLE fDNorm = _sqrt(pDRef->GetLength());

        CFieldFermion* pDdAni = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        CFieldFermion* pDdRef = dynamic_cast<CFieldFermion*>(pRef->GetCopy());
        pDdAni->Ddagger(1, 0, 0, gauges, NULL, NULL);
        pDdRef->Ddagger(1, 0, 0, gauges, NULL, NULL);
        pDdAni->AxpyMinus(pDdRef);
        const DOUBLE fDdDiff = pDdAni->GetLength();

        CFieldFermion* pDDdAni = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        CFieldFermion* pDDdRef = dynamic_cast<CFieldFermion*>(pRef->GetCopy());
        pDDdAni->ApplyOperator(EFO_F_DDdagger, 1, 0, 0, gauges, NULL, NULL);
        pDDdRef->ApplyOperator(EFO_F_DDdagger, 1, 0, 0, gauges, NULL, NULL);
        pDDdAni->AxpyMinus(pDDdRef);
        const DOUBLE fDDdDiff = pDDdAni->GetLength();

        //RHMC energy with the same pseudofermion
        pAniFermion->PrepareForHMC(1, 0, gauges, NULL);
        pAniFermion->CopyBufferTo(pRef);
        const DOUBLE fSEAni = pAniFermion->Energy(1, 0, 0, gauges, NULL, NULL);
        const DOUBLE fSERef = pRef->Energy(1, 0, 0, gauges, NULL, NULL);

        //force comparison on the same pseudofermion
        CFieldGauge* pF0Ani = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pF0Ref = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pEpsAni = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pEpsRef = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pNaikAni = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        CFieldGauge* pNaikRef = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        pF0Ani->Zero();
        pF0Ref->Zero();
        pEpsAni->Zero();
        pEpsRef->Zero();
        pNaikAni->Zero();
        pNaikRef->Zero();
        pAniFermion->CalculateF0AndNaik(pGaugeEff, pF0Ani, pEpsAni, pNaikAni);
        pRef->CalculateF0AndNaik(pGaugeEff, pF0Ref, pEpsRef, pNaikRef);
        pF0Ani->AxpyMinus(pF0Ref);
        pEpsAni->AxpyMinus(pEpsRef);
        pNaikAni->AxpyMinus(pNaikRef);

        appGeneral(_T("check4 fermion iso limit: dD=%2.6e dDdagger=%2.6e dDDdagger=%2.6e dE=%2.6e dF0=%2.6e dEps=%2.6e dNaik=%2.6e\n"),
            fDDiff, fDdDiff, fDDdDiff, abs(fSEAni - fSERef), pF0Ani->Dot(pF0Ani).x, pEpsAni->Dot(pEpsAni).x, pNaikAni->Dot(pNaikAni).x);
        if (fDDiff > fTol * (fDNorm + F(1.0)) || fDdDiff > fTol * (fDNorm + F(1.0)) || fDDdDiff > fTol * (fDNorm + F(1.0))
         || abs(fSEAni - fSERef) > fTol * (abs(fSERef) + F(1.0))
         || pF0Ani->Dot(pF0Ani).x > fTol * (pF0Ref->Dot(pF0Ref).x + F(1.0))
         || pEpsAni->Dot(pEpsAni).x > fTol * (pEpsRef->Dot(pEpsRef).x + F(1.0))
         || pNaikAni->Dot(pNaikAni).x > fTol * (pNaikRef->Dot(pNaikRef).x + F(1.0)))
        {
            LastProbem(_T("check4 fermion isotropic limit failed"));
            ++uiError;
        }

        appSafeDelete(pDAni);
        appSafeDelete(pDRef);
        appSafeDelete(pDdAni);
        appSafeDelete(pDdRef);
        appSafeDelete(pDDdAni);
        appSafeDelete(pDDdRef);
        appSafeDelete(pF0Ani);
        appSafeDelete(pF0Ref);
        appSafeDelete(pEpsAni);
        appSafeDelete(pEpsRef);
        appSafeDelete(pNaikAni);
        appSafeDelete(pNaikRef);
        appSafeDelete(pRef);
        pAniFermion->m_fXiF = fSaveXiF;
    }

    //================ check 5: temporal-only weighting on identity gauge ================
    {
        CFieldGauge* pGaugeBackup = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        pGauge->InitialField(EFIT_Identity);
        RefreshGaugeDependentState(pGauge);

        CFieldFermionHISQSU3Anisotropic* pPhi = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(pAniFermion->GetCopy());
        const Real fSaveXiF = pAniFermion->m_fXiF;
        const Real fX1 = fSaveXiF;

        //point source at the origin, on the identity gauge every value is known exactly
        BYTE* pSourceData = (BYTE*)malloc(sizeof(Real) * 6 * _HC_Volume);
        memset(pSourceData, 0, sizeof(Real) * 6 * _HC_Volume);
        ((Real*)pSourceData)[0] = F(1.0);
        pPhi->InitialWithByte(pSourceData);
        free(pSourceData);

        //the operator runs on the copies, set the anisotropy on the copies
        CFieldFermionHISQSU3Anisotropic* pD1 = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(pPhi->GetCopy());
        pD1->m_fXiF = F(1.0);
        pD1->D(1, 0, 0, (const CFieldGauge* const*)&pGaugeEff, NULL, NULL);

        CFieldFermionHISQSU3Anisotropic* pDX = dynamic_cast<CFieldFermionHISQSU3Anisotropic*>(pPhi->GetCopy());
        pDX->m_fXiF = fX1;
        pDX->D(1, 0, 0, (const CFieldGauge* const*)&pGaugeEff, NULL, NULL);

        pDX->AxpyMinus(pD1);

        //manual temporal hopping on the identity gauge, (X1 - 1) * T(phi)
        //measure the link values on the identity background: X^F, level-1 and X^L
        UINT uiGEffSize = 0;
        BYTE* pGEffData = pGaugeEff->CopyDataOut(uiGEffSize);
        const Real fEffLink = ((const Real*)pGEffData)[0];
        free(pGEffData);
        UINT uiGL1Size = 0;
        BYTE* pGL1Data = pSmearing->GetEffectiveGaugeLevel1()->CopyDataOut(uiGL1Size);
        const Real fL1Link = ((const Real*)pGL1Data)[0];
        free(pGL1Data);
        UINT uiNaikSize = 0;
        BYTE* pNaikData = pSmearing->GetNaikLink()->CopyDataOut(uiNaikSize);
        const Real fNaikLink = ((const Real*)pNaikData)[0];
        free(pNaikData);

        UINT uiSize = 0;
        BYTE* pPhiData = pPhi->CopyDataOut(uiSize);
        const UINT uiSiteCount = uiSize / sizeof(Real) / 6;
        const UINT uiLt = _HC_Lt;
        const UINT uiLz = _HC_Lz;
        const UINT uiLy = _HC_Ly;
        const Real fNaikCoef = pAniFermion->m_fNaik;
        const Real fEpsCoef = pAniFermion->m_fEpsilon * F(0.125);
        const Real fOneLinkCoef = fEffLink + fEpsCoef * fL1Link;

        //antiperiodic boundary in the temporal direction: hops wrapping the boundary change sign
        const SSmallInt4 sBC = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(pPhi->m_byFieldId);
        const UBOOL bAntiT = (sBC.w == -1);

        UINT uiDiffSize = 0;
        BYTE* pDiffData = pDX->CopyDataOut(uiDiffSize);
        DOUBLE fDiffNorm2 = 0.0;
        DOUBLE fErrNorm2 = 0.0;
        for (UINT idx = 0; idx < uiSiteCount; ++idx)
        {
            const UINT x = idx / (uiLy * uiLz * uiLt);
            const UINT y = (idx / (uiLz * uiLt)) % uiLy;
            const UINT z = (idx / uiLt) % uiLz;
            const UINT t = idx % uiLt;
            const Real fEta4 = (0 == ((x + y + z) & 1)) ? F(1.0) : F(-1.0);
            const UINT idxTp = x * uiLy * uiLz * uiLt + y * uiLz * uiLt + z * uiLt + (t + 1) % uiLt;
            const UINT idxTm = x * uiLy * uiLz * uiLt + y * uiLz * uiLt + z * uiLt + (t + uiLt - 1) % uiLt;
            const UINT idxTp3 = x * uiLy * uiLz * uiLt + y * uiLz * uiLt + z * uiLt + (t + 3) % uiLt;
            const UINT idxTm3 = x * uiLy * uiLz * uiLt + y * uiLz * uiLt + z * uiLt + (t + uiLt - 3) % uiLt;
            const Real fSignP = (bAntiT && t + 1 >= uiLt) ? F(-1.0) : F(1.0);
            const Real fSignM = (bAntiT && t < 1) ? F(-1.0) : F(1.0);
            const Real fSignP3 = (bAntiT && t + 3 >= uiLt) ? F(-1.0) : F(1.0);
            const Real fSignM3 = (bAntiT && t < 3) ? F(-1.0) : F(1.0);
            const Real* phiX = (const Real*)pPhiData;
            const Real* diffX = (const Real*)pDiffData;
            for (UINT c = 0; c < 6; ++c)
            {
                const Real fT = fEta4 * (fOneLinkCoef * (fSignP * phiX[idxTp * 6 + c] - fSignM * phiX[idxTm * 6 + c])
                    + fNaikCoef * fNaikLink * (fSignP3 * phiX[idxTp3 * 6 + c] - fSignM3 * phiX[idxTm3 * 6 + c]));
                const Real fExpect = (fX1 - F(1.0)) * fT;
                fErrNorm2 += (diffX[idx * 6 + c] - fExpect) * (diffX[idx * 6 + c] - fExpect);
                fDiffNorm2 += diffX[idx * 6 + c] * diffX[idx * 6 + c];
            }
        }

        appGeneral(_T("check5 temporal-only weighting: err=%2.6e |diff|=%2.6e (links: eff=%f l1=%f naik=%f, antiT=%d)\n"),
            _sqrt(fErrNorm2), _sqrt(fDiffNorm2), fEffLink, fL1Link, fNaikLink, bAntiT ? 1 : 0);
        free(pPhiData);
        free(pDiffData);
        if (_sqrt(fErrNorm2) > fTol * (_sqrt(fDiffNorm2) + F(1.0)))
        {
            LastProbem(_T("check5 temporal-only weighting failed"));
            ++uiError;
        }

        appSafeDelete(pPhi);
        appSafeDelete(pD1);
        appSafeDelete(pDX);
        pGaugeBackup->CopyTo(pGauge);
        RefreshGaugeDependentState(pGauge);
        appSafeDelete(pGaugeBackup);
    }

    //================ check 7: two flavors force additivity ================
    {
        TArray<CAction*> fermionActions;
        for (INT i = 0; i < appGetLattice()->m_pActionList.Num(); ++i)
        {
            CAction* pAction = appGetLattice()->m_pActionList[i];
            if (NULL != dynamic_cast<CActionFermionHISQCombined*>(pAction))
            {
                fermionActions.AddItem(pAction);
            }
        }
        if (fermionActions.Num() != 3)
        {
            appCrucial(_T("TestAnisotropicChecks: need 3 CActionFermionHISQCombined ([2], [3], [2,3]), found %d\n"), fermionActions.Num());
            ++uiError;
        }
        else
        {
            const CFieldGauge* gauges[1] = { pGauge };
            //prepare the pseudofermions once per field, NOT by the combined action again
            fermionActions[0]->PrepareForHMC(1, 0, gauges, NULL, 0);
            fermionActions[1]->PrepareForHMC(1, 0, gauges, NULL, 0);

            CFieldGauge* pF2 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            CFieldGauge* pF3 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            CFieldGauge* pF23 = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
            pF2->Zero();
            pF3->Zero();
            pF23->Zero();
            CFieldGauge* force2[1] = { pF2 };
            CFieldGauge* force3[1] = { pF3 };
            CFieldGauge* force23[1] = { pF23 };
            fermionActions[0]->CalculateForce(1, 0, gauges, NULL, force2, NULL, NULL, ESP_Once);
            fermionActions[1]->CalculateForce(1, 0, gauges, NULL, force3, NULL, NULL, ESP_Once);
            fermionActions[2]->CalculateForce(1, 0, gauges, NULL, force23, NULL, NULL, ESP_Once);

            pF23->AxpyMinus(pF2);
            pF23->AxpyMinus(pF3);
            const DOUBLE fDiff = pF23->Dot(pF23).x;
            const DOUBLE fNorm = pF2->Dot(pF2).x;
            appGeneral(_T("check7 two-flavor additivity: |F23-F2-F3|=%2.6e (|F2|=%2.6e)\n"), fDiff, fNorm);
            if (fDiff > fTol * (fNorm + F(1.0)))
            {
                LastProbem(_T("check7 two-flavor additivity failed"));
                ++uiError;
            }
            appSafeDelete(pF2);
            appSafeDelete(pF3);
            appSafeDelete(pF23);
        }
    }

    //================ check 8: solver sanity and anti-Hermiticity ================
    {
        const CFieldGauge* gauges[1] = { pGaugeEff };
        CFieldFermion* pPhi = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        pPhi->InitialField(EFIT_RandomGaussian);
        CFieldFermion* pX = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        pX->Zero();

        CSLASolver* pSolver = appGetFermionSolver(pAniFermion->m_byFieldId);
        if (NULL == pSolver)
        {
            appCrucial(_T("TestAnisotropicChecks: no solver for field 2\n"));
            ++uiError;
        }
        else
        {
            pSolver->Solve(pX, pPhi, 1, 0, 0, gauges, NULL, NULL, EFO_F_DDdagger, ESP_Once, NULL);
            CFieldFermion* pR = dynamic_cast<CFieldFermion*>(pX->GetCopy());
            pR->ApplyOperator(EFO_F_DDdagger, 1, 0, 0, gauges, NULL, NULL);
            pR->AxpyMinus(pPhi);
            const DOUBLE fRes = _sqrt(pR->GetLength()) / _sqrt(pPhi->GetLength());
            appGeneral(_T("check8 solver true residual: %2.6e\n"), fRes);
            if (fRes > F(0.0001))
            {
                LastProbem(_T("check8 solver residual failed"));
                ++uiError;
            }
            appSafeDelete(pR);
        }

        //anti-Hermiticity of the massless hopping operator: <x, D0 y> = -<D0 x, y>
        const Real fSaveMass = pAniFermion->m_f2am;
        pAniFermion->m_f2am = F(0.0);
        CFieldFermion* pXf = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        pXf->InitialField(EFIT_RandomGaussian);
        CFieldFermion* pYf = dynamic_cast<CFieldFermion*>(pAniFermion->GetCopy());
        pYf->InitialField(EFIT_RandomGaussian);
        CFieldFermion* pDY = dynamic_cast<CFieldFermion*>(pYf->GetCopy());
        pDY->D(1, 0, 0, gauges, NULL, NULL);
        CFieldFermion* pDX = dynamic_cast<CFieldFermion*>(pXf->GetCopy());
        pDX->D(1, 0, 0, gauges, NULL, NULL);
        const CLGComplex c1 = pXf->DotReal(pDY);
        const CLGComplex c2 = pDX->DotReal(pYf);
        const DOUBLE fAntiSym = _cuCabsf(_cuCaddf(c1, c2));
        const DOUBLE fScale = _cuCabsf(c1);
        pAniFermion->m_f2am = fSaveMass;
        appGeneral(_T("check8 anti-Hermiticity: |<x,D0y>+<D0x,y>|=%2.6e (|<x,D0y>|=%2.6e)\n"), fAntiSym, fScale);
        if (fAntiSym > fTol * (fScale + F(1.0)))
        {
            LastProbem(_T("check8 anti-Hermiticity failed"));
            ++uiError;
        }

        appSafeDelete(pPhi);
        appSafeDelete(pX);
        appSafeDelete(pXf);
        appSafeDelete(pYf);
        appSafeDelete(pDY);
        appSafeDelete(pDX);
    }

    return uiError;
}

__REGIST_TEST(TestAnisotropicChecks, HmcDiagnostics, TestHmcDiagnostics_AnisotropicChecks, AnisotropicChecks);

#pragma endregion


__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_SU3Gauge_KS, SU3Gauge_KS);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_TreeImproved_HISQ, TreeImproved_HISQ);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_TreeImproved_AnisotropicHISQ, TreeImproved_AnisotropicHISQ);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_OneLoop_AnisotropicHISQ, OneLoop_AnisotropicHISQ);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_RotatingHISQEvenOdd, RotatingHISQEvenOdd);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_TreeImproved_StoutCloverWilson, TreeImproved_StoutCloverWilson);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_RigidAcc, RigidAcc);
__REGIST_TEST(TestHmcDiagnostics, HmcDiagnostics, TestHmcDiagnostics_PSU3, PSU3);

#pragma endregion

//=============================================================================
// END OF FILE
//=============================================================================
