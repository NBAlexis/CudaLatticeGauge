//=============================================================================
// FILENAME : TestRotation.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [05/10/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestRotation(CParameters& sParam)
{
    Real fExpected = F(0.65);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedResDebug"), fExpected);
#endif

    //we calculate staple energy from beta = 1 - 6
    //CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    //if (NULL == pAction)
    //{
    //    return 1;
    //}
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasurePolyakovXY3D* pPolyakov3D = NULL;
    if (NULL == pMeasure)
    {
        pPolyakov3D = dynamic_cast<CMeasurePolyakovXY3D*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
        if (NULL == pPolyakov3D)
        {
            return 1;
        }
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(3, FALSE);
#endif

    //Measure
    if (NULL != pMeasure)
    {
        pMeasure->Reset();
    }
    else
    {
        pPolyakov3D->Reset();
    }
    
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(4, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(12, TRUE);
#endif

    const Real fRes = NULL == pMeasure ? pPolyakov3D->m_fLastRealResult : pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.15))
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = static_cast<Real>(appGetLattice()->m_pUpdator->GetHDiff());
    const Real fLastHDiff = appGetLattice()->m_pUpdator->GetLastHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/5) : expected >= 4. HDiff = %f : expected < 1. last HDiff = %f : expected < 0.15\n"), uiAccept, fHDiff, fLastHDiff);
#else
    appGeneral(_T("accept (%d/15) : expected >= 12. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%),  last HDiff = %f : expected < 0.15\n"), uiAccept, fHDiff, fLastHDiff);
#endif

#if _CLG_DEBUG
    if (uiAccept < 4)
#else
    if (uiAccept < 12)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fLastHDiff > F(0.15))
#else
    if (fLastHDiff > F(0.08))
#endif
    {
        ++uiError;
    }

    return uiError;
}

UINT TestAcceleration(CParameters& sParam)
{
    Real fExpected = F(0.36);

    INT iVaule = 2;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 10;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    //we calculate staple energy from beta = 1 - 6
    //CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    //if (NULL == pAction)
    //{
    //    return 1;
    //}
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(4, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/4) : expected >= 3. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 3)
#else
    if (uiAccept < 23)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}

UINT TestRigidAcceleration(CParameters& sParam)
{
    Real fExpected = F(0.55);

    INT iVaule = 2;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 10;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    //we calculate staple energy from beta = 1 - 6
    //CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    //if (NULL == pAction)
    //{
    //    return 1;
    //}
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(3, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/4) : expected >= 3. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 3)
#else
    if (uiAccept < 23)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}

UINT TestBoost(CParameters& sParam)
{
    Real fExpected = F(0.38);

    INT iVaule = 2;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 10;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    //we calculate staple energy from beta = 1 - 6
    //CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    //if (NULL == pAction)
    //{
    //    return 1;
    //}
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(4, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/4) : expected >= 3. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 3)
#else
    if (uiAccept < 23)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}

UINT TestRotationKS(CParameters& sParam)
{
    Real fExpected = F(0.65);

    INT iVaule = 3;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 15;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedResDebug"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(3, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/5) : expected >= 4. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/18) : expected >= 15. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 4)
#else
    if (uiAccept < 15)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}

UINT TestEMSimple(CParameters& sParam)
{
    Real fExpected = F(0.38);

    INT iVaule = 2;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 10;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    //we calculate staple energy from beta = 1 - 6
    //CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    //if (NULL == pAction)
    //{
    //    return 1;
    //}
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(4, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/4) : expected >= 3. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 3)
#else
    if (uiAccept < 23)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}


UINT TestRotationKSU1(CParameters& sParam)
{
    Real fExpected = F(0.65);

    INT iVaule = 3;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 15;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedResDebug"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(3, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/5) : expected >= 4. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/18) : expected >= 15. HDiff = %f : expected < 0.5 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 4)
#else
    if (uiAccept < 15)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.5))
#endif
    {
        ++uiError;
    }

    return uiError;
}


UINT TestRotationEMProjectivePlane(CParameters& sParam)
{
    Real fExpected = F(0.38);

    INT iVaule = 2;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 10;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(4, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.15))
#else
    if (appAbs(fRes - fExpected) > F(0.02))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/4) : expected >= 3. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 3)
#else
    if (uiAccept < 23)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}

UINT TestBetaGradient(CParameters& sParam)
{
    Real fExpected = F(0.048);

    INT iVaule = 2;
    sParam.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 10;
    sParam.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iAfterEquib = static_cast<UINT>(iVaule);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif

    CMeasurePolyakovXY* pMeasure = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(1, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(iBeforeEquib, FALSE);
#endif

    //Measure
    pMeasure->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(4, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(iAfterEquib, TRUE);
#endif

    const Real fRes = _cuCabsf(pMeasure->m_cAverageLoop);
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
#if _CLG_DEBUG
    if (appAbs(fRes - fExpected) > F(0.005))
#else
    if (appAbs(fRes - fExpected) > F(0.005))
#endif
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/4) : expected >= 3. HDiff = %f : expected < 1\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#else
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
#endif

#if _CLG_DEBUG
    if (uiAccept < 3)
#else
    if (uiAccept < 23)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(1.0))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestRotation, Rotation, TestRotationWilsonDiracDirichlet);
__REGIST_TEST(TestRotation, Rotation, TestRotationWilsonDiracTorus);
__REGIST_TEST(TestRotation, Rotation, TestRotationWilsonDiracProjectivePlane);

__REGIST_TEST(TestRotation, Rotation, TestRotationKSDirichlet);
__REGIST_TEST(TestRotation, Rotation, TestRotationKSTorus);
__REGIST_TEST(TestRotation, Rotation, TestRotationKSProjectivePlane);

__REGIST_TEST(TestRotation, Rotation, TestRotationEMProjectivePlane);

__REGIST_TEST(TestRotation, Rotation, TestRotationProjectivePlaneU1);

__REGIST_TEST(TestRotation, Rotation, TestRotationQuenched3D);
__REGIST_TEST(TestRotation, Rotation, TestRotationQuenchedU13D);

__REGIST_TEST(TestAcceleration, Updator, TestAcceleration);

//This has sign problem..
__REGIST_TEST(TestRigidAcceleration, Updator, TestRigidAcceleration);

__REGIST_TEST(TestBoost, Updator, TestBoost);

__REGIST_TEST(TestEMSimple, Updator, TestEMSimple);

//This is well tested
//__REGIST_TEST(TestRotationKS, Updator, TestRotationProjectivePlane);

//Special cases
__REGIST_TEST(TestBetaGradient, Updator, TestBetaGradient);

//=============================================================================
// END OF FILE
//=============================================================================
