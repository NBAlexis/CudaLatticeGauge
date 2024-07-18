//=============================================================================
// FILENAME : CSolverBiCGstab.cpp
// 
// DESCRIPTION:
// This is the class for BiCGStab solver
//
// REVISION:
//  [01/08/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CSolverBiCGstab.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverBiCGStab)

CSLASolverBiCGStab::CSLASolverBiCGStab() 
    : CSLASolver()
    , m_uiReTry(3)
    , m_uiDevationCheck(10)
    , m_uiStepCount(20)
    , m_fAccuracy(F(0.000001))
    , m_fSmallRho(F(0.00000001))
{

}

CSLASolverBiCGStab::~CSLASolverBiCGStab()
{
    CSLASolverBiCGStab::ReleaseBuffers();
}

void CSLASolverBiCGStab::Configurate(const CParameters& param)
{
    INT iValue;
    //Real fValue;
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
#if _CLG_DOUBLEFLOAT
    Real fValue = F(0.0);
    if (param.FetchValueReal(_T("Accuracy"), fValue))
    {
        m_fAccuracy = fValue;
#else
    DOUBLE dValue = 0.0;
    if (param.FetchValueDOUBLE(_T("Accuracy"), dValue))
    {
        m_fAccuracy = dValue;
#endif
        //if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        //{
        //    m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
        //    appGeneral(_T("Solver accuracy too small (%2.18f), set to be %2.18f\n"), fValue, m_fAccuracy);
        //}
    }

#if _CLG_DOUBLEFLOAT
    fValue = F(0.0);
    if (param.FetchValueReal(_T("SmallRho"), fValue))
    {
        m_fSmallRho = fValue;
#else
    dValue = 0.0;
    if (param.FetchValueDOUBLE(_T("SmallRho"), dValue))
    {
        m_fSmallRho = dValue;
#endif
        //if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        //{
        //    m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
        //    appGeneral(_T("Solver accuracy too small (%2.18f), set to be %2.18f\n"), fValue, m_fAccuracy);
        //}
    }

    //if (m_uiReTry < 2)
    //{
    //    m_fSmallRho = _CLG_FLT_MIN_;
    //}
}

void CSLASolverBiCGStab::AllocateBuffers(const CField* )
{

}

void CSLASolverBiCGStab::ReleaseBuffers()
{

}

#if 0
UBOOL CSLASolverBiCGStab::Solve1(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart) const
{
#if !_CLG_DOUBLEFLOAT
    //CField* pB = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pP = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pV = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pRh = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pS = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pT = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
    DOUBLE fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = cuCabs(pFieldB->Dot(pFieldB));
    }

    appParanoiac(_T("-- CSLASolverBiCGStab::Solve start operator: %s--\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str());

    //pFieldB->CopyTo(pB);

    //Using b as the guess, (Assuming A is near identity?)
    //r_0 = b - A x_0
    if (NULL == pStart)
    {
        pFieldB->CopyTo(pR);
        pFieldB->CopyTo(pX);
    }
    else
    {
        pStart->CopyTo(pR);
        pStart->CopyTo(pX);
    }
    pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //-A x_0
    pR->AxpyPlus(pFieldB); //pR->AxpyPlus(pB); //b - A x_0
    pR->CopyTo(pRh);
    //By Yousef Saad, we use conjugated one.
    pRh->Dagger();

    cuDoubleComplex rho = make_cuDoubleComplex(F(0.0), F(0.0));
    cuDoubleComplex last_rho = make_cuDoubleComplex(F(0.0), F(0.0));
    cuDoubleComplex alpha = make_cuDoubleComplex(F(0.0), F(0.0));
    cuDoubleComplex beta = make_cuDoubleComplex(F(0.0), F(0.0));
    cuDoubleComplex omega = make_cuDoubleComplex(F(0.0), F(0.0));

    for (UINT i = 0; i < m_uiReTry; ++i)
    {
        for (UINT j = 0; j < m_uiStepCount * m_uiDevationCheck; ++j)
        {
            //==========
            //One step
            rho = pRh->Dot(pR);//rho = rh dot r(i-1), if rho = 0, failed (assume will not)
            if (__cuCabsSqd(rho) < _CLG_FLT_MIN)
            {
                appParanoiac(_T("CSLASolverBiCGStab::rho too small:%0.18f, i: %d\n"), __cuCabsSqd(rho), i);
                if (i < m_uiReTry - 1)
                {
                    i = m_uiReTry - 1;
                }
                break;
            }

            if (0 == j) //if is the first iteration, p=r(i-1)
            {
                pR->CopyTo(pP);
            }
            else //if not the first iteration, 
            {
                //beta = last_alpha * rho /(last_omega * last_rho)
                beta = cuCdiv(cuCmul(alpha, rho), cuCmul(omega, last_rho));
                //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
                pV->ScalarMultply(_cToFloat(omega));
                pP->AxpyMinus(pV);
                pP->ScalarMultply(_cToFloat(beta));
                pP->AxpyPlus(pR);
            }

            //v(i) = A p(i)
            pP->CopyTo(pV);
            pV->ApplyOperator(uiM, pGaugeFeild);
            alpha = cuCdiv(rho, pRh->Dot(pV));//alpha = rho / (rh dot v(i))
            //s=r(i-1) - alpha v(i)
            pR->CopyTo(pS);
            pS->Axpy(make_cuComplex(-static_cast<Real>(alpha.x), -static_cast<Real>(alpha.y)), pV);

            if (0 == (j + 1) % m_uiDevationCheck)
            {
                //Normal of S is small, then stop
                const DOUBLE fDeviation = pS->Dot(pS).x / fBLength;
                appParanoiac(_T("CSLASolverBiCGStab::Solve deviation: restart:%d, iteration:%d, deviation:%8.18f\n"), i, j, fDeviation);
                if (fDeviation < m_fAccuracy)
                {
                    //pX->Axpy(alpha, pP); //This is tested a better result not to do the final step
                    pX->CopyTo(pFieldX);

                    //pB->Return();
                    pX->Return();
                    pP->Return();
                    pV->Return();
                    pR->Return();
                    pRh->Return();
                    pS->Return();
                    pT->Return();
                    return TRUE;
                }
            }

            //t=As
            pS->CopyTo(pT);
            pT->ApplyOperator(uiM, pGaugeFeild);
            omega = cuCdivf_cd_host(pS->Dot(pT), pT->Dot(pT).x);//omega = ts / tt

            //r(i)=s-omega t
            pS->CopyTo(pR);
            pR->Axpy(_make_cuComplex(-static_cast<Real>(omega.x), -static_cast<Real>(omega.y)), pT);

            //x(i)=x(i-1) + alpha p + omega s
            pX->Axpy(_cToFloat(alpha), pP);
            pX->Axpy(_cToFloat(omega), pS);

            last_rho = rho;//last_rho = rho
        }

        //we are here, means we do not converge.
        //we need to restart with a new guess, we use last X
        pX->CopyTo(pR);

        pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //-A x_0
        pR->AxpyPlus(pFieldB);  //pR->AxpyPlus(pB); //b - A x_0
        pR->CopyTo(pRh);
    }

    //The solver failed.
    appGeneral(_T("CSLASolverBiCGStab fail to solve!"));

    pX->CopyTo(pFieldX);
    //pB->Return();
    pX->Return();
    pP->Return();
    pV->Return();
    pR->Return();
    pRh->Return();
    pS->Return();
    pT->Return();
    return FALSE;
#else
    //CField* pB = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pP = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pV = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pRh = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pS = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pT = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = _cuCabsf(pFieldB->Dot(pFieldB));
    }

    appParanoiac(_T("-- CSLASolverBiCGStab::Solve start operator: %s--\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str());

    //pFieldB->CopyTo(pB);

    //Using b as the guess, (Assuming A is near identity?)
    //r_0 = b - A x_0
    if (NULL == pStart)
    {
        pFieldB->CopyTo(pR);
        pFieldB->CopyTo(pX);
    }
    else
    {
        pStart->CopyTo(pR);
        pStart->CopyTo(pX);
    }
    pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //-A x_0
    pR->AxpyPlus(pFieldB); //pR->AxpyPlus(pB); //b - A x_0
    pR->CopyTo(pRh);
    //By Yousef Saad, we use conjugated one.
    pRh->Dagger();

    CLGComplex rho = _make_cuComplex(F(0.0), F(0.0));
    CLGComplex last_rho = _make_cuComplex(F(0.0), F(0.0));
    CLGComplex alpha = _make_cuComplex(F(0.0), F(0.0));
    CLGComplex beta = _make_cuComplex(F(0.0), F(0.0));
    CLGComplex omega = _make_cuComplex(F(0.0), F(0.0));

    for (UINT i = 0; i < m_uiReTry; ++i)
    {
        for (UINT j = 0; j < m_uiStepCount * m_uiDevationCheck; ++j)
        {
            //==========
            //One step
            rho = pRh->Dot(pR);//rho = rh dot r(i-1), if rho = 0, failed (assume will not)
            if (__cuCabsSqf(rho) < _CLG_FLT_MIN)
            {
                appParanoiac(_T("CSLASolverBiCGStab::rho too small:%0.18f, i: %d\n"), __cuCabsSqf(rho), i);
                if (i < m_uiReTry - 1)
                {
                    i = m_uiReTry - 1;
                }
                break;
            }

            if (0 == j) //if is the first iteration, p=r(i-1)
            {
                pR->CopyTo(pP);
            }
            else //if not the first iteration, 
            {
                //beta = last_alpha * rho /(last_omega * last_rho)
                beta = _cuCdivf(_cuCmulf(alpha, rho), _cuCmulf(omega, last_rho));
                //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
                pV->ScalarMultply(omega);
                pP->AxpyMinus(pV);
                pP->ScalarMultply(beta);
                pP->AxpyPlus(pR);
            }

            //v(i) = A p(i)
            pP->CopyTo(pV);
            pV->ApplyOperator(uiM, pGaugeFeild);
            alpha = _cuCdivf(rho, pRh->Dot(pV));//alpha = rho / (rh dot v(i))
            //s=r(i-1) - alpha v(i)
            pR->CopyTo(pS);
            pS->Axpy(_make_cuComplex(-alpha.x, -alpha.y), pV);

            if (0 == (j + 1) % m_uiDevationCheck)
            {
                //Normal of S is small, then stop
                const Real fDeviation = pS->Dot(pS).x / fBLength;
                appParanoiac(_T("CSLASolverBiCGStab::Solve deviation: restart:%d, iteration:%d, deviation:%8.18f\n"), i, j, fDeviation);
                if (fDeviation < m_fAccuracy)
                {
                    //pX->Axpy(alpha, pP); //This is tested a better result not to do the final step
                    pX->CopyTo(pFieldX);

                    //pB->Return();
                    pX->Return();
                    pP->Return();
                    pV->Return();
                    pR->Return();
                    pRh->Return();
                    pS->Return();
                    pT->Return();
                    return TRUE;
                }
            }

            //t=As
            pS->CopyTo(pT);
            pT->ApplyOperator(uiM, pGaugeFeild);
            omega = cuCdivf_cr_host(pS->Dot(pT), pT->Dot(pT).x);//omega = ts / tt

            //r(i)=s-omega t
            pS->CopyTo(pR);
            pR->Axpy(_make_cuComplex(-omega.x, -omega.y), pT);

            //x(i)=x(i-1) + alpha p + omega s
            pX->Axpy(alpha, pP);
            pX->Axpy(omega, pS);

            last_rho = rho;//last_rho = rho
        }

        //we are here, means we do not converge.
        //we need to restart with a new guess, we use last X
        pX->CopyTo(pR);

        pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //-A x_0
        pR->AxpyPlus(pFieldB);  //pR->AxpyPlus(pB); //b - A x_0
        pR->CopyTo(pRh);
    }

    //The solver failed.
    appGeneral(_T("CSLASolverBiCGStab fail to solve!"));

    pX->CopyTo(pFieldX);
    //pB->Return();
    pX->Return();
    pP->Return();
    pV->Return();
    pR->Return();
    pRh->Return();
    pS->Return();
    pT->Return();
    return FALSE;
#endif
}
#endif

//It is tested this is better, the main difference is to let p0 = r0, and rho = r0^* by Yousef Saad.
UBOOL CSLASolverBiCGStab::Solve(CField* pFieldX, const CField* pFieldB, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    //CField* pB = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pP = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pV = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pRh = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pS = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pT = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
#if !_CLG_DOUBLEFLOAT
    DOUBLE fBLength = 1.0;
#else
    Real fBLength = F(1.0);
#endif
    if (!m_bAbsoluteAccuracy)
    {
#if !_CLG_DOUBLEFLOAT
        fBLength = cuCabs(pFieldB->Dot(pFieldB));
#else
        fBLength = _cuCabsf(pFieldB->Dot(pFieldB));
#endif
    }

    appParanoiac(_T("-- CSLASolverBiCGStab::Solve start operator: %s--\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str());

    //pFieldB->CopyTo(pB);

    //Using b as the guess, (Assuming A is near identity?)
    //r_0 = b - A x_0
    if (NULL == pStart)
    {
        pFieldB->CopyTo(pR);
        pFieldB->CopyTo(pX);
    }
    else 
    {
        pStart->CopyTo(pR);
        pStart->CopyTo(pX);
    }


    for (UINT i = 0; i < m_uiReTry; ++i)
    {
        //appGeneral(_T(" ============================ pR ===================\n"));
        //pR->DebugPrintMe();
        pR->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields, EOCT_Minus); //-A x_0
        //appGeneral(_T(" ============================ pR2 ===================\n"));
        //pR->DebugPrintMe();
        pR->AxpyPlus(pFieldB); //pR->AxpyPlus(pB); //b - A x_0
        //appGeneral(_T(" ============================ pR3 ===================\n"));
        //pR->DebugPrintMe();
        //appGeneral(_T(" ============================ ==== ===================\n"));
        pR->CopyTo(pRh);
        //By Yousef Saad, we use conjugated one.
        pRh->Dagger();
        pR->CopyTo(pP);

#if _CLG_DOUBLEFLOAT
        CLGComplex rho = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex alpha = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex beta = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex omega = _make_cuComplex(F(0.0), F(0.0));
#else
        cuDoubleComplex rho = make_cuDoubleComplex(0.0, 0.0);
        cuDoubleComplex alpha = make_cuDoubleComplex(0.0, 0.0);
        cuDoubleComplex beta = make_cuDoubleComplex(0.0, 0.0);
        cuDoubleComplex omega = make_cuDoubleComplex(0.0, 0.0);
#endif

        for (UINT j = 0; j < m_uiStepCount * m_uiDevationCheck; ++j)
        {
            //==========
            //One step
            if (0 == j)
            {
                //for j > 0, rho is calculated at the end of this loop
                rho = pRh->Dot(pR);
            }
            
#if _CLG_DOUBLEFLOAT
            if (__cuCabsSqf(rho) < ((m_uiReTry == i + 1) ? _CLG_FLT_MIN : m_fSmallRho)) //if rho = 0, restart will be better  
            {
                appParanoiac(_T("CSLASolverBiCGStab::rho too small:%0.28f\n"), __cuCabsSqf(rho));
#else
            if (__cuCabsSqfd(rho) < ((m_uiReTry == i + 1) ? _CLG_FLT_MIN : m_fSmallRho))
            {
                appParanoiac(_T("CSLASolverBiCGStab::rho too small:%0.28f\n"), __cuCabsSqfd(rho));
#endif
                i = i + 1;
                break;
            }

            //Before this step, there is precondition for p
            //v(i) = A p(i)
            pP->CopyTo(pV);
            pV->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
#if !_CLG_DOUBLEFLOAT
            alpha = cuCdiv(rho, pRh->Dot(pV));//alpha = rho / (rh dot v(i))
#else
            alpha = _cuCdivf(rho, pRh->Dot(pV));//alpha = rho / (rh dot v(i))
#endif
            
            //s=r(i-1) - alpha v(i)
            pR->CopyTo(pS);
#if !_CLG_DOUBLEFLOAT
            pS->Axpy(_make_cuComplex(-static_cast<Real>(alpha.x), -static_cast<Real>(alpha.y)), pV);
#else
            pS->Axpy(_make_cuComplex(-alpha.x, -alpha.y), pV);
#endif

            if (0 == (j + 1) % m_uiDevationCheck)
            {
                //Normal of S is small, then stop
                const DOUBLE fDeviation = pS->Dot(pS).x / fBLength;
                appParanoiac(_T("CSLASolverBiCGStab::Solve deviation: restart:%d, iteration:%d, deviation:%8.18f\n"), i, j, fDeviation);
                if (fDeviation < m_fAccuracy)
                {
                    //pX->Axpy(alpha, pP); //This is tested a better result not to do the final step
                    pX->CopyTo(pFieldX);

                    //pB->Return();
                    pX->Return();
                    pP->Return();
                    pV->Return();
                    pR->Return();
                    pRh->Return();
                    pS->Return();
                    pT->Return();
                    return TRUE;
                }
            }
            //after this step, there is precondition for s

            //t=As
            pS->CopyTo(pT);
            pT->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);

#if !_CLG_DOUBLEFLOAT
            omega = cuCdivf_cd_host(pS->Dot(pT), pT->Dot(pT).x);//omega = ts / tt
#else
            omega = cuCdivf_cr_host(pS->Dot(pT), pT->Dot(pT).x);//omega = ts / tt
#endif

            //r(i)=s-omega t
            pS->CopyTo(pR);
            
#if !_CLG_DOUBLEFLOAT
            pR->Axpy(_make_cuComplex(-static_cast<Real>(omega.x), -static_cast<Real>(omega.y)), pT);
            beta = cuCdiv(alpha, cuCmul(omega, rho));
            rho = pRh->Dot(pR);
            beta = cuCmul(beta, rho);

            //x(i)=x(i-1) + alpha p + omega s
            pX->Axpy(_cToFloat(alpha), pP);
            pX->Axpy(_cToFloat(omega), pS);

            //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
            pV->ScalarMultply(_cToFloat(omega));
            pP->AxpyMinus(pV);
            pP->ScalarMultply(_cToFloat(beta));
            pP->AxpyPlus(pR);
#else
            pR->Axpy(_make_cuComplex(-omega.x, -omega.y), pT);
            beta = _cuCdivf(alpha, _cuCmulf(omega, rho));
            rho = pRh->Dot(pR);
            beta = _cuCmulf(beta, rho);

            //x(i)=x(i-1) + alpha p + omega s
            pX->Axpy(alpha, pP);
            pX->Axpy(omega, pS);

            //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
            pV->ScalarMultply(omega);
            pP->AxpyMinus(pV);
            pP->ScalarMultply(beta);
            pP->AxpyPlus(pR);
#endif


        }

        //we are here, means we do not converge.
        //we need to restart with a new guess, we use last X
        pX->CopyTo(pR);
    }

    //The solver failed.
    appGeneral(_T("CSLASolverBiCGStab fail to solve! (If it is rho too small, maybe already done)"));

    pX->CopyTo(pFieldX);
    //pB->Return();
    pX->Return();
    pP->Return();
    pV->Return();
    pR->Return();
    pRh->Return();
    pS->Return();
    pT->Return();
    return FALSE;
}



__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================