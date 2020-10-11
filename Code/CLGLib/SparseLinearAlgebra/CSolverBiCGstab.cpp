//=============================================================================
// FILENAME : CSolverBiCGStab.cpp
// 
// DESCRIPTION:
// This is the class for BiCGStab solver
//
// REVISION:
//  [01/08/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverBiCGStab)

CSLASolverBiCGStab::CSLASolverBiCGStab() 
    : CSLASolver()
    , m_uiReTry(1)
    , m_uiDevationCheck(10)
    , m_uiStepCount(20)
    , m_fAccuracy(F(0.000001))
{

}

CSLASolverBiCGStab::~CSLASolverBiCGStab()
{
    CSLASolverBiCGStab::ReleaseBuffers();
}

void CSLASolverBiCGStab::Configurate(const CParameters& param)
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

void CSLASolverBiCGStab::AllocateBuffers(const CField* )
{

}

void CSLASolverBiCGStab::ReleaseBuffers()
{

}

UBOOL CSLASolverBiCGStab::Solve1(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart) const
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
    Real fBLength = F(1.0);
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
#if !_CLG_DOUBLEFLOAT
            rho = _cToFloat(pRh->Dot(pR));
#else
            rho = pRh->Dot(pR);//rho = rh dot r(i-1), if rho = 0, failed (assume will not)
#endif
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
#if !_CLG_DOUBLEFLOAT
            alpha = _cuCdivf(rho, _cToFloat(pRh->Dot(pV)));//alpha = rho / (rh dot v(i))
#else
            alpha = _cuCdivf(rho, pRh->Dot(pV));//alpha = rho / (rh dot v(i))
#endif
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
#if !_CLG_DOUBLEFLOAT
            omega = cuCdivf_cr_host(_cToFloat(pS->Dot(pT)), static_cast<Real>(pT->Dot(pT).x));//omega = ts / tt
#else
            omega = cuCdivf_cr_host(pS->Dot(pT), pT->Dot(pT).x);//omega = ts / tt
#endif

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
}

//It is tested this is better, the main difference is to let p0 = r0, and rho = r0^* by Yousef Saad.
UBOOL CSLASolverBiCGStab::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
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
    Real fBLength = F(1.0);
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
        pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //-A x_0
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

        CLGComplex rho = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex alpha = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex beta = _make_cuComplex(F(0.0), F(0.0));
        CLGComplex omega = _make_cuComplex(F(0.0), F(0.0));

        for (UINT j = 0; j < m_uiStepCount * m_uiDevationCheck; ++j)
        {
            //==========
            //One step
            if (0 == j)
            {
                //for j > 0, rho is calculated at the end of this loop
#if !_CLG_DOUBLEFLOAT
                rho = _cToFloat(pRh->Dot(pR));
#else
                rho = pRh->Dot(pR);
#endif
            }
            
            if (__cuCabsSqf(rho) < _CLG_FLT_MIN) //if rho = 0, failed (assume will not)    
            {
                appParanoiac(_T("CSLASolverBiCGStab::rho too small:%0.28f\n"), __cuCabsSqf(rho));
                if (i < m_uiReTry - 1)
                {
                    i = m_uiReTry - 1;
                }
                break;
            }

            //Before this step, there is precondition for p
            //v(i) = A p(i)
            pP->CopyTo(pV);
            pV->ApplyOperator(uiM, pGaugeFeild);
#if !_CLG_DOUBLEFLOAT
            alpha = _cuCdivf(rho, _cToFloat(pRh->Dot(pV)));//alpha = rho / (rh dot v(i))
#else
            alpha = _cuCdivf(rho, pRh->Dot(pV));//alpha = rho / (rh dot v(i))
#endif
            
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
            //after this step, there is precondition for s

            //t=As
            pS->CopyTo(pT);
            pT->ApplyOperator(uiM, pGaugeFeild);

#if !_CLG_DOUBLEFLOAT
            omega = cuCdivf_cr_host(_cToFloat(pS->Dot(pT)), static_cast<Real>(pT->Dot(pT).x));//omega = ts / tt
#else
            omega = cuCdivf_cr_host(pS->Dot(pT), pT->Dot(pT).x);//omega = ts / tt
#endif

            //r(i)=s-omega t
            pS->CopyTo(pR);
            pR->Axpy(_make_cuComplex(-omega.x, -omega.y), pT);
            beta = _cuCdivf(alpha, _cuCmulf(omega, rho));
#if !_CLG_DOUBLEFLOAT
            rho = _cToFloat(pRh->Dot(pR));
#else
            rho = pRh->Dot(pR);
#endif
            beta = _cuCmulf(beta, rho);

            //x(i)=x(i-1) + alpha p + omega s
            pX->Axpy(alpha, pP);
            pX->Axpy(omega, pS);

            //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
            pV->ScalarMultply(omega);
            pP->AxpyMinus(pV);
            pP->ScalarMultply(beta);
            pP->AxpyPlus(pR);
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