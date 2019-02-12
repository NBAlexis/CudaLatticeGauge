//=============================================================================
// FILENAME : CSolverBiCGStab.cpp
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [1/8/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverBiCGStab)

CSLASolverBiCGStab::CSLASolverBiCGStab() 
    : m_pOwner(NULL) 
    , m_uiReTry(1)
    , m_uiDevationCheck(10)
    , m_uiStepCount(20)
    , m_fAccuracy(F(0.000001))
    , m_bAbsoluteAccuracy(FALSE)
{

}

CSLASolverBiCGStab::~CSLASolverBiCGStab()
{
    ReleaseBuffers();
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
    }
}

void CSLASolverBiCGStab::AllocateBuffers(const CField* )
{

}

void CSLASolverBiCGStab::ReleaseBuffers()
{

}

UBOOL CSLASolverBiCGStab::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart)
{
    CField* pB = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
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

    pFieldB->CopyTo(pB);

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
    pR->ApplyOperator(uiM, pGaugeFeild); //A x_0
    pR->ScalarMultply(F(-1.0)); //-A x_0
    pR->AxpyPlus(pB); //b - A x_0
    pR->CopyTo(pRh);

    Real rho = 0;
    Real last_rho = 0;
    Real alpha = 0;
    Real beta = 0;
    Real omega = 0;

    for (UINT i = 0; i < m_uiReTry; ++i)
    {
        for (UINT j = 0; j < m_uiStepCount * m_uiDevationCheck; ++j)
        {
            //==========
            //One step
            rho = _cuCabsf(pRh->Dot(pR));//rho = rh dot r(i-1), if rho = 0, failed (assume will not)
            if (appAbs(rho) < FLT_MIN)
            {
                appParanoiac(_T("CSLASolverBiCGStab::rho too small:%0.18f\n"), rho);
                break;
            }

            if (0 == j) //if is the first iteration, p=r(i-1)
            {
                pR->CopyTo(pP);
            }
            else //if not the first iteration, 
            {
                //beta = last_alpha * rho /(last_omega * last_rho)
                beta = alpha * rho / (omega * last_rho);
                //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
                pV->ScalarMultply(omega);
                pP->AxpyMinus(pV);
                pP->ScalarMultply(beta);
                pP->AxpyPlus(pR);
            }

            //v(i) = A p(i)
            pP->CopyTo(pV);
            pV->ApplyOperator(uiM, pGaugeFeild);

            alpha = rho / (_cuCabsf(pRh->Dot(pV)));//alpha = rho / (rh dot v(i))

            //s=r(i-1) - alpha v(i)
            pR->CopyTo(pS);
            pS->Axpy(-alpha, pV);

            if (0 != j && (0 == j % m_uiDevationCheck))
            {
                //Normal of S is small, then stop
                Real fDeviation = _cuCabsf(pS->Dot(pS)) / fBLength;
                appParanoiac(_T("CSLASolverBiCGStab::Solve deviation: restart:%d, iteration:%d, deviation:%8.18f\n"), i, j, fDeviation);
                if (fDeviation < m_fAccuracy)
                {
                    //pX->Axpy(alpha, pP); //This is tested a better result not to do the final step
                    pX->CopyTo(pFieldX);

                    pB->Return();
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

            omega = _cuCabsf(pT->Dot(pS)) / _cuCabsf(pT->Dot(pT));//omega = ts / tt

            //r(i)=s-omega t
            pS->CopyTo(pR);
            pR->Axpy(-omega, pT);

            //x(i)=x(i-1) + alpha p + omega s
            pX->Axpy(alpha, pP);
            pX->Axpy(omega, pS);

            last_rho = rho;//last_rho = rho
        }

        //we are here, means we do not converge.
        //we need to restart with a new guess, we use last X
        pX->CopyTo(pR);

        pR->ApplyOperator(uiM, pGaugeFeild); //A x_0
        pR->ScalarMultply(-1); //-A x_0
        pR->AxpyPlus(pB); //b - A x_0
        pR->CopyTo(pRh);
    }

    //The solver failed.
    appCrucial(_T("CSLASolverBiCGStab fail to solve!"));

    pX->CopyTo(pFieldX);
    pB->Return();
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