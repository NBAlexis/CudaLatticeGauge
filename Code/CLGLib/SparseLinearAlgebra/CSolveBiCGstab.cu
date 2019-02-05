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
    , m_pS(NULL)
    , m_pT(NULL)
    , m_pR(NULL)
    , m_pX(NULL)
    , m_pRh(NULL)
    , m_pP(NULL)
    , m_pV(NULL)
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

void CSLASolverBiCGStab::AllocateBuffers(const CField* pField)
{
    ReleaseBuffers();
    m_pS = pField->GetZero();
    m_pT = pField->GetZero();
    m_pR = pField->GetZero();
    m_pX = pField->GetZero();
    m_pRh = pField->GetZero();
    m_pP = pField->GetZero();
    m_pV = pField->GetZero();
}

void CSLASolverBiCGStab::ReleaseBuffers()
{
    appSafeDelete(m_pS);
    appSafeDelete(m_pT);
    appSafeDelete(m_pR);
    appSafeDelete(m_pX);
    appSafeDelete(m_pRh);
    appSafeDelete(m_pP);
    appSafeDelete(m_pV);
}

UBOOL CSLASolverBiCGStab::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM)
{
    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = _cuCabsf(pFieldB->Dot(pFieldB));
    }

    appParanoiac(_T("-- CSLASolverBiCGStab::Solve start --\n"));

    //Using b as the guess, (Assuming M is near identity?)
    pFieldB->CopyTo(m_pX);

    //r_0 = b - A x_0
    pFieldB->CopyTo(m_pR); 
    m_pR->ApplyOperator(uiM, pGaugeFeild); //A x_0
    //appGeneral(_T("==================== kai = %f ==================\n"), m_pR->m_)
    m_pR->ScalarMultply(F(-1.0)); //-A x_0
    m_pR->AxpyPlus(m_pX); //b - A x_0
    m_pR->CopyTo(m_pRh);

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
            rho = _cuCabsf(m_pRh->Dot(m_pR));//rho = rh dot r(i-1), if rho = 0, failed (assume will not)
            if (appAbs(rho) < F(0.00000001))
            {
                break;
            }

            if (0 == j) //if is the first iteration, p=r(i-1)
            {
                m_pR->CopyTo(m_pP);
            }
            else //if not the first iteration, 
            {
                //beta = last_alpha * rho /(last_omega * last_rho)
                beta = alpha * rho / (omega * last_rho);
                //p(i) = r(i-1)+beta( p(i-1) - last_omega v(i-1) )
                m_pV->ScalarMultply(omega);
                m_pP->AxpyMinus(m_pV);
                m_pP->ScalarMultply(beta);
                m_pP->AxpyPlus(m_pR);
            }

            //v(i) = A p(i)
            m_pP->CopyTo(m_pV);
            m_pV->ApplyOperator(uiM, pGaugeFeild);

            alpha = rho / (_cuCabsf(m_pRh->Dot(m_pV)));//alpha = rho / (rh dot v(i))

            //s=r(i-1) - alpha v(i)
            m_pR->CopyTo(m_pS);
            m_pS->Axpy(-alpha, m_pV);

            if (0 == (j - 1) % m_uiDevationCheck)
            {
                //Normal of S is small, then stop
                Real fDeviation = _cuCabsf(m_pS->Dot(m_pS)) / fBLength;
                appParanoiac(_T("CSLASolverBiCGStab::Solve deviation: restart:%d, iteration:%d, deviation:%8.18f\n"), i, j, fDeviation);
                if (fDeviation < m_fAccuracy)
                {
                    //m_pX->Axpy(alpha, m_pP);
                    m_pX->CopyTo(pFieldX);
                    return TRUE;
                }
            }

            //t=As
            m_pS->CopyTo(m_pT);
            m_pT->ApplyOperator(uiM, pGaugeFeild);

            omega = _cuCabsf(m_pT->Dot(m_pS)) / _cuCabsf(m_pT->Dot(m_pT));//omega = ts / tt

            //r(i)=s-omega t
            m_pS->CopyTo(m_pR);
            m_pR->Axpy(-omega, m_pT);

            //x(i)=x(i-1) + alpha p + omega s
            m_pX->Axpy(alpha, m_pP);
            m_pX->Axpy(omega, m_pS);

            last_rho = rho;//last_rho = rho
        }

        //we are here, means we do not converge.
        //we need to restart with a new guess, we use last X
        m_pX->CopyTo(m_pR);

        m_pR->ApplyOperator(uiM, pGaugeFeild); //A x_0
        m_pR->ScalarMultply(-1); //-A x_0
        m_pR->AxpyPlus(m_pX); //b - A x_0
        m_pR->CopyTo(m_pRh);
    }

    //The solver failed.
    appCrucial(_T("CSLASolverBiCGStab fail to solve!"));
    return FALSE;
}



__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================