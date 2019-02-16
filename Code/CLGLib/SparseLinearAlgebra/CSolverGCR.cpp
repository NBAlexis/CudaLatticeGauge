//=============================================================================
// FILENAME : CSolverGCR.cpp
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [02/16/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverGCR)

CSLASolverGCR::CSLASolverGCR()
    : CSLASolver()
    , m_uiReStart(3)
    , m_uiMaxDim(20)
    , m_uiIterateNumber(50)
    , m_uiCheckError(5)
    , m_fAccuracy(F(0.000001))
    , m_bAbsoluteAccuracy(FALSE)
{

}

CSLASolverGCR::~CSLASolverGCR()
{
    ReleaseBuffers();
}

void CSLASolverGCR::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;

    if (param.FetchValueINT(_T("MaxDim"), iValue))
    {
        m_uiMaxDim = static_cast<UINT>(iValue);
    }

    if (m_uiMaxDim < 5 || m_uiMaxDim > _kMaxStep)
    {
        appCrucial(_T("Max Dim must >= 5 and <= 100, set to default (20)"));
        m_uiMaxDim = 20;
    }

    if (param.FetchValueINT(_T("Restart"), iValue))
    {
        m_uiReStart = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("Iterate"), iValue))
    {
        m_uiIterateNumber = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("DiviationStep"), iValue))
    {
        m_uiCheckError = static_cast<UINT>(iValue);
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

void CSLASolverGCR::AllocateBuffers(const CField* )
{

}

void CSLASolverGCR::ReleaseBuffers()
{

}

UBOOL CSLASolverGCR::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart)
{
    TArray<CField*> pP;
    TArray<CField*> pAP;
    TArray<Real> length_AP;
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        pP.AddItem(pVectors);
        pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        pAP.AddItem(pVectors);
        length_AP.AddItem(F(0.0));
    }

    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pAAP = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = pFieldB->Dot(pFieldB).x;
    }

    appParanoiac(_T("-- CSLASolverGCR::Solve start operator: %s-- fLength = %f --\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str(), fBLength);

    //set initial gauss x0 = b or pStart
    if (NULL == pStart)
    {
        pFieldB->CopyTo(pX);
    }
    else
    {
        pStart->CopyTo(pX);
    }

    Real fLastDiavation = F(0.0);
    for (UINT i = 0; i < m_uiReStart; ++i)
    {
        //r = b - A x0, p0 = r
        pX->CopyTo(pR); 
        pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //r = -A x0
        pR->AxpyPlus(pFieldB); //r = b-Ax0
        pR->CopyTo(pP[0]);

        for (UINT jj = 0; jj < m_uiIterateNumber; ++jj)
        {
            UINT j = jj % m_uiMaxDim;

            pP[j]->CopyTo(pAP[j]);
            pAP[j]->ApplyOperator(uiM, pGaugeFeild);
            length_AP[j] = pAP[j]->Dot(pAP[j]).x;
            //appParanoiac(_T("length p = %f ap = %f r = %f\n"), pP[j]->Dot(pP[j]).x, length_AP[j], pR->Dot(pR).x);
            _Complex al = cuCdivf_cr(pAP[j]->Dot(pR), length_AP[j]);

            pX->Axpy(al, pP[j]);
            pR->Axpy(_make_cuComplex(-al.x, -al.y), pAP[j]);

            if (0 == ((jj + 1) % m_uiCheckError))
            {
                fLastDiavation = pR->Dot(pR).x;
                appParanoiac(_T("CSLASolverGCR::Solve deviation: ---- diviation ----. restart = %d itera = %d divation = %f\n"), i, jj, fLastDiavation);
                if (fLastDiavation < m_fAccuracy * fBLength)
                {
                    pX->CopyTo(pFieldX);

                    pX->Return();
                    pR->Return();
                    pAAP->Return();
                    for (UINT k = 0; k < m_uiMaxDim; ++k)
                    {
                        pP[k]->Return();
                        pAP[k]->Return();
                    }
                    return TRUE;
                }
            }
            
            if (jj != m_uiIterateNumber - 1) //otherwise, just restart
            {
                //p(j+1) = r
                UINT nextPjIndex = (jj + 1) % m_uiMaxDim;
                pAP[j]->CopyTo(pP[nextPjIndex]);
                pAP[j]->CopyTo(pAAP);
                pAAP->ApplyOperator(uiM, pGaugeFeild);

                for (UINT k = 0; k < appMin(jj, m_uiMaxDim); ++k)
                {
                    if (k != nextPjIndex)
                    {
                        _Complex beta = cuCdivf_cr(pAP[k]->Dot(pAAP), -length_AP[j]);
                        pP[nextPjIndex]->Axpy(beta, pP[k]);
                    }
                }
            }
        }
        appParanoiac(_T("CSLASolverGCR::Solve deviation: ---- restart ----. last divation = %f\n"), fLastDiavation);
    }

    appParanoiac(_T("CSLASolverGCR::Solve failed: last divation = %f\n"), fLastDiavation);

    pX->CopyTo(pFieldX);

    pX->Return();
    pR->Return();
    pAAP->Return();
    for (UINT k = 0; k < m_uiMaxDim; ++k)
    {
        pP[k]->Return();
        pAP[k]->Return();
    }
    return FALSE;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================