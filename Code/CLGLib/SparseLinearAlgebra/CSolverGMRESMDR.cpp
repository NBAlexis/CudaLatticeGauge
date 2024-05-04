//=============================================================================
// FILENAME : CSLASolverGMRESMDR.cpp
// 
// DESCRIPTION:
// This is the class for GMRES Solver
//
// REVISION:
//  [03/24/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverGMRESMDR)

void CSLASolverGMRESMDR::GenerateCUFirstTime(CField* pX, CField* pR, const CField* pFieldB, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    EFieldOperator uiM)
{
    QRFactorizationOfUk();

    //Change Ck to AQk
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        m_lstU[i]->CopyTo(m_lstC[i]);
        m_lstC[i]->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields);
    }

    //m_pHostHmGm
    //m_pHostB
    //m_pHostPk
    switch (m_eDeflationType)
    {
    case EEDT_REV:
    {
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            for (UINT j = 0; j < m_uiKDim; ++j)
            {
                //Q^+ AQ
#if !_CLG_DOUBLEFLOAT
                m_pHostHmGm[i * m_uiKDim + j] = _cToFloat(m_lstU[i]->Dot(m_lstC[j]));
#else
                m_pHostHmGm[i * m_uiKDim + j] = m_lstU[i]->Dot(m_lstC[j]);
#endif
            }
        }

        checkCudaErrors(cudaMemcpy(m_pDeviceHm, m_pHostHmGm, sizeof(CLGComplex) * m_uiKDim * m_uiKDim, cudaMemcpyHostToDevice));
        m_pHelper->EigenValueProblem(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiKDim, m_uiKDim);
    }
    break;
    case EEDT_HEV:
    {
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            for (UINT j = 0; j < m_uiKDim; ++j)
            {
                //Q^+ A^+ Q
#if !_CLG_DOUBLEFLOAT
                m_pHostB[i * m_uiKDim + j] = _cToFloat(m_lstC[i]->Dot(m_lstU[j]));
                if (j >= i)
                {
                    //Q^+A^+ AQ
                    m_pHostHmGm[i * m_uiKDim + j] = _cToFloat(m_lstC[i]->Dot(m_lstC[j]));
            }
#else
                m_pHostB[i * m_uiKDim + j] = m_lstC[i]->Dot(m_lstU[j]);
                if (j >= i)
                {
                    //Q^+A^+ AQ
                    m_pHostHmGm[i * m_uiKDim + j] = m_lstC[i]->Dot(m_lstC[j]);
                }
#endif

            }
        }
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            for (UINT j = 0; j < i; ++j)
            {
                m_pHostHmGm[i * m_uiKDim + j] = _cuConjf(m_pHostHmGm[j * m_uiKDim + i]);
            }
        }

        checkCudaErrors(cudaMemcpy(m_pDeviceA, m_pHostHmGm, sizeof(CLGComplex) * m_uiKDim * m_uiKDim, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(m_pDeviceB, m_pHostB, sizeof(CLGComplex) * m_uiKDim * m_uiKDim, cudaMemcpyHostToDevice));

        m_pHelper->GeneralizedEigenValueProblem(m_pDeviceA, m_pDeviceB, m_pDeviceEigenValue, m_pDevicePk, m_uiKDim, m_uiKDim, FALSE);
    }
    break;
    case EEDT_SVD:
    {
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            for (UINT j = i; j < m_uiKDim; ++j)
            {
                //Q^+A^+ AQ
#if !_CLG_DOUBLEFLOAT
                m_pHostHmGm[i * m_uiKDim + j] = _cToFloat(m_lstC[i]->Dot(m_lstC[j]));
#else
                m_pHostHmGm[i * m_uiKDim + j] = m_lstC[i]->Dot(m_lstC[j]);
#endif
            }
        }
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            for (UINT j = 0; j < i; ++j)
            {
                m_pHostHmGm[i * m_uiKDim + j] = _cuConjf(m_pHostHmGm[j * m_uiKDim + i]);
            }
        }

        checkCudaErrors(cudaMemcpy(m_pDeviceHm, m_pHostHmGm, sizeof(CLGComplex) * m_uiKDim * m_uiKDim, cudaMemcpyHostToDevice));
        m_pHelper->EigenValueProblem(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiKDim, m_uiKDim);
    }
    break;
    }

    //Uk = Q Pk
    m_pFieldMatrix->VectorMultiplyMatrix(m_lstU, m_lstV, m_pDevicePk, m_uiKDim, m_uiKDim);

    //Initial X and R
    pX->CopyTo(pR);
    pR->ApplyOperator(uiM, gaugeNum, bosonNum, gaugeFields, bosonFields, EOCT_Minus); //r0 = -A x0
    pR->AxpyPlus(pFieldB); //r0 = b-Ax0

#if !_CLG_DOUBLEFLOAT
    m_cLastDiviation = _cToFloat(pR->Dot(pR));
#else
    m_cLastDiviation = pR->Dot(pR);
#endif
    m_fDiviation = _hostsqrt(__cuCabsSqf(m_cLastDiviation));
}

/**
* [U,R]=U
*/
void CSLASolverGMRESMDR::QRFactorizationOfUk()
{
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        for (UINT j = 0; j < m_uiKDim; ++j)
        {
            m_pHostTmpR[i * m_uiKDim + j] = _make_cuComplex(F(0.0), F(0.0));
        }
    }

    //QR of AY
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
#if !_CLG_DOUBLEFLOAT
        const Real fLength = static_cast<Real>(_hostsqrtd(m_lstU[i]->Dot(m_lstU[i]).x));
#else
        const Real fLength = _hostsqrt(m_lstU[i]->Dot(m_lstU[i]).x);
#endif
        m_pHostTmpR[i * m_uiKDim + i] = _make_cuComplex(fLength, F(0.0));
        m_lstU[i]->ScalarMultply(F(1.0) / fLength);
        for (UINT j = i + 1; j < m_uiKDim; ++j)
        {
#if !_CLG_DOUBLEFLOAT
            m_pHostTmpR[i * m_uiKDim + j] = _cToFloat(m_lstU[i]->Dot(m_lstU[j]));
#else
            m_pHostTmpR[i * m_uiKDim + j] = m_lstU[i]->Dot(m_lstU[j]);
#endif
            m_lstU[j]->Axpy(
                _make_cuComplex(
                    -m_pHostTmpR[i * m_uiKDim + j].x,
                    -m_pHostTmpR[i * m_uiKDim + j].y),
                m_lstU[i]);
        }
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================