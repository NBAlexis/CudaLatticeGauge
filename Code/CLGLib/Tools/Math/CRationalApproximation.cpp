//=============================================================================
// FILENAME : CRationalApproximation.cpp
// 
// DESCRIPTION:
//
//
//
// REVISION:
//  [mm/dd/yy]
//  [01/29/2025 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/CRationalApproximation.h"

__BEGIN_NAMESPACE

CRatinalApproximation::CRatinalApproximation(const TArray<Real>& parameters)
{
    m_uiDegree = parameters.Num() / 2;
    appAssert(static_cast<INT>(m_uiDegree * 2 + 1) == parameters.Num());

    m_fC = parameters[0];
    for (UINT i = 0; i < m_uiDegree; ++i)
    {
        m_lstA.AddItem(parameters[1 + i]);
        m_lstB.AddItem(parameters[1 + m_uiDegree + i]);
    }
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(Real) * m_uiDegree));
    appSimpleCopyHD(m_pDeviceData, m_lstA.GetData(), sizeof(Real) * m_uiDegree);
}

void CRatinalApproximation::Initial(const TArray<Real>& parameters)
{
    m_lstA.RemoveAll();
    m_lstB.RemoveAll();
    m_uiDegree = parameters.Num() / 2;
    appAssert(static_cast<INT>(m_uiDegree * 2 + 1) == parameters.Num());
    m_fC = parameters[0];
    for (UINT i = 0; i < m_uiDegree; ++i)
    {
        m_lstA.AddItem(parameters[1 + i]);
        m_lstB.AddItem(parameters[1 + m_uiDegree + i]);
    }
    if (NULL != m_pDeviceData)
    {
        checkCudaErrors(__cudaFree(m_pDeviceData));
    }
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(Real) * m_uiDegree));
    appSimpleCopyHD(m_pDeviceData, m_lstA.GetData(), sizeof(Real) * m_uiDegree);
}

CRatinalApproximation::~CRatinalApproximation()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

CLGAPI CRatinalApproximationSet GRASet;

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
