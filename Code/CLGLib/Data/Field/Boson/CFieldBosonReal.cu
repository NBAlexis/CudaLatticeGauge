//=============================================================================
// FILENAME : CFieldBosonReal.cu
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [mm/dd/yy]
//  [11/04/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldBosonReal.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldBosonReal)

__CLGIMPLEMENT_CLASS(CFieldBosonRealD_ZSlice)

__global__ void _CLG_LAUNCH_BOUND
_kernelSetZSlice(Real* pDeviceData, BYTE zslice, Real fValue)
{
    intokernalInt4;
    if (sSite4.z == zslice)
    {
        intokernal;
        pDeviceData[uiSiteIndex] = fValue;
    }
}

void CFieldBosonReal::SetZSlice(BYTE zIndex, Real fValue)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelSetZSlice, block, threads, m_pDeviceData, zIndex, fValue);
}

void CFieldBosonRealD_ZSlice::InitialOtherParameters(CParameters& param)
{
    CFieldBosonReal::InitialOtherParameters(param);

    param.FetchValueArrayBYTE(_T("ZSliceList"), m_lstFixedZSlce);
    param.FetchValueArrayReal(_T("ValueList"), m_lstFixedValue);
}

void CFieldBosonRealD_ZSlice::FixBoundary(EFixBoundary eType)
{
    //appGeneral(_T("fix boundary %d\n"), eType);
    if (EFB_Field == eType)
    {
        for (INT i = 0; i < m_lstFixedZSlce.Num() && i < m_lstFixedValue.Num(); ++i)
        {
            SetZSlice(m_lstFixedZSlce[i], m_lstFixedValue[i]);
        }
    }
    else
    {
        for (INT i = 0; i < m_lstFixedZSlce.Num() && i < m_lstFixedValue.Num(); ++i)
        {
            SetZSlice(m_lstFixedZSlce[i], F(0.0));
        }
    }
}

void CFieldBosonRealD_ZSlice::CopyParamTo(CField* U) const
{
    CFieldBosonReal::CopyParamTo(U);
    CFieldBosonRealD_ZSlice* target = dynamic_cast<CFieldBosonRealD_ZSlice*>(U);
    if (NULL != target)
    {
        target->m_lstFixedZSlce = m_lstFixedZSlce;
        target->m_lstFixedValue = m_lstFixedValue;
    }
}

CCString CFieldBosonRealD_ZSlice::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldBosonReal::GetInfos(tab);
    sRet = sRet + tab + _T("ZSliceList : ") + appToString(m_lstFixedZSlce) + _T("\n");
    sRet = sRet + tab + _T("ValueList : ") + appToString(m_lstFixedValue) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================