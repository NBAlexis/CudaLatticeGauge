//=============================================================================
// FILENAME : CGaugeSmearing.cu
// 
// DESCRIPTION:
// put some common functions here
//
// REVISION:
//  [mm/dd/yy]
//  [10/04/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleWithoutT(
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceSU3* staple)
{
    intokernal;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;     //3
    UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //18
#endif

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        //plaqCount = 6, they are:
        //j=0,1,2,3 remove mu, for example, when linkIndex is for y, 1, it caches 10,10-,12,12-,13,13-

        //If only spatial, we only use the first 4
        #pragma unroll
        for (int i = 0; i < 4; ++i)
        {
            const SIndex& first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAllLink];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAllLink];
                deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }

        staple[linkIndex] = res;
    }
    staple[_deviceGetLinkIndex(uiSiteIndex, 3)] = deviceSU3::makeSU3Id();
}

#pragma endregion

void CGaugeSmearing::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    INT iValue = 0;
    params.FetchValueINT(_T("HasT"), iValue);
    m_bHasT = (0 != iValue);

    iValue = 1;
    params.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    if (!params.FetchValueINT(_T("Iterate"), iValue))
    {
        appWarning(_T("CGaugeSmearing: Iterate not set, set to 1 by defualt."));
    }
    else if (iValue > 0)
    {
        m_uiIterate = static_cast<UINT>(iValue);
    }
    else
    {
        appWarning(_T("CGaugeSmearing: Iterate not correct, set to 1 by defualt."));
    }
    pOwner->m_pGaugeSmearing[m_byFieldId] = this;

    iValue = 0;
    params.FetchValueINT(_T("Update"), iValue);
    m_bCalledWhenUpdate = (0 != iValue);
}

CCString CGaugeSmearing::GetInfos(const CCString& tab) const
{
    CCString sRet = CBase::GetInfos(tab);
    sRet = sRet + tab + _T("HasT : ") + appToString(m_bHasT) + _T("\n");
    sRet = sRet + tab + _T("Iterate : ") + appToString(m_uiIterate) + _T("\n");
    sRet = sRet + tab + _T("FieldId : ") + appToString(m_byFieldId) + _T("\n");
    sRet = sRet + tab + _T("Called when Update : ") + appToString(m_bCalledWhenUpdate) + _T("\n");
    return sRet;
}

void CGaugeSmearing::CalculateSpatialFatLink(const CFieldGauge* pGauge, CFieldGauge* pFatlink) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    if (NULL == pFatlink || EFT_GaugeSU3 != pFatlink->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pStapleSU3 = dynamic_cast<CFieldGaugeSU3*>(pFatlink);
    preparethread;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleWithoutT, block, threads,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStapleSU3->m_pDeviceData);
#else
    _LAUNCH_KERNEL(_kernelStapleWithoutT, block, threads, 
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
        pStapleSU3->m_pDeviceData);
#endif
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================