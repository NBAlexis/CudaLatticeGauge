//=============================================================================
// FILENAME : CGaugeSmearing.cu
// 
// DESCRIPTION:
// put some common functions here
//
// REVISION:
//  [10/04/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleWithoutT(
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3* staple)
{
    intokernal;

    const UINT plaqLengthm1 = plaqLength - 1;     //3
    UINT plaqCountAll = plaqCount * plaqLengthm1; //18

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
            const SIndex& first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
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
    if (!params.FetchValueINT(_T("Iterate"), iValue))
    {
        appCrucial(_T("CGaugeSmearing: Iterate not set, set to 1 by defualt."));
    }
    else if (iValue > 0)
    {
        m_uiIterate = static_cast<UINT>(iValue);
    }
    else
    {
        appCrucial(_T("CGaugeSmearing: Iterate not correct, set to 1 by defualt."));
    }

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
    _kernelStapleWithoutT << <block, threads >> > (
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStapleSU3->m_pDeviceData);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================