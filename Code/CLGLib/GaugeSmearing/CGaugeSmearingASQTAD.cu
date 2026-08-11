//=============================================================================
// FILENAME : CGaugeSmearingASQTAD.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [11/14/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CGaugeSmearingASQTAD.h"

#define _CLG_ASQTAD_INCLUDE_DEBUG 0

__BEGIN_NAMESPACE

#if _CLG_ASQTAD_INCLUDE_DEBUG

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkTTwoType(
    const deviceGauge* __restrict__ pDeviceData1,
    const deviceGauge* __restrict__ pDeviceData2,
    BYTE deviceData2Index,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceGauge sRet = _makeId<deviceGauge>();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        const deviceGauge* pDeviceData = (deviceData2Index == i) ? pDeviceData2 : pDeviceData1;

        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _dagger(sRet);
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkTThreeType(
    const deviceGauge* __restrict__ pDeviceData1,
    const deviceGauge* __restrict__ pDeviceData2,
    const deviceGauge* __restrict__ pDeviceData3,
    BYTE deviceData2Index,
    BYTE deviceData3Index,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceGauge sRet = _makeId<deviceGauge>();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        const deviceGauge* pDeviceData = (deviceData2Index == i) ? pDeviceData2 :
            ((deviceData3Index == i) ? pDeviceData3 : pDeviceData1);

        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _dagger(sRet);
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

/**
*
* -- f0 -->
*         |
*         |
*         v
* \       \
*  \       \
*   \_______\
* 3-2-3-1-2
*
* -- f0 -->
* |
* |
* |
* \       \
*  \       \
*   \_______\
* 2-1-3-2-3
*
* ---f0--->
* |       |
* |       |
* |       v
* \
*  \
*   \_______
* 3-2-3-1-2
*
* ---f0--->
* |       |
* |       |
* |       v
* \       \
*  \       \
*   \       \
* 2-3-1-3-2
*
* ---f0--->
* |       |
* |       |
* |       v
*         \
*          \
*    _______\
* 2-1-3-2-3
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5Force_Naive0(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalInt4dirC;

    for (BYTE idir1m = 0; idir1m < uiDir; ++idir1m)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir1m);
        forceres[linkIndex] = _makeZero<deviceGauge>();

        for (BYTE idir2m = 0; idir2m < uiDir; ++idir2m)
        {
            if (idir2m == idir1m)
            {
                continue;
            }

            for (BYTE idir3m = 0; idir3m < uiDir; ++idir3m)
            {
                if (idir3m == idir1m || idir3m == idir2m)
                {
                    continue;
                }

                //3-2-3-1-2 -- 5 f0
                //2-1-3-2-3 -- 1 f0
                //3-2-3-1-2 -- 2 f0
                //2-3-1-3-2 -- 3 f0
                //2-1-3-2-3 -- 4 f0
                const INT idir1 = idir1m + 1;
                const INT idir2 = idir2m + 1;
                const INT idir3 = idir3m + 1;

                const INT path1ff[5] = { idir3, idir2, -idir3, idir1, -idir2 };
                const INT path1fb[5] = { idir3, -idir2, -idir3, idir1, idir2 };
                const INT path1bf[5] = { -idir3, idir2, idir3, idir1, -idir2 };
                const INT path1bb[5] = { -idir3, -idir2, idir3, idir1, idir2 };

                const INT path2ff[5] = { idir2, idir1, idir3, -idir2, -idir3 };
                const INT path2fb[5] = { idir2, idir1, -idir3, -idir2, idir3 };
                const INT path2bf[5] = { -idir2, idir1, idir3, idir2, -idir3 };
                const INT path2bb[5] = { -idir2, idir1, -idir3, idir2, idir3 };

                const INT path3ff[5] = { idir3, idir2, -idir3, idir1, -idir2 };
                const INT path3fb[5] = { idir3, -idir2, -idir3, idir1, idir2 };
                const INT path3bf[5] = { -idir3, idir2, idir3, idir1, -idir2 };
                const INT path3bb[5] = { -idir3, -idir2, idir3, idir1, idir2 };

                const INT path4ff[5] = { idir2, idir3, idir1, -idir3, -idir2 };
                const INT path4fb[5] = { idir2, -idir3, idir1, idir3, -idir2 };
                const INT path4bf[5] = { -idir2, idir3, idir1, -idir3, idir2 };
                const INT path4bb[5] = { -idir2, -idir3, idir1, idir3, idir2 };

                const INT path5ff[5] = { idir2, idir1, idir3, -idir2, -idir3 };
                const INT path5fb[5] = { idir2, idir1, -idir3, -idir2, idir3 };
                const INT path5bf[5] = { -idir2, idir1, idir3, idir2, -idir3 };
                const INT path5bb[5] = { -idir2, idir1, -idir3, idir2, idir3 };

                //===================== first loop =======================
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 4, sSite4, 5, byFieldId, path1ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 4, sSite4, 5, byFieldId, path1fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 4, sSite4, 5, byFieldId, path1bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 4, sSite4, 5, byFieldId, path1bb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 0, sSite4, 5, byFieldId, path2ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 0, sSite4, 5, byFieldId, path2fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 0, sSite4, 5, byFieldId, path2bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 0, sSite4, 5, byFieldId, path2bb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3bb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4bb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5bb));
            }
        }
    }
}

/**
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
*
* -- f0 -->
*         |
*         |
*         v
* \       \
*  \       \
*   \_______\
* 3-2-3-1-2
*
*    f0  0      1       2       3
* me
* 0      -    zy ty   yz tz   yt zt     1 2   1 2   1 2
* 1   zx tx     -     xz tz   xt zt     1 2   0 2   0 2
* 2   yx tx   xy ty     -     xt yt     0 2   0 2   0 1
* 3   yx zx   xy zy   xz yz     -       0 1   0 1   0 1
*
* -- f0 -->
* |
* |
* |
* \       \
*  \       \
*   \_______\
* 2-1-3-2-3
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5Force_Naive1(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalInt4dirC;

    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const BYTE t1idx[4][6] = {
        {1, 2, 1, 2, 1, 2},
        {1, 2, 0, 2, 0, 2},
        {0, 2, 0, 2, 0, 1},
        {0, 1, 0, 1, 0, 1}
    };

    for (BYTE idir1m = 0; idir1m < uiDir; ++idir1m)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir1m);
        forceres[linkIndex] = _makeZero<deviceGauge>();
        BYTE p5id = 0;
        for (BYTE idir2m = 0; idir2m < uiDir; ++idir2m)
        {
            if (idir2m == idir1m)
            {
                continue;
            }

            for (BYTE idir3m = 0; idir3m < uiDir; ++idir3m)
            {
                if (idir3m == idir1m || idir3m == idir2m)
                {
                    continue;
                }

                const deviceGauge* p3 = pointersp3[t1idx[idir1m][p5id]];

                //3-2-3-1-2 -- 5 f0
                //2-1-3-2-3 -- 1 f0
                //3-2-3-1-2 -- 2 f0
                //2-3-1-3-2 -- 3 f0
                //2-1-3-2-3 -- 4 f0
                const INT idir1 = idir1m + 1;
                const INT idir2 = idir2m + 1;
                const INT idir3 = idir3m + 1;

                const INT path1f[3] = { idir2, idir1, -idir2 };
                const INT path1b[3] = { -idir2, idir1, idir2 };

                const INT path2f[3] = { idir2, idir1, -idir2 };
                const INT path2b[3] = { -idir2, idir1, idir2 };

                const INT path3ff[5] = { idir3, idir2, -idir3, idir1, -idir2 };
                const INT path3fb[5] = { idir3, -idir2, -idir3, idir1, idir2 };
                const INT path3bf[5] = { -idir3, idir2, idir3, idir1, -idir2 };
                const INT path3bb[5] = { -idir3, -idir2, idir3, idir1, idir2 };

                const INT path4ff[5] = { idir2, idir3, idir1, -idir3, -idir2 };
                const INT path4fb[5] = { idir2, -idir3, idir1, idir3, -idir2 };
                const INT path4bf[5] = { -idir2, idir3, idir1, -idir3, idir2 };
                const INT path4bb[5] = { -idir2, -idir3, idir1, idir3, idir2 };

                const INT path5ff[5] = { idir2, idir1, idir3, -idir2, -idir3 };
                const INT path5fb[5] = { idir2, idir1, -idir3, -idir2, idir3 };
                const INT path5bf[5] = { -idir2, idir1, idir3, idir2, -idir3 };
                const INT path5bb[5] = { -idir2, idir1, -idir3, idir2, idir3 };

                //===================== first loop =======================
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 2, 0, sSite4, 3, byFieldId, path1f));
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 2, 0, sSite4, 3, byFieldId, path1b));

                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 0, 2, sSite4, 3, byFieldId, path2f));
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 0, 2, sSite4, 3, byFieldId, path2b));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 1, sSite4, 5, byFieldId, path3bb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4bb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 3, sSite4, 5, byFieldId, path5bb));

                ++p5id;
            }
        }
    }
}

/**
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
*
* ---f0--->
* |       |
* |       |
* |       v
* \
*  \
*   \_______
* 3-2-3-1-2
*
*    f0  0      1       2       3
* me
* 0      -    zy ty   yz tz   yt zt     1 2   1 2   1 2
* 1   zx tx     -     xz tz   xt zt     1 2   0 2   0 2
* 2   yx tx   xy ty     -     xt yt     0 2   0 2   0 1
* 3   yx zx   xy zy   xz yz     -       0 1   0 1   0 1
*
*
* ---f0--->
* |       |
* |       |
* |       v
*         \
*          \
*    _______\
* 2-1-3-2-3
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5Force_Naive2(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    const deviceGauge* __restrict__ p3f_1,
    const deviceGauge* __restrict__ p3f_2,
    const deviceGauge* __restrict__ p3f_3,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalInt4dirC;

    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };

    const BYTE t1idx[4][6] = {
        {1, 2, 1, 2, 1, 2},
        {1, 2, 0, 2, 0, 2},
        {0, 2, 0, 2, 0, 1},
        {0, 1, 0, 1, 0, 1}
    };

    for (BYTE idir1m = 0; idir1m < uiDir; ++idir1m)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir1m);
        forceres[linkIndex] = _makeZero<deviceGauge>();
        BYTE p5id = 0;
        for (BYTE idir2m = 0; idir2m < uiDir; ++idir2m)
        {
            if (idir2m == idir1m)
            {
                continue;
            }

            for (BYTE idir3m = 0; idir3m < uiDir; ++idir3m)
            {
                if (idir3m == idir1m || idir3m == idir2m)
                {
                    continue;
                }

                const deviceGauge* p3 = pointersp3[t1idx[idir1m][p5id]];
                const deviceGauge* pf3 = pointerspf3[t1idx[idir1m][p5id]];

                const INT idir1 = idir1m + 1;
                const INT idir2 = idir2m + 1;
                const INT idir3 = idir3m + 1;

                const INT pathf[3] = { idir2, idir1, -idir2 };
                const INT pathb[3] = { -idir2, idir1, idir2 };

                const INT path4ff[5] = { idir2, idir3, idir1, -idir3, -idir2 };
                const INT path4fb[5] = { idir2, -idir3, idir1, idir3, -idir2 };
                const INT path4bf[5] = { -idir2, idir3, idir1, -idir3, idir2 };
                const INT path4bb[5] = { -idir2, -idir3, idir1, idir3, idir2 };

                //===================== first loop =======================
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 2, 0, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 2, 0, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 0, 2, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 0, 2, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 0, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 0, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 2, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 2, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4ff));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4fb));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4bf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf0, 2, sSite4, 5, byFieldId, path4bb));

                ++p5id;
            }
        }
    }
}

/**
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
*
* ---f0--->
* |       |
* |       |
* |       v
* \       \
*  \       \ dir2
*   \       \
* 2-3-1-3-2
* me is dir1, f0 is dir1
*
*  dir2  0      1       2       3
* me
* 0      -    zx tx   yx tx   yx zx     1 2   0 2   0 1
* 1   zy ty     -     xy ty   xy zy     1 2   0 2   0 1
* 2   yz tz   xz tz     -     xz yz     1 2   0 2   0 1
* 3   yt zt   xt zt   xt yt     -       1 2   0 2   0 1
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5Force_Naive3(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    const deviceGauge* __restrict__ p3f_1,
    const deviceGauge* __restrict__ p3f_2,
    const deviceGauge* __restrict__ p3f_3,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalInt4dirC;

    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };

    const BYTE t1idx[4][6] = {
        {1, 2, 1, 2, 1, 2},
        {1, 2, 0, 2, 0, 2},
        {0, 2, 0, 2, 0, 1},
        {0, 1, 0, 1, 0, 1}
    };
    const BYTE t2idx[6] = { 1, 2, 0, 2, 0, 1 };

    for (BYTE idir1m = 0; idir1m < uiDir; ++idir1m)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir1m);
        forceres[linkIndex] = _makeZero<deviceGauge>();
        BYTE p5id = 0;
        for (BYTE idir2m = 0; idir2m < uiDir; ++idir2m)
        {
            if (idir2m == idir1m)
            {
                continue;
            }

            for (BYTE idir3m = 0; idir3m < uiDir; ++idir3m)
            {
                if (idir3m == idir1m || idir3m == idir2m)
                {
                    continue;
                }

                const deviceGauge* p3 = pointersp3[t1idx[idir1m][p5id]];
                const deviceGauge* pf3 = pointerspf3[t1idx[idir1m][p5id]];
                const deviceGauge* pf3_2 = pointerspf3[t2idx[p5id]];

                const INT idir1 = idir1m + 1;
                const INT idir2 = idir2m + 1;
                const INT idir3 = idir3m + 1;

                const INT pathf[3] = { idir2, idir1, -idir2 };
                const INT pathb[3] = { -idir2, idir1, idir2 };

                //===================== first loop =======================
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 2, 0, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 2, 0, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 0, 2, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTThreeType(pDeviceData, pf0, p3, 0, 2, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 0, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 0, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 2, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3, 2, sSite4, 3, byFieldId, pathb));

                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3_2, 1, sSite4, 3, byFieldId, pathf));
                _add(forceres[linkIndex], _deviceLinkTTwoType(pDeviceData, pf3_2, 1, sSite4, 3, byFieldId, pathb));

                ++p5id;
            }
        }
    }
}

#endif

//Means to optimize register usage, however, it looks that the PTX has already optimize the register usage
template<typename deviceGauge>
static __device__ __inline__ deviceGauge __deviceOneStaple0(const SIndex* __restrict__ s, const deviceGauge* __restrict__ p1, const deviceGauge* __restrict__ p2, const deviceGauge* __restrict__ p3, BYTE byFieldId)
{
    SIndex site = s[0];
    deviceGauge res(_deviceGetGaugeBCT(byFieldId, p1, site));
    if (site.NeedToDagger())
    {
        _dagger(res);
    }
    site = s[1];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p2, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p2, site));
    }
    site = s[2];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p3, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p3, site));
    }
    return res;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge __deviceOneStaple1(const SIndex * __restrict__ s, const deviceGauge* __restrict__ p1, const deviceGauge* __restrict__ p2, const deviceGauge* __restrict__ p3, BYTE byFieldId)
{
    SIndex site = s[0];
    deviceGauge res(_deviceGetGaugeBCDirForceSIndexT(p1, site));
    if (site.NeedToDagger())
    {
        _dagger(res);
    }
    site = s[1];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p2, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p2, site));
    }
    site = s[2];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p3, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p3, site));
    }
    return res;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge __deviceOneStaple2(const SIndex* __restrict__ s, const deviceGauge* __restrict__ p1, const deviceGauge* __restrict__ p2, const deviceGauge* __restrict__ p3, BYTE byFieldId)
{
    SIndex site = s[0];
    deviceGauge res(_deviceGetGaugeBCT(byFieldId, p1, site));
    if (site.NeedToDagger())
    {
        _dagger(res);
    }
    site = s[1];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCDirForceSIndexT(p2, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCDirForceSIndexT(p2, site));
    }
    site = s[2];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p3, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p3, site));
    }
    return res;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge __deviceOneStaple3(const SIndex* __restrict__ s, const deviceGauge* __restrict__ p1, const deviceGauge* __restrict__ p2, const deviceGauge* __restrict__ p3, BYTE byFieldId)
{
    SIndex site = s[0];
    deviceGauge res(_deviceGetGaugeBCT(byFieldId, p1, site));
    if (site.NeedToDagger())
    {
        _dagger(res);
    }
    site = s[1];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p2, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p2, site));
    }
    site = s[2];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCDirForceSIndexT(p3, site));
    }
    else
    {
        _mul(res, _deviceGetGaugeBCDirForceSIndexT(p3, site));
    }
    return res;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge __deviceOneStapleBuff(const SIndex* __restrict__ s, const deviceGauge* __restrict__ p1, const deviceGauge* __restrict__ p2, const deviceGauge* __restrict__ p3, BYTE byFieldId, CLGComplex* buff)
{
    SIndex site = s[0];
    deviceGauge res(_deviceGetGaugeBCT(byFieldId, p1, site));
    if (site.NeedToDagger())
    {
        _dagger(res);
    }
    site = s[1];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p2, site), buff);
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p2, site), buff);
    }
    site = s[2];
    if (site.NeedToDagger())
    {
        _muldag(res, _deviceGetGaugeBCT(byFieldId, p3, site), buff);
    }
    else
    {
        _mul(res, _deviceGetGaugeBCT(byFieldId, p3, site), buff);
    }
    return res;
}

#pragma region kernel



/**
* assert 6 == plaqCount
* assert 4 == plaqLength
* 
* for dir = mu,
* for nu from 0 to 3 and skip mu
* p3[i] stores staple_{mu,nu} with forward (y,x,-y) + backward (-y,x,y)
*   link  staple   staple  staple
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat3(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* p3_1,
    deviceGauge* p3_2,
    deviceGauge* p3_3,
    BYTE byFieldId)
{
    intokernalEDir_NoDir(3);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    deviceGauge* pointers[3] = { p3_1 , p3_2, p3_3 };

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //we store staple_i-forward + staple_i-backward

    //if (uiLinkIndex < 4)
    //{
    //    UINT byIdxStart0 = 6U * elementIdx + uiLinkIndex * plaqCountAllLink;
    //    SSmallInt4 site1 = __deviceSiteIndexToInt4(pCachedStapleIndex[byIdxStart0 + 0].m_uiSiteIndex);
    //    SSmallInt4 site2 = __deviceSiteIndexToInt4(pCachedStapleIndex[byIdxStart0 + 1].m_uiSiteIndex);
    //    SSmallInt4 site3 = __deviceSiteIndexToInt4(pCachedStapleIndex[byIdxStart0 + 2].m_uiSiteIndex);
    //    SSmallInt4 site4 = __deviceSiteIndexToInt4(pCachedStapleIndex[byIdxStart0 + 3].m_uiSiteIndex);
    //    SSmallInt4 site5 = __deviceSiteIndexToInt4(pCachedStapleIndex[byIdxStart0 + 4].m_uiSiteIndex);
    //    SSmallInt4 site6 = __deviceSiteIndexToInt4(pCachedStapleIndex[byIdxStart0 + 5].m_uiSiteIndex);

    //    printf("dir:%d forward element:%d link1(%d,%d,%d,%d)_%d_%d, link2(%d,%d,%d,%d)_%d_%d, link3(%d,%d,%d,%d)_%d_%d\n",
    //        uiLinkIndex,
    //        elementIdx,
    //        site1.x, site1.y, site1.z, site1.w,
    //        pCachedStapleIndex[byIdxStart0 + 0].m_byDir,
    //        pCachedStapleIndex[byIdxStart0 + 0].NeedToDagger(),
    //        site2.x, site2.y, site2.z, site2.w,
    //        pCachedStapleIndex[byIdxStart0 + 1].m_byDir,
    //        pCachedStapleIndex[byIdxStart0 + 1].NeedToDagger(),
    //        site3.x, site3.y, site3.z, site3.w,
    //        pCachedStapleIndex[byIdxStart0 + 2].m_byDir,
    //        pCachedStapleIndex[byIdxStart0 + 2].NeedToDagger()
    //        );

    //    printf("dir:%d backward element:%d link1(%d,%d,%d,%d)_%d_%d, link2(%d,%d,%d,%d)_%d_%d, link3(%d,%d,%d,%d)_%d_%d\n",
    //        uiLinkIndex,
    //        elementIdx,
    //        site4.x, site4.y, site4.z, site4.w,
    //        pCachedStapleIndex[byIdxStart0 + 3].m_byDir,
    //        pCachedStapleIndex[byIdxStart0 + 3].NeedToDagger(),
    //        site5.x, site5.y, site5.z, site5.w,
    //        pCachedStapleIndex[byIdxStart0 + 4].m_byDir,
    //        pCachedStapleIndex[byIdxStart0 + 4].NeedToDagger(),
    //        site6.x, site6.y, site6.z, site6.w,
    //        pCachedStapleIndex[byIdxStart0 + 5].m_byDir,
    //        pCachedStapleIndex[byIdxStart0 + 5].NeedToDagger()
    //    );
    //}

    //elementidx is staple
    //======================= forward ============================
    UINT byIdxStart = 6U * elementIdx + uiLinkIndex * plaqCountAllLink;
    const SIndex& first = pCachedStapleIndex[byIdxStart];
    deviceGauge res(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

    if (first.NeedToDagger())
    {
        _dagger(res);
    }
    #pragma unroll
    for (BYTE j = 1U; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul = _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

        if (nextlink.NeedToDagger())
        {
            _muldag(res, toMul);
        }
        else
        {
            _mul(res, toMul);
        }
    }
    pointers[elementIdx][uiLinkIndex] = res;

    //======================= backward ============================
    byIdxStart = byIdxStart + 3;
    const SIndex& second = pCachedStapleIndex[byIdxStart];
    res = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);

    if (second.NeedToDagger())
    {
        _dagger(res);
    }
    #pragma unroll
    for (BYTE j = 1U; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul = _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

        if (nextlink.NeedToDagger())
        {
            _muldag(res, toMul);
        }
        else
        {
            _mul(res, toMul);
        }
    }
    _add(pointers[elementIdx][uiLinkIndex], res);
}

/**
* information on p3[link]:
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx  
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
* 
* for p5: (six for each link, 6x4 if backward-forward are counted)
*  link   p5_1_1  p5_1_2  p5_2_1  p5_2_2  p5_3_1  p5_3_2
*   0      y-zx    y-tx    z-yx    z-tx    t-yx    t-zx     
*   1      x-zy    x-ty    z-xy    z-ty    t-xy    t-zy
*   2      x-yz    x-tz    y-xz    y-tz    t-xz    t-yz
*   3      x-yt    x-zt    y-xt    y-zt    z-xt    z-yt
* 
* staple for p5: just ommit the middle direction, so for p5_i_j, it is the i-th staple
*  link   p5_1_1  p5_1_2  p5_2_1  p5_2_2  p5_3_1  p5_1_2
*   0       yx      yx      zx      zx      tx    tx
*   1       xy      xy      zy      zy      ty    ty
* 
* p3 for p5: just skip the link
*  link   p5_1_1  p5_1_2  p5_2_1  p5_2_2  p5_3_1  p5_3_2
*   0        1       2       0       2       0       1
*   1        1       2       0       2       0       1
*   ...
* 
* for example, p5_1_1 was:
* 
* -- y -->
*         |
* x       | p3_zx(-z) + tx(-t)
*         v
* <- -y --
*/
__device__ __constant__ constexpr BYTE _fat5_byp3Index[6] = { 1, 2, 0, 2, 0, 1 };

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* p3_1,
    const deviceGauge* p3_2,
    const deviceGauge* p3_3,
    deviceGauge* p5_1_1,
    deviceGauge* p5_1_2,
    deviceGauge* p5_2_1,
    deviceGauge* p5_2_2,
    deviceGauge* p5_3_1,
    deviceGauge* p5_3_2,
    BYTE byFieldId)
{
    intokernalEDir_NoDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    deviceGauge* pointersp5[6] = { p5_1_1 , p5_1_2, p5_2_1, p5_2_2, p5_3_1, p5_3_2 };
    //const BYTE uiStapleIndex[6] = {0, 0, 1, 1, 2, 2};
    //const BYTE byp3Index[6] = {1, 2, 0, 2, 0 , 1};

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //we store staple_i-forward + staple_i-backward
    //p5idx = elementidx
    //use p3_1 for p5_

    const BYTE stapleidx = (elementIdx >> 1U);
    //======================= forward ============================
    UINT byIdxStart = 6U * stapleidx + uiLinkIndex * plaqCountAllLink;
    const SIndex& first = pCachedStapleIndex[byIdxStart];
    deviceGauge res(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

    if (first.NeedToDagger())
    {
        _dagger(res);
    }
    #pragma unroll
    for (BYTE j = 1U; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul = (1U == j) ? 
            pointersp3[_fat5_byp3Index[elementIdx]][_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]
          : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

        if (nextlink.NeedToDagger())
        {
            _muldag(res, toMul);
        }
        else
        {
            _mul(res, toMul);
        }
    }
    pointersp5[elementIdx][uiLinkIndex] = res;

    //======================= backward ============================
    byIdxStart = byIdxStart + 3U;
    const SIndex& second = pCachedStapleIndex[byIdxStart];
    res = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);

    if (second.NeedToDagger())
    {
        _dagger(res);
    }
    #pragma unroll
    for (BYTE j = 1U; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul = (1U == j) ? 
            pointersp3[_fat5_byp3Index[elementIdx]][_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]
          : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

        if (nextlink.NeedToDagger())
        {
            _muldag(res, toMul);
        }
        else
        {
            _mul(res, toMul);
        }
    }
    _add(pointersp5[elementIdx][uiLinkIndex], res);
}

/**
* information on p3:
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
*
* -- y -->
*         |
* x       | p3_yx(-y)
*         v
* <- -y --
* 
* Note: we shall not use this directly, because yx(-y) + (-y)xy is cached instead of yx(-y) alone, but
* 
* ---  y  -->
* <-- -y ---
* |
* x
* |
* v
* ---  y  -->
* <-- -y ---
* 
* = x-link
* so, we can still use the cached staple, but need to substract 6*lepageCoefficient for x it was (y,z,t,-y,-z,-t)
* 
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADLepage(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* p3_1,
    const deviceGauge* p3_2,
    const deviceGauge* p3_3,
    deviceGauge* resoult,
    Real fLepage,
    BYTE byFieldId)
{
    intokernalDir_NoDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };

    deviceGauge lepageres = _makeZero<deviceGauge>();
    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //we store staple_i-forward + staple_i-backward
    #pragma unroll
    for (BYTE p5idx = 0U; p5idx < 3U; ++p5idx)
    {
        //======================= forward ============================
        UINT byIdxStart = 6U * p5idx + uiLinkIndex * plaqCountAllLink;
        const SIndex& first = pCachedStapleIndex[byIdxStart];
        deviceGauge res(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(res);
        }
        #pragma unroll
        for (BYTE j = 1U; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
            const deviceGauge& toMul = (1U == j) ? 
                pointersp3[p5idx][_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]
              : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

            if (nextlink.NeedToDagger())
            {
                _muldag(res, toMul);
            }
            else
            {
                _mul(res, toMul);
            }
        }
        _add(lepageres, res);

        //======================= backward ============================
        byIdxStart = byIdxStart + 3U;
        const SIndex& second = pCachedStapleIndex[byIdxStart];
        res = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);

        if (second.NeedToDagger())
        {
            _dagger(res);
        }
        #pragma unroll
        for (BYTE j = 1U; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
            const deviceGauge& toMul = (1U == j) ? 
                pointersp3[p5idx][_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]
              : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

            if (nextlink.NeedToDagger())
            {
                _muldag(res, toMul);
            }
            else
            {
                _mul(res, toMul);
            }
        }
        _add(lepageres, res);
    }
    _mul(lepageres, fLepage);
    _add(resoult[uiLinkIndex], lepageres);
}

/**
* 
*  link   p5_1_1  p5_1_2  p5_2_1  p5_2_2  p5_3_1  p5_3_2
*   0      y-zx    y-tx    z-yx    z-tx    t-yx    t-zx
*   1      x-zy    x-ty    z-xy    z-ty    t-xy    t-zy
*   2      x-yz    x-tz    y-xz    y-tz    t-xz    t-yz
*   3      x-yt    x-zt    y-xt    y-zt    z-xt    z-yt
* 
* -- y -->
*         |
* x       | p5(ztx(-t)(-z) + p5(tzx(-z)(-t)
*         v
* <- -y --
* 
*  six for each link, 6x8 if backward-forward are counted
*  link     0            1            2            3
*   0       -        y-ztx tzx    z-ytx tyx    t-yzx zyx    3 5   1 4   0 2
*   1   x-zty tzy        -        z-xty txy    t-xzy zxy    3 5   1 4   0 2
*   2   x-ytz tyz    y-xtz txz        -        t-xyz yxz    3 5   1 4   0 2
*   3   x-yzt zyt    y-xzt zxt    z-xyt yxt        -        3 5   1 4   0 2
* 
*/
__device__ __constant__ constexpr BYTE _fat7_byp5Index[6] = { 3, 5, 1, 4, 0, 2 };

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat7(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* p5_1_1,
    const deviceGauge* p5_1_2,
    const deviceGauge* p5_2_1,
    const deviceGauge* p5_2_2,
    const deviceGauge* p5_3_1,
    const deviceGauge* p5_3_2,
    deviceGauge* resoult,
    Real fFat7,
    BYTE byFieldId)
{
    intokernalDir_NoDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp5[6] = { p5_1_1, p5_1_2, p5_2_1, p5_2_2, p5_3_1, p5_3_2 };
    //const BYTE uiStapleIndex[6] = { 0, 0, 1, 1, 2, 2 };
    //const BYTE byp5Index[6] = { 3, 5, 1, 4, 0 , 2 };

    deviceGauge resoult_for_thislink = _makeZero<deviceGauge>();
    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //we store staple_i-forward + staple_i-backward

    #pragma unroll
    for (BYTE p7idx = 0U; p7idx < 6U; ++p7idx)
    {
        const BYTE stapleidx = (p7idx >> 1U);
        //======================= forward ============================
        UINT byIdxStart = 6U * stapleidx + uiLinkIndex * plaqCountAllLink;
        const SIndex& first = pCachedStapleIndex[byIdxStart];
        deviceGauge res(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(res);
        }
        #pragma unroll
        for (BYTE j = 1U; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
            const deviceGauge& toMul = (1U == j) ? 
                pointersp5[_fat7_byp5Index[p7idx]][_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]
              : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

            if (nextlink.NeedToDagger())
            {
                _muldag(res, toMul);
            }
            else
            {
                _mul(res, toMul);
            }
        }
        _add(resoult_for_thislink, res);

        //======================= backward ============================
        byIdxStart = byIdxStart + 3U;
        const SIndex& second = pCachedStapleIndex[byIdxStart];
        res = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);

        if (second.NeedToDagger())
        {
            _dagger(res);
        }
        #pragma unroll
        for (BYTE j = 1U; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
            const deviceGauge& toMul = (1U == j) ? 
                pointersp5[_fat7_byp5Index[p7idx]][_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]
              : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

            if (nextlink.NeedToDagger())
            {
                _muldag(res, toMul);
            }
            else
            {
                _mul(res, toMul);
            }
        }
        _add(resoult_for_thislink, res);
    }

    _mul(resoult_for_thislink, fFat7);
    _add(resoult[uiLinkIndex], resoult_for_thislink);
}

/**
* calculate force also using staple, but replace one of the links with f0
*/
//template<typename deviceGauge>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelASQTADFat3Force(
//    const deviceGauge* __restrict__ pDeviceData,
//    const deviceGauge* __restrict__ pf0,
//    const SIndex* __restrict__ pCachedStapleIndex,
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    BYTE plaqLength, BYTE plaqCountPerLink,
//#endif
//    deviceGauge* forceres,
//    deviceGauge* p3f_1,
//    deviceGauge* p3f_2,
//    deviceGauge* p3f_3,
//    BYTE byFieldId)
//{
//    intokernalEDir_NoDir(3);
//
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    const UINT plaqLengthm1 = plaqLength - 1; //3
//    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
//#endif
//    deviceGauge* pointers[3] = { p3f_1 , p3f_2, p3f_3 };
//
//    if (0 == elementIdx)
//    {
//        forceres[uiLinkIndex] = _makeZero<deviceGauge>();
//    }
//    __syncthreads();
//    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
//        //replace them for 3 times and add them
//        //======================= forward ============================
//    UINT byIdxStart = 6 * elementIdx + uiLinkIndex * plaqCountAllLink;
//    const SIndex& first = pCachedStapleIndex[byIdxStart];
//    deviceGauge res1(_deviceGetGaugeBCDirForceSIndexT(pf0, first));
//    deviceGauge res2(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
//    deviceGauge res3(res2);
//
//    if (first.NeedToDagger())
//    {
//        _dagger(res1);
//        _dagger(res2);
//        _dagger(res3);
//    }
//    #pragma unroll
//    for (BYTE j = 1; j < plaqLengthm1; ++j)
//    {
//        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
//        const deviceGauge& toMul1 = _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);
//        const deviceGauge& toMul2 = (1 == j) ? _deviceGetGaugeBCDirForceSIndexT(pf0, nextlink) : toMul1;
//        const deviceGauge& toMul3 = (1 != j) ? _deviceGetGaugeBCDirForceSIndexT(pf0, nextlink) : toMul1;
//
//        if (nextlink.NeedToDagger())
//        {
//            _muldag(res1, toMul1);
//            _muldag(res2, toMul2);
//            _muldag(res3, toMul3);
//        }
//        else
//        {
//            _mul(res1, toMul1);
//            _mul(res2, toMul2);
//            _mul(res3, toMul3);
//        }
//    }
//
//    _add(res1, res2);
//    _add(res1, res3);
//    #pragma unroll
//    for (BYTE i = 0; i < 3; ++i)
//    {
//        if (i == elementIdx)
//        {
//            _add(forceres[uiLinkIndex], res1);
//        }
//        __syncthreads();
//    }
//    pointers[elementIdx][uiLinkIndex] = res2;
//    __syncthreads();
//    //======================= backward ============================
//    byIdxStart = byIdxStart + 3;
//    const SIndex& second = pCachedStapleIndex[byIdxStart];
//    res1 = _deviceGetGaugeBCDirForceSIndexT(pf0, second);
//    res2 = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);
//    res3 = res2;
//
//    if (second.NeedToDagger())
//    {
//        _dagger(res1);
//        _dagger(res2);
//        _dagger(res3);
//    }
//    #pragma unroll
//    for (BYTE j = 1; j < plaqLengthm1; ++j)
//    {
//        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
//        const deviceGauge& toMul1 = _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);
//        const deviceGauge& toMul2 = (1 == j) ? _deviceGetGaugeBCDirForceSIndexT(pf0, nextlink) : toMul1;
//        const deviceGauge& toMul3 = (1 != j) ? _deviceGetGaugeBCDirForceSIndexT(pf0, nextlink) : toMul1;
//
//        if (nextlink.NeedToDagger())
//        {
//            _muldag(res1, toMul1);
//            _muldag(res2, toMul2);
//            _muldag(res3, toMul3);
//        }
//        else
//        {
//            _mul(res1, toMul1);
//            _mul(res2, toMul2);
//            _mul(res3, toMul3);
//        }
//    }
//
//    _add(res1, res2);
//    _add(res1, res3);
//    #pragma unroll
//    for (BYTE i = 0; i < 3; ++i)
//    {
//        if (i == elementIdx)
//        {
//            _add(forceres[uiLinkIndex], res1);
//        }
//        __syncthreads();
//    }
//    _add(pointers[elementIdx][uiLinkIndex], res2);
//}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat3Force_Optimize(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* forceres,
    deviceGauge* p3f_1,
    deviceGauge* p3f_2,
    deviceGauge* p3f_3,
    BYTE byFieldId)
{
    intokernalEDir_NoDir(3);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    deviceGauge* pointers[3] = { p3f_1 , p3f_2, p3f_3 };

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //replace them for 3 times and add them
    //======================= forward ============================

    const SIndex* forward = pCachedStapleIndex + 6U * elementIdx + uiLinkIndex * plaqCountAllLink;
    deviceGauge res = __deviceOneStaple1(forward, pf0, pDeviceData, pDeviceData, byFieldId);
    deviceGauge res2 = __deviceOneStaple2(forward, pDeviceData, pf0, pDeviceData, byFieldId);
    _add(res, __deviceOneStaple3(forward, pDeviceData, pDeviceData, pf0, byFieldId));

    const SIndex* backward = forward + 3U;
    _add(res, __deviceOneStaple1(backward, pf0, pDeviceData, pDeviceData, byFieldId));
    _add(res2, __deviceOneStaple2(backward, pDeviceData, pf0, pDeviceData, byFieldId));
    _add(res, __deviceOneStaple3(backward, pDeviceData, pDeviceData, pf0, byFieldId));

    pointers[elementIdx][uiLinkIndex] = res2;
    _add(res, res2);

    if (0 == elementIdx)
    {
        forceres[uiLinkIndex] = res;
    }
    #pragma unroll
    for (BYTE i = 1U; i < 3U; ++i)
    {
        __syncthreads();
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res);
        }
    }
}

/**
* for a five staple, it replace each link with f0
* 
* 
* -- f0 -->
*         |
*         |
*         v
* \       \
*  \       \
*   \_______\
* 
* -- f0 -->
* |       
* |       
* |       
* \       \
*  \       \
*   \_______\
* 
* 
* 
* ---f0--->
* |       |
* |       |
* |       v
* \        
*  \        
*   \_______ 
* 
* ---f0--->
* |       |
* |       |
* |       v
* \       \
*  \       \
*   \       \
* 
* ---f0--->
* |       |
* |       |
* |       v
*         \
*          \
*    _______\
* 
* The first two is
* 
* -- f0 -->
*         |
*         |
*         v
* <--p3---

*
* -- f0 -->
* |
* |
* |
* <--p3----
* 
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
* so
* link  f0-1   f0-2   f0-3
*  0     y      z      t
*       1 2    1 2    1 2    y-x-zyz y-x-tyt    z-x-yzy z-x-tzt    t-x-yty t-x-ztz
*        0      0      0     zyz-x-y tyt-x-y    ...
* 
*  1     x      z      t
*       1 2    0 2    0 2    x-y-zxz x-y-txt    z-y-xzx z-y-tzt    t-y-xtx t-y-ztz
*        0      1      1      
* 
*  2     x      y      t
*       0 2    0 2    0 1    x-z-yxy x-z-txt    y-z-xyx y-z-tyt    t-z-xtx t-z-yty
*        1      1      2
* 
*  3     x      y      z
*       0 1    0 1    0 1    x-t-yxy x-t-zxz    y-t-xyx y-t-zyz    z-t-xzx z-t-yzy
*        2      2      2
*
* the last three is
* 
* ---pf3-->
* \
*  \
*   \_______
*
* ---pf3-->
* \       \
*  \       \
*   \       \
*
*
* ---pf3-->
*          \
*           \
*     _______\
* 
* The first one and the third one is as same as above, let's consider the second case
* 
*   link   p3_1    p3_2    p3_3
*    0      yx      zx      tx
*    1      xy      zy      ty
*    2      xz      yz      tz
*    3      xt      yt      zt
* so
* link  f0-1   f0-2   f0-3
*  0     y      z      t
*       1 2    0 2    0 1    y-zxz-y y-txt-y    z-yxy-z z-txt-z    t-yxy-t t-zxz-t
*        0      1      2
*
*  1     x      z      t
*       1 2    0 2    0 1    x-zyz-x x-tyt-x    z-xyx-z z-tyt-z    t-xyx-t t-zyz-t
*        0      1      2
*
*  2     x      y      t
*       1 2    0 2    0 1    x-yzy-x x-tzt-x    y-xzx-y y-tzt-y    t-xzx-t t-yzy-t
*        0      1      2
*
*  3     x      y      z
*       1 2    0 2    0 1    x-yty-x x-ztz-x    y-xtx-y y-ztz-y    z-xtx-z z-yty-z
*        0      1      2
* 
* t1 for x,y,z,t is [000],[011],[112],[222]
* t2 for x,y,z,t is [012]
* 
* p3[t1] - U -  f0
* f0 - U -  p3[t1]
* pf3[t1]  - U - U
* U  - pf3[t2] - U
* U  - U -  pf3[t1]
*/

__device__ __constant__ constexpr BYTE _fat5_t1idx[4][6] = {
    {1, 2, 1, 2, 1, 2},
    {1, 2, 0, 2, 0, 2},
    {0, 2, 0, 2, 0, 1},
    {0, 1, 0, 1, 0, 1}
};
#define _fat5_t2idx _fat5_byp3Index

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5Force(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    const deviceGauge* __restrict__ p3f_1,
    const deviceGauge* __restrict__ p3f_2,
    const deviceGauge* __restrict__ p3f_3,
    deviceGauge* forceres,
    deviceGauge* p5f_1_1,
    deviceGauge* p5f_1_2,
    deviceGauge* p5f_2_1,
    deviceGauge* p5f_2_2,
    deviceGauge* p5f_3_1,
    deviceGauge* p5f_3_2,
    BYTE byFieldId)
{
    intokernalEDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };
    deviceGauge* pointerspf5[6] = { p5f_1_1, p5f_1_2, p5f_2_1, p5f_2_2, p5f_3_1, p5f_3_2 };
    //const BYTE uiStapleIndex[6] = {0, 0, 1, 1, 2, 2};
    //const BYTE byp3Index[6] = {1, 2, 0, 2, 0 , 1};

    //const BYTE t1idx[4][6] = {
    //    {1, 2, 1, 2, 1, 2},
    //    {1, 2, 0, 2, 0, 2},
    //    {0, 2, 0, 2, 0, 1},
    //    {0, 1, 0, 1, 0, 1}
    //};
    //const BYTE t2idx[6] = {1, 2, 0, 2, 0, 1};

    if (0 == elementIdx)
    {
        forceres[uiLinkIndex] = _makeZero<deviceGauge>();
    }
    __syncthreads();

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //we store staple_i-forward + staple_i-backward
    const deviceGauge* p3t1 = pointersp3[_fat5_t1idx[dir][elementIdx]];
    const deviceGauge* pf3t1 = pointerspf3[_fat5_t1idx[dir][elementIdx]];
    const deviceGauge* pf3t2 = pointerspf3[_fat5_t2idx[elementIdx]];
    const BYTE stapleidx = (elementIdx >> 1U);
    //======================= forward ============================
    UINT byIdxStart = 6 * stapleidx + uiLinkIndex * plaqCountAllLink;
    const SIndex& first = pCachedStapleIndex[byIdxStart];
    deviceGauge res1(_deviceGetGaugeBCT(byFieldId, pf0, first));
    deviceGauge res2(_deviceGetGaugeBCT(byFieldId, p3t1, first));
    deviceGauge res3(_deviceGetGaugeBCT(byFieldId, pf3t1, first));
    deviceGauge res4(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
    deviceGauge res5(res4);

    if (first.NeedToDagger())
    {
        _dagger(res1);
        _dagger(res2);
        _dagger(res3);
        _dagger(res4);
        _dagger(res5);
    }
    #pragma unroll
    for (BYTE j = 1; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul1 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink) : _deviceGetGaugeBCT(byFieldId, p3t1, nextlink);
        const deviceGauge& toMul2 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf0, nextlink);
        const deviceGauge& toMul3 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);
        const deviceGauge& toMul4 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pf3t2, nextlink) : toMul3;
        const deviceGauge& toMul5 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf3t1, nextlink);

        if (nextlink.NeedToDagger())
        {
            _muldag(res1, toMul1);
            _muldag(res2, toMul2);
            _muldag(res3, toMul3);
            _muldag(res4, toMul4);
            _muldag(res5, toMul5);
        }
        else
        {
            _mul(res1, toMul1);
            _mul(res2, toMul2);
            _mul(res3, toMul3);
            _mul(res4, toMul4);
            _mul(res5, toMul5);
        }
    }

    _add(res1, res2);
    _add(res1, res3);
    _add(res1, res4);
    _add(res1, res5);
    #pragma unroll
    for (BYTE i = 0; i < 6; ++i)
    {
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res1);
        }
        __syncthreads();
    }
    pointerspf5[elementIdx][uiLinkIndex] = res4;
    __syncthreads();

    //======================= backward ============================
    byIdxStart = byIdxStart + 3;
    const SIndex& second = pCachedStapleIndex[byIdxStart];
    res1 = _deviceGetGaugeBCT(byFieldId, pf0, second);
    res2 = _deviceGetGaugeBCT(byFieldId, p3t1, second);
    res3 = _deviceGetGaugeBCT(byFieldId, pf3t1, second);
    res4 = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);
    res5 = res4;

    if (second.NeedToDagger())
    {
        _dagger(res1);
        _dagger(res2);
        _dagger(res3);
        _dagger(res4);
        _dagger(res5);
    }
    #pragma unroll
    for (BYTE j = 1; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul1 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink) : _deviceGetGaugeBCT(byFieldId, p3t1, nextlink);
        const deviceGauge& toMul2 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf0, nextlink);
        const deviceGauge& toMul3 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);
        const deviceGauge& toMul4 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pf3t2, nextlink) : toMul3;
        const deviceGauge& toMul5 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf3t1, nextlink);

        if (nextlink.NeedToDagger())
        {
            _muldag(res1, toMul1);
            _muldag(res2, toMul2);
            _muldag(res3, toMul3);
            _muldag(res4, toMul4);
            _muldag(res5, toMul5);
        }
        else
        {
            _mul(res1, toMul1);
            _mul(res2, toMul2);
            _mul(res3, toMul3);
            _mul(res4, toMul4);
            _mul(res5, toMul5);
        }
    }

    _add(res1, res2);
    _add(res1, res3);
    _add(res1, res4);
    _add(res1, res5);
    #pragma unroll
    for (BYTE i = 0; i < 6; ++i)
    {
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res1);
        }
        __syncthreads();
    }
    _add(pointerspf5[elementIdx][uiLinkIndex], res4);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat5Force_Optimize(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    const deviceGauge* __restrict__ p3f_1,
    const deviceGauge* __restrict__ p3f_2,
    const deviceGauge* __restrict__ p3f_3,
    deviceGauge* forceres,
    deviceGauge* p5f_1_1,
    deviceGauge* p5f_1_2,
    deviceGauge* p5f_2_1,
    deviceGauge* p5f_2_2,
    deviceGauge* p5f_3_1,
    deviceGauge* p5f_3_2,
    BYTE byFieldId)
{
    intokernalEDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };
    deviceGauge* pointerspf5[6] = { p5f_1_1, p5f_1_2, p5f_2_1, p5f_2_2, p5f_3_1, p5f_3_2 };

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    //we store staple_i-forward + staple_i-backward
    const deviceGauge* p3t1 = pointersp3[_fat5_t1idx[dir][elementIdx]];
    const deviceGauge* pf3t1 = pointerspf3[_fat5_t1idx[dir][elementIdx]];
    const deviceGauge* pf3t2 = pointerspf3[_fat5_t2idx[elementIdx]];
    const BYTE stapleidx = (elementIdx >> 1U);

    //======================= forward ============================
    const SIndex* forward = pCachedStapleIndex + 6U * stapleidx + uiLinkIndex * plaqCountAllLink;

    deviceGauge res = __deviceOneStaple1(forward, pf0, pDeviceData, p3t1, byFieldId);
    _add(res, __deviceOneStaple3(forward, p3t1, pDeviceData, pf0, byFieldId));
    _add(res, __deviceOneStaple1(forward, pf3t1, pDeviceData, pDeviceData, byFieldId));
    deviceGauge res4 = __deviceOneStaple2(forward, pDeviceData, pf3t2, pDeviceData, byFieldId);
    _add(res, __deviceOneStaple3(forward, pDeviceData, pDeviceData, pf3t1, byFieldId));

    const SIndex* backward = forward + 3U;
    _add(res, __deviceOneStaple1(backward, pf0, pDeviceData, p3t1, byFieldId));
    _add(res, __deviceOneStaple3(backward, p3t1, pDeviceData, pf0, byFieldId));
    _add(res, __deviceOneStaple1(backward, pf3t1, pDeviceData, pDeviceData, byFieldId));
    _add(res4, __deviceOneStaple2(backward, pDeviceData, pf3t2, pDeviceData, byFieldId));
    _add(res, __deviceOneStaple3(backward, pDeviceData, pDeviceData, pf3t1, byFieldId));

    pointerspf5[elementIdx][uiLinkIndex] = res4;
    _add(res, res4);
    if (0 == elementIdx)
    {
        forceres[uiLinkIndex] = res;
    }
    #pragma unroll
    for (BYTE i = 1U; i < 6U; ++i)
    {
        __syncthreads();
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res);
        }
    }
    
}

/**
* 
* 30 for each, if me is x, y can be pm
* 
* --f0--
* |    |
* |    |
* |    |
* 
* ----------
* |        |
* |        f0
* |    ____|
* 
* ---------
* |        |
* |        f0
* |____    |
* 
* ----------
* |        |
* f0       |
* |    ____|
*
* ---------
* |        |
* f0       |
* |____    |
* 
* there is no way to use the cache (except for the first one), unless we cache plus and minus p3 seperately ...
* for cases 1-4, the dir of p3,pf3 is dir1-dir2, for the last case, is dir2-dir1
* 
* we do the last four first
* it is either
* x y (x) x y with - + + + -    1,    4
* y x (x) y x with + + + - - 0,    3,
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADLepageForceA(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalEDir(12U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    BYTE nu = elementIdx >> 2U;
    nu = (dir + nu + 1U) & 3U;

    //0 for yxxyx, 1 for xyxxy
    const BYTE caseidx = (elementIdx&1U);
    //0,1,3,4
    const BYTE f0idx = caseidx + 3U * ((elementIdx&2U) >> 1U);
    //change sign
    //when caseidx = 0, change 0 and 3
    //when caseidx = 1, change 1 and 4
    //so
    // caseidx  0  1
    //    0    -1  1 
    //    1     1 -1
    //    2     1  1
    //    3    -1  1
    //    4     1 -1

    const SCHAR sa = static_cast<SCHAR>((caseidx << 1) - 1);
    const SCHAR sb = -sa;

    const SCHAR sign1 = static_cast<SCHAR>(1 - (static_cast<INT>(caseidx) << 1));
    const SCHAR sign4 = -sign1;

    SCHAR munu[2] = { static_cast<SCHAR>(dir + 1), static_cast<SCHAR>(nu + 1) };

    //+y
    SCHAR path[5] = { 
        static_cast<SCHAR>(sign1 * munu[1 - caseidx]),
        static_cast<SCHAR>(munu[caseidx]),
        static_cast<SCHAR>(munu[0]),
        static_cast<SCHAR>(sign4 * munu[1 - caseidx]),
        static_cast<SCHAR>(-munu[caseidx])
    };

    //if (0 == uiSiteIndex)
    //{
    //    printf("elementidx=%d, dir=%d, nu=%d, caseidx=%d, f0idx=%d, path=[%d,%d,%d,%d,%d]\n", elementIdx, dir, nu, caseidx, f0idx, path[0], path[1], path[2], path[3], path[4]);
    //}

    deviceGauge force = _deviceLinkTTwoField(pDeviceData, pf0, sSite4, 5, byFieldId, f0idx, path);

    //-y
    path[0] = sa * path[0];
    path[1] = sb * path[1];
    path[3] = sa * path[3];
    path[4] = sb * path[4];

    //if (0 == uiSiteIndex)
    //{
    //    printf("elementidx=%d, dir=%d, nu=%d, caseidx=%d, f0idx=%d, path=[%d,%d,%d,%d,%d]\n", elementIdx, dir, nu, caseidx, f0idx, path[0], path[1], path[2], path[3], path[4]);
    //}

    _add(force, _deviceLinkTTwoField(pDeviceData, pf0, sSite4, 5, byFieldId, f0idx, path));

    if (0U == elementIdx)
    {
        forceres[uiLinkIndex] = force;
    }
    __syncthreads();
    #pragma unroll
    for (BYTE i = 1U; i < 12U; ++i)
    {
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], force);
        }
        __syncthreads();
    }
}

/**
* --f0--
* |    |
* |    |
* |    |
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADLepageForceB(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalEDir(3U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    const SCHAR mu = dir + 1;
    const SCHAR nu = ((dir + elementIdx + 1)&3) + 1;

    SCHAR path[5] = {
        static_cast<SCHAR>(nu),
        static_cast<SCHAR>(nu),
        static_cast<SCHAR>(mu),
        static_cast<SCHAR>(-nu),
        static_cast<SCHAR>(-nu)
    };

    deviceGauge force = _deviceLinkTTwoField(pDeviceData, pf0, sSite4, 5, byFieldId, 2, path);
    path[0] = -path[0];
    path[1] = -path[1];
    path[3] = -path[3];
    path[4] = -path[4];
    _add(force, _deviceLinkTTwoField(pDeviceData, pf0, sSite4, 5, byFieldId, 2, path));

    #pragma unroll
    for (BYTE i = 0U; i < 3U; ++i)
    {
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], force);
        }
        __syncthreads();
    }
}

//template<typename deviceGauge>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelASQTADLepageForce(
//    const deviceGauge* __restrict__ pDeviceData,
//    const deviceGauge* __restrict__ pf0,
//    const SIndex* __restrict__ pCachedStapleIndex,
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    BYTE plaqLength, BYTE plaqCountPerLink,
//#endif
//    deviceGauge* forceres,
//    BYTE byFieldId)
//{
//    intokernalEDir(3);
//
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    const UINT plaqLengthm1 = plaqLength - 1; //3
//    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
//#endif
//    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
//    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };
//
//    //const BYTE t1idx[4][3] = {
//    //    {0, 0, 0},
//    //    {0, 1, 1},
//    //    {1, 1, 2},
//    //    {2, 2, 2}
//    //};
//
//    if (0 == elementIdx)
//    {
//        forceres[uiLinkIndex] = _makeZero<deviceGauge>();
//    }
//    __syncthreads();
//
//    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
//    //we store staple_i-forward + staple_i-backward
//    const deviceGauge* p3t1 = pointersp3[_lapage_t1idx[dir][elementIdx]];
//    const deviceGauge* pf3t1 = pointerspf3[_lapage_t1idx[dir][elementIdx]];
//    const deviceGauge* pf3t2 = pointerspf3[elementIdx];
//
//    //======================= forward ============================
//    UINT byIdxStart = 6 * elementIdx + uiLinkIndex * plaqCountAllLink;
//    const SIndex& first = pCachedStapleIndex[byIdxStart];
//    deviceGauge res1(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
//    deviceGauge res2(_deviceGetGaugeBCT(byFieldId, pf3t1, first));
//    deviceGauge res3(_deviceGetGaugeBCT(byFieldId, pf0, first));
//    deviceGauge res4(_deviceGetGaugeBCT(byFieldId, p3t1, first));
//    deviceGauge res5(res1);
//
//    if (first.NeedToDagger())
//    {
//        _dagger(res1);
//        _dagger(res2);
//        _dagger(res3);
//        _dagger(res4);
//        _dagger(res5);
//    }
//    #pragma unroll
//    for (BYTE j = 1; j < plaqLengthm1; ++j)
//    {
//        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
//        const deviceGauge& toMul1 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink) : _deviceGetGaugeBCT(byFieldId, pf3t1, nextlink);
//        const deviceGauge& toMul2 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);
//        const deviceGauge& toMul3 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p3t1, nextlink);
//        const deviceGauge& toMul4 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf0, nextlink);
//        const deviceGauge& toMul5 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pf3t2, nextlink) : toMul2;
//
//        if (nextlink.NeedToDagger())
//        {
//            _muldag(res1, toMul1);
//            _muldag(res2, toMul2);
//            _muldag(res3, toMul3);
//            _muldag(res4, toMul4);
//            _muldag(res5, toMul5);
//        }
//        else
//        {
//            _mul(res1, toMul1);
//            _mul(res2, toMul2);
//            _mul(res3, toMul3);
//            _mul(res4, toMul4);
//            _mul(res5, toMul5);
//        }
//    }
//
//    _add(res1, res2);
//    _add(res1, res3);
//    _add(res1, res4);
//    _add(res1, res5);
//    #pragma unroll
//    for (BYTE i = 0; i < 3; ++i)
//    {
//        if (i == elementIdx)
//        {
//            _add(forceres[uiLinkIndex], res1);
//        }
//        __syncthreads();
//    }
//
//    //======================= backward ============================
//    byIdxStart = byIdxStart + 3;
//    const SIndex& second = pCachedStapleIndex[byIdxStart];
//    res1 = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);
//    res2 = _deviceGetGaugeBCT(byFieldId, pf3t1, second);
//    res3 = _deviceGetGaugeBCT(byFieldId, pf0, second);
//    res4 = _deviceGetGaugeBCT(byFieldId, p3t1, second);
//    res5 = res1;
//
//    if (second.NeedToDagger())
//    {
//        _dagger(res1);
//        _dagger(res2);
//        _dagger(res3);
//        _dagger(res4);
//        _dagger(res5);
//    }
//    #pragma unroll
//    for (BYTE j = 1; j < plaqLengthm1; ++j)
//    {
//        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
//        const deviceGauge& toMul1 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink) : _deviceGetGaugeBCT(byFieldId, pf3t1, nextlink);
//        const deviceGauge& toMul2 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);
//        const deviceGauge& toMul3 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p3t1, nextlink);
//        const deviceGauge& toMul4 = (1 == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf0, nextlink);
//        const deviceGauge& toMul5 = (1 == j) ? _deviceGetGaugeBCT(byFieldId, pf3t2, nextlink) : toMul2;
//
//        if (nextlink.NeedToDagger())
//        {
//            _muldag(res1, toMul1);
//            _muldag(res2, toMul2);
//            _muldag(res3, toMul3);
//            _muldag(res4, toMul4);
//            _muldag(res5, toMul5);
//        }
//        else
//        {
//            _mul(res1, toMul1);
//            _mul(res2, toMul2);
//            _mul(res3, toMul3);
//            _mul(res4, toMul4);
//            _mul(res5, toMul5);
//        }
//    }
//    _add(res1, res2);
//    _add(res1, res3);
//    _add(res1, res4);
//    _add(res1, res5);
//    #pragma unroll
//    for (BYTE i = 0; i < 3; ++i)
//    {
//        if (i == elementIdx)
//        {
//            _add(forceres[uiLinkIndex], res1);
//        }
//        __syncthreads();
//    }
//}


//template<typename deviceGauge>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelASQTADLepageForce_Optimize(
//    const deviceGauge* __restrict__ pDeviceData,
//    const deviceGauge* __restrict__ pf0,
//    const SIndex* __restrict__ pCachedStapleIndex,
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    BYTE plaqLength, BYTE plaqCountPerLink,
//#endif
//    const deviceGauge* __restrict__ p3_1,
//    const deviceGauge* __restrict__ p3_2,
//    const deviceGauge* __restrict__ p3_3,
//    const deviceGauge* __restrict__ p3f_1,
//    const deviceGauge* __restrict__ p3f_2,
//    const deviceGauge* __restrict__ p3f_3,
//    deviceGauge* forceres,
//    BYTE byFieldId)
//{
//    intokernalEDir(3);
//
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    const UINT plaqLengthm1 = plaqLength - 1; //3
//    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
//#endif
//    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
//    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };
//
//    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
//    //we store staple_i-forward + staple_i-backward
//    const deviceGauge* p3t1 = pointersp3[_lapage_t1idx[dir][elementIdx]];
//    const deviceGauge* pf3t1 = pointerspf3[_lapage_t1idx[dir][elementIdx]];
//    const deviceGauge* pf3t2 = pointerspf3[elementIdx];
//
//    //======================= forward ============================
//
//    const SIndex* forward = pCachedStapleIndex + 6U * elementIdx + uiLinkIndex * plaqCountAllLink;
//    deviceGauge res = __deviceOneStaple3(forward, pDeviceData, pDeviceData, pf3t1, byFieldId);
//    _add(res, __deviceOneStaple1(forward, pf3t1, pDeviceData, pDeviceData, byFieldId));
//    _add(res, __deviceOneStaple1(forward, pf0, pDeviceData, p3t1, byFieldId));
//    _add(res, __deviceOneStaple3(forward, p3t1, pDeviceData, pf0, byFieldId));
//    _add(res, __deviceOneStaple2(forward, pDeviceData, pf3t2, pDeviceData, byFieldId));
//
//    const SIndex* backward = forward + 3U;
//    _add(res, __deviceOneStaple3(backward, pDeviceData, pDeviceData, pf3t1, byFieldId));
//    _add(res, __deviceOneStaple1(backward, pf3t1, pDeviceData, pDeviceData, byFieldId));
//    _add(res, __deviceOneStaple1(backward, pf0, pDeviceData, p3t1, byFieldId));
//    _add(res, __deviceOneStaple3(backward, p3t1, pDeviceData, pf0, byFieldId));
//    _add(res, __deviceOneStaple2(backward, pDeviceData, pf3t2, pDeviceData, byFieldId));
//
//    if (0 == elementIdx)
//    {
//        forceres[uiLinkIndex] = res;
//    }
//    __syncthreads();
//    #pragma unroll
//    for (BYTE i = 1U; i < 3U; ++i)
//    {
//        __syncthreads();
//        if (i == elementIdx)
//        {
//            _add(forceres[uiLinkIndex], res);
//        }
//    }
//}

/**
* 
* ---f0--->
*         |
*         |
*         |
*  \      \
*   \      \
*    \      \
*     |      |
*     |      |
*     |______|
* 
* ---f0--->
* |       
* |        
* |        
*  \      \
*   \      \
*    \      \
*     |      |
*     |      |
*     |______|
*      
*
* ---f0--->
* |       |
* |       |
* |       |
*  \       
*   \       
*    \       
*     |      |
*     |      |
*     |______|
*
* ---f0--->
* |       |
* |       |
* |       |
*         \
*          \
*           \
*     |      |
*     |      |
*     |______|
* 
* ---f0--->
* |       |
* |       |
* |       |
*  \      \
*   \      \
*    \      \
*     |       
*     |       
*     |______ 
* 
* ---f0--->
* |       |
* |       |
* |       |
*  \      \
*   \      \
*    \      \
*            |
*            |
*      ______|
* 
* 
* ---f0--->
* |       |
* |       |
* |       |
*  \      \
*   \      \
*    \      \
*     |      |
*     |      |
*     |      |
* 
* The first two is:
* 
* ---f0--->   ---f0--->
*         |   |
*         |   |
*  <--p5--|   |<--p5---
* 
* me: dir1
* f0: dir2
* p5: dir3-dir4-dir2, dir4-dir3-dir2
* 
*  link   p5_1_1  p5_1_2  p5_2_1  p5_2_2  p5_3_1  p5_3_2
*   0      y-zx    y-tx    z-yx    z-tx    t-yx    t-zx
*   1      x-zy    x-ty    z-xy    z-ty    t-xy    t-zy
*   2      x-yz    x-tz    y-xz    y-tz    t-xz    t-yz
*   3      x-yt    x-zt    y-xt    y-zt    z-xt    z-yt
* 
*  d1 d2    0             1             2             3
*   0       -         xy-zty tzy    xz-ytz tyz    xt-yzt zyt     3 5   3 5   3 5
*   1   yx-ztx tzx        -         yz-xtz txz    yt-xzt zxt     3 5   1 4   1 4
*   2   zx-ytx tyx    zy-xty txy        -         zt-xyt yxt     
*   3   tx-yzx zyx    ty-xzy zxy    tz-xyz yxz        -     
* 
* The second two is:
* 
* ---pf3-->    ---pf3-->
* \                     \
*  \                     \
*   \__p3____       __p3__\
*
*  me:  dir1
*  p3:  dir3-dir2   dir4-dir2
*  pf3: dir4-dir2   dir3-dir2
* 
* then is
* 
* --pf5-->   --pf5-->
* |                 |
* |                 |
* |<------   <-------
*
* me: dir1
* pf5: dir3-dir4-dir2, same as the first case
* 
* The last one is:
* 
* ---pf5-->
* |       |
* |       |
* |       |
* 
* me: dir1
* pf5: dir4-dir3-dir1 dir3-dir4-dir1
* 
* same as fat7,  3 5   1 4   0 2
* 
*/

#define _fat7_t1idx _fat5_t1idx

__device__ __constant__ constexpr BYTE _fat7_t2idx[4][6] = {
    {2, 1, 2, 1, 2, 1},
    {2, 1, 2, 0, 2, 0},
    {2, 0, 2, 0, 1, 0},
    {1, 0, 1, 0, 1, 0}
};

__device__ __constant__ constexpr BYTE _fat7_t3idx[4][6] = {
    {3, 5, 3, 5, 3, 5},
    {3, 5, 1, 4, 1, 4},
    {1, 4, 1, 4, 0, 2},
    {0, 2, 0, 2, 0, 2}
};

#define _fat7_t4idx _fat7_byp5Index

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat7Force(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    const deviceGauge* __restrict__ p3f_1,
    const deviceGauge* __restrict__ p3f_2,
    const deviceGauge* __restrict__ p3f_3,
    const deviceGauge* __restrict__ p5_1_1,
    const deviceGauge* __restrict__ p5_1_2,
    const deviceGauge* __restrict__ p5_2_1,
    const deviceGauge* __restrict__ p5_2_2,
    const deviceGauge* __restrict__ p5_3_1,
    const deviceGauge* __restrict__ p5_3_2,
    const deviceGauge* __restrict__ p5f_1_1,
    const deviceGauge* __restrict__ p5f_1_2,
    const deviceGauge* __restrict__ p5f_2_1,
    const deviceGauge* __restrict__ p5f_2_2,
    const deviceGauge* __restrict__ p5f_3_1,
    const deviceGauge* __restrict__ p5f_3_2,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalEDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };
    const deviceGauge* pointersp5[6] = { p5_1_1, p5_1_2, p5_2_1, p5_2_2, p5_3_1, p5_3_2 };
    const deviceGauge* pointerspf5[6] = { p5f_1_1, p5f_1_2, p5f_2_1, p5f_2_2, p5f_3_1, p5f_3_2 };

    //constexpr BYTE uiStapleIndex[6] = { 0, 0, 1, 1, 2, 2 };

    //constexpr BYTE t1idx[4][6] = {
    //    {1, 2, 1, 2, 1, 2},
    //    {1, 2, 0, 2, 0, 2},
    //    {0, 2, 0, 2, 0, 1},
    //    {0, 1, 0, 1, 0, 1}
    //};
    //constexpr BYTE t2idx[4][6] = {
    //    {2, 1, 2, 1, 2, 1},
    //    {2, 1, 2, 0, 2, 0},
    //    {2, 0, 2, 0, 1, 0},
    //    {1, 0, 1, 0, 1, 0}
    //};

    //constexpr BYTE t3idx[4][6] = {
    //    {3, 5, 3, 5, 3, 5},
    //    {3, 5, 1, 4, 1, 4},
    //    {1, 4, 1, 4, 0, 2},
    //    {0, 2, 0, 2, 0, 2}
    //};
    //constexpr BYTE t4idx[6] = { 3, 5, 1, 4, 0, 2 };

    if (0 == elementIdx)
    {
        forceres[uiLinkIndex] = _makeZero<deviceGauge>();
    }
    __syncthreads();

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    const deviceGauge* p3t1 = pointersp3[_fat7_t1idx[dir][elementIdx]];
    const deviceGauge* p3ft2 = pointerspf3[_fat7_t2idx[dir][elementIdx]];

    const deviceGauge* p5t3 = pointersp5[_fat7_t3idx[dir][elementIdx]];
    const deviceGauge* p5ft3 = pointerspf5[_fat7_t3idx[dir][elementIdx]];
    const deviceGauge* p5ft4 = pointerspf5[_fat7_t4idx[elementIdx]];
    const BYTE stapleidx = (elementIdx >> 1U);

    //======================= forward ============================
    UINT byIdxStart = 6U * stapleidx + uiLinkIndex * plaqCountAllLink;
    const SIndex& first = pCachedStapleIndex[byIdxStart];
    deviceGauge res1(_deviceGetGaugeBCT(byFieldId, pf0, first));
    deviceGauge res2(_deviceGetGaugeBCT(byFieldId, p5t3, first));

    deviceGauge res3(_deviceGetGaugeBCT(byFieldId, p3t1, first));
    deviceGauge res4(_deviceGetGaugeBCT(byFieldId, p3ft2, first));

    deviceGauge res5(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
    deviceGauge res6(_deviceGetGaugeBCT(byFieldId, p5ft3, first));

    deviceGauge res7(res5);

    if (first.NeedToDagger())
    {
        _dagger(res1);
        _dagger(res2);
        _dagger(res3);
        _dagger(res4);
        _dagger(res5);
        _dagger(res6);
        _dagger(res7);
    }
    #pragma unroll
    for (BYTE j = 1U; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul1 = (1U == j) ? _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink) : _deviceGetGaugeBCT(byFieldId, p5t3, nextlink);
        const deviceGauge& toMul2 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf0, nextlink);

        const deviceGauge& toMul3 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p3ft2, nextlink);
        const deviceGauge& toMul4 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p3t1, nextlink);

        const deviceGauge& toMul5 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p5ft3, nextlink);
        const deviceGauge& toMul6 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

        const deviceGauge& toMul7 = (1U == j) ? _deviceGetGaugeBCT(byFieldId, p5ft4, nextlink) : toMul6;

        if (nextlink.NeedToDagger())
        {
            _muldag(res1, toMul1);
            _muldag(res2, toMul2);
            _muldag(res3, toMul3);
            _muldag(res4, toMul4);
            _muldag(res5, toMul5);
            _muldag(res6, toMul6);
            _muldag(res7, toMul7);
        }
        else
        {
            _mul(res1, toMul1);
            _mul(res2, toMul2);
            _mul(res3, toMul3);
            _mul(res4, toMul4);
            _mul(res5, toMul5);
            _mul(res6, toMul6);
            _mul(res7, toMul7);
        }
    }
    _add(res1, res2);
    _add(res1, res3);
    _add(res1, res4);
    _add(res1, res5);
    _add(res1, res6);
    _add(res1, res7);
    #pragma unroll
    for (BYTE i = 0U; i < 6U; ++i)
    {
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res1);
        }
        __syncthreads();
    }

    //======================= backward ============================
    byIdxStart = byIdxStart + 3U;
    const SIndex& second = pCachedStapleIndex[byIdxStart];
    res1 = _deviceGetGaugeBCT(byFieldId, pf0, second);
    res2 = _deviceGetGaugeBCT(byFieldId, p5t3, second);

    res3 = _deviceGetGaugeBCT(byFieldId, p3t1, second);
    res4 = _deviceGetGaugeBCT(byFieldId, p3ft2, second);

    res5 = _deviceGetGaugeBCT(byFieldId, pDeviceData, second);
    res6 = _deviceGetGaugeBCT(byFieldId, p5ft3, second);

    res7 = res5;

    if (second.NeedToDagger())
    {
        _dagger(res1);
        _dagger(res2);
        _dagger(res3);
        _dagger(res4);
        _dagger(res5);
        _dagger(res6);
        _dagger(res7);
    }
    #pragma unroll
    for (BYTE j = 1U; j < plaqLengthm1; ++j)
    {
        const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
        const deviceGauge& toMul1 = (1U == j) ? _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink) : _deviceGetGaugeBCT(byFieldId, p5t3, nextlink);
        const deviceGauge& toMul2 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pf0, nextlink);

        const deviceGauge& toMul3 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p3ft2, nextlink);
        const deviceGauge& toMul4 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p3t1, nextlink);

        const deviceGauge& toMul5 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, p5ft3, nextlink);
        const deviceGauge& toMul6 = (1U == j) ? toMul1 : _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

        const deviceGauge& toMul7 = (1U == j) ? _deviceGetGaugeBCT(byFieldId, p5ft4, nextlink) : toMul6;

        if (nextlink.NeedToDagger())
        {
            _muldag(res1, toMul1);
            _muldag(res2, toMul2);
            _muldag(res3, toMul3);
            _muldag(res4, toMul4);
            _muldag(res5, toMul5);
            _muldag(res6, toMul6);
            _muldag(res7, toMul7);
        }
        else
        {
            _mul(res1, toMul1);
            _mul(res2, toMul2);
            _mul(res3, toMul3);
            _mul(res4, toMul4);
            _mul(res5, toMul5);
            _mul(res6, toMul6);
            _mul(res7, toMul7);
        }
    }
    _add(res1, res2);
    _add(res1, res3);
    _add(res1, res4);
    _add(res1, res5);
    _add(res1, res6);
    _add(res1, res7);
    #pragma unroll
    for (BYTE i = 0U; i < 6U; ++i)
    {
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res1);
        }
        __syncthreads();
    }
}


template<typename deviceGauge, INT N>
__global__ void _CLG_LAUNCH_BOUND
_kernelASQTADFat7Force_Optimize(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pf0,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    const deviceGauge* __restrict__ p3_1,
    const deviceGauge* __restrict__ p3_2,
    const deviceGauge* __restrict__ p3_3,
    const deviceGauge* __restrict__ p3f_1,
    const deviceGauge* __restrict__ p3f_2,
    const deviceGauge* __restrict__ p3f_3,
    const deviceGauge* __restrict__ p5_1_1,
    const deviceGauge* __restrict__ p5_1_2,
    const deviceGauge* __restrict__ p5_2_1,
    const deviceGauge* __restrict__ p5_2_2,
    const deviceGauge* __restrict__ p5_3_1,
    const deviceGauge* __restrict__ p5_3_2,
    const deviceGauge* __restrict__ p5f_1_1,
    const deviceGauge* __restrict__ p5f_1_2,
    const deviceGauge* __restrict__ p5f_2_1,
    const deviceGauge* __restrict__ p5f_2_2,
    const deviceGauge* __restrict__ p5f_3_1,
    const deviceGauge* __restrict__ p5f_3_2,
    deviceGauge* forceres,
    BYTE byFieldId)
{
    intokernalEDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const deviceGauge* pointersp3[3] = { p3_1, p3_2, p3_3 };
    const deviceGauge* pointerspf3[3] = { p3f_1, p3f_2, p3f_3 };
    const deviceGauge* pointersp5[6] = { p5_1_1, p5_1_2, p5_2_1, p5_2_2, p5_3_1, p5_3_2 };
    const deviceGauge* pointerspf5[6] = { p5f_1_1, p5f_1_2, p5f_2_1, p5f_2_2, p5f_3_1, p5f_3_2 };

    //for each dir, 6 staples are cached, for dir2 = 0 to 3, skip dir1, and forward for one, backward for the other
    const deviceGauge* p3t1 = pointersp3[_fat7_t1idx[dir][elementIdx]];
    const deviceGauge* p3ft2 = pointerspf3[_fat7_t2idx[dir][elementIdx]];

    const deviceGauge* p5t3 = pointersp5[_fat7_t3idx[dir][elementIdx]];
    const deviceGauge* p5ft3 = pointerspf5[_fat7_t3idx[dir][elementIdx]];
    const deviceGauge* p5ft4 = pointerspf5[_fat7_t4idx[elementIdx]];
    const BYTE stapleidx = (elementIdx >> 1U);

    //======================= forward ============================
    const SIndex* forward = pCachedStapleIndex + 6U * stapleidx + uiLinkIndex * plaqCountAllLink;
    CLGComplex buff[N];
    deviceGauge res = __deviceOneStapleBuff(forward, pf0, pDeviceData, p5t3, byFieldId, buff);
    _add(res, __deviceOneStapleBuff(forward, p5t3, pDeviceData, pf0, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(forward, p3t1, pDeviceData, p3ft2, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(forward, p3ft2, pDeviceData, p3t1, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(forward, pDeviceData, pDeviceData, p5ft3, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(forward, p5ft3, pDeviceData, pDeviceData, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(forward, pDeviceData, p5ft4, pDeviceData, byFieldId, buff));

    const SIndex* backward = forward + 3U;
    _add(res, __deviceOneStapleBuff(backward, pf0, pDeviceData, p5t3, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(backward, p5t3, pDeviceData, pf0, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(backward, p3t1, pDeviceData, p3ft2, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(backward, p3ft2, pDeviceData, p3t1, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(backward, pDeviceData, pDeviceData, p5ft3, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(backward, p5ft3, pDeviceData, pDeviceData, byFieldId, buff));
    _add(res, __deviceOneStapleBuff(backward, pDeviceData, p5ft4, pDeviceData, byFieldId, buff));

    if (0U == elementIdx)
    {
        forceres[uiLinkIndex] = res;
    }
    #pragma unroll
    for (BYTE i = 1U; i < 6U; ++i)
    {
        __syncthreads();
        if (i == elementIdx)
        {
            _add(forceres[uiLinkIndex], res);
        }
    }
}

//template<typename deviceGauge>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelASQTADFatForce_Optimize(
//    const SIndex* __restrict__ pCachedStapleIndex,
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    UINT plaqLength, UINT plaqCountPerLink,
//#endif
//    BYTE dir, BYTE idxoffset,
//    const deviceGauge* __restrict__ p1,
//    const deviceGauge* __restrict__ p2,
//    const deviceGauge* __restrict__ p3,
//    deviceGauge* forceres,
//    BYTE byFieldId)
//{
//    intokernal;
//#if !_CLG_ASSUME_SQUARE_LATTICE
//    const UINT plaqLengthm1 = plaqLength - 1; //3
//    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
//#endif
//    const UINT uiLinkIndex = (uiSiteIndex << 2U) | dir;
//    const SIndex* forward = pCachedStapleIndex + uiLinkIndex * plaqCountAllLink + idxoffset;
//    _add(forceres[uiLinkIndex], __deviceOneStaple(forward, p1, p2, p3, byFieldId));
//}


/**
* 
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectRationalApproximation(
    deviceGauge* pDeviceData,
    deviceGauge** pForce,
    CLGComplex* detres,
    const Real* __restrict__ pCoefficients,
    UINT order,
    UBOOL bSU3,
    Real fOneOverNc)
{
    intokernalDir_NoDir;
    //c0 = pCoefficients[0];
    //a[i] = pCoefficients[1 + i];
    //b[i] = pCoefficients[1 + i + order];

    
    const deviceGauge udu = _dagmulC(pDeviceData[uiLinkIndex], pDeviceData[uiLinkIndex]); //Ud U
    deviceGauge res = _makeZero<deviceGauge>();
    for (UINT i = 0; i < order; ++i)
    {
        deviceGauge udub = _addC(udu, pCoefficients[1 + i + order]); //Ud U + bi
        udub = _rcpC(udub); // udub = 1/(Ud U+bi)
        pForce[i][uiLinkIndex] = udub; // udub = 1/(Ud U+bi)
        _mul(udub, pCoefficients[1 + i]); //udub = ai/(Ud U+bi)
        _add(res, udub);
    }
    
    //res = sum ai/(Ud U+bi)
    _add(res, pCoefficients[0]); //res = c0 + sum ai/(Ud U+bi)
    pForce[order][uiLinkIndex] = res;
    res = _mulC(pDeviceData[uiLinkIndex], res); //U [c0 + sum ai/(Ud U+bi)]
    
    //deviceGauge testres = res;
    //_muldag(testres, res);
    //_sub(testres, F(1.0));
    //Real fDelta = _lensq(testres);
    //
    //if (fDelta > F(0.0001))
    //{
    //    printf("delta %f\n", fDelta);
    //    _print(pDeviceData[uiLinkIndex]);
    //}

    if (bSU3)
    {
        //if (0 == linkIndex)
        //{
        //    printf("==== raw u =====");
        //    _print(pDeviceData[linkIndex]);
        //}
        CLGComplex det = _detv(res);

        //atan2 is used in power, so it is automatically -pi/3 to pi/3
        det = __cuCpowerf(det, -fOneOverNc);
        detres[uiLinkIndex] = det;
        //if (2 == uiLinkIndex)
        //{
        //    printf("==== u3 =====");
        //    _print(res);
        //    printf("==== det =====");
        //    _print(det);
        //}
        _mul(res, det);
        //if (2 == uiLinkIndex)
        //{
        //    printf("==== su3 =====");
        //    _print(res);
        //}
    }
    
    pDeviceData[uiLinkIndex] = res;
    
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectRationalApproximationForce(
    const deviceGauge* __restrict__ pGaugeNotProjected,
    deviceGauge* pf0,
    const deviceGauge* const * __restrict__ pForce,
    const Real* __restrict__ pCoefficients,
    UINT order)
{
    intokernalDir_NoDir;
    //c0 = pCoefficients[0];
    //a[i] = pCoefficients[1 + i];
    //b[i] = pCoefficients[1 + i + order];

    deviceGauge udf1 = _dagmulC(pGaugeNotProjected[uiLinkIndex], pf0[uiLinkIndex]); //Ud f0
    _add(udf1, _daggerC(udf1));
    deviceGauge res1 = _makeZero<deviceGauge>();
    for (UINT i = 0; i < order; ++i)
    {
        //const deviceGauge& udub = pForce[i][uiLinkIndex]; //udub = 1/(Ud U+bi)
        //deviceGauge udub1 = _mulC(udub, _mulC(udf1, udub)); //udub = 1/(Ud U+bi) Ud f0 1/(Ud U+bi)
        //_mul(udub1, pCoefficients[1 + i]); //udub = ai/(Ud U+bi) Ud f0 1/(Ud U+bi)
        //_add(res1, udub1);

        deviceGauge udub = pForce[i][uiLinkIndex];
        _mul(udub, udf1);
        _mul(udub, pForce[i][uiLinkIndex]);
        _mul(udub, pCoefficients[1 + i]);
        _add(res1, udub);
    }

    res1 = _mulC(pGaugeNotProjected[uiLinkIndex], res1);
    _mul(pf0[uiLinkIndex], pForce[order][uiLinkIndex]);
    _sub(pf0[uiLinkIndex], res1); //f0 (c0 + sum ai/(Ud U+bi)) - U sum ai/(Ud U+bi) Ud f0 1/(Ud U+bi)
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectRationalApproximationForceWithProjection(
    const deviceGauge* __restrict__ pGaugeNotProjected,
    const deviceGauge* __restrict__ pGaugeProjected,
    deviceGauge* pf0,
    const deviceGauge* const* __restrict__ pForce,
    const CLGComplex* __restrict__ detres,
    const Real* __restrict__ pCoefficients,
    UINT order,
    Real fOneOverNc)
{
    intokernalDir_NoDir;
    //c0 = pCoefficients[0];
    //a[i] = pCoefficients[1 + i];
    //b[i] = pCoefficients[1 + i + order];

    if (0 == uiLinkIndex)
    {
        _print(pf0[uiLinkIndex]);
    }

    _mul(pf0[uiLinkIndex], _cuConjf(detres[uiLinkIndex])); //J0
    deviceGauge r = _rcpC(pGaugeProjected[uiLinkIndex]); //Mh^-1
    deviceGauge udf1 = _muldagC(pGaugeProjected[uiLinkIndex], pf0[uiLinkIndex]); //Mh J0d
    _mul(r, _tr(udf1));
    _mul(r, fOneOverNc); //Now this is R
    _dagger(r);
    _sub(pf0[uiLinkIndex], r); //K0

    udf1 = _dagmulC(pGaugeNotProjected[uiLinkIndex], pf0[uiLinkIndex]); //Ad K0
    _add(udf1, _daggerC(udf1));
    deviceGauge res1 = _makeZero<deviceGauge>();
    //deviceGauge res2 = _makeZero<deviceGauge>();
    if (0 == uiLinkIndex)
    {
        _print(pGaugeNotProjected[uiLinkIndex]);
        _print(pGaugeProjected[uiLinkIndex]);
        _print(_rcpC(pGaugeProjected[uiLinkIndex]));
    }

    for (UINT i = 0; i < order; ++i)
    {
        const deviceGauge& udub = pForce[i][uiLinkIndex]; //udub = 1/(Ud U+bi)
        deviceGauge udub1 = _mulC(udub, _mulC(udf1, udub)); //udub = 1/(Ud U+bi) Ud f0 1/(Ud U+bi)
        //deviceGauge udub2 = _mulC(udub, _mulC(udf2, udub)); //udub = 1/(Ud U+bi) Ud f0 1/(Ud U+bi)
        _mul(udub1, pCoefficients[1 + i]); //udub = ai/(Ud U+bi) Ud f0 1/(Ud U+bi)
        //_mul(udub2, pCoefficients[1 + i]); //udub = ai/(Ud U+bi) Ud f0 1/(Ud U+bi)
        _add(res1, udub1);
        //_add(res2, udub2);
    }

    //_sub(res1, res2);
    res1 = _mulC(pGaugeNotProjected[uiLinkIndex], res1);

    _mul(pf0[uiLinkIndex], pForce[order][uiLinkIndex]);
    _sub(pf0[uiLinkIndex], res1);

    if (0 == uiLinkIndex)
    {
        _print(pf0[uiLinkIndex]);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectRationalApproximationForceWithProjection2(
    const deviceGauge* __restrict__ pGaugeProjected,
    deviceGauge* pf0,
    const CLGComplex* __restrict__ detres,
    Real fOneOverNc)
{
    intokernalDir_NoDir;
    //c0 = pCoefficients[0];
    //a[i] = pCoefficients[1 + i];
    //b[i] = pCoefficients[1 + i + order];

    //if (2 == uiLinkIndex)
    //{
    //    printf("==== raw f0 =====");
    //    _print(pf0[uiLinkIndex]);

    //    printf("==== det =====");
    //    _print(detres[uiLinkIndex]);
    //}
    //_mul(pf0[uiLinkIndex], detres[uiLinkIndex]);
    //_dagger(pf0[uiLinkIndex]);
    _mul(pf0[uiLinkIndex], _cuConjf(detres[uiLinkIndex])); //J0
    deviceGauge udf1 = _dagmulC(pGaugeProjected[uiLinkIndex], pf0[uiLinkIndex]); 
    deviceGauge r = pGaugeProjected[uiLinkIndex];
    r = _rcpC(r);
    _dagger(r);

    _mul(r, _mulC(_tr(udf1), fOneOverNc));
    _sub(pf0[uiLinkIndex], r); //K0
    //_dagger(pf0[uiLinkIndex]);

    //_dagger(r);
    //_add(pf0[linkIndex], r);
    //if (2 == uiLinkIndex)
    //{
    //    printf("==== after dressed =====");
    //    _print(pf0[uiLinkIndex]);
    //}
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectSU3Force(
    const deviceGauge* __restrict__ pGaugeProjected,
    deviceGauge* pf0,
    const DOUBLE* __restrict__ detphase,
    DOUBLE fOneOverNc)
{
    intokernalDir_NoDir;
    //c0 = pCoefficients[0];
    //a[i] = pCoefficients[1 + i];
    //b[i] = pCoefficients[1 + i + order];

    //if (2 == uiLinkIndex)
    //{
    //    printf("==== raw f0 =====");
    //    _print(pf0[uiLinkIndex]);

    //    printf("==== det =====");
    //    _print(detres[uiLinkIndex]);
    //}
    //_mul(pf0[uiLinkIndex], detres[uiLinkIndex]);
    //_dagger(pf0[uiLinkIndex]);
    _mul(pf0[uiLinkIndex], _make_cuComplex(
        static_cast<Real>(cos(detphase[uiLinkIndex])),
        static_cast<Real>(sin(-detphase[uiLinkIndex]))
    )); //J0
    //deviceGauge udf1 = _dagmulC(pGaugeProjected[uiLinkIndex], pf0[uiLinkIndex]);
    deviceGauge udf1 = pf0[uiLinkIndex];
    _muldag(udf1, pGaugeProjected[uiLinkIndex]);
    deviceGauge r = pGaugeProjected[uiLinkIndex];

    _mul(r, _mulC(_tr(udf1), fOneOverNc));
    _sub(pf0[uiLinkIndex], r); //K0
    //_dagger(pf0[uiLinkIndex]);

    //_dagger(r);
    //_add(pf0[linkIndex], r);
    //if (2 == uiLinkIndex)
    //{
    //    printf("==== after dressed =====");
    //    _print(pf0[uiLinkIndex]);
    //}
}

//================ for SU3 we can use Cayley Hamilton method to calculate sqrt(M) =================
//Although use deviceGauge, it only works for SU3

//Q = U^+U, so Q is Hermitian
template<typename deviceGauge>
__device__ __inline__ void _deviceUVW(const deviceGauge& Q, const deviceGauge& Q2, DOUBLE& u, DOUBLE& v, DOUBLE& w) = delete;

template<>
__device__ __inline__ void _deviceUVW<deviceSU3>(const deviceSU3& Q, const deviceSU3& Q2, DOUBLE& u, DOUBLE& v, DOUBLE& w)
{
    //trace of Hermitian matrix is real
    //DOUBLE c0 = _tr(Q).x / 3.0;
    const DOUBLE a11 = Q.m_me[0].x;
    const DOUBLE a22 = Q.m_me[4].x;
    const DOUBLE a33 = Q.m_me[8].x;

    const DOUBLE a12r = Q.m_me[1].x;
    const DOUBLE a12i = Q.m_me[1].y;
    const DOUBLE a13r = Q.m_me[2].x;
    const DOUBLE a13i = Q.m_me[2].y;
    const DOUBLE a23r = Q.m_me[5].x;
    const DOUBLE a23i = Q.m_me[5].y;

    DOUBLE c0 = (a11 + a22 + a33) / 3.0;
    //c0 = c0 + a22;
    //c0 = c0 + a33;
    //c0 = c0 / 3.0;
    //const DOUBLE c1 = _tr(Q2).x / 2.0;
    //DOUBLE c1 = Q2.m_me[0].x;
    //c1 = c1 + Q2.m_me[4].x;
    //c1 = c1 + Q2.m_me[8].x;
    //c1 = c1 / 2.0;
    DOUBLE c1 = (a11 * a11 + a22 * a22 + a33 * a33) * 0.5;
    const DOUBLE abssq12 = a12r * a12r + a12i * a12i;
    const DOUBLE abssq13 = a13r * a13r + a13i * a13i;
    const DOUBLE abssq23 = a23r * a23r + a23i * a23i;
    c1 += abssq12 + abssq13 + abssq23;

    //deviceSU3 Q3 = Q;
    //_mul(Q3, Q2);
    //const DOUBLE c2 = _tr(Q3).x / 3.0;

    //a11^3 + a22^3 + a33^3 + 3 (a11 + a22) a12sq + 3 (a11 + a33) a13sq + 3 (a22 + a33) a23sq
    //+ 6 (a12r a13i a23i - a12i a13r a23i + a12i a13i a23r + a12r a13r a23r)
    DOUBLE c2 = (a11 * a11 * a11 + a22 * a22 * a22 + a33 * a33 * a33) / 3.0;
    c2 += (a11 + a22) * abssq12;
    c2 += (a11 + a33) * abssq13;
    c2 += (a22 + a33) * abssq23;
    DOUBLE s = a12r * a13i * a23i;
    s -= a12i * a13r * a23i;
    s += a12i * a13i * a23r;
    s += a12r * a13r * a23r;
    c2 += 2.0 * s;

    //Eq.~(31) of hep-lat/0702028 here c0 is c0/3
    s = c1 / 3.0 - c0 * c0 / 2.0;
    DOUBLE r;
    if (s < _CLG_FLT_MIN_)
    {
        //printf("are we here 1? s=%f\n", s);
        s = 0.0;
        //when s = 0, r is irrelevant
        r = 0.0;
    }
    else
    {
        s = sqrt(s);
        r = c2 / 2.0 - c0 * c1 + c0 * c0 * c0;
        r = r / (s * s * s);
        if (r > 1.0 - _CLG_FLT_MIN_)
        {
            //printf("are we here 2? r = %f\n", r);
            r = 0.0;
        }
        else if (r < -1.0 + _CLG_FLT_MIN_)
        {
            //pi / 3
            //printf("are we here 3? r = %f\n", r);
            r = 1.04719755119659774615421446109317;
        }
        else
        {
            r = acos(r) / 3.0;
        }
        
        s = 2.0 * s;
    }

    //printf("c0=%f, c1=%f, c2=%f, s=%f, r=%f\n", c0, c1, c2, s, r);

    //Eq.~(32) note: g0,1,2 are eignvalues of Q, which is none-negative (must be positive, because if it was 0, it cannot be invertable)
    DOUBLE g0 = c0 + s * cos(r - 2.094395102393195492308428922186);
    DOUBLE g1 = c0 + s * cos(r);
    DOUBLE g2 = c0 + s * cos(r + 2.094395102393195492308428922186);
    //printf("c0=%f, s=%f, r=%f, g0=%f, g1=%f, g2=%f\n", g0, g1, g2);

    //Eq.~(29)
    u = sqrt(g0) + sqrt(g1) + sqrt(g2);
    v = sqrt(g0 * g1) + sqrt(g1 * g2) + sqrt(g0 * g2);
    w = sqrt(g0 * g1 * g2);

    //if (w < 1.0e-12)
    //{
    //    printf("warning, w is very small! w = %f\n", w);
    //}
}

__device__ __inline__ void _deviceSqrtf012(DOUBLE &f0, DOUBLE& f1, DOUBLE& f2, const DOUBLE& u, const DOUBLE& v, const DOUBLE& w)
{
    const DOUBLE uv = u * v;
    const DOUBLE u2 = u * u;
    const DOUBLE d = 1.0 / (w * (uv - w));
    f0 = (-w * (u2 + v) + uv * v) * d;
    f1 = (-w - u * u2 + 2.0 * uv) * d;
    f2 = u * d;

    //if (abs(w * (uv - w)) < 1.0e-5)
    //{
    //    printf("warning, w(uv-w) is very small! w(uv-w) = %.15f\n", uv - w);
    //}
}

template<typename deviceGauge>
__device__ __inline__ cuDoubleComplex _deviceDoubleDet(const deviceGauge& g) = delete;

template<>
__device__ __inline__ cuDoubleComplex _deviceDoubleDet<>(const deviceSU3& g)
{
    cuDoubleComplex m0 = _cToDouble(g.m_me[0]);
    cuDoubleComplex m1 = _cToDouble(g.m_me[1]);
    cuDoubleComplex m2 = _cToDouble(g.m_me[2]);
    cuDoubleComplex m3 = _cToDouble(g.m_me[3]);
    cuDoubleComplex m4 = _cToDouble(g.m_me[4]);
    cuDoubleComplex m5 = _cToDouble(g.m_me[5]);
    cuDoubleComplex m6 = _cToDouble(g.m_me[6]);
    cuDoubleComplex m7 = _cToDouble(g.m_me[7]);
    cuDoubleComplex m8 = _cToDouble(g.m_me[8]);

    return cuCsub(
        cuCadd(
            cuCadd(
                cuCmul(cuCmul(m0, m4), m8),
                cuCmul(cuCmul(m1, m5), m6)
            ),
            cuCmul(cuCmul(m3, m7), m2)
        ),
        cuCadd(
            cuCadd(
                cuCmul(cuCmul(m2, m4), m6),
                cuCmul(cuCmul(m1, m3), m8)
            ),
            cuCmul(cuCmul(m0, m5), m7)
        )
    );
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectCaylayHamilton(
    deviceGauge* pDeviceData,
    DOUBLE* coefficients,
    deviceGauge* pQ,
    deviceGauge* pQ2,
    deviceGauge* pInverseSqrtQ2,
    DOUBLE* detphase,
    UBOOL bSU3,
    DOUBLE fOneOverNc)
{
    intokernalDir_NoDir;

    const UINT idxstart = 5 * uiLinkIndex;
    DOUBLE f0 = 0.0;
    DOUBLE& f1 = coefficients[idxstart];
    DOUBLE& f2 = coefficients[idxstart + 1];
    DOUBLE& u = coefficients[idxstart + 2];
    DOUBLE& v = coefficients[idxstart + 3];
    DOUBLE& w = coefficients[idxstart + 4];
    deviceGauge& Q = pQ[uiLinkIndex];
    deviceGauge& Q2 = pQ2[uiLinkIndex];
    deviceGauge& InverseSqrtQ2 = pInverseSqrtQ2[uiLinkIndex];

    Q = pDeviceData[uiLinkIndex];
    _dagmul(Q, pDeviceData[uiLinkIndex]);
    Q2 = Q;
    _mul(Q2, Q);

    deviceSU3::deviceUVW(Q, Q2, u, v, w);

    //if (0 == uiLinkIndex)
    //{
    //    printf("%f, %f, %f\n", u, v, w);
    //}
    deviceSU3::deviceSqrtf012(f0, f1, f2, u, v, w);

    InverseSqrtQ2 = Q;
    _mul(InverseSqrtQ2, f1);
    _add(InverseSqrtQ2, _mulC(Q2, f2));
    _add(InverseSqrtQ2, f0);

    _mul(pDeviceData[uiLinkIndex], InverseSqrtQ2);

    //deviceGauge testres = pDeviceData[uiLinkIndex];
    //_muldag(testres, pDeviceData[uiLinkIndex]);
    //_sub(testres, F(1.0));
    //Real fDelta = _lensq(testres);
    //if (fDelta > F(0.00000001))
    //{
    //    printf("delta %.12f\n", fDelta);
    //    _print(pDeviceData[uiLinkIndex]);
    //}

    if (bSU3)
    {
        //if (0 == linkIndex)
        //{
        //    printf("==== raw u =====");
        //    _print(pDeviceData[linkIndex]);
        //}
        cuDoubleComplex det = _deviceDoubleDet(pDeviceData[uiLinkIndex]);
        //if (abs(_cuCabsf(det) - 1.0) > 0.000001)
        //{
        //    printf("warning, det is not 1, det = ");
        //    _print(det);
        //}

        //atan2 is used in power, so it is automatically -pi/3 to pi/3
        //det = __cuCpowerf(det, -fOneOverNc);
        detphase[uiLinkIndex] = -fOneOverNc * atan2(det.y, det.x);
        //if (2 == uiLinkIndex)
        //{
        //    printf("==== u3 =====");
        //    _print(res);
        //    printf("==== det =====");
        //    _print(det);
        //}
        _mul(pDeviceData[uiLinkIndex], _make_cuComplex(
            static_cast<Real>(cos(detphase[uiLinkIndex])),
            static_cast<Real>(sin(detphase[uiLinkIndex]))
            ));
        //if (2 == uiLinkIndex)
        //{
        //    printf("==== su3 =====");
        //    _print(res);
        //}

        //cuDoubleComplex detcheck = _deviceDoubleDet(pDeviceData[uiLinkIndex]);
        //if (abs(detcheck.x - F(1.0)) > F(0.0000001) || abs(detcheck.y) > F(0.0000001))
        //{
        //    printf("after projection, det is not 1, det = %.12f + %.12f i\n", detcheck.x, detcheck.y);
        //}
    }
}

__device__ __inline__ void _deviceCMatrixElements(DOUBLE& c00, DOUBLE& c01, DOUBLE& c02, DOUBLE& c11, DOUBLE& c12, DOUBLE& c22, const DOUBLE& u, const DOUBLE& v, const DOUBLE& w)
{
    const DOUBLE u2 = u * u;
    const DOUBLE u3 = u2 * u;
    const DOUBLE u4 = u3 * u;
    //const DOUBLE u5 = u4 * u;
    //const DOUBLE u6 = u5 * u;
    //const DOUBLE u7 = u6 * u;
    //const DOUBLE u8 = u7 * u;
    const DOUBLE u6 = u4 * u2;

    const DOUBLE v2 = v * v;
    const DOUBLE v3 = v2 * v;
    const DOUBLE v4 = v3 * v;
    //const DOUBLE v5 = v4 * v;
    //const DOUBLE v6 = v5 * v;

    const DOUBLE w2 = w * w;
    const DOUBLE w3 = w2 * w;
    const DOUBLE w4 = w3 * w;
    //const DOUBLE w5 = w4 * w;

    const DOUBLE d1 = 1.0 / (u * v - w);
    DOUBLE d = d1 / w;
    d = 0.5 * d * d * d;

    //c00 = -w3 * u6;
    //c00 += 3 * v * w3 * u4;
    //c00 += 3 * v4 * w * u4;
    //c00 -= v6 * u3;
    //c00 -= 4 * w4 * u3;
    //c00 -= 12 * v3 * w2 * u3;
    //c00 += 16 * v2 * w3 * u2;
    //c00 += 3 * v5 * w * u2;
    //c00 -= 8 * v * w4 * u;
    //c00 -= 3 * v4 * w2 * u;
    //c00 += w5; 
    //c00 += v3 * w3;
    //c00 *= d;

    // 1/2 (-(v^3/w^3) + (3 u v)/w^2 - 3/w + (u^4 (u^2 + v))/(-u v + w)^3 + (4 u^3)/(-u v + w)^2 + 2/(-u v + w))
    c00 = -v3 / w3;
    c00 += 3 * u * v / w2;
    c00 -= 3 / w;
    c00 -= u4 * (u2 + v) * d1 * d1 * d1;
    c00 += 4 * u3 * d1 * d1;
    c00 -= 2 * d1;
    c00 *= 0.5;

    //c01 = -w2 * u7;
    //c01 -= v2 * w * u6;
    //c01 += v4 * u5;
    //c01 += 6 * v * w2 * u5;
    //c01 -= 5 * w3 * u4;
    //c01 -= v3 * w * u4;
    //c01 -= 2 * v5 * u3;
    //c01 -= 6 * v2 * w2 * u3;
    //c01 += 10 * v * w3 * u2;
    //c01 += 6 * v4 * w * u2;
    //c01 -= 3 * w4 * u;
    //c01 -= 6 * v3 * w2 * u;
    //c01 += 2 * v2 * w3;
    //c01 *= d;

    //u^3 (u^2 - 2 v) v^4 - u^2 v^2 (u^4 + u^2 v - 6 v^2) w - u(u ^ 6 - 6 u ^ 4 v + 6 u ^ 2 v ^ 2 + 6 v ^ 3) w ^ 2 + (-5 u ^ 4 + 10 u ^ 2 v + 2 v ^ 2) w ^ 3 - 3 u w ^ 4
    const DOUBLE u2m2v = u2 - 2 * v;
    c01 = u3 * u2m2v * v4;
    c01 -= u2 * v2 * (u4 + u2 * v - 6 * v2) * w;
    c01 -= u * (u6 - 6 * u4 * v + 6 * u2 * v2 + 6 * v3) * w2;
    c01 += (-5 * u4 + 10 * u2 * v + 2 * v2) * w3;
    c01 -= 3 * u * w4;
    c01 *= d;

    //c02 = w2 * u5;
    //c02 += v2 * w * u4;
    //c02 -= v4 * u3;
    //c02 -= 4 * v * w2 * u3;
    //c02 += 4 * w3 * u2;
    //c02 += 3 * v3 * w * u2;
    //c02 -= 3 * v2 * w2 * u;
    //c02 += v * w3;
    //c02 *= d;

    //-u^3 v^4 + u^2 v^2 (u^2 + 3 v) w + u(u ^ 4 - 4 u ^ 2 v - 3 v ^ 2) w ^ 2 + (4 u ^ 2 + v) w ^ 3
    c02 = -u3 * v4;
    c02 += u2 * v2 * (u2 + 3 * v) * w;
    c02 += u * (u4 - 4 * u2 * v - 3 * v2) * w2;
    c02 += (4 * u2 + v) * w3;
    c02 *= d;

    //c11 = -w * u8;
    //c11 -= v2 * u7;
    //c11 += 7 * v * w * u6;
    //c11 += 4 * v3 * u5;
    //c11 -= 5 * w2 * u5;
    //c11 -= 16 * v2 * w * u4;
    //c11 -= 4 * v4 * u3;
    //c11 += 16 * v * w2 * u3;
    //c11 -= 3 * w3 * u2;
    //c11 += 12 * v3 * w * u2;
    //c11 -= 12 * v2 * w2 * u;
    //c11 += 3 * v * w3;
    //c11 *= d;

    //- ( (u^3 - 2 u v)^2 (u v^2 + u^2 w - 3 v w) + u(5 u ^ 2 - 6 v) (u ^ 2 - 2 v) w ^ 2 + 3 (u ^ 2 - v) w ^ 3)
    c11 = u * u2m2v;
    c11 =  c11 * c11 * (u * v2 + u2 * w - 3 * v * w);
    c11 += u * (5 * u2 - 6 * v) * u2m2v * w2;
    c11 += 3 * (u2 - v) * w3;
    c11 *= -d;

    //c12 = w * u6;
    //c12 += v2 * u5;
    //c12 -= 5 * v * w * u4;
    //c12 -= 2 * v3 * u3;
    //c12 += 4 * w2 * u3;
    //c12 += 6 * v2 * w * u2;
    //c12 -= 6 * v * w2 * u;
    //c12 += w3;
    //c12 *= d;

    //u^3 (u^2 - 2 v) v^2 + u^2 (u^4 - 5 u^2 v + 6 v^2) w + w^2 (4 u^3 - 6 u v + w)
    c12 = u3 * u2m2v * v2;
    c12 += u2 * (u4 - 5 * u2 * v + 6 * v2) * w;
    c12 += w2 * (4 * u3 - 6 * u * v + w);
    c12 *= d;

    //c22 = -w * u4;
    //c22 -= v2 * u3;
    //c22 += 3 * v * w * u2;
    //c22 -= 3 * w2 * u;
    //c22 *= d;

    c22 = -u * (u2 * v2 + u3 * w - 3 * u * v * w + 3 * w2) * d;
}

//b0 + b1 Q + b2 Q^2
template<typename deviceGauge>
__device__ __inline__ void _deviceB(const deviceGauge& Q, const deviceGauge& Q2, deviceGauge& B, const DOUBLE& b0, const DOUBLE& b1, const DOUBLE& b2)
{
    B = Q;
    _mul(B, b1);
    _add(B, _mulC(Q2, b2));
    _add(B, b0);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelProjectCaylayHamiltonForce(
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForce,
    const DOUBLE* __restrict__ coefficients,
    const deviceGauge* __restrict__ pQ,
    const deviceGauge* __restrict__ pQ2,
    const deviceGauge* __restrict__ pInverseSqQ)
{
    intokernalDir_NoDir;

    const UINT idxstart = 5 * uiLinkIndex;
    const DOUBLE& f1 = coefficients[idxstart];
    const DOUBLE& f2 = coefficients[idxstart + 1];
    const DOUBLE& u = coefficients[idxstart + 2];
    const DOUBLE& v = coefficients[idxstart + 3];
    const DOUBLE& w = coefficients[idxstart + 4];
    const deviceGauge& Q = pQ[uiLinkIndex];
    const deviceGauge& Q2 = pQ2[uiLinkIndex];
    const deviceGauge& inverseSqQ = pInverseSqQ[uiLinkIndex];

    deviceGauge& Sigma = pForce[uiLinkIndex];
    const deviceGauge& U = pDeviceData[uiLinkIndex];

    DOUBLE c00, c01, c02, c11, c12, c22;
    _deviceCMatrixElements(c00, c01, c02, c11, c12, c22, u, v, w);

    deviceGauge A = U;
     _dagmul(A, Sigma);
    //M = force.U + Ud.forced
    _re2(A); //A = M

    _muldag(Sigma, inverseSqQ); //Sigma = f0 U12+

    deviceGauge B;
    DOUBLE a = _retr(A); //a0
    //printf("trace1.y=%f", _tr(A).y);
    _deviceB(Q, Q2, B, c00, c01, c02); //Sigma = B0
    _mul(B, a); // B = a0 B0

    deviceGauge C = A;
    _mul(C, Q);
    a = _retr(C); //a1
    //printf("trace2.y=%f", _tr(C).y);
    _deviceB(Q, Q2, C, c01, c11, c12);
    _mul(C, a);
    _add(B, C); // B = a0 B0 + a1 B1

    C = A;
    _mul(C, Q2);
    a = _retr(C); //a2
    //printf("trace3.y=%f", _tr(C).y);
    _deviceB(Q, Q2, C, c02, c12, c22);
    _mul(C, a);
    _add(B, C); //B = sum an Bn

    //C = U;
    //_mul(C, B); //C = U0 (sum an Bn)
    //_add(Sigma, C); //Sigma = f0 U12+ + U0 (sum an Bn)

    C = Q;
    _mul(C, A); //C = QM
    _re2(C); //C = (Q M + M Q)
    
    _mul(C, f2); //C = F2 (Q M + M Q)
    _mul(A, f1); //A = F1 M
    _add(A, C); //A = F1 M + F2 (Q M + M Q)
    _add(A, B); //A = sum an Bn + F1 M + F2 (Q M + M Q)

    C = U;
    _mul(C, A); //C = U0.(sum an Bn) + F1 U0 M + F2 U0 (Q M + M Q)
    _add(Sigma, C); //Sigma= f0 U12+ + U0 (sum an Bn) + F1 U0 M + F2 U0 (Q M + M Q)
}

#pragma endregion

template<typename gaugetype, INT matrixN>
CGaugeSmearingASQTAD<gaugetype, matrixN>::~CGaugeSmearingASQTAD()
{
    //no matter whether it is returned, it will be destroyed
    //if (NULL != m_pOriginalGauge)
    //{
    //    m_pOriginalGauge->Return();
    //    m_pOriginalGauge = NULL;
    //}
    if (NULL != m_pDeviceRationalApproximation)
    {
        checkCudaErrors(__cudaFree(m_pDeviceRationalApproximation));
        m_pDeviceRationalApproximation = NULL;
    }
    if (NULL != m_pProjForcePtr)
    {
        checkCudaErrors(__cudaFree(m_pProjForcePtr));
        m_pProjForcePtr = NULL;
    }
    if (NULL != m_pDeviceDet)
    {
        checkCudaErrors(__cudaFree(m_pDeviceDet));
        m_pDeviceDet = NULL;
    }
    if (NULL != m_pDeviceDetPhase)
    {
        checkCudaErrors(__cudaFree(m_pDeviceDetPhase));
        m_pDeviceDetPhase = NULL;
    }
    if (NULL != m_pCaylayHamiltonCache)
    {
        checkCudaErrors(__cudaFree(m_pCaylayHamiltonCache));
        m_pCaylayHamiltonCache = NULL;
    }
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    CGaugeSmearing::Initial(pOwner, params);

    //m_bProj = TRUE;
    //m_bProjDet = TRUE;
    m_bCalledWhenUpdate = TRUE;
    //Real fra[25] = { F(0.0), F(0.00943108618345698), F(0.0122499930158508), F(0.0187308029056777), F(0.0308130330025528),
    //                 F(0.0521206555919226), F(0.0890870585774984), F(0.153090120000215), F(0.26493803350899), F(0.466760251501358),
    //                 F(0.866223656646014), F(1.8819154073627), F(6.96033769739192), F(5.23045292201785e-05), F(0.000569214182255549), 
    //                 F(0.00226724207135389), F(0.00732861083302471), F(0.0222608882919378), F(0.0662886891030569), F(0.196319420401789), 
    //                 F(0.582378159903323), F(1.74664271771668), F(5.42569216297222), F(18.850085313508), F(99.6213166072174) };
    //checkCudaErrors(__cudaMalloc((void**)&m_pDeviceRationalApproximation, sizeof(Real) * 25));
    //checkCudaErrors(cudaMemcpy(m_pDeviceRationalApproximation, fra, sizeof(Real) * 25, cudaMemcpyHostToDevice));
    //m_uiRationalApproximationOrder = 12;

    //for (INT i = 0; i < 25; ++i)
    //{
    //    appGeneral(_T("%f\n"), fra[i]);
    //}

    //when b is small, single float cannot correctly calculate the force

    INT iVaule = 1;
    params.FetchValueINT(_T("Proj"), iVaule);
    m_bProj = (0 != iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("ProjSUN"), iVaule);
    m_bProjDet = (0 != iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("UseCaylayHamilton"), iVaule);
    m_bUseCaylayHamilton = (0 != iVaule);

    if (!UseCaylayHamilton())
    {
        Real fra[11] = { F(0.28017824602695207), F(0.019466137466992547), F(0.03561358805516765),
                          F(0.08215727655939538), F(0.21113622525403777), F(0.7946025292571568), F(0.00020653817365416855),
                          F(0.00302707751043529), F(0.020073267806506062), F(0.1251758627154607), F(1.0029328744648631) };
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceRationalApproximation, sizeof(Real) * 11));
        checkCudaErrors(cudaMemcpy(m_pDeviceRationalApproximation, fra, sizeof(Real) * 11, cudaMemcpyHostToDevice));
        m_uiRationalApproximationOrder = 5;

        //Real fra[9] = { F(0.0954113732807719), F(0.20212594413375087), F(0.30595220754840236),
        //            F(0.6336901399796908), F(2.269185448392556), F(0.02335218756717348),
        //            F(0.2862905885952218), F(1.4973064554346185), F(9.531358573058725) };
        //checkCudaErrors(__cudaMalloc((void**)&m_pDeviceRationalApproximation, sizeof(Real) * 9));
        //checkCudaErrors(cudaMemcpy(m_pDeviceRationalApproximation, fra, sizeof(Real) * 9, cudaMemcpyHostToDevice));
        //m_uiRationalApproximationOrder = 4;
    }

    params.FetchValueReal(_T("Origin"), m_fOriginal);
    params.FetchValueReal(_T("Fat3"), m_fFat3);
    params.FetchValueReal(_T("Fat5"), m_fFat5);
    params.FetchValueReal(_T("Fat7"), m_fFat7);
    params.FetchValueReal(_T("Lepage"), m_fLepage);

    //I don't know why, but see the line-144, 145 of quda/tests/utils/llfat_utils.cpp
    //m_fOriginal = m_fOriginal - F(6.0) * m_fLepage;
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::GaugeSmearing(CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject)
{
    //if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    //{
    //    appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
    //    return;
    //}

    //test staggered phase
    //CFieldGauge* pAfterPhase = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
    //pOrignal->CopyTo(pAfterPhase);
    //pAfterPhase->ApplyStaggeredPhase();
    //pAfterPhase->CopyTo(pGauge);

    appParanoiac(_T("CGaugeSmearingASQTAD::GaugeSmearing\n"));

    if (NULL == m_pP3_1)
    {
        m_pP3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP3_3 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        m_pP5_1_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_1_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_2_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_2_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_3_1 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
        m_pP5_3_2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

        if (m_bProj)
        {
            m_pGaugeNotProjected = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));

            if (UseCaylayHamilton())
            {
                m_pQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pQ2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                m_pInverseSqrtQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
                checkCudaErrors(__cudaMalloc((void**)&m_pCaylayHamiltonCache, sizeof(DOUBLE) * _HC_Volume * _HC_Dir * 5));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDetPhase, sizeof(DOUBLE) * _HC_Dir * _HC_Volume));
                }
            }
            else
            {
                TArray<gaugetype*> ptrs;
                for (UINT i = 0; i < m_uiRationalApproximationOrder + 1; ++i)
                {
                    m_pProjForce.AddItem(dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__)));
                    CFieldGaugeLink<gaugetype, matrixN>* pProjForce = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(m_pProjForce[i]);
                    ptrs.AddItem(pProjForce->m_pDeviceData);
                }
                checkCudaErrors(__cudaMalloc((void**)&m_pProjForcePtr, sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1)));
                checkCudaErrors(cudaMemcpy(m_pProjForcePtr, ptrs.GetData(), sizeof(gaugetype*) * (m_uiRationalApproximationOrder + 1), cudaMemcpyHostToDevice));

                if (m_bProjDet)
                {
                    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceDet, sizeof(CLGComplex) * _HC_Dir * _HC_Volume));
                }
            }
        }
    }

    Fat357Lepage(pGauge,
        NULL,
        m_pP3_1,
        m_pP3_2,
        m_pP3_3,
        m_pP5_1_1,
        m_pP5_1_2,
        m_pP5_2_1,
        m_pP5_2_2,
        m_pP5_3_1,
        m_pP5_3_2,
        m_fOriginal,
        m_fFat3,
        m_fFat5,
        m_fFat7,
        m_fLepage,
        pOrignal
        //pAfterPhase
    );

    if (m_bProj)
    {
        //pGauge->ApplyStaggeredPhase();
        if (UseCaylayHamilton())
        {
            ProjectCaylayHamilton(pGauge,
                m_pGaugeNotProjected,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            ProjectRationalApproximation(pGauge,
                m_pGaugeNotProjected,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
        //pGauge->ApplyStaggeredPhase();
    }

    //pGauge->ApplyStaggeredPhase(1);

    //m_pP3_1->ScalarMultply(F(-1.0));
    //m_pP3_2->ScalarMultply(F(-1.0));
    //m_pP3_3->ScalarMultply(F(-1.0));
    //m_pP3_1->ApplyStaggeredPhase();
    //m_pP3_2->ApplyStaggeredPhase();
    //m_pP3_3->ApplyStaggeredPhase();
    //m_pP5_1_1->ApplyStaggeredPhase();
    //m_pP5_1_2->ApplyStaggeredPhase();
    //m_pP5_2_1->ApplyStaggeredPhase();
    //m_pP5_2_2->ApplyStaggeredPhase();
    //m_pP5_3_1->ApplyStaggeredPhase();
    //m_pP5_3_2->ApplyStaggeredPhase();
    //pGauge->ApplyStaggeredPhase();
    //pAfterPhase->Return();
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::Fat357Lepage(class CFieldGauge* pGauge, class CFieldGauge* preserve, class CFieldGauge* pp3_1, class CFieldGauge* pp3_2, class CFieldGauge* pp3_3,
    class CFieldGauge* pp5_1_1, class CFieldGauge* pp5_1_2, class CFieldGauge* pp5_2_1, class CFieldGauge* pp5_2_2, class CFieldGauge* pp5_3_1, class CFieldGauge* pp5_3_2,
    Real fOrignal, Real fFat3, Real fFat5, Real fFat7, Real fLepage, const class CFieldGauge* pGaugeOrignal)
{
    if (NULL != preserve)
    {
        pGauge->CopyTo(preserve);
    }

    //3 fields to store 3-staple
    const CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3Original = NULL != preserve ? dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(preserve) : dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pGaugeOrignal);
    CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pGauge);
    const BYTE byFieldId = pGauge->m_byFieldId;
    CFieldGaugeLink<gaugetype, matrixN>* p3_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp3_1);
    CFieldGaugeLink<gaugetype, matrixN>* p3_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp3_2);
    CFieldGaugeLink<gaugetype, matrixN>* p3_3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp3_3);
    CFieldGaugeLink<gaugetype, matrixN>* p5_1_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp5_1_1);
    CFieldGaugeLink<gaugetype, matrixN>* p5_1_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp5_1_2);
    CFieldGaugeLink<gaugetype, matrixN>* p5_2_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp5_2_1);
    CFieldGaugeLink<gaugetype, matrixN>* p5_2_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp5_2_2);
    CFieldGaugeLink<gaugetype, matrixN>* p5_3_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp5_3_1);
    CFieldGaugeLink<gaugetype, matrixN>* p5_3_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pp5_3_2);

    preparethreadEVar(_HC_Dir * 3, blockfat3, threadsfat3);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelASQTADFat3<gaugetype>, blockfat3, threadsfat3, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        p3_1->m_pDeviceData,
        p3_2->m_pDeviceData,
        p3_3->m_pDeviceData,
        byFieldId
        );
#else
    _LAUNCH_KERNEL(_kernelASQTADFat3<gaugetype>, blockfat3, threadsfat3, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        p3_1->m_pDeviceData,
        p3_2->m_pDeviceData,
        p3_3->m_pDeviceData,
        byFieldId
    );
#endif

    preparethreadEVar(_HC_Dir * 6, blockfat5, threadsfat5);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelASQTADFat5<gaugetype>, blockfat5, threadsfat5, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        p3_1->m_pDeviceData,
        p3_2->m_pDeviceData,
        p3_3->m_pDeviceData,
        p5_1_1->m_pDeviceData,
        p5_1_2->m_pDeviceData,
        p5_2_1->m_pDeviceData,
        p5_2_2->m_pDeviceData,
        p5_3_1->m_pDeviceData,
        p5_3_2->m_pDeviceData,
        byFieldId
        );
#else
    _LAUNCH_KERNEL(_kernelASQTADFat5<gaugetype>, blockfat5, threadsfat5, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        p3_1->m_pDeviceData,
        p3_2->m_pDeviceData,
        p3_3->m_pDeviceData,
        p5_1_1->m_pDeviceData,
        p5_1_2->m_pDeviceData,
        p5_2_1->m_pDeviceData,
        p5_2_2->m_pDeviceData,
        p5_3_1->m_pDeviceData,
        p5_3_2->m_pDeviceData,
        byFieldId
    );
#endif

    //why this is one-link - 6 * lepage, see comments in _kernelASQTADLepage
    pGauge->ScalarMultply(fOrignal - F(6.0) * fLepage);
    pGauge->Axpy(fFat3, p3_1);
    pGauge->Axpy(fFat3, p3_2);
    pGauge->Axpy(fFat3, p3_3);

    pGauge->Axpy(fFat5, p5_1_1);
    pGauge->Axpy(fFat5, p5_1_2);
    pGauge->Axpy(fFat5, p5_2_1);
    pGauge->Axpy(fFat5, p5_2_2);
    pGauge->Axpy(fFat5, p5_3_1);
    pGauge->Axpy(fFat5, p5_3_2);

    //pGauge->DebugPrintMe();

    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelASQTADLepage<gaugetype>, block, threads, pGaugeSU3Original->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        p3_1->m_pDeviceData,
        p3_2->m_pDeviceData,
        p3_3->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        fLepage,
        byFieldId
        );
    _LAUNCH_KERNEL(_kernelASQTADFat7<gaugetype>, block, threads, pGaugeSU3Original->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        p5_1_1->m_pDeviceData,
        p5_1_2->m_pDeviceData,
        p5_2_1->m_pDeviceData,
        p5_2_2->m_pDeviceData,
        p5_3_1->m_pDeviceData,
        p5_3_2->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        fFat7,
        byFieldId
    );
#else
    if (abs(fLepage) > _CLG_FLT_EPSILON)
    {
        _LAUNCH_KERNEL(_kernelASQTADLepage<gaugetype>, block, threads, pGaugeSU3Original->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            p3_1->m_pDeviceData,
            p3_2->m_pDeviceData,
            p3_3->m_pDeviceData,
            pGaugeSU3->m_pDeviceData,
            fLepage,
            byFieldId
        );
    }

    _LAUNCH_KERNEL(_kernelASQTADFat7<gaugetype>, block, threads, pGaugeSU3Original->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        p5_1_1->m_pDeviceData,
        p5_1_2->m_pDeviceData,
        p5_2_1->m_pDeviceData,
        p5_2_2->m_pDeviceData,
        p5_3_1->m_pDeviceData,
        p5_3_2->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        fFat7,
        byFieldId
    );
#endif
    //pGauge->DebugPrintMe();
    pGauge->FixBoundary(EFB_Field);
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat3Only(
    const CFieldGauge* pGauge,
    CFieldGauge* p3_1, CFieldGauge* p3_2, CFieldGauge* p3_3)
{
    const CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pGauge);
    const BYTE byFieldId = pGauge->m_byFieldId;
    CFieldGaugeLink<gaugetype, matrixN>* p3_1_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p3_1);
    CFieldGaugeLink<gaugetype, matrixN>* p3_2_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p3_2);
    CFieldGaugeLink<gaugetype, matrixN>* p3_3_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p3_3);

    preparethreadEVar(_HC_Dir * 3, blockfat3, threadsfat3);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelASQTADFat3<gaugetype>, blockfat3, threadsfat3, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        p3_1_link->m_pDeviceData,
        p3_2_link->m_pDeviceData,
        p3_3_link->m_pDeviceData,
        byFieldId
        );
#else
    _LAUNCH_KERNEL(_kernelASQTADFat3<gaugetype>, blockfat3, threadsfat3, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        p3_1_link->m_pDeviceData,
        p3_2_link->m_pDeviceData,
        p3_3_link->m_pDeviceData,
        byFieldId
    );
#endif
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::CalculateFat5Only(
    const CFieldGauge* pGauge,
    const CFieldGauge* p3_1, const CFieldGauge* p3_2, const CFieldGauge* p3_3,
    CFieldGauge* p5_1_1, CFieldGauge* p5_1_2, CFieldGauge* p5_2_1,
    CFieldGauge* p5_2_2, CFieldGauge* p5_3_1, CFieldGauge* p5_3_2)
{
    const CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pGauge);
    const BYTE byFieldId = pGauge->m_byFieldId;
    const CFieldGaugeLink<gaugetype, matrixN>* p3_1_link = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(p3_1);
    const CFieldGaugeLink<gaugetype, matrixN>* p3_2_link = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(p3_2);
    const CFieldGaugeLink<gaugetype, matrixN>* p3_3_link = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(p3_3);

    CFieldGaugeLink<gaugetype, matrixN>* p5_1_1_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p5_1_1);
    CFieldGaugeLink<gaugetype, matrixN>* p5_1_2_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p5_1_2);
    CFieldGaugeLink<gaugetype, matrixN>* p5_2_1_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p5_2_1);
    CFieldGaugeLink<gaugetype, matrixN>* p5_2_2_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p5_2_2);
    CFieldGaugeLink<gaugetype, matrixN>* p5_3_1_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p5_3_1);
    CFieldGaugeLink<gaugetype, matrixN>* p5_3_2_link = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(p5_3_2);

    preparethreadEVar(_HC_Dir * 6, blockfat5, threadsfat5);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelASQTADFat5<gaugetype>, blockfat5, threadsfat5, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        p3_1_link->m_pDeviceData,
        p3_2_link->m_pDeviceData,
        p3_3_link->m_pDeviceData,
        p5_1_1_link->m_pDeviceData,
        p5_1_2_link->m_pDeviceData,
        p5_2_1_link->m_pDeviceData,
        p5_2_2_link->m_pDeviceData,
        p5_3_1_link->m_pDeviceData,
        p5_3_2_link->m_pDeviceData,
        byFieldId
        );
#else
    _LAUNCH_KERNEL(_kernelASQTADFat5<gaugetype>, blockfat5, threadsfat5, pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        p3_1_link->m_pDeviceData,
        p3_2_link->m_pDeviceData,
        p3_3_link->m_pDeviceData,
        p5_1_1_link->m_pDeviceData,
        p5_1_2_link->m_pDeviceData,
        p5_2_1_link->m_pDeviceData,
        p5_2_2_link->m_pDeviceData,
        p5_3_1_link->m_pDeviceData,
        p5_3_2_link->m_pDeviceData,
        byFieldId
    );
#endif
}

template<typename gaugetype, INT matrixN>
void  CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximation(class CFieldGauge* pGauge, class CFieldGauge* pGaugeBeforeProj, gaugetype** forcePointers,
    const Real* rationalCoeffs, UINT uiOrder, UBOOL bToSU3, CLGComplex* det)
{
    CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pGauge);
    //const BYTE byFieldId = pGauge->m_byFieldId;

    preparethreadDir;
    pGauge->CopyTo(pGaugeBeforeProj);
    _LAUNCH_KERNEL(_kernelProjectRationalApproximation<gaugetype>, block, threads, pGaugeSU3->m_pDeviceData,
        forcePointers,
        det,
        rationalCoeffs,
        uiOrder,
        bToSU3,
        static_cast<Real>(F(1.0) / matrixN)
        );
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamilton(class CFieldGauge* pGauge, class CFieldGauge* pGaugeBeforeProj, class CFieldGauge* pQ, class CFieldGauge* pQ2, 
    class CFieldGauge* pInverseSqrtQ, DOUBLE* constants, UBOOL bToSU3, DOUBLE* detphase)
{
    CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pGauge);
    CFieldGaugeLink<gaugetype, matrixN>* pQSU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pQ);
    CFieldGaugeLink<gaugetype, matrixN>* pQ2SU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pQ2);
    CFieldGaugeLink<gaugetype, matrixN>* pInverseSqrtQSU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pInverseSqrtQ);

    preparethreadDir;
    pGauge->CopyTo(pGaugeBeforeProj);
    _LAUNCH_KERNEL(_kernelProjectCaylayHamilton<gaugetype>, block, threads, 
        pGaugeSU3->m_pDeviceData,
        constants,
        pQSU3->m_pDeviceData,
        pQ2SU3->m_pDeviceData,
        pInverseSqrtQSU3->m_pDeviceData,
        detphase,
        bToSU3,
        1.0 / matrixN
    );
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::DerivateOnU(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    pf0->FixBoundary(EFB_Force);

    if (m_bProj)
    {
        //pf0->ApplyStaggeredPhase();
        if (UseCaylayHamilton())
        {
            ProjectCaylayHamiltonForce(pEffectiveGauge,
                m_pGaugeNotProjected,
                pf0,
                m_pQ,
                m_pQ2,
                m_pInverseSqrtQ,
                m_pCaylayHamiltonCache,
                m_bProjDet,
                m_pDeviceDetPhase);
        }
        else
        {
            ProjectRationalApproximationForce(pEffectiveGauge,
                m_pGaugeNotProjected,
                pf0,
                m_pProjForcePtr,
                m_pDeviceRationalApproximation,
                m_uiRationalApproximationOrder,
                m_bProjDet,
                m_pDeviceDet);
        }
        //pf0->ApplyStaggeredPhase();
    }

    SmearingForce(pOrignalGauge,
        pf0,
        m_pP3_1,
        m_pP3_2,
        m_pP3_3,
        m_pP5_1_1,
        m_pP5_1_2,
        m_pP5_2_1,
        m_pP5_2_2,
        m_pP5_3_1,
        m_pP5_3_2,
        m_fOriginal,
        m_fFat3,
        m_fFat5,
        m_fFat7,
        m_fLepage
    );
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::SmearingForce(const class CFieldGauge* pOrignalGauge, class CFieldGauge* pf0,
    const class CFieldGauge* pp3_1, const class CFieldGauge* pp3_2, const class CFieldGauge* pp3_3,
    const class CFieldGauge* pp5_1_1, const class CFieldGauge* pp5_1_2, const class CFieldGauge* pp5_2_1, const class CFieldGauge* pp5_2_2, const class CFieldGauge* pp5_3_1, const class CFieldGauge* pp5_3_2,
    Real fOrignal, Real fFat3, Real fFat5, Real fFat7, Real fLepage)
{
    const CFieldGaugeLink<gaugetype, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pOrignalGauge);
    CFieldGaugeLink<gaugetype, matrixN>* pf0SU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pf0);
    const BYTE byFieldId = pGaugeSU3->m_byFieldId;
    const CFieldGaugeLink<gaugetype, matrixN>* p3_1 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp3_1);
    const CFieldGaugeLink<gaugetype, matrixN>* p3_2 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp3_2);
    const CFieldGaugeLink<gaugetype, matrixN>* p3_3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp3_3);

    const CFieldGaugeLink<gaugetype, matrixN>* p5_1_1 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp5_1_1);
    const CFieldGaugeLink<gaugetype, matrixN>* p5_1_2 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp5_1_2);
    const CFieldGaugeLink<gaugetype, matrixN>* p5_2_1 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp5_2_1);
    const CFieldGaugeLink<gaugetype, matrixN>* p5_2_2 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp5_2_2);
    const CFieldGaugeLink<gaugetype, matrixN>* p5_3_1 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp5_3_1);
    const CFieldGaugeLink<gaugetype, matrixN>* p5_3_2 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pp5_3_2);

    CFieldGaugeLink<gaugetype, matrixN>* pf_res_all = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* pf_res = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));

    CFieldGaugeLink<gaugetype, matrixN>* p3f_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p3f_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p3f_3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));

    CFieldGaugeLink<gaugetype, matrixN>* p5f_1_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p5f_1_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p5f_2_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p5f_2_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p5f_3_1 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
    CFieldGaugeLink<gaugetype, matrixN>* p5f_3_2 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));

    //preparethreadDir;
    preparethreadEVar(_HC_Dir * 3, block3, thread3);
    preparethreadEVar(_HC_Dir * 6, block6, thread6);
    {
        _RECORD2(CGaugeSmearingASQTAD::SmearingForce::_kernelASQTADFat3Force_Optimize, a);
        
        //_LAUNCH_KERNEL(_kernelASQTADFat3Force, block3, thread3, pGaugeSU3->m_pDeviceData,
        //    pf0SU3->m_pDeviceData,
        //    appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        //    appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        //    appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        //    pf_res_all->m_pDeviceData,
        //    p3f_1->m_pDeviceData,
        //    p3f_2->m_pDeviceData,
        //    p3f_3->m_pDeviceData,
        //    byFieldId
        //    );
#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelASQTADFat3Force_Optimize<gaugetype>, block3, thread3, pGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            pf_res_all->m_pDeviceData,
            p3f_1->m_pDeviceData,
            p3f_2->m_pDeviceData,
            p3f_3->m_pDeviceData,
            byFieldId
            );
#else
        _LAUNCH_KERNEL(_kernelASQTADFat3Force_Optimize<gaugetype>, block3, thread3, pGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            pf_res_all->m_pDeviceData,
            p3f_1->m_pDeviceData,
            p3f_2->m_pDeviceData,
            p3f_3->m_pDeviceData,
            byFieldId
        );
#endif
        pf_res_all->ScalarMultply(fFat3);
    }

    {
        _RECORD2(CGaugeSmearingASQTAD::SmearingForce::_kernelASQTADFat5Force, a);
#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelASQTADFat5Force_Optimize<gaugetype>, block6, thread6, pGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            p3_1->m_pDeviceData,
            p3_2->m_pDeviceData,
            p3_3->m_pDeviceData,
            p3f_1->m_pDeviceData,
            p3f_2->m_pDeviceData,
            p3f_3->m_pDeviceData,
            pf_res->m_pDeviceData,
            p5f_1_1->m_pDeviceData,
            p5f_1_2->m_pDeviceData,
            p5f_2_1->m_pDeviceData,
            p5f_2_2->m_pDeviceData,
            p5f_3_1->m_pDeviceData,
            p5f_3_2->m_pDeviceData,
            byFieldId
            );
#else
        _LAUNCH_KERNEL(_kernelASQTADFat5Force_Optimize<gaugetype>, block6, thread6, pGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            p3_1->m_pDeviceData,
            p3_2->m_pDeviceData,
            p3_3->m_pDeviceData,
            p3f_1->m_pDeviceData,
            p3f_2->m_pDeviceData,
            p3f_3->m_pDeviceData,
            pf_res->m_pDeviceData,
            p5f_1_1->m_pDeviceData,
            p5f_1_2->m_pDeviceData,
            p5f_2_1->m_pDeviceData,
            p5f_2_2->m_pDeviceData,
            p5f_3_1->m_pDeviceData,
            p5f_3_2->m_pDeviceData,
            byFieldId
        );
#endif
        pf_res_all->Axpy(fFat5, pf_res);
    }

    if (abs(fLepage) > _CLG_FLT_EPSILON)
    {
        _RECORD2(CGaugeSmearingASQTAD::SmearingForce::_kernelASQTADLepageForce, a);
//      preparethreadEVar(_HC_Dir * 3, blocklapage, thread);
//#if !_CLG_ASSUME_SQUARE_LATTICE
//        _LAUNCH_KERNEL(_kernelASQTADLepageForce_Optimize<gaugetype>, block3, thread3, pGaugeSU3->m_pDeviceData,
//            pf0SU3->m_pDeviceData,
//            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
//            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
//            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
//            p3_1->m_pDeviceData,
//            p3_2->m_pDeviceData,
//            p3_3->m_pDeviceData,
//            p3f_1->m_pDeviceData,
//            p3f_2->m_pDeviceData,
//            p3f_3->m_pDeviceData,
//            pf_res->m_pDeviceData,
//            byFieldId
//            );
//#else
//        _LAUNCH_KERNEL(_kernelASQTADLepageForce_Optimize<gaugetype>, block3, thread3, pGaugeSU3->m_pDeviceData,
//            pf0SU3->m_pDeviceData,
//            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
//            p3_1->m_pDeviceData,
//            p3_2->m_pDeviceData,
//            p3_3->m_pDeviceData,
//            p3f_1->m_pDeviceData,
//            p3f_2->m_pDeviceData,
//            p3f_3->m_pDeviceData,
//            pf_res->m_pDeviceData,
//            byFieldId
//        );
//#endif

        {
            preparethreadEDir(12);
            _LAUNCH_KERNEL(_kernelASQTADLepageForceA<gaugetype>, block, threads, 
                pGaugeSU3->m_pDeviceData,
                pf0SU3->m_pDeviceData,
                pf_res->m_pDeviceData,
                byFieldId
            );
        }
        {
            preparethreadEDir(3);
            _LAUNCH_KERNEL(_kernelASQTADLepageForceB<gaugetype>, block, threads,
                pGaugeSU3->m_pDeviceData,
                pf0SU3->m_pDeviceData,
                pf_res->m_pDeviceData,
                byFieldId
            );
        }
        

        pf_res_all->Axpy(fLepage, pf_res);
    }

    {
        _RECORD2(CGaugeSmearingASQTAD::SmearingForce::_kernelASQTADFat7Force, a);

        //preparethread;
        //pf_res->Zero();
        //const gaugetype* fat7_pointersp3[3] = { p3_1->m_pDeviceData, p3_2->m_pDeviceData, p3_3->m_pDeviceData };
        //const gaugetype* fat7_pointerspf3[3] = { p3f_1->m_pDeviceData, p3f_2->m_pDeviceData, p3f_3->m_pDeviceData };
        //const gaugetype* fat7_pointersp5[6] = { p5_1_1->m_pDeviceData, p5_1_2->m_pDeviceData, p5_2_1->m_pDeviceData, p5_2_2->m_pDeviceData, p5_3_1->m_pDeviceData, p5_3_2->m_pDeviceData };
        //const gaugetype* fat7_pointerspf5[6] = { p5f_1_1->m_pDeviceData, p5f_1_2->m_pDeviceData, p5f_2_1->m_pDeviceData, p5f_2_2->m_pDeviceData, p5f_3_1->m_pDeviceData, p5f_3_2->m_pDeviceData };
        //constexpr BYTE fat7_uiStapleIndex[6] = { 0, 0, 1, 1, 2, 2 };
        //constexpr BYTE fat7_t1idx[4][6] = {
        //    {1, 2, 1, 2, 1, 2},
        //    {1, 2, 0, 2, 0, 2},
        //    {0, 2, 0, 2, 0, 1},
        //    {0, 1, 0, 1, 0, 1}
        //};
        //constexpr BYTE fat7_t2idx[4][6] = {
        //    {2, 1, 2, 1, 2, 1},
        //    {2, 1, 2, 0, 2, 0},
        //    {2, 0, 2, 0, 1, 0},
        //    {1, 0, 1, 0, 1, 0}
        //};
        //constexpr BYTE fat7_t3idx[4][6] = {
        //    {3, 5, 3, 5, 3, 5},
        //    {3, 5, 1, 4, 1, 4},
        //    {1, 4, 1, 4, 0, 2},
        //    {0, 2, 0, 2, 0, 2}
        //};
        //constexpr BYTE fat7_t4idx[6] = { 3, 5, 1, 4, 0, 2 };
        //for (BYTE dir = 0; dir < 4; ++dir)
        //{
        //    for (BYTE elementIdx = 0; elementIdx < 6; ++elementIdx)
        //    {
        //        const gaugetype* p3t1 = fat7_pointersp3[fat7_t1idx[dir][elementIdx]];
        //        const gaugetype* p3ft2 = fat7_pointerspf3[fat7_t2idx[dir][elementIdx]];
        //        const gaugetype* p5t3 = fat7_pointersp5[fat7_t3idx[dir][elementIdx]];
        //        const gaugetype* p5ft3 = fat7_pointerspf5[fat7_t3idx[dir][elementIdx]];
        //        const gaugetype* p5ft4 = fat7_pointerspf5[fat7_t4idx[elementIdx]];
        //        const gaugetype* pf0 = pf0SU3->m_pDeviceData;
        //        const gaugetype* p0 = pGaugeSU3->m_pDeviceData;
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            pf0, p0, p5t3, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            p5t3, p0, pf0, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            p3t1, p0, p3ft2, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            p3ft2, p0, p3t1, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            p0, p0, p5ft3, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            p5ft3, p0, p0, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx],
        //            p0, p5ft4, p0, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            pf0, p0, p5t3, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            p5t3, p0, pf0, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            p3t1, p0, p3ft2, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            p3ft2, p0, p3t1, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            p0, p0, p5ft3, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            p5ft3, p0, p0, pf_res->m_pDeviceData, byFieldId);
        //        _LAUNCH_KERNEL(_kernelASQTADFatForce_Optimize, block, threads, appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId], dir, 6U * fat7_uiStapleIndex[elementIdx] + 3U,
        //            p0, p5ft4, p0, pf_res->m_pDeviceData, byFieldId);
        //    }
        //}


#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelASQTADFat7Force_Optimize TMPARG(gaugetype, matrixN), block6, thread6, pGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            p3_1->m_pDeviceData,
            p3_2->m_pDeviceData,
            p3_3->m_pDeviceData,
            p3f_1->m_pDeviceData,
            p3f_2->m_pDeviceData,
            p3f_3->m_pDeviceData,
            p5_1_1->m_pDeviceData,
            p5_1_2->m_pDeviceData,
            p5_2_1->m_pDeviceData,
            p5_2_2->m_pDeviceData,
            p5_3_1->m_pDeviceData,
            p5_3_2->m_pDeviceData,
            p5f_1_1->m_pDeviceData,
            p5f_1_2->m_pDeviceData,
            p5f_2_1->m_pDeviceData,
            p5f_2_2->m_pDeviceData,
            p5f_3_1->m_pDeviceData,
            p5f_3_2->m_pDeviceData,
            pf_res->m_pDeviceData,
            byFieldId
            );
#else
        _LAUNCH_KERNEL(_kernelASQTADFat7Force_Optimize TMPARG(gaugetype, matrixN), block6, thread6, pGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
            p3_1->m_pDeviceData,
            p3_2->m_pDeviceData,
            p3_3->m_pDeviceData,
            p3f_1->m_pDeviceData,
            p3f_2->m_pDeviceData,
            p3f_3->m_pDeviceData,
            p5_1_1->m_pDeviceData,
            p5_1_2->m_pDeviceData,
            p5_2_1->m_pDeviceData,
            p5_2_2->m_pDeviceData,
            p5_3_1->m_pDeviceData,
            p5_3_2->m_pDeviceData,
            p5f_1_1->m_pDeviceData,
            p5f_1_2->m_pDeviceData,
            p5f_2_1->m_pDeviceData,
            p5f_2_2->m_pDeviceData,
            p5f_3_1->m_pDeviceData,
            p5f_3_2->m_pDeviceData,
            pf_res->m_pDeviceData,
            byFieldId
        );
#endif
        pf_res_all->Axpy(fFat7, pf_res);
    }

    pf0->ScalarMultply(fOrignal);
    pf0->AxpyPlus(pf_res_all);

    pf_res_all->Return();
    pf_res->Return();

    p3f_1->Return();
    p3f_2->Return();
    p3f_3->Return();

    p5f_1_1->Return();
    p5f_1_2->Return();
    p5f_2_1->Return();
    p5f_2_2->Return();
    p5f_3_1->Return();
    p5f_3_2->Return();
    _CHECKCUDA;
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectRationalApproximationForce(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pGaugeBeforeProj, class CFieldGauge* pf0, 
    const gaugetype* const* forcePointers, const Real* rationalCoeffs, UINT uiOrder, UBOOL bToSU3, const CLGComplex* det)
{
    preparethreadDir;
    const CFieldGaugeLink<gaugetype, matrixN>* pEffectiveGaugeSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pEffectiveGauge);
    const CFieldGaugeLink<gaugetype, matrixN>* pGaugeNotProjected = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pGaugeBeforeProj);
    CFieldGaugeLink<gaugetype, matrixN>* pf0SU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pf0);

    if (bToSU3)
    {
        _LAUNCH_KERNEL(_kernelProjectRationalApproximationForceWithProjection2<gaugetype>, block, threads,
            pEffectiveGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            det,
            static_cast<Real>(F(1.0) / matrixN)
            );
    }
    _LAUNCH_KERNEL(_kernelProjectRationalApproximationForce<gaugetype>, block, threads,
        pGaugeNotProjected->m_pDeviceData,
        pf0SU3->m_pDeviceData,
        forcePointers,
        rationalCoeffs,
        uiOrder);
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingASQTAD<gaugetype, matrixN>::ProjectCaylayHamiltonForce(const class CFieldGauge* pEffectiveGauge, const class CFieldGauge* pGaugeBeforeProj, class CFieldGauge* pf0,
    const class CFieldGauge* pQ, const class CFieldGauge* pQ2, const class CFieldGauge* pInverseSqrtQ, const DOUBLE* constants, UBOOL bToSU3, const DOUBLE* detphase)
{
    preparethreadDir;
    const CFieldGaugeLink<gaugetype, matrixN>* pEffectiveGaugeSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pEffectiveGauge);
    const CFieldGaugeLink<gaugetype, matrixN>* pGaugeNotProjected = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pGaugeBeforeProj);

    const CFieldGaugeLink<gaugetype, matrixN>* pQSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pQ);
    const CFieldGaugeLink<gaugetype, matrixN>* pQ2SU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pQ2);
    const CFieldGaugeLink<gaugetype, matrixN>* pInverseSqrtQSU3 = dynamic_cast<const CFieldGaugeLink<gaugetype, matrixN>*>(pInverseSqrtQ);
    CFieldGaugeLink<gaugetype, matrixN>* pf0SU3 = dynamic_cast<CFieldGaugeLink<gaugetype, matrixN>*>(pf0);

    if (bToSU3)
    {
        _LAUNCH_KERNEL(_kernelProjectSU3Force<gaugetype>, block, threads,
            pEffectiveGaugeSU3->m_pDeviceData,
            pf0SU3->m_pDeviceData,
            detphase,
            1.0 / matrixN
        );

        //_LAUNCH_KERNEL(_kernelProjectRationalApproximationForceWithProjection2<gaugetype>, block, threads,
        //    pEffectiveGaugeSU3->m_pDeviceData,
        //    pf0SU3->m_pDeviceData,
        //    det,
        //    static_cast<Real>(F(1.0) / matrixN)
        //);
    }
    _LAUNCH_KERNEL(_kernelProjectCaylayHamiltonForce<gaugetype>, block, threads,
        pGaugeNotProjected->m_pDeviceData,
        pf0SU3->m_pDeviceData,
        constants,
        pQSU3->m_pDeviceData,
        pQ2SU3->m_pDeviceData,
        pInverseSqrtQSU3->m_pDeviceData);
    //pf0SU3->Zero();
}

template<typename gaugetype, INT matrixN>
CCString CGaugeSmearingASQTAD<gaugetype, matrixN>::GetInfos(const CCString &tab) const
{
    CCString sRet = CGaugeSmearing::GetInfos(tab);
    sRet = sRet + tab + _T("Original : ") + appToString(m_fOriginal) + _T("\n");
    sRet = sRet + tab + _T("Fat3     : ") + appToString(m_fFat3) + _T("\n");
    sRet = sRet + tab + _T("Fat5     : ") + appToString(m_fFat5) + _T("\n");
    sRet = sRet + tab + _T("Fat7     : ") + appToString(m_fFat7) + _T("\n");
    sRet = sRet + tab + _T("Lepage   : ") + appToString(m_fLepage) + _T("\n");
    sRet = sRet + tab + _T("ProjectU : ") + appToString(m_bProj) + _T("\n");
    return sRet;
}

template class CGaugeSmearingASQTAD<deviceSU3, 3>;
__CLGIMPLEMENT_CLASS(CGaugeSmearingASQTADSU3)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================