//=============================================================================
// FILENAME : DeviceInlineStaggeredGamma.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [07/06/2024 nbale]
//=============================================================================

#ifndef _DEVICEINLINE_STAGGEREDGAMMA_H_
#define _DEVICEINLINE_STAGGEREDGAMMA_H_

__BEGIN_NAMESPACE

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _devicePlaneDiagonalT(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2)
{
    INT dir1[2];

    dir1[0] = dim1; dir1[1] = dim2;
    deviceGauge sRet(_deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    _mul(sRet, F(0.5));
    return sRet;
}

/**
 * dim1, 2, 3 =
 * 1: x, -1: -x
 * 2: y, -2: -y
 * 3: z, -3: -z
 * 4: t, -4: -t
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceCubicDiagonalT(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2, SBYTE dim3)
{
    INT dir1[3];

    dir1[0] = dim1; dir1[1] = dim2; dir1[2] = dim3;
    deviceGauge sRet(_deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim1; dir1[1] = dim3; dir1[2] = dim2;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1; dir1[2] = dim3;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim3; dir1[2] = dim1;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim1; dir1[2] = dim2;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim2; dir1[2] = dim1;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    _mul(sRet, OneOver6);
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceHyperCubicDiagonalT(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2, SBYTE dim3, SBYTE dim4)
{
    deviceGauge sRet = _makeZero<deviceGauge>();
    const SBYTE dim1234[4] = { dim1, dim2, dim3, dim4 };
    INT dir1[4];
    SBYTE dim234[3];
    for (BYTE k = 0; k < 4; ++k)
    {
        dir1[0] = dim1234[k];
        for (BYTE k2 = 0; k2 < 3; ++k2)
        {
            BYTE idx = k2 + 1 + k;
            idx = idx > 3 ? (idx - 4) : idx;
            dim234[k2] = dim1234[idx];
        }

        dir1[1] = dim234[0]; dir1[2] = dim234[1]; dir1[3] = dim234[2];
        _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[0]; dir1[2] = dim234[2]; dir1[3] = dim234[1];
        _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[0]; dir1[3] = dim234[2];
        _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[2]; dir1[3] = dim234[0];
        _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[0]; dir1[3] = dim234[1];
        _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[1]; dir1[3] = dim234[0];
        _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 4, byFieldId, dir1));
    }

    _mul(sRet, OneOver24);
    return sRet;
}

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINE_STAGGEREDGAMMA_H_

//=============================================================================
// END OF FILE
//=============================================================================
