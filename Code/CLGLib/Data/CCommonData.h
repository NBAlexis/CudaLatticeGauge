//=============================================================================
// FILENAME : CCommonData.h
// 
// DESCRIPTION:
// This is the class for the common data
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CCOMMONDATA_H_
#define _CCOMMONDATA_H_

#pragma region constants

//DC for device constant
//HC for host constant

#define _DC_Dim (_constIntegers[ECI_Dim])
#define _HC_Dim (appGetCudaHelper()->m_ConstIntegers[ECI_Dim])

#define _DC_Dir (_constIntegers[ECI_Dir])
#define _HC_Dir (appGetCudaHelper()->m_ConstIntegers[ECI_Dir])

#define _DC_Lx (_constIntegers[ECI_Lx])
#define _HC_Lx (appGetCudaHelper()->m_ConstIntegers[ECI_Lx])
#define _DC_Ly (_constIntegers[ECI_Ly])
#define _HC_Ly (appGetCudaHelper()->m_ConstIntegers[ECI_Ly])
#define _DC_Lz (_constIntegers[ECI_Lz])
#define _HC_Lz (appGetCudaHelper()->m_ConstIntegers[ECI_Lz])
#define _DC_Lt (_constIntegers[ECI_Lt])
#define _HC_Lt (appGetCudaHelper()->m_ConstIntegers[ECI_Lt])

#define _DC_Volumn (_constIntegers[ECI_Volumn])
#define _HC_Volumn (appGetCudaHelper()->m_ConstIntegers[ECI_Volumn])

#define _DC_MultX (_constIntegers[ECI_MultX])
#define _HC_MultX (appGetCudaHelper()->m_ConstIntegers[ECI_MultX])
#define _DC_MultY (_constIntegers[ECI_MultY])
#define _HC_MultY (appGetCudaHelper()->m_ConstIntegers[ECI_MultY])
#define _DC_MultZ (_constIntegers[ECI_MultZ])
#define _HC_MultZ (appGetCudaHelper()->m_ConstIntegers[ECI_MultZ])

#define _DC_GridDimZT (_constIntegers[ECI_GridDimZT])
#define _HC_GridDimZT (appGetCudaHelper()->m_ConstIntegers[ECI_GridDimZT])

#define _DC_DecompX (_constIntegers[ECI_DecompX])
#define _HC_DecompX (appGetCudaHelper()->m_ConstIntegers[ECI_DecompX])
#define _DC_DecompY (_constIntegers[ECI_DecompY])
#define _HC_DecompY (appGetCudaHelper()->m_ConstIntegers[ECI_DecompY])
#define _DC_DecompZ (_constIntegers[ECI_DecompZ])
#define _HC_DecompZ (appGetCudaHelper()->m_ConstIntegers[ECI_DecompZ])
#define _DC_DecompLx (_constIntegers[ECI_DecompLx])
#define _HC_DecompLx (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLx])
#define _DC_DecompLy (_constIntegers[ECI_DecompLy])
#define _HC_DecompLy (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLy])
#define _DC_DecompLz (_constIntegers[ECI_DecompLz])
#define _HC_DecompLz (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLz])
#define _DC_ThreadCountPerBlock (_constIntegers[ECI_ThreadCountPerBlock])
#define _HC_ThreadCountPerBlock (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadCountPerBlock])

#define _HC_SUN (appGetCudaHelper()->m_ConstIntegers[ECI_SUN])
#define _HC_PlaqutteCount (appGetCudaHelper()->m_ConstIntegers[ECI_PlaqutteCount])
#define _HC_LinkCount (appGetCudaHelper()->m_ConstIntegers[ECI_LinkCount])
#define _HC_ThreadConstraint (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadConstaint])
#define _HC_SummationDecompose (appGetCudaHelper()->m_ConstIntegers[ECI_SummationDecompose])

#define _DC_Seed (_constIntegers[ECI_RandomSeed])
#define _HC_Seed (appGetCudaHelper()->m_ConstIntegers[ECI_RandomSeed])
#define _DC_ExpPrecision (_constIntegers[ECI_ExponentPrecision])

#define _DC_ActionListL (_constIntegers[ECI_ActionListLength])
#define _HC_ActionListL (appGetCudaHelper()->m_ConstIntegers[ECI_ActionListLength])

#define _D_ComplexThreadBuffer (appGetCudaHelper()->m_pComplexBufferThreadCount)
#define _D_RealThreadBuffer (appGetCudaHelper()->m_pRealBufferThreadCount)
//#define _D_IndexBuffer (appGetCudaHelper()->m_pIndexBuffer)


#define _DC_InverseSqrtLink16 (_constFloats[ECF_InverseSqrtLink16])

#pragma endregion

__BEGIN_NAMESPACE

enum { kMaxFieldCount = 8, kMaxActionCount = 8 };

__DEFINE_ENUM(EFieldType,

    EFT_GaugeSU3,
    EFT_FermionWilsonSquareSU3,
    EFT_Max,
    EFT_ForceDword = 0x7fffffff,

    )



//====================================================
// Some common structures
//====================================================
#pragma region index functions

//========================================
// implement after CLatticeData is known

__device__ __inline__ static UINT _deviceGetSiteIndex(const UINT* coord)
{
    return coord[0] * _DC_MultX + coord[1] * _DC_MultY + coord[2] * _DC_MultZ + coord[3];
}
__device__ __inline__ static UINT _deviceGetLinkIndex(UINT siteIndex, UINT dir)
{
    return siteIndex * _DC_Dir + dir;
}
__device__ __inline__ static UINT _deviceGetLinkIndex(const UINT* coord, UINT dir)
{
    return _deviceGetSiteIndex(coord) * _DC_Dir + dir;
}

/**
* for site, dir_plus_one = 0
* for link, dir_plus_one = dir + 1
*/
__device__ __inline__ static UINT _deviceGetFatIndex(const UINT* coord, UINT dir_plus_one)
{
    return _deviceGetSiteIndex(coord) * (_DC_Dir + 1) + dir_plus_one;
}

__device__ __inline__ static UINT _deviceGetFatIndex(UINT uiSiteIndex, UINT dir_plus_one)
{
    return uiSiteIndex * (_DC_Dir + 1) + dir_plus_one;
}

/**
* int4.xyzw = x, y, z, t
*/
__device__ __inline__ static int4 __deviceSiteIndexToInt4(UINT siteIndex)
{
    int4 xyzt;
    xyzt.x = siteIndex / _DC_MultX;
    xyzt.y = (siteIndex % _DC_MultX) / _DC_MultY;
    xyzt.z = (siteIndex % _DC_MultY) / _DC_MultZ;
    xyzt.w = (siteIndex % _DC_MultZ);
    return xyzt;
}
__device__ __inline__ static int4 __deviceLinkIndexToInt4(UINT linkIndex)
{
    return __deviceSiteIndexToInt4(linkIndex / _DC_Dir);
}
__device__ __inline__ static int4 __deviceFatIndexToInt4(UINT fatIndex)
{
    return __deviceSiteIndexToInt4(fatIndex / (_DC_Dir + 1));
}

#pragma endregion

enum
{
    _kDagger = 0x01,
};

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    struct alignas(8) SIndex
    {
        __device__ SIndex()
            : m_uiSiteIndex(0)
            , m_byDir(0)
            , m_byTag(0)
            , m_byBoundaryFieldId(0)
        {

        }

        __device__ SIndex(UINT uiIndex, BYTE dir = 0, BYTE indexTag = 0, BYTE bcField = 0)
            : m_uiSiteIndex(uiIndex)
            , m_byDir(dir)
            , m_byTag(indexTag)
            , m_byBoundaryFieldId(bcField)
        {

        }

        __device__ SIndex(const SIndex& other)
            : m_uiSiteIndex(other.m_uiSiteIndex)
            , m_byDir(other.m_byDir)
            , m_byTag(other.m_byTag)
            , m_byBoundaryFieldId(other.m_byBoundaryFieldId)
        {

        }

        __device__ __inline__ void DebugPrint() const
        {
            int4 xyzt = __deviceSiteIndexToInt4(m_uiSiteIndex);
            printf("(xyzt:%d,%d,%d,%d)_(%x)%s\n", xyzt.x, xyzt.y, xyzt.z, xyzt.w, m_byDir, NeedToDagger() ? "^-1" : "");
        }

        __device__ __inline__ UBOOL NeedToDagger() const { return 0 != (_kDagger & m_byTag); }
        __device__ __inline__ UBOOL NeedBoundaryField() const { return 0 != m_byBoundaryFieldId; }

        UINT m_uiSiteIndex;
        BYTE m_byDir;
        BYTE m_byTag;
        BYTE m_byBoundaryFieldId;
        BYTE m_byUnused;
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__device__ __inline__ static SIndex __deviceLinkIndexToSIndex(UINT linkIndex)
{
    return SIndex(linkIndex / _DC_Dir, static_cast<BYTE>(linkIndex % _DC_Dir));
}

__device__ __inline__ static SIndex __deviceSiteIndexToSIndex(UINT siteIndex)
{
    return SIndex(siteIndex);
}

//Those are not constants, but commonly used parameters
class CLGAPI CCommonData
{
public:

    static Real m_fBeta;
    static Real m_fKai;
};

__END_NAMESPACE

#endif //#ifndef _CCOMMONDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================