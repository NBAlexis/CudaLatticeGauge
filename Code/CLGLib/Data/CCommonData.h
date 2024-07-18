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
#define _DC_Diri (static_cast<INT>(_constIntegers[ECI_Dir]))
#define _HC_Dir (appGetCudaHelper()->m_ConstIntegers[ECI_Dir])
#define _HC_Diri (static_cast<INT>(appGetCudaHelper()->m_ConstIntegers[ECI_Dir]))

#define _DC_Lx (_constIntegers[ECI_Lx])
#define _HC_Lx (appGetCudaHelper()->m_ConstIntegers[ECI_Lx])
#define _DC_Ly (_constIntegers[ECI_Ly])
#define _HC_Ly (appGetCudaHelper()->m_ConstIntegers[ECI_Ly])
#define _DC_Lz (_constIntegers[ECI_Lz])
#define _HC_Lz (appGetCudaHelper()->m_ConstIntegers[ECI_Lz])
#define _DC_Lt (_constIntegers[ECI_Lt])
#define _HC_Lt (appGetCudaHelper()->m_ConstIntegers[ECI_Lt])

#define _DC_Lxi (static_cast<INT>(_constIntegers[ECI_Lx]))
#define _DC_Lyi (static_cast<INT>(_constIntegers[ECI_Ly]))
#define _DC_Lzi (static_cast<INT>(_constIntegers[ECI_Lz]))
#define _DC_Lti (static_cast<INT>(_constIntegers[ECI_Lt]))
#define _HC_Lxi (static_cast<INT>(appGetCudaHelper()->m_ConstIntegers[ECI_Lx]))
#define _HC_Lyi (static_cast<INT>(appGetCudaHelper()->m_ConstIntegers[ECI_Ly]))
#define _HC_Lzi (static_cast<INT>(appGetCudaHelper()->m_ConstIntegers[ECI_Lz]))
#define _HC_Lti (static_cast<INT>(appGetCudaHelper()->m_ConstIntegers[ECI_Lt]))

#define _DC_Center SSmallInt4(_constIntegers[ECI_Center])
#define _DC_Centerx (_constSignedIntegers[ECSI_CenterX])
#define _DC_Centery (_constSignedIntegers[ECSI_CenterY])
#define _DC_Centerz (_constSignedIntegers[ECSI_CenterZ])
#define _DC_Centert (_constSignedIntegers[ECSI_CenterT])
#define _HC_Center SSmallInt4(appGetCudaHelper()->m_ConstIntegers[ECI_Center])
#define _HC_Centerx (appGetCudaHelper()->m_ConstSignedIntegers[ECSI_CenterX])
#define _HC_Centery (appGetCudaHelper()->m_ConstSignedIntegers[ECSI_CenterY])
#define _HC_Centerz (appGetCudaHelper()->m_ConstSignedIntegers[ECSI_CenterZ])
#define _HC_Centert (appGetCudaHelper()->m_ConstSignedIntegers[ECSI_CenterT])

#define _DC_Volume (_constIntegers[ECI_Volume])
#define _HC_Volume (appGetCudaHelper()->m_ConstIntegers[ECI_Volume])
#define _DC_Volume_xyz (_constIntegers[ECI_Volume_xyz])
#define _HC_Volume_xyz (appGetCudaHelper()->m_ConstIntegers[ECI_Volume_xyz])

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

//Remember not to make plaqutte count device, 
//Because the plaqutte count can be figured out after boundary condition (so not a constant)
#define _HC_PlaqutteCount (appGetCudaHelper()->m_ConstIntegers[ECI_PlaqutteCount])

#define _HC_LinkCount (appGetCudaHelper()->m_ConstIntegers[ECI_LinkCount])
#define _HC_ThreadConstraint (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadConstaint])
#define _HC_ThreadConstraintX (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadConstaintX])
#define _HC_ThreadConstraintY (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadConstaintY)
#define _HC_ThreadConstraintZ (appGetCudaHelper()->m_ConstIntegers[ECI_ThreadConstaintZ])
#define _HC_SummationDecompose (appGetCudaHelper()->m_ConstIntegers[ECI_SummationDecompose])

#define _DC_Seed (_constIntegers[ECI_RandomSeed])
#define _HC_Seed (appGetCudaHelper()->m_ConstIntegers[ECI_RandomSeed])
#define _DC_ExpPrecision (_constIntegers[ECI_ExponentPrecision])
#define _HC_ExpPrecision (appGetCudaHelper()->m_ConstIntegers[ECI_ExponentPrecision])

#define _DC_ActionListL (_constIntegers[ECI_ActionListLength])
#define _HC_ActionListL (appGetCudaHelper()->m_ConstIntegers[ECI_ActionListLength])

#define _D_ComplexThreadBuffer (appGetCudaHelper()->m_pComplexBufferThreadCount)
#define _D_RealThreadBuffer (appGetCudaHelper()->m_pRealBufferThreadCount)

//#define _D_IndexBuffer (appGetCudaHelper()->m_pIndexBuffer)


#define _DC_InverseSqrtLink16 (_constFloats[ECF_InverseSqrtLink16])

#define _DC_ALog (_constIntegers[ECI_UseLogADefinition])
#define _HC_ALog (appGetCudaHelper()->m_ConstIntegers[ECI_UseLogADefinition])

#pragma endregion

__BEGIN_NAMESPACE

inline class CCudaHelper* appGetCudaHelper();

enum 
{
    kMaxActionCount = 16, 
    kMaxPlaqutteCache = 32,
};

__DEFINE_ENUM(EFieldType,

    EFT_GaugeSU3,

    EFT_GaugeU1,
    EFT_GaugeReal,

    EFT_GaugeSU2,
    EFT_GaugeSUN,
    EFT_GaugeSU4,
    EFT_GaugeSU5,
    EFT_GaugeSU6,
    EFT_GaugeSU7,
    EFT_GaugeSU8,

    EFT_FermionWilsonSquareSU3,
    EFT_FermionStaggeredSU3,
    EFT_FermionStaggeredU1,

    EFT_BosonComplex,
    EFT_BosonComplexVector2,
    EFT_BosonComplexVector3,
    EFT_BosonComplexVector4,
    EFT_BosonComplexVector5,
    EFT_BosonComplexVector6,
    EFT_BosonComplexVector7,
    EFT_BosonComplexVector8,

    EFT_Max,
    EFT_ForceDword = 0x7fffffff,

    )


__DEFINE_ENUM(ESolverPhase,
    ESP_StartTrajectory,
    ESP_InTrajectory,
    ESP_EndTrajectory,
    ESP_Once,

    ESP_ForceDWORD = 0x7fffffff,
    )


#define _SSMALLINT4(intd) ((SSmallInt4*)(&intd))

#define _CSSMALLINT4(intd) ((const SSmallInt4*)(&intd))

#if defined(__cplusplus)
    extern "C" {
#endif /* __cplusplus */
    //instead of int4
    struct CLGAPI SSmallInt4
    {
        __device__ __host__ SSmallInt4() {}
        __device__ __host__ SSmallInt4(const SSmallInt4& other) : x(other.x), y(other.y), z(other.z), w(other.w) {}
        __device__ __host__ SSmallInt4(SBYTE inx, SBYTE iny, SBYTE inz, SBYTE inw) : x(inx), y(iny), z(inz), w(inw) {}
        __device__ __host__ SSmallInt4(UINT uiData) : m_uiData(uiData) {}

        union
        {
            UINT m_uiData;
            SBYTE m_byData4[4];
            struct 
            {
                SBYTE x, y, z, w;
            };
        };

        __device__ __inline__ INT X() const
        {
            return static_cast<INT>(x);
        }
        __device__ __inline__ INT Y() const
        {
            return static_cast<INT>(y);
        }
        __device__ __inline__ INT Z() const
        {
            return static_cast<INT>(z);
        }
        __device__ __inline__ INT T() const
        {
            return static_cast<INT>(w);
        }

        __device__ __inline__ UBOOL IsOdd() const
        {
            return (x + y + z + w) & 1;
        }

        /**
         * eta_{mu}(n) = (-1)^{sum (nu<mu)}
         */
        __device__ __inline__ UBOOL EtaOdd(BYTE nu) const
        {
            SWORD sSum = 0;
            for (BYTE byIdx = 0; byIdx < nu && byIdx < 4; ++byIdx)
            {
                sSum += m_byData4[byIdx];
            }
            return sSum & 1;
        }

        __device__ __inline__ void Add(const SSmallInt4& other)
        {
            x = x + other.x;
            y = y + other.y;
            z = z + other.z;
            w = w + other.w;
        }

        __device__ __inline__ SSmallInt4 AddC(const SSmallInt4& other) const
        {
            SSmallInt4 ret;
            ret.x = x + other.x;
            ret.y = y + other.y;
            ret.z = z + other.z;
            ret.w = w + other.w;
            return ret;
        }

        __device__ __inline__ void Sub(const SSmallInt4& other)
        {
            x = x - other.x;
            y = y - other.y;
            z = z - other.z;
            w = w - other.w;
        }

        __device__ __inline__ SSmallInt4 SubC(const SSmallInt4& other) const
        {
            SSmallInt4 ret;
            ret.x = x - other.x;
            ret.y = y - other.y;
            ret.z = z - other.z;
            ret.w = w - other.w;
            return ret;
        }

        __device__ __inline__ UINT _deviceToSiteIndex() const
        {
            return static_cast<UINT>(((x * _DC_Ly + y) * _DC_Lz + z) * _DC_Lt + w);
        }

        __host__ __inline__ UINT _hostToSiteIndex() const
        {
            return static_cast<UINT>(((x * _HC_Ly + y) * _HC_Lz + z) * _HC_Lt + w);
        }

        __device__ __inline__ UBOOL Out() const
        {
            return x < 0 || x >= _DC_Lxi
                || y < 0 || y >= _DC_Lyi
                || z < 0 || z >= _DC_Lzi
                || w < 0 || w >= _DC_Lti;
        }
    };
#if defined(__cplusplus)
}
#endif /* __cplusplus */

template<>
inline CCString appToString(const SSmallInt4& content)
{
    CCString sret;
    sret.Format(_T("[%d, %d, %d, %d]"),
        static_cast<INT>(content.x),
        static_cast<INT>(content.y),
        static_cast<INT>(content.z),
        static_cast<INT>(content.w));
    return sret;
}

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
__device__ __inline__ static UINT _deviceGetSiteIndex(const SSmallInt4& coord)
{
    return static_cast<UINT>(coord.x * _DC_MultX + coord.y * _DC_MultY + coord.z * _DC_MultZ + coord.w);
}
__device__ __inline__ static UINT _deviceGetLinkIndex(UINT siteIndex, BYTE dir)
{
    return siteIndex * _DC_Dir + dir;
}
__device__ __inline__ static UINT _deviceGetLinkIndex(const UINT* coord, BYTE dir)
{
    return _deviceGetSiteIndex(coord) * _DC_Dir + dir;
}

/**
* for site, dir_plus_one = 0
* for link, dir_plus_one = dir + 1
*/
__device__ __inline__ static UINT _deviceGetFatIndex(const UINT* coord, BYTE dir_plus_one)
{
    return _deviceGetSiteIndex(coord) * (_DC_Dir + 1) + dir_plus_one;
}

__device__ __inline__ static UINT _deviceGetFatIndex(UINT uiSiteIndex, BYTE dir_plus_one)
{
    return uiSiteIndex * (_DC_Dir + 1) + dir_plus_one;
}

/**
* SSmallInt4.xyzw = x, y, z, t
*/
__device__ __inline__ static SSmallInt4 __deviceSiteIndexToInt4(UINT siteIndex)
{
    SSmallInt4 xyzt;
    xyzt.x = static_cast<SBYTE>(siteIndex / _DC_MultX);
    xyzt.y = static_cast<SBYTE>((siteIndex % _DC_MultX) / _DC_MultY);
    xyzt.z = static_cast<SBYTE>((siteIndex % _DC_MultY) / _DC_MultZ);
    xyzt.w = static_cast<SBYTE>((siteIndex % _DC_MultZ));
    return xyzt;
}

__device__ __inline__ static SSmallInt4 __deviceLinkIndexToInt4(UINT linkIndex)
{
    return __deviceSiteIndexToInt4(linkIndex / _DC_Dir);
}

__device__ __inline__ static SSmallInt4 __deviceFatIndexToInt4(UINT fatIndex)
{
    return __deviceSiteIndexToInt4(fatIndex / (_DC_Dir + 1));
}

#pragma endregion

#pragma region Host functions

inline static SSmallInt4 __hostSiteIndexToInt4(UINT siteIndex)
{
    SSmallInt4 xyzt;
    xyzt.x = static_cast<SBYTE>(siteIndex / _HC_MultX);
    xyzt.y = static_cast<SBYTE>((siteIndex % _HC_MultX) / _HC_MultY);
    xyzt.z = static_cast<SBYTE>((siteIndex % _HC_MultY) / _HC_MultZ);
    xyzt.w = static_cast<SBYTE>((siteIndex % _HC_MultZ));
    return xyzt;
}

inline static UINT _hostGetSiteIndex(const SSmallInt4& coord)
{
    return static_cast<UINT>(coord.x * _HC_MultX + coord.y * _HC_MultY + coord.z * _HC_MultZ + coord.w);
}

#pragma endregion

//at most 8 tags
enum
{
    _kDaggerOrOpposite  = 0x01,
    _kDirichlet         = 0x02,
    _kGlue              = 0x04,
};

#define _SINDEX(longlongdata) ((SIndex*)&longlongdata)
#define _CSINDEX(longlongdata) ((const SIndex*)&longlongdata)

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    typedef struct
    {
        ULONGLONG m_data[kMaxPlaqutteCache];
    } SCachedIndexArray;

    struct alignas(8) SIndex
    {
        __device__ SIndex()
            : m_uiSiteIndex(0)
            , m_byDir(0)
            , m_byTag(0)
            , m_byBoundaryFieldId(0)
            , m_byReginId(0)
        {

        }

        __device__ SIndex(UINT uiIndex, BYTE dir = 0, BYTE indexTag = 0, BYTE bcField = 0, BYTE byRegionId = 0)
            : m_uiSiteIndex(uiIndex)
            , m_byDir(dir)
            , m_byTag(indexTag)
            , m_byBoundaryFieldId(bcField)
            , m_byReginId(byRegionId)
        {

        }

        __device__ SIndex(const SIndex& other)
            : m_uiSiteIndex(other.m_uiSiteIndex)
            , m_byDir(other.m_byDir)
            , m_byTag(other.m_byTag)
            , m_byBoundaryFieldId(other.m_byBoundaryFieldId)
            , m_byReginId(other.m_byReginId)
        {

        }

        __device__ __inline__ void DebugPrint() const
        {
            const SSmallInt4 xyzt = __deviceSiteIndexToInt4(m_uiSiteIndex);
            printf("%s(xyzt:%d,%d,%d,%d)_(%x)%s\n", NeedToOpposite() ? "-" : "", xyzt.x, xyzt.y, xyzt.z, xyzt.w, m_byDir, NeedToDagger() ? "^-1" : "");
        }

        __device__ __inline__ UBOOL NeedToDagger() const { return 0 != (_kDaggerOrOpposite & m_byTag); }
        __device__ __inline__ UBOOL NeedToOpposite() const { return 0 != (_kDaggerOrOpposite & m_byTag); }
        __device__ __inline__ UBOOL IsDirichlet() const { return 0 != (_kDirichlet & m_byTag); }

        __device__ __inline__ SIndex DaggerC() const
        {
            SIndex ret = *this;
            ret.m_byTag = ret.m_byTag ^ _kDaggerOrOpposite;
            return ret;
        }

        union 
        {
            ULONGLONG m_ullData;

            struct
            {
                UINT m_uiSiteIndex;
                BYTE m_byDir;
                BYTE m_byTag;

                //NOTE, THIS IS NOT USING (AND SHOULD NOT BE USED)
                BYTE m_byBoundaryFieldId;

                /**
                * For miscellaneous usage
                */
                BYTE m_byReginId;
            };
        };
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

#if !_CLG_DOUBLEFLOAT
    static DOUBLE m_fBeta;
    static DOUBLE m_fOmega;
#else
    static Real m_fBeta;
    static Real m_fOmega;
#endif
    static Real m_fKai;

    static UBOOL m_bStoreStaple;
    static UBOOL m_bStoreLastSolution;
    static UBOOL m_bStochasticGaussian;

    //Used in rotating frame. Since the fermion fields are copied,
    //it is convinient to set all parameters at one place
    
    //static SSmallInt4 m_sCenter;

    //Use for acceleration, Since the fermion fields are copied,
    static Real m_fG;

    static UINT m_uiMaxThreadPerBlock;

    //No room for other paramters in solvers, so put it here, maybe move to somewhere else later
    static Real m_fShiftedMass;

    //External Electric - Magnetic field
    static Real m_fBz;
    static Real m_fEz;
};

__END_NAMESPACE

#endif //#ifndef _CCOMMONDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================