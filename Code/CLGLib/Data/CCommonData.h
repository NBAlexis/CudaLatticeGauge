//=============================================================================
// FILENAME : CCommonData.h
// 
// DESCRIPTION:
// This is the class for the common data
//
// REVISION:
//  [mm/dd/yy]
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CCOMMONDATA_H_
#define _CCOMMONDATA_H_

#pragma region constants

//DC for device constant
//HC for host constant
#if _CLG_ASSUME_SQUARE_LATTICE
#define _DC_Dim (4)
#define _HC_Dim (4)
#define _DC_Dir (4)
#define _DC_Diri (4)
#define _HC_Dir (4)
#define _HC_Diri (4)
#else
#define _DC_Dim (_constIntegers[ECI_Dim])
#define _HC_Dim (appGetCudaHelper()->m_ConstIntegers[ECI_Dim])
#define _DC_Dir (_constIntegers[ECI_Dir])
#define _DC_Diri (static_cast<INT>(_constIntegers[ECI_Dir]))
#define _HC_Dir (appGetCudaHelper()->m_ConstIntegers[ECI_Dir])
#define _HC_Diri (static_cast<INT>(appGetCudaHelper()->m_ConstIntegers[ECI_Dir]))
#endif

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
#define _DC_VolumeHalf (_constIntegers[ECI_VolumeHalf])
#define _HC_VolumeHalf (appGetCudaHelper()->m_ConstIntegers[ECI_VolumeHalf])
#define _DC_Volume_xyz (_constIntegers[ECI_Volume_xyz])
#define _HC_Volume_xyz (appGetCudaHelper()->m_ConstIntegers[ECI_Volume_xyz])
#define _DC_Volume_xyt (_constIntegers[ECI_Volume_xyt])
#define _HC_Volume_xyt (appGetCudaHelper()->m_ConstIntegers[ECI_Volume_xyt])
#define _DC_Volume_xzt (_constIntegers[ECI_Volume_xzt])
#define _HC_Volume_xzt (appGetCudaHelper()->m_ConstIntegers[ECI_Volume_xzt])
#define _DC_Volume_yzt (_constIntegers[ECI_Volume_yzt])
#define _HC_Volume_yzt (appGetCudaHelper()->m_ConstIntegers[ECI_Volume_yzt])

#define _DC_MultX (_constIntegers[ECI_MultX])
#define _HC_MultX (appGetCudaHelper()->m_ConstIntegers[ECI_MultX])
#define _DC_MultY (_constIntegers[ECI_MultY])
#define _HC_MultY (appGetCudaHelper()->m_ConstIntegers[ECI_MultY])
#define _DC_MultZ (_constIntegers[ECI_MultZ])
#define _HC_MultZ (appGetCudaHelper()->m_ConstIntegers[ECI_MultZ])

#define _DC_GridDimZT (_constIntegers[ECI_GridDimZT])
#define _HC_GridDimZT (appGetCudaHelper()->m_ConstIntegers[ECI_GridDimZT])

//Multi-GPU (Phase 1): global lattice, process grid, per-rank global offset, halo width.
//On single-GPU builds these are set to the full lattice / [1,1,1,1] / 0 / 2.
#define _DC_GlobalLx (_constIntegers[ECI_GlobalLx])
#define _HC_GlobalLx (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalLx])
#define _DC_GlobalLy (_constIntegers[ECI_GlobalLy])
#define _HC_GlobalLy (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalLy])
#define _DC_GlobalLz (_constIntegers[ECI_GlobalLz])
#define _HC_GlobalLz (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalLz])
#define _DC_GlobalLt (_constIntegers[ECI_GlobalLt])
#define _HC_GlobalLt (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalLt])

#define _DC_OffsetX (_constIntegers[ECI_GlobalOffsetX])
#define _HC_OffsetX (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalOffsetX])
#define _DC_OffsetY (_constIntegers[ECI_GlobalOffsetY])
#define _HC_OffsetY (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalOffsetY])
#define _DC_OffsetZ (_constIntegers[ECI_GlobalOffsetZ])
#define _HC_OffsetZ (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalOffsetZ])
#define _DC_OffsetT (_constIntegers[ECI_GlobalOffsetT])
#define _HC_OffsetT (appGetCudaHelper()->m_ConstIntegers[ECI_GlobalOffsetT])

#define _DC_HaloWidth (_constIntegers[ECI_HaloWidth])
#define _HC_HaloWidth (appGetCudaHelper()->m_ConstIntegers[ECI_HaloWidth])

//Process grid factors per direction. >1 means this direction is split across
//ranks (its out-of-lattice neighbours live on another rank and must be fetched
//from halo storage rather than periodic-wrapped within this rank). On single-GPU
//builds every factor is 1, so the halo path is never taken.
#define _DC_GpuGridX (_constIntegers[ECI_GpuGridX])
#define _HC_GpuGridX (appGetCudaHelper()->m_ConstIntegers[ECI_GpuGridX])
#define _DC_GpuGridY (_constIntegers[ECI_GpuGridY])
#define _HC_GpuGridY (appGetCudaHelper()->m_ConstIntegers[ECI_GpuGridY])
#define _DC_GpuGridZ (_constIntegers[ECI_GpuGridZ])
#define _HC_GpuGridZ (appGetCudaHelper()->m_ConstIntegers[ECI_GpuGridZ])
#define _DC_GpuGridT (_constIntegers[ECI_GpuGridT])
#define _HC_GpuGridT (appGetCudaHelper()->m_ConstIntegers[ECI_GpuGridT])

//NOTE: the local->global coordinate helper lives in Data/Lattice/CIndexData.h
//(_deviceSIndexToGlobalInt4, on the baked 32-bit global coordinate table).

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

#define _DC_DecompBlock (_constIntegers[ECI_DecompAllBlock])
#define _DC_DecompThread (_constIntegers[ECI_DecompAllThread])
#define _DC_DecompBlockDir (_constIntegers[ECI_DecompAllBlockDir])
#define _DC_DecompThreadDir (_constIntegers[ECI_DecompAllThreadDir])

#define _HC_DecompBlock (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllBlock])
#define _HC_DecompThread (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllThread])
#define _HC_DecompBlockDir (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllBlockDir])
#define _HC_DecompThreadDir (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllThreadDir])

#define _DC_DecompBlockHalf (_constIntegers[ECI_DecompAllBlockHalf])
#define _DC_DecompThreadHalf (_constIntegers[ECI_DecompAllThreadHalf])
#define _DC_DecompBlockDirHalf (_constIntegers[ECI_DecompAllBlockDirHalf])
#define _DC_DecompThreadDirHalf (_constIntegers[ECI_DecompAllThreadDirHalf])

#define _HC_DecompBlockHalf (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllBlockHalf])
#define _HC_DecompThreadHalf (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllThreadHalf])
#define _HC_DecompBlockDirHalf (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllBlockDirHalf])
#define _HC_DecompThreadDirHalf (appGetCudaHelper()->m_ConstIntegers[ECI_DecompAllThreadDirHalf])

#define _DC_DecompX3D (_constIntegers[ECI_DecompX3D])
#define _HC_DecompX3D (appGetCudaHelper()->m_ConstIntegers[ECI_DecompX3D])
#define _DC_DecompY3D (_constIntegers[ECI_DecompY3D])
#define _HC_DecompY3D (appGetCudaHelper()->m_ConstIntegers[ECI_DecompY3D])
#define _DC_DecompZ3D (_constIntegers[ECI_DecompZ3D])
#define _HC_DecompZ3D (appGetCudaHelper()->m_ConstIntegers[ECI_DecompZ3D])
#define _DC_DecompLx3D (_constIntegers[ECI_DecompLx3D])
#define _HC_DecompLx3D (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLx3D])
#define _DC_DecompLy3D (_constIntegers[ECI_DecompLy3D])
#define _HC_DecompLy3D (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLy3D])
#define _DC_DecompLz3D (_constIntegers[ECI_DecompLz3D])
#define _HC_DecompLz3D (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLz3D])

#define _DC_DecompX3DXYT (_constIntegers[ECI_DecompX3DXYT])
#define _HC_DecompX3DXYT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompX3DXYT])
#define _DC_DecompY3DXYT (_constIntegers[ECI_DecompY3DXYT])
#define _HC_DecompY3DXYT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompY3DXYT])
#define _DC_DecompZ3DXYT (_constIntegers[ECI_DecompZ3DXYT])
#define _HC_DecompZ3DXYT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompZ3DXYT])
#define _DC_DecompLx3DXYT (_constIntegers[ECI_DecompLx3DXYT])
#define _HC_DecompLx3DXYT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLx3DXYT])
#define _DC_DecompLy3DXYT (_constIntegers[ECI_DecompLy3DXYT])
#define _HC_DecompLy3DXYT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLy3DXYT])
#define _DC_DecompLz3DXYT (_constIntegers[ECI_DecompLz3DXYT])
#define _HC_DecompLz3DXYT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLz3DXYT])

#define _DC_DecompX3DXZT (_constIntegers[ECI_DecompX3DXZT])
#define _HC_DecompX3DXZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompX3DXZT])
#define _DC_DecompY3DXZT (_constIntegers[ECI_DecompY3DXZT])
#define _HC_DecompY3DXZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompY3DXZT])
#define _DC_DecompZ3DXZT (_constIntegers[ECI_DecompZ3DXZT])
#define _HC_DecompZ3DXZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompZ3DXZT])
#define _DC_DecompLx3DXZT (_constIntegers[ECI_DecompLx3DXZT])
#define _HC_DecompLx3DXZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLx3DXZT])
#define _DC_DecompLy3DXZT (_constIntegers[ECI_DecompLy3DXZT])
#define _HC_DecompLy3DXZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLy3DXZT])
#define _DC_DecompLz3DXZT (_constIntegers[ECI_DecompLz3DXZT])
#define _HC_DecompLz3DXZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLz3DXZT])

#define _DC_DecompX3DYZT (_constIntegers[ECI_DecompX3DYZT])
#define _HC_DecompX3DYZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompX3DYZT])
#define _DC_DecompY3DYZT (_constIntegers[ECI_DecompY3DYZT])
#define _HC_DecompY3DYZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompY3DYZT])
#define _DC_DecompZ3DYZT (_constIntegers[ECI_DecompZ3DYZT])
#define _HC_DecompZ3DYZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompZ3DYZT])
#define _DC_DecompLx3DYZT (_constIntegers[ECI_DecompLx3DYZT])
#define _HC_DecompLx3DYZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLx3DYZT])
#define _DC_DecompLy3DYZT (_constIntegers[ECI_DecompLy3DYZT])
#define _HC_DecompLy3DYZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLy3DYZT])
#define _DC_DecompLz3DYZT (_constIntegers[ECI_DecompLz3DYZT])
#define _HC_DecompLz3DYZT (appGetCudaHelper()->m_ConstIntegers[ECI_DecompLz3DYZT])

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


#define _HC_GaugeMomentumFactor (appGetCudaHelper()->m_ConstFloats[ECF_GaugeMomentumFactor])

#define _DC_ALog (_constIntegers[ECI_UseLogADefinition])
#define _HC_ALog (appGetCudaHelper()->m_ConstIntegers[ECI_UseLogADefinition])

#define _HC_Profiler (appGetCudaHelper()->m_ConstIntegers[ECI_Profiler])

#define _DC_MILCSTAGGEREDPHASE (_constIntegers[ECI_MILC_StaggeredPhase])
#define _HC_MILCSTAGGEREDPHASE (appGetCudaHelper()->m_ConstIntegers[ECI_MILC_StaggeredPhase])

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
    EFT_GaugeSU3_12,

    EFT_GaugeU1,
    EFT_GaugeReal,

    EFT_GaugeZ2,
    EFT_GaugeZ3,
    EFT_GaugeZ4,
    EFT_GaugeZ5,
    EFT_GaugeZ6,

    EFT_GaugeD3,
    EFT_GaugeD4,
    EFT_GaugeD8,

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
    EFT_FermionStaggeredSU2,
    EFT_FermionStaggeredSU4,
    EFT_FermionStaggeredSU5,
    EFT_FermionStaggeredSU6,
    EFT_FermionStaggeredSU7,
    EFT_FermionStaggeredSU8,

    EFT_BosonComplex,
    EFT_BosonComplexVector2,
    EFT_BosonComplexVector3,
    EFT_BosonComplexVector4,
    EFT_BosonComplexVector5,
    EFT_BosonComplexVector6,
    EFT_BosonComplexVector7,
    EFT_BosonComplexVector8,

    EFT_BosonReal,

    EFT_Tensor2Real,
    EFT_Tensor2Complex,
    EFT_Tensor2SU2,
    EFT_Tensor2SU3,
    EFT_Tensor2Z3,

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


__DEFINE_ENUM(ECacheCall,
    ECC_BeforeGaugeUpdate,
    ECC_BeforeFermionUpdateBeforeSmearing,
    ECC_BeforeFermionUpdateAfterSmearing,
    ECC_BeforeAllUpdateBeforeSmearing,
    ECC_BeforeAllUpdateAfterSmearing,
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
        __device__ __host__ SSmallInt4(SCHAR inx, SCHAR iny, SCHAR inz, SCHAR inw) : x(inx), y(iny), z(inz), w(inw) {}
        __device__ __host__ SSmallInt4(UINT uiData) : m_uiData(uiData) {}

        union
        {
            UINT m_uiData;
            SCHAR m_byData4[4];
            struct 
            {
                SCHAR x, y, z, w;
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
         * Only use when baking!
         * 
         * baking eta after constants are set
         */
        __device__ __inline__ UBOOL EtaOdd(BYTE nu) const
        {
            SWORD sSum = 0;
            if (_DC_MILCSTAGGEREDPHASE)
            {
                switch (nu)
                {
                    case 0:
                        sSum = m_byData4[3];
                        break;
                    case 1:
                        sSum = m_byData4[3] + m_byData4[0];
                        break;
                    case 2:
                        sSum = m_byData4[3] + m_byData4[0] + m_byData4[1];
                        break;
                    case 3:
                        sSum = 0;
                        break;
                    default:
                        sSum = m_byData4[0] + m_byData4[1] + m_byData4[2] + m_byData4[3];
                        break;
                }

            }
            else
            {
                for (BYTE byIdx = 0; byIdx < nu && byIdx < 4; ++byIdx)
                {
                    sSum += m_byData4[byIdx];
                }
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

        __host__ __inline__ UBOOL operator==(const SSmallInt4& Other) const
        {
            return x == Other.x && y == Other.y && z == Other.z && w == Other.w;
        }
    };

    //32-bit int4, used for global coordinates which may exceed SCHAR range
    struct CLGAPI SInt4
    {
        __device__ __host__ SInt4() {}
        __device__ __host__ SInt4(const SInt4& other) : x(other.x), y(other.y), z(other.z), w(other.w) {}
        __device__ __host__ SInt4(INT inx, INT iny, INT inz, INT inw) : x(inx), y(iny), z(inz), w(inw) {}

        INT x, y, z, w;

        __host__ __inline__ UBOOL operator==(const SInt4& Other) const
        {
            return x == Other.x && y == Other.y && z == Other.z && w == Other.w;
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
#if _CLG_ASSUME_SQUARE_LATTICE
    return (siteIndex << 2U) | dir;
#else
    return siteIndex * _DC_Dir + dir;
#endif
}
__device__ __inline__ static UINT _deviceGetLinkIndex(const UINT* coord, BYTE dir)
{
#if _CLG_ASSUME_SQUARE_LATTICE
    return (_deviceGetSiteIndex(coord) << 2U) | dir;
#else
    return _deviceGetSiteIndex(coord) * _DC_Dir + dir;
#endif
}

__device__ __inline__ static UINT _deviceGetSiteFromLink(UINT linkIndex)
{
#if _CLG_ASSUME_SQUARE_LATTICE
    return linkIndex >> 2U;
#else
    return linkIndex / _DC_Dir;
#endif
}

__device__ __inline__ static UINT _deviceGetDirFromLink(UINT linkIndex)
{
#if _CLG_ASSUME_SQUARE_LATTICE
    return linkIndex | 3U;
#else
    return linkIndex % _DC_Dir;
#endif
}



/**
* SSmallInt4.xyzw = x, y, z, t
*/
__device__ __inline__ static SSmallInt4 __deviceSiteIndexToInt4(UINT siteIndex);
//{
//    SSmallInt4 xyzt;
//    xyzt.x = static_cast<SCHAR>(siteIndex / _DC_MultX);
//    xyzt.y = static_cast<SCHAR>((siteIndex % _DC_MultX) / _DC_MultY);
//    xyzt.z = static_cast<SCHAR>((siteIndex % _DC_MultY) / _DC_MultZ);
//    xyzt.w = static_cast<SCHAR>((siteIndex % _DC_MultZ));
//    return xyzt;
//}

__device__ __inline__ static SSmallInt4 __deviceSiteIndexToInt4Baking(UINT siteIndex)
{
    SSmallInt4 xyzt;
    xyzt.x = static_cast<SCHAR>(siteIndex / _DC_MultX);
    xyzt.y = static_cast<SCHAR>((siteIndex % _DC_MultX) / _DC_MultY);
    xyzt.z = static_cast<SCHAR>((siteIndex % _DC_MultY) / _DC_MultZ);
    xyzt.w = static_cast<SCHAR>((siteIndex % _DC_MultZ));
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
    xyzt.x = static_cast<SCHAR>(siteIndex / _HC_MultX);
    xyzt.y = static_cast<SCHAR>((siteIndex % _HC_MultX) / _HC_MultY);
    xyzt.z = static_cast<SCHAR>((siteIndex % _HC_MultY) / _HC_MultZ);
    xyzt.w = static_cast<SCHAR>((siteIndex % _HC_MultZ));
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
    _kOutside           = 0x08,
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
        __device__ __inline__ Real NeedToOppositeCoeff() const { return static_cast<Real>(1 - 2 * static_cast<INT>(NeedToOpposite())); }
        __device__ __inline__ UBOOL IsDirichlet() const { return 0 != (_kDirichlet & m_byTag); }
        __device__ __inline__ UBOOL IsOutside() const { return 0 != (_kOutside & m_byTag); }

        /**
        * Improve-1 (multi-GPU-improve1.md 3.7, I8): the unique named invalid
        * sentinel. An invalid SIndex has NO legal global coordinate and must
        * never be dereferenced; never rely on "index > local volume" to detect
        * invalid entries.
        */
        static const UINT _kInvalidSiteIndex = 0xFFFFFFFFU;

        __device__ __inline__ UBOOL IsInvalid() const { return _kInvalidSiteIndex == m_uiSiteIndex; }

        /**
        * Improve-1 (3.7): fixed IsHalo semantics -- not invalid AND the index
        * falls in [localVolume, localVolume + haloCount). Do NOT equate _kGlue
        * with "cross-rank halo"; _kGlue also carries boundary-topology
        * semantics.
        */
        __device__ __inline__ UBOOL IsHalo(UINT uiLocalVolume, UINT uiHaloCount) const
        {
            return !IsInvalid() && m_uiSiteIndex >= uiLocalVolume && m_uiSiteIndex < uiLocalVolume + uiHaloCount;
        }

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

/**
* Multi-GPU: resolve an SIndex (as returned by the baked neighbour tables) to
* the GLOBAL coordinate of the site it refers to (Plan 1.4-R1/P4-1.1).
*
* Position-dependent physics -- rotation / acceleration / cylinder actions and
* the angular-momentum measurements -- feeds site coordinates straight into the
* physics formula (e.g. (sSite4.x - _DC_Centerx) * fOmega). After the lattice is
* split those coordinates are local while the centre is global, which would
* silently produce wrong numbers. Such kernels MUST wrap their site coordinate
* in this call.
*
* Improve-1 (3.7): implemented in CIndexData.h as a plain lookup into the baked
* 32-bit table CIndexData::m_pGlobalCoordinateTable (local slots
* [0, _DC_Volume), then halo slots), so it never narrows through SCHAR and never
* reverse-scans the face/edge/corner halo layout per call. An invalid SIndex has
* no legal global coordinate and fails an assert in debug builds instead of
* returning (0,0,0,0).
*
* On single-GPU builds every offset is 0 and the halo part is empty, so the
* baked table is the identity and behaviour is bit-identical to before.
*/
__device__ __inline__ static SInt4 _deviceSIndexToGlobalInt4(const SIndex& sIndex);

//Those are not constants, but commonly used parameters
class CLGAPI CCommonData
{
public:

//#if !_CLG_DOUBLEFLOAT
//    static DOUBLE m_fBeta;
//#else
//    static Real m_fBeta;
//#endif
    //to be removed
    //static Real m_fKai;

    //static UBOOL m_bStoreStaple;
    static UBOOL m_bStoreLastSolution;
    static UBOOL m_bStochasticGaussian;

    //Used in rotating frame. Since the fermion fields are copied,
    //it is convinient to set all parameters at one place
    
    //static SSmallInt4 m_sCenter;

    //Use for acceleration, Since the fermion fields are copied,
    //to be removed
    static Real m_fG;

    static UINT m_uiMaxThreadPerBlock;

    //No room for other paramters in solvers, so put it here, maybe move to somewhere else later
    static Real m_fShiftedMass;

    //to be removed
    //External Electric - Magnetic field
    //static Real m_fBz;
    //static Real m_fEz;
};

static inline void appBlockThreads(UINT threads, UINT& blockvar, UINT& threadvar)
{
    blockvar = threads > CCommonData::m_uiMaxThreadPerBlock ? appCeil(threads, CCommonData::m_uiMaxThreadPerBlock) : 1;
    threadvar = threads > CCommonData::m_uiMaxThreadPerBlock ? appCeil(threads, blockvar) : threads;
}

static inline void appBlockThreadsE(UINT threads, UINT elements, UINT& blockvar, UINT& threadvar)
{
    threadvar = elements * (CCommonData::m_uiMaxThreadPerBlock / elements);
    blockvar = (threads * elements + threadvar - 1) / threadvar;
}

__END_NAMESPACE

#endif //#ifndef _CCOMMONDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================