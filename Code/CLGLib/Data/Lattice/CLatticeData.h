//=============================================================================
// FILENAME : CLatticeData.h
// 
// DESCRIPTION:
// This is the class for the lattce data
// NOTE:: We only have 4D case, 3D = 1xLxLxL, and 2D= 1x1xLxL
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CLATTICEDATA_H_
#define _CLATTICEDATA_H_

__BEGIN_NAMESPACE

struct CLGAPI CDeviceLattice
{
    //One can only create a device lattice on device
    __device__ CDeviceLattice() {}

    __device__ __inline__ CRandomSchrage* GetDeviceRandom() const { return m_pDeviceRandom; }

    UINT m_uiVolumn;
    UINT m_uiDim;
    UINT m_uiDir;
    UINT m_uiTLength;
    FLOAT m_fBeta;
    UINT m_uiLatticeLength[CCommonData::kMaxDim];
    UINT m_uiLatticeDecompose[CCommonData::kLatticeDecompose * 2];
    UINT m_uiLatticeMultipy[CCommonData::kMaxDim - 1];

    class CIndex* m_pIndex;
    CRandomSchrage* m_pDeviceRandom;
};

class CLGAPI CLatticeData
{
    static CLatticeData* m_pInstance;
    CDeviceLattice* m_pDeviceInstance;

public:
    static void Create() { if (NULL == m_pInstance) { m_pInstance = new CLatticeData(); } }

    static __host__ CLatticeData* GetInstance() { Create(); return m_pInstance; }
    __device__ __inline__ CDeviceLattice* GetDeviceInstance() const  { return m_pDeviceInstance; }

    /**
    * One thread deal with only data[x / m_pLatticeDecompose[0], y / m_pLatticeDecompose[1], z / m_pLatticeDecompose[2]]
    * with m_pLatticeDecompose[3] * m_pLatticeDecompose[4] * m_pLatticeDecompose[5] blocks
    * For 1D, m_pLatticeDecompose[1] = m_pLatticeDecompose[2] = 1
    * For 2D, m_pLatticeDecompose[2] = 1
    */
    UINT m_uiVolumn;
    UINT m_uiDim;
    UINT m_uiDir;
    UINT m_uiLatticeLength[CCommonData::kMaxDim];
    UINT m_uiTLength; //this is special because T dir is not decomposed to thread blocks
    UINT m_uiLatticeDecompose[CCommonData::kLatticeDecompose * 2];
    FLOAT m_fBeta;

    /*
    * SU3(x=(x,y,z,t))_{n=a*3+b}=
    * m_pData[(
        (x*m_uiLatticeLength[1]*m_uiLatticeLength[2]*m_uiLatticeLength[3] 
       + y*m_uiLatticeLength[2]*m_uiLatticeLength[3] 
       + z*m_uiLatticeLength[3] 
       + t)
       * m_uiDir + dir) * elementCount + n]
    *
    * m_uiLatticeMultipy[0] = m_uiLatticeLength[1]*m_uiLatticeLength[2]*m_uiLatticeLength[3]
    * m_uiLatticeMultipy[1] = m_uiLatticeLength[2]*m_uiLatticeLength[3]
    * m_uiLatticeMultipy[2] = m_uiLatticeLength[3]
    * for field on set, dir = 1
    */
    UINT m_uiLatticeMultipy[CCommonData::kMaxDim - 1];

    STRING m_sFields[CCommonData::kMaxFieldCount];
    class CField* m_pFields[CCommonData::kMaxFieldCount];

    class CFieldGauge* m_pGaugeField;
    class CFieldGauge* m_pGaugeFieldStaple;
    class CFieldGauge* m_pGaugeFieldMomentum;
    class CFieldGauge* m_pGaugeFieldForce;

private:

    /**
    * Initial with CommonData
    */
    CLatticeData();
    ~CLatticeData();
};

#pragma region Index Calculation

__device__ __inline__ 
static UINT _deviceGetSiteIndex(const CDeviceLattice * pLattice, const UINT* coord)
{
    return coord[0] * pLattice->m_uiLatticeMultipy[0] + coord[1] * pLattice->m_uiLatticeMultipy[1] + coord[2] * pLattice->m_uiLatticeMultipy[2];
}

__device__ __inline__ 
static UINT _deviceGetLinkIndex(const class CDeviceLattice * pLattice, UINT siteIndex, UINT dir)
{
    return siteIndex * pLattice->m_uiDir + dir;
}

__device__ __inline__ 
static UINT _deviceGetLinkIndex(const CDeviceLattice * pLattice, const UINT* coord, UINT dir)
{
    return _deviceGetSiteIndex(pLattice, coord) * pLattice->m_uiDir + dir;
}

__device__ __inline__ 
static UINT _deviceGetFatIndex(const CDeviceLattice * pLattice, const UINT* coord, UINT dir_plus_one)
{
    return _deviceGetSiteIndex(pLattice, coord) * (pLattice->m_uiDir + 1) + dir_plus_one;
}

__device__ __inline__ 
static int4 __deviceSiteIndexToInt4(const CDeviceLattice * pLattice, UINT siteIndex)
{
    int4 xyzt;
    xyzt.x = siteIndex / pLattice->m_uiLatticeMultipy[0];
    xyzt.y = (siteIndex % pLattice->m_uiLatticeMultipy[0]) / pLattice->m_uiLatticeMultipy[1];
    xyzt.z = (siteIndex % pLattice->m_uiLatticeMultipy[1]) / pLattice->m_uiLatticeMultipy[2];
    xyzt.w = (siteIndex % pLattice->m_uiLatticeMultipy[2]);
    return xyzt;
}

__device__ __inline__ 
static int4 __deviceLinkIndexToInt4(const CDeviceLattice * pLattice, UINT linkIndex)
{
    return __deviceSiteIndexToInt4(pLattice, linkIndex / pLattice->m_uiDir);
}

__device__ __inline__ 
static int4 __deviceFatIndexToInt4(const CDeviceLattice * pLattice, UINT fatIndex)
{
    return __deviceSiteIndexToInt4(pLattice, fatIndex / (pLattice->m_uiDir + 1));
}

#pragma endregion Index Calculation

__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================