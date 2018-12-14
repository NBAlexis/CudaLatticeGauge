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

//class CLGAPI CDeviceLattice
//{
//public:
//    //One can only create a device lattice on device
//    CDeviceLattice() {}
//
//
//    //UINT m_uiVolumn;
//    //UINT m_uiDim;
//    //UINT m_uiDir;
//    //UINT m_uiTLength;
//    //FLOAT m_fBeta;
//    //UINT m_uiLatticeLength[CCommonData::kMaxDim];
//    //UINT m_uiLatticeDecompose[CCommonData::kLatticeDecompose * 2];
//    //UINT m_uiLatticeMultipy[CCommonData::kMaxDim - 1];
//
//    class CIndex* m_pIndex;
//};

class CLGAPI CLatticeData
{
    //static CLatticeData* m_pInstance;
    //CDeviceLattice* m_pDeviceInstance;

public:
    /**
    * Initial with CommonData
    */
    CLatticeData();
    ~CLatticeData();

    //static void Create() { if (NULL == m_pInstance) { m_pInstance = new CLatticeData(); } }
    //static void Release() { appSafeDelete(m_pInstance); }
    //~CLatticeData() { ; }

    //static __host__ CLatticeData* GetInstance() { Create(); return m_pInstance; }
    //__device__ __inline__ CDeviceLattice* GetDeviceInstance() const  { return m_pDeviceInstance; }

    /**
    * One thread deal with only data[x / m_pLatticeDecompose[0], y / m_pLatticeDecompose[1], z / m_pLatticeDecompose[2]]
    * with m_pLatticeDecompose[3] * m_pLatticeDecompose[4] * m_pLatticeDecompose[5] blocks
    * For 1D, m_pLatticeDecompose[1] = m_pLatticeDecompose[2] = 1
    * For 2D, m_pLatticeDecompose[2] = 1
    */
    //UINT m_uiVolumn;
    //UINT m_uiDim;
    //UINT m_uiDir;
    //UINT m_uiLatticeLength[CCommonData::kMaxDim];
    //UINT m_uiTLength; //this is special because T dir is not decomposed to thread blocks
    //UINT m_uiLatticeDecompose[CCommonData::kLatticeDecompose * 2];
    //FLOAT m_fBeta;

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
    //UINT m_uiLatticeMultipy[CCommonData::kMaxDim - 1];

    //CCString m_sFields[CCommonData::kMaxFieldCount];
    //class CField* m_pFields[CCommonData::kMaxFieldCount];
    class CRandom* m_pRandom;
    class CRandomSchrage* m_pRandomSchrage;
    class CIndex* m_pIndex;

    class CFieldGauge* m_pGaugeField;
    class CFieldGauge* m_pGaugeFieldStaple;

    //this is a device copy
    //see:
    //https://stackoverflow.com/questions/53781421/cuda-the-member-field-with-device-ptr-and-device-member-function-to-visit-it-i
    class CRandom* m_pDeviceRandom;
    class CRandomSchrage* m_pDeviceRandomSchrage;
    class CIndex* m_pDeviceIndex;

    class CFieldGauge* m_pDeviceGaugeField;
    class CFieldGauge* m_pDeviceGaugeFieldStaple;
};


__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================