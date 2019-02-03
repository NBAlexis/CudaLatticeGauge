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
//    //Real m_fBeta;
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

    void CreateFermionSolver(const CCString& sSolver, const CParameters& param, const class CField* pFermionField);
    void OnUpdatorConfigurationAccepted();
    void OnUpdatorFinished();
    void GetPlaquetteLengthCount(UINT& plaqLength, UINT& countPerSite, UINT& countPerLink);

    class CRandom* m_pRandom;

    class CFieldGauge* m_pGaugeField;
    class CFieldGauge* m_pGaugeFieldStaple;
    THashMap<BYTE, class CField*> m_pFieldMap;

    /**
    * \note we assume the first action is gauge action
    */
    TArray<class CAction*> m_pActionList;
    THashMap<BYTE, class CAction*> m_pActionMap;

    class CUpdator* m_pUpdator;
    class CMeasurementManager* m_pMeasurements;

    //this is a device copy
    //see:
    //https://stackoverflow.com/questions/53781421/cuda-the-member-field-with-device-ptr-and-device-member-function-to-visit-it-i
    class CRandom* m_pDeviceRandom;
    class CIndex* m_pDeviceIndex;

    class CFieldGauge* m_pDeviceGaugeField;
    class CFieldGauge* m_pDeviceGaugeFieldStaple;

    class CSLASolver* m_pFermionSolver;
    

    class CField* GetFieldById(BYTE byId) const { return m_pFieldMap.GetAt(byId); }
    class CAction* GetActionById(BYTE byId) const { return m_pActionMap.GetAt(byId); }
};

inline class CSLASolver* appGetFermionSolver();

__END_NAMESPACE

#endif //#ifndef _CLATTICEDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================