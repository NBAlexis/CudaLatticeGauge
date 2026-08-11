//=============================================================================
// FILENAME : CMeasure.cpp
// 
// DESCRIPTION:
// Some common functions
//
//
// REVISION:
//  [09/26/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
* Initial as zero
*/

template <class T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialZero_XYPlane(T* pBuffer)
{
    _Zero<T>(pBuffer[threadIdx.x + blockIdx.x * blockDim.x]);
}

/**
* Average over z and t
*/
template <class T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageOverZT_XYPlane(T* pBuffer)
{
    const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    _div(pBuffer[_ixy], _DC_Lz * _DC_Lt);
}

template <class T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSlice(T* resZ, UINT uiMax)
{
    const UINT idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < uiMax)
    {
        _Zero<T>(resZ[idx]);
    }
}

//__global__ void _CLG_LAUNCH_BOUND
//_kernelXY_To_RAverage_C(const UINT* __restrict__ pCount, CLGComplex* pValue, UINT uiMax)
//{
//    const UINT uiIdx = threadIdx.x + blockIdx.x * blockDim.x;
//    if (uiIdx < uiMax && pCount[uiIdx] > 0)
//    {
//        pValue[uiIdx].x = pValue[uiIdx].x / static_cast<Real>(pCount[uiIdx]);
//        pValue[uiIdx].y = pValue[uiIdx].y / static_cast<Real>(pCount[uiIdx]);
//    }
//}

template<class T>
__global__ void _CLG_LAUNCH_BOUND
_kernelXY_To_R(
    const T* __restrict__ jgsXY,
    UINT uiMax,
    BYTE byFieldId,
    T* result,
    UINT* pCount,
    UBOOL bShiftCenter)
{
    //Note that, generally, uiXY = (threadIdx.x + blockIdx.x * blockDim.x) does not hold anymore
    //Here, we use a special decomposition
    const UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    const INT iX = static_cast<INT>(uiXY / _DC_Ly);
    const INT iY = static_cast<INT>(uiXY % _DC_Ly);
    INT iC;
    const INT iCenterX = _DC_Centerx;
    const INT iCenterY = _DC_Centery;
    if (bShiftCenter)
    {
        iC = (((iCenterX - iX) * 2) - 1) * (((iCenterX - iX) * 2) - 1)
           + (((iCenterY - iY) * 2) - 1) * (((iCenterY - iY) * 2) - 1);
    }
    else
    {
        iC = (iCenterX - iX) * (iCenterX - iX)
           + (iCenterY - iY) * (iCenterY - iY);
    }
    //printf("Center: %d, %d iXY: %d iC: %d\n", iCenterX, iCenterY, uiXY, iC);
    SSmallInt4 sSite4;
    sSite4.z = _DC_Centerz;
    sSite4.w = _DC_Centert;
    sSite4.x = static_cast<SCHAR>(iX);
    sSite4.y = static_cast<SCHAR>(iY);
    if (iC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        if (NULL != pCount)
        {
            atomicAdd(&pCount[iC], 1);
        }
        _atomicAdd(&result[iC], jgsXY[uiXY]);
    }
}

template<class T>
__global__ void _CLG_LAUNCH_BOUND
_kernelXY_To_RAverage(const UINT* __restrict__ pCount, T* pValue, UINT uiMax)
{
    const UINT uiIdx = threadIdx.x + blockIdx.x * blockDim.x;
    if (uiIdx < uiMax && pCount[uiIdx] > 0)
    {
        _div(pValue[uiIdx], pCount[uiIdx]);
    }
}

template<class T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialDist(UINT* pCount, T* pValue, UINT uiMaxR)
{
    const UINT uiIdx = (threadIdx.x + blockIdx.x * blockDim.x);
    if (uiIdx < uiMaxR)
    {
        if (NULL != pCount)
        {
            pCount[uiIdx] = 0;
        }
        if (NULL != pValue)
        {
            _Zero<T>(pValue[uiIdx]);
        }
    }
}

#pragma endregion

void CMeasure::Average()
{
    appPushLogDate(FALSE);
    if (m_lstRealResults.Num() > 0)
    {
        appAssert(m_uiConfigurationCount == static_cast<UINT>(m_lstRealResults.Num()));
        m_fAverageRealRes = F(0.0);
        for (INT i = 0; i < m_lstRealResults.Num(); ++i)
        {
            m_fAverageRealRes += m_lstRealResults[i];
        }
        m_fAverageRealRes = m_fAverageRealRes / m_uiConfigurationCount;
        appParanoiac(_T(" === Averaged (%d measures) === %f\n"), m_uiConfigurationCount, m_fAverageRealRes);
    }

    if (m_lstComplexResults.Num() > 0)
    {
        appAssert(m_uiConfigurationCount == static_cast<UINT>(m_lstComplexResults.Num()));
        m_cAverageCmpRes = _zeroc;
        for (INT i = 0; i < m_lstComplexResults.Num(); ++i)
        {
            m_cAverageCmpRes.x += m_lstComplexResults[i].x;
            m_cAverageCmpRes.y += m_lstComplexResults[i].y;
        }
        m_cAverageCmpRes.x = m_cAverageCmpRes.x / m_uiConfigurationCount;
        m_cAverageCmpRes.y = m_cAverageCmpRes.y / m_uiConfigurationCount;
        appParanoiac(_T(" === Averaged (%d measures) === %f + %f \n"), m_uiConfigurationCount, m_cAverageCmpRes.x, m_cAverageCmpRes.y);
    }
    appPopLogDate();
}

void CMeasure::GlobalSumReal(DOUBLE& fValue) const
{
#if _CLG_MULTI_GPU
    //P4-3.1: global sum of a locally-measured value so every rank holds the
    //GLOBAL result (measurements accumulate per-rank partial sums on the local
    //sub-lattice). No-op on single-GPU / unsplit builds.
    if (NULL != appGetComm())
    {
        appGetComm()->AllreduceSum(fValue);
    }
#endif
}

void CMeasure::GlobalSumRealArray(DOUBLE* pValues, UINT uiCount) const
{
#if _CLG_MULTI_GPU
    //P4-3.5: global sum of a locally-measured DOUBLE array (e.g. Lt slices) so
    //every rank holds the GLOBAL profile. No-op on single-GPU / unsplit builds.
    if (NULL != appGetComm())
    {
        appGetComm()->AllreduceSum(pValues, uiCount);
    }
#endif
}

DOUBLE CMeasure::GlobalPlaqutteCount() const
{
#if _CLG_MULTI_GPU
    //P4-3.2: the local count is the per-rank share on a decomposed lattice;
    //derive the global count from the global lattice (Dir*(Dir-1)/2 plaquettes
    //per site) so normalized energies are rank-count independent.
    if (NULL != appGetComm())
    {
        const UINT* pG = appGetComm()->GlobalLattice();
        const DOUBLE uiVol = static_cast<DOUBLE>(pG[0]) * pG[1] * pG[2] * pG[3];
        const DOUBLE uiDir = static_cast<DOUBLE>(_HC_Dir);
        return uiVol * uiDir * (uiDir - 1.0) / 2.0;
    }
#endif
    return static_cast<DOUBLE>(_HC_PlaqutteCount);
}

void CMeasure::GlobalSumComplexArray(CLGComplex* pValues, UINT uiCount) const
{
#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && NULL != pValues && uiCount > 0)
    {
#if _CLG_DOUBLEFLOAT
        appGetComm()->AllreduceSum(reinterpret_cast<cuDoubleComplex*>(pValues), uiCount);
#else
        //Stage through doubles: cuComplex is {float x, float y} contiguous.
        TArray<DOUBLE> staging;
        staging.AddItem(static_cast<DOUBLE>(pValues[0].x));
        staging.AddItem(static_cast<DOUBLE>(pValues[0].y));
        for (UINT i = 1; i < uiCount; ++i)
        {
            staging.AddItem(static_cast<DOUBLE>(pValues[i].x));
            staging.AddItem(static_cast<DOUBLE>(pValues[i].y));
        }
        appGetComm()->AllreduceSum(staging.GetData(), uiCount * 2);
        for (UINT i = 0; i < uiCount; ++i)
        {
            pValues[i].x = static_cast<Real>(staging[i * 2]);
            pValues[i].y = static_cast<Real>(staging[i * 2 + 1]);
        }
#endif
    }
#endif
}

void CMeasure::GlobalSumComplex(CLGComplex& fValue) const
{
    //I9: scalar convenience wrapper over the array variant.
    GlobalSumComplexArray(&fValue, 1);
}

#if !_CLG_DOUBLEFLOAT
void CMeasure::GlobalSumRealArray(Real* pValues, UINT uiCount) const
{
#if _CLG_MULTI_GPU
    //I9: float-build Real array; CLGComm stages through DOUBLE internally.
    if (NULL != appGetComm() && NULL != pValues && uiCount > 0)
    {
        appGetComm()->AllreduceSum(pValues, uiCount);
    }
#endif
}
#endif

DOUBLE CMeasure::GlobalL(UINT uiDir) const
{
#if _CLG_MULTI_GPU
    if (NULL != appGetComm())
    {
        const UINT* pG = appGetComm()->GlobalLattice();
        return static_cast<DOUBLE>(pG[uiDir]);
    }
#endif
    switch (uiDir)
    {
    case 0: return static_cast<DOUBLE>(_HC_Lx);
    case 1: return static_cast<DOUBLE>(_HC_Ly);
    case 2: return static_cast<DOUBLE>(_HC_Lz);
    default: return static_cast<DOUBLE>(_HC_Lt);
    }
}

void CMeasure::WriteRealListToFile(const CCString& sFileName) const
{
    WriteRealArray(sFileName, m_lstRealResults);
}

void CMeasure::WriteCmpListToFile(const CCString& sFileName) const
{
    WriteComplexArray(sFileName, m_lstComplexResults);
}

void CMeasure::Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_pLatticeData = pLatticeData;
    m_byId = byId;

    INT iNeedGaugeSmearing = 0;
    param.FetchValueINT(_T("GaugeSmearing"), iNeedGaugeSmearing);
    m_bNeedSmearing = 0 != iNeedGaugeSmearing;

    INT iValue = 0;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFermionFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    param.FetchValueArrayBYTE(_T("GaugeFields"), m_lstGaugeFieldIds);
    param.FetchValueArrayBYTE(_T("BosonFields"), m_lstBosonFieldIds);

    if (0 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        m_lstGaugeFieldIds.AddItem(1);
    }
}

//UINT CMeasure::GetDefaultMatrixN() const
//{
//    if (m_lstGaugeFieldIds.Num() > 0)
//    {
//        const CFieldGauge* gauge = dynamic_cast<const CFieldGauge*>(appGetLattice()->GetFieldById(m_lstGaugeFieldIds[0]));
//        if (NULL != gauge)
//        {
//            return gauge->MatrixN();
//        }
//    }
//
//    return appGetLattice()->GetDefaultSUN();
//}

void CMeasure::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple)
{
    if (1 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pAcceptGauge, m_lstGaugeFieldIds[0]);
        OnConfigurationAcceptedSingleField(pAcceptGauge[idx], (NULL == pCorrespondingStaple) ? NULL : pCorrespondingStaple[idx]);
        return;
    }
    appCrucial(_T("OnConfigurationAccepted not implemented!\n"));
}

void CMeasure::SourceSanning(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site)
{
    if (1 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pAcceptGauge, m_lstGaugeFieldIds[0]);
        SourceSanningSingleField(pAcceptGauge[idx], (NULL == pCorrespondingStaple) ? NULL : pCorrespondingStaple[idx], sources, site);
        return;
    }
    appCrucial(_T("SourceSanning not implemented!\n"));
}

void CMeasure::OnConfigurationAcceptedZ4(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd)
{
    if (1 == m_lstGaugeFieldIds.Num() && 0 == m_lstBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pAcceptGauge, m_lstGaugeFieldIds[0]);
        OnConfigurationAcceptedZ4SingleField(pAcceptGauge[idx], (NULL == pCorrespondingStaple) ? NULL : pCorrespondingStaple[idx], pZ4, pInverseZ4, bStart, bEnd);
        return;
    }
    appCrucial(_T("OnConfigurationAcceptedZ4 not implemented!\n"));
}


/**
* array[x, y] = 0
*/
template<class T>
void CMeasure::_ZeroXYPlane(T* pDeviceRes)
{
    const dim3 block(_HC_DecompX, 1, 1);
    const dim3 threads(_HC_DecompLx, 1, 1);
    _LAUNCH_KERNEL(_kernelInitialZero_XYPlane<T>, block, threads, pDeviceRes);
}


#define _IMPLEMENT_CMEASURE_ZeroXYPlane(t) \
template void CMeasure::_ZeroXYPlane<t>(t* pDeviceRes);


template<class T>
void CMeasure::_ZeroSlice(T* pDeviceRes, BYTE byDir)
{
    if (0 == byDir)
    {
        UINT uib;
        UINT uit;
        appBlockThreads(_HC_Lx, uib, uit);
        _LAUNCH_KERNEL(_kernelInitialSlice<T>, uib, uit, pDeviceRes, _HC_Lx);
    }
    else if (1 == byDir)
    {
        UINT uib;
        UINT uit;
        appBlockThreads(_HC_Ly, uib, uit);
        _LAUNCH_KERNEL(_kernelInitialSlice<T>, uib, uit, pDeviceRes, _HC_Ly);
    }
    else if (2 == byDir)
    {
        dim3 blockz(_HC_DecompY, 1, 1);
        dim3 threadz(_HC_DecompLy, 1, 1);
        _LAUNCH_KERNEL(_kernelInitialSlice<T>, blockz, threadz, pDeviceRes, _HC_Lz);
    }
    else
    {
        UINT uib;
        UINT uit;
        appBlockThreads(_HC_Lt, uib, uit);
        _LAUNCH_KERNEL(_kernelInitialSlice<T>, uib, uit, pDeviceRes, _HC_Lt);
    }
}

#define _IMPLEMENT_CMEASURE_ZeroSlice(t) \
template void CMeasure::_ZeroSlice<t>(t* pDeviceRes, BYTE byDir);

/**
* array[x, y] = array[x, y] / (lz * lt)
*/
template<class T>
void CMeasure::_AverageXYPlane(T* pDeviceRes)
{
    const dim3 block(_HC_DecompX, 1, 1);
    const dim3 threads(_HC_DecompLx, 1, 1);
    _LAUNCH_KERNEL(_kernelInitialZero_XYPlane<T>, block, threads, pDeviceRes);
}

#define _IMPLEMENT_CMEASURE_AverageXYPlane(t) \
template void CMeasure::_AverageXYPlane<t>(t* pDeviceRes);

template<class T>
void CMeasure::XYDataToRdistri(
    UBOOL bShiftCenter,
    const T* source,
    UINT* count,
    T* result,
    UINT uiMaxR,
    UBOOL bCalculateCounter,
    BYTE byFieldId)
{
    const dim3 block2(_HC_DecompX, 1, 1);
    const dim3 threads2(_HC_DecompLx, 1, 1);

    __SIMPLEDECOMPOSE(uiMaxR + 1);
    //const dim3 block3(1, 1, 1);
    //const dim3 threads3(uiMaxR + 1, 1, 1);

    _LAUNCH_KERNEL(_kernelInitialDist<T>, block, thread, 
        bCalculateCounter ? count : NULL,
        result,
        uiMaxR + 1
        );

    _LAUNCH_KERNEL(_kernelXY_To_R<T>, block2, threads2,
        source,
        uiMaxR,
        byFieldId,
        result,
        bCalculateCounter ? count : NULL,
        bShiftCenter
        );

    _LAUNCH_KERNEL(_kernelXY_To_RAverage<T>, block, thread, count, result, uiMaxR + 1);
}

#define _IMPLEMENT_CMEASURE_XYDataToRdistri(t) \
template void CMeasure::XYDataToRdistri<t>(UBOOL bShiftCenter, const t* source, UINT* count, t* result, UINT uiMaxR, UBOOL bCalculateCounter, BYTE byFieldId);


_IMPLEMENT_CMEASURE_ZeroXYPlane(Real)
_IMPLEMENT_CMEASURE_ZeroSlice(Real)
_IMPLEMENT_CMEASURE_AverageXYPlane(Real)
_IMPLEMENT_CMEASURE_XYDataToRdistri(Real)
_IMPLEMENT_CMEASURE_ZeroXYPlane(CLGComplex)
_IMPLEMENT_CMEASURE_ZeroSlice(CLGComplex)
_IMPLEMENT_CMEASURE_AverageXYPlane(CLGComplex)
_IMPLEMENT_CMEASURE_XYDataToRdistri(CLGComplex)
#if _CLG_DOUBLEFLOAT
_IMPLEMENT_CMEASURE_ZeroXYPlane(FLOAT)
_IMPLEMENT_CMEASURE_ZeroSlice(FLOAT)
_IMPLEMENT_CMEASURE_AverageXYPlane(FLOAT)
_IMPLEMENT_CMEASURE_XYDataToRdistri(FLOAT)
_IMPLEMENT_CMEASURE_ZeroXYPlane(cuComplex)
_IMPLEMENT_CMEASURE_ZeroSlice(cuComplex)
_IMPLEMENT_CMEASURE_AverageXYPlane(cuComplex)
_IMPLEMENT_CMEASURE_XYDataToRdistri(cuComplex)
#else
_IMPLEMENT_CMEASURE_ZeroXYPlane(DOUBLE)
_IMPLEMENT_CMEASURE_ZeroSlice(DOUBLE)
_IMPLEMENT_CMEASURE_AverageXYPlane(DOUBLE)
_IMPLEMENT_CMEASURE_XYDataToRdistri(DOUBLE)
_IMPLEMENT_CMEASURE_ZeroXYPlane(cuDoubleComplex)
_IMPLEMENT_CMEASURE_ZeroSlice(cuDoubleComplex)
_IMPLEMENT_CMEASURE_AverageXYPlane(cuDoubleComplex)
_IMPLEMENT_CMEASURE_XYDataToRdistri(cuDoubleComplex)
#endif

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================