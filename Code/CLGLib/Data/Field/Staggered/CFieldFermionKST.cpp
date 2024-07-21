//=============================================================================
// FILENAME : CFieldFermionKST.cpp
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [12/08/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldFermionKST.h"

__BEGIN_NAMESPACE



#pragma region Kernel

//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelAxpyPlusFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther)
//{
//    intokernal;
//    _add(pMe[uiSiteIndex], pOther[uiSiteIndex]);
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelAxpyMinusFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther)
//{
//    intokernal;
//    _sub(pMe[uiSiteIndex], pOther[uiSiteIndex]);
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelAxpyComplexFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther, CLGComplex a)
//{
//    intokernal;
//    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelMulFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther, UBOOL bConj)
//{
//    intokernal;
//    if (bConj)
//    {
//        _dagger(pMe[uiSiteIndex]);
//    }
//    _mul(pMe[uiSiteIndex], pOther[uiSiteIndex]);
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelAxpyRealFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther, Real a)
//{
//    intokernal;
//    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelDotFermionKST(const deviceVector* __restrict__ pMe, const deviceVector* __restrict__ pOther, cuDoubleComplex* result
//)
//{
//    intokernal;
//    result[uiSiteIndex] = _cToDouble(_dot(pMe[uiSiteIndex], pOther[uiSiteIndex]));
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelScalarMultiplyComplexKST(deviceVector* pMe, CLGComplex a)
//{
//    intokernal;
//    _mul(pMe[uiSiteIndex], a);
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelScalarMultiplyRealKST(deviceVector* pMe, Real a)
//{
//    intokernal;
//    _mul(pMe[uiSiteIndex], a);
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelFermionKSConjugateT(deviceVector* pDeviceData)
//{
//    intokernal;
//    _dagger(pDeviceData[uiSiteIndex]);
//}
//
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelMakePointSourceKST(deviceVector* pDeviceData, UINT uiDesiredSite, BYTE byColor)
//{
//    intokernal;
//    if (uiSiteIndex == uiDesiredSite)
//    {
//        pDeviceData[uiSiteIndex] = _makeColorVector<deviceVector>(byColor);
//    }
//    else
//    {
//        pDeviceData[uiSiteIndex] = _makeZero<deviceVector>();
//    }
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelMakePointSourceKSOneT(deviceVector* pDeviceData, UINT uiDesiredSite)
//{
//    intokernal;
//    if (uiSiteIndex == uiDesiredSite)
//    {
//        pDeviceData[uiSiteIndex] = _makeId<deviceVector>();
//    }
//    else
//    {
//        pDeviceData[uiSiteIndex] = _makeZero<deviceVector>();
//    }
//}
//
//template<typename deviceVector>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelMakeWallSourceKST(deviceVector* pDeviceData, 
//    INT uiDesiredT, UINT uiShift, BYTE color, BYTE byFieldID)
//{
//    intokernalOnlyInt4;
//
//    //pDeviceData[uiSiteIndex] = _makeZero<deviceVector>_makeZero();
//    //We shall not set zero here! [2024/7/21/ why not: because offset site will be set to non-zero!!]
//
//    if ( (0 == (sSite4.x & 1))
//      && (0 == (sSite4.y & 1))
//      && (0 == (sSite4.z & 1))
//      && (uiDesiredT < 0 || uiDesiredT == sSite4.w))
//    {
//        //sSite4 is no longer used
//        sSite4.x = sSite4.x + static_cast<SBYTE>(uiShift & 1);
//        sSite4.y = sSite4.y + static_cast<SBYTE>((uiShift >> 1) & 1);
//        sSite4.z = sSite4.z + static_cast<SBYTE>((uiShift >> 2) & 1);
//        const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][__bi(sSite4)];
//        if (!sIdx.IsDirichlet())
//        {
//            pDeviceData[sIdx.m_uiSiteIndex] = _makeColorVector<deviceVector>(color);
//        }
//    }
//}

#pragma endregion




//template<typename deviceVector, typename deviceGauge, INT vectorN>
//CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CFieldFermionKST()
//    : CFieldFermionKS()
//    , m_pRationalFieldPointers(NULL)
//{
//    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount));
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//CFieldFermionKST<deviceVector, deviceGauge, vectorN>::~CFieldFermionKST()
//{
//    checkCudaErrors(__cudaFree(m_pDeviceData));
//    if (NULL != m_pRationalFieldPointers)
//    {
//        checkCudaErrors(cudaFree(m_pRationalFieldPointers));
//    }
//    if (NULL != m_pMDNumerator)
//    {
//        checkCudaErrors(cudaFree(m_pMDNumerator));
//    }
//}

/**
*
*/
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialField(EFieldInitialType eInitialType)
//{
//    preparethread;
//    _kernelInitialFermionKST << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Dagger()
//{
//    preparethread;
//    _kernelFermionKSConjugateT << <block, threads >> > (m_pDeviceData);
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
//{
//    if (eFieldType != EFFT_CLGBin)
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialFieldWithFile: Only support CLG Bin File\n"));
//        return;
//    }
//
//    UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * vectorN * m_uiSiteCount);
//    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
//    if (NULL == data)
//    {
//        appCrucial(_T("File not found: %s\n"), sFileName.c_str());
//        _FAIL_EXIT;
//    }
//    if (uiSize != static_cast<UINT>(sizeof(Real) * 2 * vectorN * m_uiSiteCount))
//    {
//        appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), static_cast<UINT>(sizeof(Real) * 2 * vectorN * m_uiSiteCount), uiSize);
//        _FAIL_EXIT;
//    }
//    InitialWithByte(data);
//    free(data);
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialWithByte(BYTE* byData)
//{
//    deviceVector* readData = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
//    for (UINT i = 0; i < m_uiSiteCount; ++i)
//    {
//        Real thisSite[2 * vectorN];
//        memcpy(thisSite, byData + i * sizeof(Real) * 2 * vectorN, sizeof(Real) * 2 * vectorN);
//        for (UINT k = 0; k < 2 * vectorN; ++k)
//        {
//            _setelement(readData[i], k, thisSite[k]);
//        }
//    }
//    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyHostToDevice));
//    free(readData);
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialOtherParameters(CParameters& params)
//{
//    CFieldFermionKS::InitialOtherParameters(params);
//
//    if (NULL != m_pRationalFieldPointers)
//    {
//        checkCudaErrors(cudaFree(m_pRationalFieldPointers));
//    }
//    checkCudaErrors(cudaMalloc((void**)&m_pRationalFieldPointers, sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree));
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DebugPrintMe() const
//{
//    deviceVector* toprint = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
//    checkCudaErrors(cudaMemcpy(toprint, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
//    for (UINT uiSite = 0; uiSite < m_uiSiteCount; ++uiSite)
//    {
//        const SSmallInt4 site4 = __hostSiteIndexToInt4(uiSite);
//        appGeneral(_T(" --- %d,%d,%d,%d --- %s\n"),
//            site4.x, site4.y, site4.z, site4.w, appToString(toprint[uiSite]).c_str()
//        );
//    }
//
//    appSafeFree(toprint);
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyTo(CField* U) const
//{
//    if (NULL == U || GetFieldType() != U->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return;
//    }
//
//    CFieldFermionKS::CopyTo(U);
//
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(U);
//    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//
//    if (NULL != pField->m_pMDNumerator)
//    {
//        checkCudaErrors(cudaFree(pField->m_pMDNumerator));
//    }
//    if (NULL != pField->m_pRationalFieldPointers)
//    {
//        checkCudaErrors(cudaFree(pField->m_pRationalFieldPointers));
//    }
//
//    checkCudaErrors(cudaMalloc((void**)&pField->m_pRationalFieldPointers, sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree));
//    checkCudaErrors(cudaMalloc((void**)&pField->m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree));
//    checkCudaErrors(cudaMemcpy(pField->m_pMDNumerator, m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree, cudaMemcpyDeviceToDevice));
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::AxpyPlus(const CField* x)
//{
//    if (NULL == x || GetFieldType() != x->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return;
//    }
//    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
//
//    preparethread;
//    _kernelAxpyPlusFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::AxpyMinus(const CField* x)
//{
//    if (NULL == x || GetFieldType() != x->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return;
//    }
//    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
//
//    preparethread;
//    _kernelAxpyMinusFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Axpy(Real a, const CField* x)
//{
//    if (NULL == x || GetFieldType() != x->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return;
//    }
//    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
//
//    preparethread;
//    _kernelAxpyRealFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Axpy(const CLGComplex& a, const CField* x)
//{
//    if (NULL == x || GetFieldType() != x->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return;
//    }
//    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
//
//    preparethread;
//    _kernelAxpyComplexFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Mul(const CField* other, UBOOL bDagger)
//{
//    if (NULL == other || GetFieldType() != other->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return;
//    }
//    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(other);
//
//    preparethread;
//    _kernelMulFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, bDagger);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//cuDoubleComplex CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Dot(const CField* x) const
//{
//    if (NULL == x || GetFieldType() != x->GetFieldType())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
//        return make_cuDoubleComplex(0, 0);
//    }
//    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
//    preparethread;
//    _kernelDotFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);
//
//    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ScalarMultply(const CLGComplex& a)
//{
//    preparethread;
//    _kernelScalarMultiplyComplexKST << <block, threads >> > (m_pDeviceData, a);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ScalarMultply(Real a)
//{
//    preparethread;
//    _kernelScalarMultiplyRealKST << <block, threads >> > (m_pDeviceData, a);
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ApplyGamma(EGammaMatrix eGamma)
//{
//    appCrucial(_T("Not implemented yet...\n"));
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
//{
//    if (NULL == pGauge || vectorN != pGauge->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//    InitialField(EFIT_Zero);
//
//    //If it was gamma_mu or gamma_5 or sigmaij, it is i gamma mu and i gamma 5, therefore multiply -i
//    UBOOL bImag = (GAMMA1 == eGamma) || (GAMMA2 == eGamma) || (GAMMA3 == eGamma) || (GAMMA4 == eGamma) || (GAMMA5 == eGamma)
//        || (SIGMA12 == eGamma) || (SIGMA31 == eGamma) || (SIGMA41 == eGamma) || (SIGMA23 == eGamma) || (SIGMA42 == eGamma) || (SIGMA43 == eGamma);
//
//    appApplyGammaKS(
//        m_pDeviceData,
//        pPooled->m_pDeviceData,
//        pFieldSU3->m_pDeviceData,
//        eGamma,
//        m_bEachSiteEta,
//        FALSE,
//        F(0.5),
//        bImag ? EOCT_Complex : EOCT_None,
//        F(1.0),
//        bImag ? _make_cuComplex(F(0.0), -F(1.0)) : _onec,
//        m_byFieldId,
//        pGauge->m_byFieldId
//    );
//
//    pPooled->Return();
//}


//Kai should be part of D operator
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//
//    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
//        FALSE, eCoeffType, fRealCoeff, cCompCoeff);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    if (m_bDiagonalMass)
//    {
//        appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
//    }
//
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//
//    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
//        FALSE, eCoeffType, fRealCoeff, cCompCoeff);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::D0S(const CField* pGauge)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    if (m_bDiagonalMass)
//    {
//        appCrucial(_T("In the cass mass is not a number, should not in here except for check anti-hermiticity!\n"));
//    }
//
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//
//    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, F(0.0),
//        FALSE, EOCT_None, F(1.0), _onec);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//UINT CFieldFermionKST<deviceVector, deviceGauge, vectorN>::TestAntiHermitianS(const CFieldGauge* pGauge) const
//{
//    const UINT uiVolume = _HC_Volume;
//    const UINT uiRealVolume = vectorN * uiVolume;
//    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
//    deviceVector* hostData = (deviceVector*)malloc(sizeof(deviceVector) * uiVolume);
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* v = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    for (UINT i = 0; i < uiVolume; ++i)
//    {
//        const SSmallInt4 point = __hostSiteIndexToInt4(i);
//        for (UINT j = 0; j < vectorN; ++j)
//        {
//            SFermionBosonSource source;
//            source.m_byColorIndex = static_cast<BYTE>(j);
//            source.m_eSourceType = EFS_Point;
//            source.m_sSourcePoint = point;
//            v->InitialAsSource(source);
//            v->D0S(pGauge);
//
//            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceVector) * uiVolume, cudaMemcpyDeviceToHost));
//
//            const UINT x = i * vectorN + j;
//            for (UINT k = 0; k < uiVolume; ++k)
//            {
//                for (UINT kcolor = 0; kcolor < vectorN; ++kcolor)
//                {
//                    matrixElement[(vectorN * k + kcolor) * uiRealVolume + x] = _make_cuComplex(_element(hostData[k], 2 * kcolor), _element(hostData[k], 2 * kcolor + 1));
//                }
//            }
//            appGeneral(_T("%d / %d have been done\n"), x, uiRealVolume);
//        }
//    }
//
//    UINT uiE = 0;
//    UINT uiWrong = 0;
//    //List all results
//    for (UINT i = 0; i < uiRealVolume * uiRealVolume; ++i)
//    {
//        const UINT x = i / uiRealVolume;
//        const UINT y = i % uiRealVolume;
//        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / vectorN);
//        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / vectorN);
//        const UINT daggerIdx = y * uiRealVolume + x;
//        const BYTE cx = x % vectorN;
//        const BYTE cy = y % vectorN;
//
//        if (_cuCabsf(matrixElement[i]) > F(0.0000001))
//        {
//            ++uiE;
//            if (appAbs(matrixElement[i].x + matrixElement[daggerIdx].x) > F(0.0000001)
//             || appAbs(matrixElement[i].y - matrixElement[daggerIdx].y) > F(0.0000001))
//            {
//                ++uiWrong;
//                appGeneral(_T("[(%d, %d, %d, %d)_(%d)-(%d, %d, %d, %d)_(%d)]: D = %f + %f I   Ddagger = %f + %f I\n"),
//                    xSite.x, xSite.y, xSite.z, xSite.w, cx, 
//                    ySite.x, ySite.y, ySite.z, ySite.w, cy, 
//                    matrixElement[i].x, matrixElement[i].y,
//                    matrixElement[daggerIdx].x, matrixElement[daggerIdx].y);
//            }
//        }
//    }
//    v->Return();
//    appSafeFree(matrixElement);
//    appSafeFree(hostData);
//    appGeneral(_T("%d none zero element checked, %d wrong found...\n"), uiE, uiWrong);
//    return uiWrong;
//}

//Kai should be part of D operator
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//
//    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
//        TRUE, eCoeffType, fRealCoeff, cCompCoeff);
//
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    if (m_bDiagonalMass)
//    {
//        appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
//    }
//
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//
//    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
//        TRUE, eCoeffType, fRealCoeff, cCompCoeff);
//
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
//        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
//    //why only apply coeff in the next step?
//    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
//        FALSE, eCoeffType, fRealCoeff, cCompCoeff);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
//        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
//    //why only apply coeff in the next step?
//    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
//        FALSE, eCoeffType, fRealCoeff, cCompCoeff);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//    if (m_bDiagonalMass)
//    {
//        appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
//    }
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
//        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
//    //why only apply coeff in the next step?
//    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
//        FALSE, eCoeffType, fRealCoeff, cCompCoeff);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
//{
//    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
//    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
//    {
//        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
//        return;
//    }
//
//    Real fRealCoeff = fCoeffReal;
//    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
//    if (EOCT_Minus == eCoeffType)
//    {
//        eCoeffType = EOCT_Real;
//        fRealCoeff = F(-1.0);
//    }
//    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
//
//    DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
//        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
//    //why only apply coeff in the next step?
//    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
//        FALSE, eCoeffType, fRealCoeff, cCompCoeff);
//
//    pPooled->Return();
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialAsSource(const SFermionBosonSource& sourceData)
//{
//    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
//    switch (sourceData.m_eSourceType)
//    {
//    case EFS_Point:
//    {
//        preparethread;
//        if (sourceData.m_byColorIndex < vectorN)
//        {
//            _kernelMakePointSourceKST << <block, threads >> > (m_pDeviceData, uiSiteIndex, sourceData.m_byColorIndex);
//        }
//        else
//        {
//            _kernelMakePointSourceKSOneT << <block, threads >> > (m_pDeviceData, uiSiteIndex);
//        }
//    }
//    break;
//    case EFS_Wall:
//    {
//        preparethread;
//        _kernelInitialFermionKST << <block, threads >> > (m_pDeviceData, m_byFieldId,EFIT_Zero);
//        _kernelMakeWallSourceKST << <block, threads >> > (
//            m_pDeviceData,
//            static_cast<INT>(sourceData.m_sSourcePoint.w),
//            sourceData.m_bySpinIndex,
//            sourceData.m_byColorIndex,
//            m_byFieldId);
//    }
//    break;
//    default:
//        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionBosonSource, sourceData.m_eSourceType).c_str());
//        break;
//    }
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//BYTE* CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyDataOut(UINT& uiSize) const
//{
//    deviceVector* toSave = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
//    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 2 * vectorN);
//    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
//    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
//    for (UINT i = 0; i < m_uiSiteCount; ++i)
//    {
//        Real oneSite[2 * vectorN];
//        for (UINT k = 0; k < 2 * vectorN; ++k)
//        {
//            oneSite[k] = _element(toSave[i], k);
//        }
//        memcpy(saveData + sizeof(Real) * i * 2 * vectorN, oneSite, sizeof(Real) * 2 * vectorN);
//    }
//
//    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
//    //free(saveData);
//    free(toSave);
//    return saveData;
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//BYTE* CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyDataOutFloat(UINT& uiSize) const
//{
//    deviceVector* toSave = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
//    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * 2 * vectorN);
//    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
//    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
//    for (UINT i = 0; i < m_uiSiteCount; ++i)
//    {
//        FLOAT oneSite[2 * vectorN];
//        for (UINT k = 0; k < 2 * vectorN; ++k)
//        {
//            oneSite[k] = static_cast<FLOAT>(_element(toSave[i], k));
//        }
//        memcpy(saveData + sizeof(FLOAT) * i * 2 * vectorN, oneSite, sizeof(FLOAT) * 2 * vectorN);
//    }
//
//    free(toSave);
//    return saveData;
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//BYTE* CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyDataOutDouble(UINT& uiSize) const
//{
//    deviceVector* toSave = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
//    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiSiteCount * 2 * vectorN);
//    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
//    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
//    for (UINT i = 0; i < m_uiSiteCount; ++i)
//    {
//        DOUBLE oneSite[2 * vectorN];
//        for (UINT k = 0; k < 2 * vectorN; ++k)
//        {
//            oneSite[k] = static_cast<DOUBLE>(_element(toSave[i], k));
//        }
//        memcpy(saveData + sizeof(DOUBLE) * i * 2 * vectorN, oneSite, sizeof(DOUBLE) * 2 * vectorN);
//    }
//
//    free(toSave);
//    return saveData;
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//TArray<CFieldFermion*> CFieldFermionKST<deviceVector, deviceGauge, vectorN>::GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const
//{
//    TArray<CFieldFermion*> ret;
//    for (UINT j = 0; j < vectorN; ++j)
//    {
//        ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));
//        if (NULL == ret[j])
//        {
//            appCrucial(_T("GetSourcesAtSiteFromPool failed!\n"));
//            _FAIL_EXIT;
//        }
//    }
//
//    for (BYTE c = 0; c < vectorN; ++c)
//    {
//        SFermionBosonSource sourceData;
//        sourceData.m_eSourceType = EFS_Point;
//        sourceData.m_sSourcePoint = site;
//        sourceData.m_byColorIndex = c;
//        sourceData.m_bySpinIndex = 0;
//
//        ret[c]->InitialAsSource(sourceData);
//
//        if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
//        {
//            ret[c]->m_fLength = ret[c]->Dot(ret[c]).x;
//        }
//
//        ret[c]->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
//    }
//    return ret;
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::PrepareForHMCOnlyRandomize()
//{
//    preparethread;
//    _kernelInitialFermionKST << <block, threads >> > (
//        m_pDeviceData,
//        m_byFieldId,
//        EFIT_RandomGaussian);
//}
//
//template<typename deviceVector, typename deviceGauge, INT vectorN>
//void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
//{
//    D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
//}


#pragma region Field Matrix Operation

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//CFieldMatrixOperationKST<deviceVector, deviceGauge, vectorN>::CFieldMatrixOperationKST()
//{
//    m_pHostResBuffer = (deviceVector**)malloc(sizeof(deviceVector*) * _kFieldMatrixMaxDim);
//    m_pHostLeftBuffer = (deviceVector**)malloc(sizeof(deviceVector*) * _kFieldMatrixMaxDim);
//    checkCudaErrors(cudaMalloc((void**)&m_pResBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim));
//    checkCudaErrors(cudaMalloc((void**)&m_pLeftBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim));
//}

//template<typename deviceVector, typename deviceGauge, INT vectorN>
//CFieldMatrixOperationKST<deviceVector, deviceGauge, vectorN>::~CFieldMatrixOperationKST()
//{
//    free(m_pHostResBuffer);
//    free(m_pHostLeftBuffer);
//
//    checkCudaErrors(cudaFree(m_pResBuffer));
//    checkCudaErrors(cudaFree(m_pLeftBuffer));
//}



//__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldMatrixOperationKST, U1, CLGComplex, CLGComplex, 1)

__CLGIMPLEMENT_CLASS(CFieldFermionKSU1)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================