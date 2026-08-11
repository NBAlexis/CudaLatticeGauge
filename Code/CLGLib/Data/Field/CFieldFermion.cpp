//=============================================================================
// FILENAME : CFieldFermion.cpp
// 
// DESCRIPTION:
// There are functions for common fermions
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "Staggered/CFieldFermionKST.h"
#include "Measurement/CMeasureAngularMomentumKS.h"

__BEGIN_NAMESPACE

void ExportDiagnalWilsonSU3(const CCString& sFileName, EMeasureDiagnal eType, 
    INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields,
    const CFieldFermionWilsonSquareSU3* pFermion)
{
    UBOOL bOnlyRed = TRUE;
    TArray <TArray<CLGComplex>> rets;
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(pFermion->m_byFieldId, _T(__FILE__), __LINE__));
    if (NULL == pF1)
    {
        appCrucial(_T("CMeasureDiagnal::ExportDiagnalStaggeredSU3 only work with CFieldFermionKSSU3 and CFieldGaugeSU3"));
        return;
    }

    UINT uiSiteCount = pF1->GetSiteCount();
    deviceWilsonVectorSU3* hostv = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * uiSiteCount);

    for (UINT x = 0; x < uiSiteCount; ++x)
    {
        for (BYTE spinor = 0; spinor < 4; ++spinor)
        {
            BYTE maxC = bOnlyRed ? 1 : 3;
            for (BYTE c = 0; c < maxC; ++c)
            {
                TArray<CLGComplex> ret;
                SFermionBosonSource source;
                source.m_eSourceType = EFS_Point;
                source.m_byColorIndex = c;
                source.m_bySpinIndex = spinor;
                source.m_sSourcePoint = __hostSiteIndexToInt4(x);
                pF1->InitialAsSource(source);
                Real scale = F(1.0);

                switch (eType)
                {
                case EMD_D:
                    {
                        pF1->D(gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields);
                    }
                    break;
                case EMD_InverseD:
                    {
                        pF1->InverseD(gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields);
                    }
                    break;
                case EMD_Gamma1:
                case EMD_Gamma2:
                case EMD_Gamma3:
                case EMD_Gamma4:
                case EMD_Gamma5:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, static_cast<EGammaMatrix>(static_cast<INT>(GAMMA1) + static_cast<INT>(eType - EMD_Gamma1)));
                    }
                    break;
                case EMD_Sigma12:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, SIGMA12);
                    }
                    break;
                case EMD_Sigma13:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, SIGMA31);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Sigma14:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, SIGMA41);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Sigma23:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, SIGMA23);
                    }
                    break;
                case EMD_Sigma24:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, SIGMA42);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Sigma34:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, SIGMA43);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Gamma51:
                case EMD_Gamma52:
                case EMD_Gamma53:
                case EMD_Gamma54:
                    {
                        pF1->ApplyGamma(0, 0, NULL, NULL, static_cast<EGammaMatrix>(static_cast<INT>(GAMMA51) + static_cast<INT>(eType - EMD_Gamma51)));
                    }
                    break;
                default:
                    {
                        appCrucial(_T("eType not implemented for ExportDiagnalWilsonSU3: %s"), __ENUM_TO_STRING(EMeasureDiagnal, eType).c_str());
                    }
                    break;
                }

                checkCudaErrors(cudaMemcpy(hostv, pF1->m_pDeviceData, sizeof(deviceSU3Vector) * uiSiteCount, cudaMemcpyDeviceToHost));
                for (UINT y = 0; y < uiSiteCount; ++y)
                {
                    for (BYTE spinor2 = 0; spinor2 < 4; ++spinor2)
                    {
                        for (BYTE c2 = 0; c2 < maxC; ++c2)
                        {
                            ret.AddItem(cuCmulf_cr(hostv[y].m_d[spinor2].m_ve[c2], scale));
                        }
                    }
                }

                rets.AddItem(ret);
            }
        }
    }

    WriteComplexArray2(sFileName, rets);
    appSafeFree(hostv);
}

void ExportDiagnalStaggeredSU3(const CCString& sFileName, EMeasureDiagnal eType, 
    INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields,
    const CFieldFermionKSSU3* pFermion)
{
    UBOOL bOnlyRed = TRUE;
    TArray <TArray<CLGComplex>> rets;
    CFieldFermionKSSU3* pF1 = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(pFermion->m_byFieldId, _T(__FILE__), __LINE__));
    if (NULL == pF1)
    {
        appCrucial(_T("CMeasureDiagnal::ExportDiagnalStaggeredSU3 only work with CFieldFermionKSSU3 and CFieldGaugeSU3"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(gaugeFields[0]);
    CFieldFermionKSSU3* pF2 = dynamic_cast<CFieldFermionKSSU3*>(pF1->GetCopy());

    UINT uiSiteCount = pF1->GetSiteCount();
    deviceSU3Vector* hostv = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * uiSiteCount);

    for (UINT x = 0; x < uiSiteCount; ++x)
    {
        BYTE maxC = bOnlyRed ? 1 : 3;
        for (BYTE c = 0; c < maxC; ++c)
        {
            TArray<CLGComplex> ret;

            SFermionBosonSource source;
            source.m_eSourceType = EFS_Point;
            source.m_byColorIndex = c;
            source.m_sSourcePoint = __hostSiteIndexToInt4(x);
            pF1->InitialAsSource(source);

            switch (eType)
            {
            case EMD_D:
                {
                    pF1->D(gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields);
                }
                break;
            case EMD_InverseD:
                {
                    pF1->InverseD(gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields);
                }
                break;
            case EMD_Gamma1:
            case EMD_Gamma2:
            case EMD_Gamma3:
            case EMD_Gamma4:
            case EMD_Gamma5:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, static_cast<EGammaMatrix>(static_cast<INT>(GAMMA1) + static_cast<INT>(eType - EMD_Gamma1)));
                }
                break;
            case EMD_Sigma12:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA12);
                }
                break;
            case EMD_Sigma13:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA31);
                }
                break;
            case EMD_Sigma14:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA41);
                }
                break;
            case EMD_Sigma23:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA23);
                }
                break;
            case EMD_Sigma24:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA42);
                }
                break;
            case EMD_Sigma34:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA43);
                }
                break;
            case EMD_Gamma51:
            case EMD_Gamma52:
            case EMD_Gamma53:
            case EMD_Gamma54:
                {
                    pF1->ApplyGamma(gaugeNum, bosonNum, gaugeFields, pBoson, static_cast<EGammaMatrix>(static_cast<INT>(GAMMA51) + static_cast<INT>(eType - EMD_Gamma51)));
                }
                break;
            case EMD_Oribital:
                {
                    pF1->CopyTo(pF2);
                    CMeasureAngularMomentumKS::ApplyOrbitalMatrix(pF1->m_pDeviceData, pF2->m_pDeviceData, pGaugeSU3->m_pDeviceData, pF1->m_byFieldId, pGaugeSU3->m_byFieldId);
                }
                break;
            case EMD_Spin:
                {
                    pF1->CopyTo(pF2);
                    CMeasureAngularMomentumKS::ApplySpinMatrix(pF1->m_pDeviceData, pF2->m_pDeviceData, pGaugeSU3->m_pDeviceData, pF1->m_byFieldId, pGaugeSU3->m_byFieldId);
                }
                break;
            default:
                {
                    appCrucial(_T("eType not implemented for ExportDiagnalWilsonSU3: %s"), __ENUM_TO_STRING(EMeasureDiagnal, eType).c_str());
                }
                break;
            }
            checkCudaErrors(cudaMemcpy(hostv, pF1->m_pDeviceData, sizeof(deviceSU3Vector) * uiSiteCount, cudaMemcpyDeviceToHost));

            for (UINT y = 0; y < uiSiteCount; ++y)
            {
                for (BYTE c2 = 0; c2 < maxC; ++c2)
                {
                    ret.AddItem(hostv[y].m_ve[c2]);
                }
            }

            rets.AddItem(ret);
        }
    }

    WriteComplexArray2Simple(sFileName, rets);
    appSafeFree(hostv);
    appSafeDelete(pF2);
}


CFieldFermion::CFieldFermion()
    : CField()
    , m_uiSiteCount(_HC_Volume)
    //Set for real at allocation time by the concrete subclass. 0 keeps single-GPU
    /// unsplit builds unchanged (mirrors CFieldGauge::m_uiHaloLinkCount).
    , m_uiHaloSiteCount(0)
    , m_bEvenPseudofermion(FALSE)
    , m_iMCIndex(-1)
    , m_iMDIndex(-1)
    , m_eRational(ER_NoRational)
    // , m_bDoperatorUseEffectiveGauge(FALSE)
{

}

CFieldMatrixOperation* CFieldMatrixOperation::Create(EFieldType ef)
{
    if (ef == EFT_FermionWilsonSquareSU3)
    {
        return new CFieldMatrixOperationWilsonSquareSU3();
    }

    if (ef == EFT_FermionStaggeredSU3)
    {
        return new CFieldMatrixOperationKSSU3();
    }

    appCrucial(_T("Matrix operation for field type %s not implemented!\n"), __ENUM_TO_STRING(EFieldType, ef).c_str());
    return NULL;
}

UBOOL CFieldFermion::RationalApproximation(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, INT iRationalIndex, UBOOL bAction)
{
    CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
    if (NULL == solver)
    {
        return FALSE;
    }
    _RECORD(CFieldFermion::RationalApproximation);
    TArray<CField*> solutions;
    TArray<CLGComplex> shifts;
    for (UINT i = 0; i < GRASet.m_pRASet[iRationalIndex]->m_uiDegree; ++i)
    {
        CField* pPooled = appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__);
        solutions.AddItem(pPooled);
        shifts.AddItem(_make_cuComplex(GRASet.m_pRASet[iRationalIndex]->m_lstB[i], F(0.0)));
    }

    solver->Solve(solutions, shifts, this, gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields, op);

    UINT uiStart = 0;
    if (bAction)
    {
        //Zero();
        solutions[0]->CopyBufferTo(this);
        ScalarMultply(GRASet.m_pRASet[iRationalIndex]->m_lstA[0]);
        uiStart = 1;
        solutions[0]->Return();
    }
    else
    {
        ScalarMultply(GRASet.m_pRASet[iRationalIndex]->m_fC);
    }

    for (UINT i = uiStart; i < GRASet.m_pRASet[iRationalIndex]->m_uiDegree; ++i)
    {
        Axpy(GRASet.m_pRASet[iRationalIndex]->m_lstA[i], solutions[i]);
        solutions[i]->Return();
    }
    return TRUE;
}

UBOOL CFieldFermion::RationalApproximationPooled(EFieldOperator op, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, INT iRationalIndex, TArray<CField*>& solutions) const
{
    CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
    if (NULL == solver)
    {
        return FALSE;
    }
    TArray<CLGComplex> shifts;
    for (UINT i = 0; i < GRASet.m_pRASet[iRationalIndex]->m_uiDegree; ++i)
    {
        CField* pPooled = appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__);
        solutions.AddItem(pPooled);
        shifts.AddItem(_make_cuComplex(GRASet.m_pRASet[iRationalIndex]->m_lstB[i], F(0.0)));
    }

    solver->Solve(solutions, shifts, this, gaugeNum, bosonNum, tensor2Num, gaugeFields, pBoson, tensor2Fields, op);
    return TRUE;
}

UBOOL CFieldFermion::InverseD(INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_D);
}

UBOOL CFieldFermion::InverseDdagger(INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_Ddagger);
}

UBOOL CFieldFermion::InverseDDdagger(INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_DDdagger);
}

UBOOL CFieldFermion::InverseDD(INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_DD);
}

void CFieldFermion::InitialOtherParameters(CParameters& params)
{
    CField::InitialOtherParameters(params);

    INT iEven = 0;
    params.FetchValueINT(_T("Even"), iEven);
    m_bEvenPseudofermion = (0 != iEven);

    TArray<Real> coeffs;
    params.FetchValueArrayReal(_T("MD"), coeffs);
    if (0 == coeffs.Num())
    {
        appGeneral(_T("no rational approximation configured for fermin MC!\n"));
        coeffs.AddItem(F(1.0));
    }
    m_iMDIndex = GRASet.Add(coeffs);

    params.FetchValueArrayReal(_T("MC"), coeffs);
    if (0 == coeffs.Num())
    {
        appGeneral(_T("no rational approximation configured for fermin MD!\n"));
        coeffs.AddItem(F(1.0));
    }
    m_iMCIndex = GRASet.Add(coeffs);

    CCString sEnumValue = _T("ER_NoRational");
    params.FetchStringValue(_T("Rational"), sEnumValue);
    m_eRational = __STRING_TO_ENUM(ERational, sEnumValue);
}

CCString CFieldFermion::GetInfos(const CCString& tab) const
{
    CCString sRet = CField::GetInfos(tab);
    sRet = sRet + tab + _T("Even Pseudofermion : ") + appToString(m_bEvenPseudofermion) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (energy and force) : ") + appToString(GRASet.m_pRASet[m_iMDIndex]->m_fC) + _T("Num:") + appToString(GRASet.m_pRASet[m_iMDIndex]->m_lstA) + _T("Don:") + appToString(GRASet.m_pRASet[m_iMDIndex]->m_lstB) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (Gaussian noise momentum) : ") + appToString(GRASet.m_pRASet[m_iMCIndex]->m_fC) + _T("Num:") + appToString(GRASet.m_pRASet[m_iMCIndex]->m_lstA) + _T("Don:") + appToString(GRASet.m_pRASet[m_iMCIndex]->m_lstB) + _T("\n");
    sRet = sRet + tab + _T("Rational : ") + __ENUM_TO_STRING(ERational, m_eRational) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================