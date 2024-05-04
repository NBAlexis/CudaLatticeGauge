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

__BEGIN_NAMESPACE

void ExportDiagnalWilsonSU3(const CCString& sFileName, EMeasureDiagnal eType, 
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson,
    const CFieldFermionWilsonSquareSU3* pFermion)
{
    UBOOL bOnlyRed = TRUE;
    TArray <TArray<CLGComplex>> rets;
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(pFermion->m_byFieldId));
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
                SFermionSource source;
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
                        pF1->D(gaugeNum, bosonNum, gaugeFields, pBoson);
                    }
                    break;
                case EMD_InverseD:
                    {
                        pF1->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
                    }
                    break;
                case EMD_Gamma1:
                case EMD_Gamma2:
                case EMD_Gamma3:
                case EMD_Gamma4:
                case EMD_Gamma5:
                    {
                        pF1->ApplyGamma(static_cast<EGammaMatrix>(static_cast<INT>(GAMMA1) + static_cast<INT>(eType - EMD_Gamma1)));
                    }
                    break;
                case EMD_Sigma12:
                    {
                        pF1->ApplyGamma(SIGMA12);
                    }
                    break;
                case EMD_Sigma13:
                    {
                        pF1->ApplyGamma(SIGMA31);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Sigma14:
                    {
                        pF1->ApplyGamma(SIGMA41);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Sigma23:
                    {
                        pF1->ApplyGamma(SIGMA23);
                    }
                    break;
                case EMD_Sigma24:
                    {
                        pF1->ApplyGamma(SIGMA42);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Sigma34:
                    {
                        pF1->ApplyGamma(SIGMA43);
                        scale = F(-1.0);
                    }
                    break;
                case EMD_Gamma51:
                case EMD_Gamma52:
                case EMD_Gamma53:
                case EMD_Gamma54:
                    {
                        pF1->ApplyGamma(static_cast<EGammaMatrix>(static_cast<INT>(GAMMA51) + static_cast<INT>(eType - EMD_Gamma51)));
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
    INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson,
    const CFieldFermionKSSU3* pFermion)
{
    UBOOL bOnlyRed = TRUE;
    TArray <TArray<CLGComplex>> rets;
    CFieldFermionKSSU3* pF1 = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(pFermion->m_byFieldId));
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

            SFermionSource source;
            source.m_eSourceType = EFS_Point;
            source.m_byColorIndex = c;
            source.m_sSourcePoint = __hostSiteIndexToInt4(x);
            pF1->InitialAsSource(source);

            switch (eType)
            {
            case EMD_D:
                {
                    pF1->D(gaugeNum, bosonNum, gaugeFields, pBoson);
                }
                break;
            case EMD_InverseD:
                {
                    pF1->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
                }
                break;
            case EMD_Gamma1:
            case EMD_Gamma2:
            case EMD_Gamma3:
            case EMD_Gamma4:
            case EMD_Gamma5:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, static_cast<EGammaMatrix>(static_cast<INT>(GAMMA1) + static_cast<INT>(eType - EMD_Gamma1)));
                }
                break;
            case EMD_Sigma12:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA12);
                }
                break;
            case EMD_Sigma13:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA31);
                }
                break;
            case EMD_Sigma14:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA41);
                }
                break;
            case EMD_Sigma23:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA23);
                }
                break;
            case EMD_Sigma24:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA42);
                }
                break;
            case EMD_Sigma34:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, SIGMA43);
                }
                break;
            case EMD_Gamma51:
            case EMD_Gamma52:
            case EMD_Gamma53:
            case EMD_Gamma54:
                {
                    pF1->ApplyGammaKS(gaugeNum, bosonNum, gaugeFields, pBoson, static_cast<EGammaMatrix>(static_cast<INT>(GAMMA51) + static_cast<INT>(eType - EMD_Gamma51)));
                }
                break;
            case EMD_Oribital:
                {
                    pF1->CopyTo(pF2);
                    CMeasureAngularMomentumKS::ApplyOrbitalMatrix(pF1->m_pDeviceData, pF2->m_pDeviceData, pGaugeSU3->m_pDeviceData, pF1->m_byFieldId);
                }
                break;
            case EMD_Spin:
                {
                    pF1->CopyTo(pF2);
                    CMeasureAngularMomentumKS::ApplySpinMatrix(pF1->m_pDeviceData, pF2->m_pDeviceData, pGaugeSU3->m_pDeviceData, pF1->m_byFieldId);
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

UBOOL CFieldFermion::InverseD(INT gaugeNum, INT bosonNum,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, gaugeFields, bosonFields, EFO_F_D);
}

UBOOL CFieldFermion::InverseDdagger(INT gaugeNum, INT bosonNum,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, gaugeFields, bosonFields, EFO_F_Ddagger);
}

UBOOL CFieldFermion::InverseDDdagger(INT gaugeNum, INT bosonNum,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, gaugeFields, bosonFields, EFO_F_DDdagger);
}

UBOOL CFieldFermion::InverseDD(INT gaugeNum, INT bosonNum,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields)
{
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, gaugeNum, bosonNum, gaugeFields, bosonFields, EFO_F_DD);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================