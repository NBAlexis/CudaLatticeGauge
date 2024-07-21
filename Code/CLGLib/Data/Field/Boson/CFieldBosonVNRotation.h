//=============================================================================
// FILENAME : CFieldBosonVNRotation.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [07/13/2024 nbale]
//=============================================================================
#include "CFieldBosonVN.h"

#define __DEFINE_ROTATION_BOSON(N) \
__CLG_REGISTER_HELPER_HEADER(CFieldBosonVNRotationSU##N) \
class CLGAPI CFieldBosonVNRotationSU##N : public CFieldBosonVNRotation<deviceSU##N##Vector, deviceSU##N> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldBosonVNRotationSU##N) \
public: \
    EFieldType GetFieldType() const override { return EFT_BosonComplexVector##N; } \
    UINT VectorN() const override { return N; } \
    UINT FloatN() const override { return 2 * N; } \
};

#ifndef _CFIELDBOSONVNROTATION_H_
#define _CFIELDBOSONVNROTATION_H_

__BEGIN_NAMESPACE



template<typename deviceDataBoson, typename deviceDataGauge>
class __DLL_EXPORT CFieldBosonVNRotation : public CFieldBosonVN<deviceDataBoson, deviceDataGauge>
{
public:
    CFieldBosonVNRotation()
        : CFieldBosonVN<deviceDataBoson, deviceDataGauge>()
        , m_fOmega(0.0)
        , m_pDevicePath(NULL)
        , m_bShiftCenter(FALSE)
        , m_pfCx(NULL)
        , m_pfCy(NULL)
        , m_pfCxy(NULL)
        , m_pfCxShift(NULL)
        , m_pfCyShift(NULL)
        , m_pfCxyShiftXYPP(NULL)
        , m_pfCxyShiftYXPP(NULL)
        , m_pfCxyShiftYXPM(NULL)
        , m_pfCxyShiftXYMP(NULL)
        , m_pfCxCySq(NULL)
        , m_pfCxCySqShift(NULL)
        , m_pfCxSq(NULL)
        , m_pfCySq(NULL)
        , m_pfCxSqShift(NULL)
        , m_pfCySqShift(NULL)
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::AllocatePathBuffer(&m_pDevicePath);

        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCx(&m_pfCx);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCy(&m_pfCy);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxy(&m_pfCxy);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxShift(&m_pfCxShift);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCyShift(&m_pfCyShift);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftXYPP(&m_pfCxyShiftXYPP);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftYXPP(&m_pfCxyShiftYXPP);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftYXPM(&m_pfCxyShiftYXPM);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxyShiftXYMP(&m_pfCxyShiftXYMP);

        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxCySq(&m_pfCxCySq);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxCySqShift(&m_pfCxCySqShift);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxSq(&m_pfCxSq);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCySq(&m_pfCySq);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCxSqShift(&m_pfCxSqShift);
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyFunctionPointCySqShift(&m_pfCySqShift);
    }

    ~CFieldBosonVNRotation()
    {
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::FreePathBuffer(m_pDevicePath);
    }

    void CopyTo(CField* U) const override
    {
        CFieldBosonVN<deviceDataBoson, deviceDataGauge>::CopyTo(U);
        CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>* pU = static_cast<CFieldBosonVNRotation<deviceDataBoson, deviceDataGauge>*>(U);
        if (NULL != pU)
        {
            pU->m_fOmega = m_fOmega;
            pU->m_bShiftCenter = m_bShiftCenter;
        }
    }

    void InitialOtherParameters(CParameters& param) override
    {
        CFieldBosonVN<deviceDataBoson, deviceDataGauge>::InitialOtherParameters(param);

        DOUBLE fValue = 0.1;
        if (param.FetchValueDOUBLE(_T("Omega"), fValue))
        {
            m_fOmega = fValue;
        }

        INT iValue = 0;
        if (param.FetchValueINT(_T("ShiftCenter"), iValue))
        {
            m_bShiftCenter = (0 != iValue);
        }
    }

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CFieldBosonVN<deviceDataBoson, deviceDataGauge>::GetInfos(tab);
        sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
        sRet = sRet + tab + _T("ShiftCenter : ") + appToString(m_bShiftCenter) + _T("\n");
        return sRet;
    }

    DOUBLE m_fOmega;
    INT* m_pDevicePath;
    UBOOL m_bShiftCenter;
    _deviceCoeffFunctionPointerTwoSites m_pfCx;
    _deviceCoeffFunctionPointerTwoSites m_pfCy;
    _deviceCoeffFunctionPointerTwoSites m_pfCxy;
    _deviceCoeffFunctionPointerTwoSites m_pfCxShift;
    _deviceCoeffFunctionPointerTwoSites m_pfCyShift;
    _deviceCoeffFunctionPointerTwoSites m_pfCxyShiftXYPP;
    _deviceCoeffFunctionPointerTwoSites m_pfCxyShiftYXPP;
    _deviceCoeffFunctionPointerTwoSites m_pfCxyShiftYXPM;
    _deviceCoeffFunctionPointerTwoSites m_pfCxyShiftXYMP;
    _deviceCoeffFunctionPointer m_pfCxCySq;
    _deviceCoeffFunctionPointer m_pfCxCySqShift;

    _deviceCoeffFunctionPointer m_pfCxSq;
    _deviceCoeffFunctionPointer m_pfCySq;
    _deviceCoeffFunctionPointer m_pfCxSqShift;
    _deviceCoeffFunctionPointer m_pfCySqShift;

    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const override
    {
        CFieldBosonVN<deviceDataBoson, deviceDataGauge>::ForceOnGauge(gaugeNum, bosonNum, pGauge, pGaugeForce, pBoson);

        if (this->m_byGaugeFieldIds.Num() < 1)
        {
            appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
            return;
        }
        INT gaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, pGauge, this->m_byGaugeFieldIds[0]);
        if (gaugeidx < 0 || gaugeidx >= this->m_byGaugeFieldIds.Num())
        {
            appCrucial(_T("CFieldBosonUN ForceOnGauge: there is no gauge!"));
            return;
        }

        const CFieldGauge* gauge = pGauge[gaugeidx];
        CFieldGauge* gaugeforce = pGaugeForce[gaugeidx];

        if (NULL == gauge || this->VectorN() != gauge->MatrixN())
        {
            appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
            return;
        }

        if (NULL == gaugeforce || this->VectorN() != gaugeforce->MatrixN())
        {
            appCrucial(_T("CFieldBosonUN can only play with gauge UN!"));
            return;
        }

#pragma region +2 Re[Omega y pt px ]

        INT hostPath[2] = { 1, 4 };
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            -m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2
        );

        hostPath[0] = 4; hostPath[1] = 1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            -m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2
        );

        hostPath[0] = 4; hostPath[1] = -1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2
        );

        hostPath[0] = -1; hostPath[1] = 4;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2
        );

#pragma endregion

#pragma region -2 Re[Omega x pt py ]

        hostPath[0] = 2; hostPath[1] = 4;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2
        );

        hostPath[0] = 4; hostPath[1] = 2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2
        );

        hostPath[0] = 4; hostPath[1] = -2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            -m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2
        );

        hostPath[0] = -2; hostPath[1] = 4;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            -m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2
        );

#pragma endregion

#pragma region -2 Re[Omega x y px py ]

        hostPath[0] = 1; hostPath[1] = 2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftXYPP : m_pfCxy,
            m_pDevicePath,
            2
        );

        hostPath[0] = 2; hostPath[1] = 1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftYXPP : m_pfCxy,
            m_pDevicePath,
            2
        );

        hostPath[0] = 2; hostPath[1] = -1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            -m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftYXPM : m_pfCxy,
            m_pDevicePath,
            2
        );

        hostPath[0] = -1; hostPath[1] = 2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLinkForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            -m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftXYMP : m_pfCxy,
            m_pDevicePath,
            2
        );

#pragma endregion

#pragma region +Omega^2 y^2 px^2 + Omega^2 x^2 py^2, shift term

        this->PartialSqForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCySqShift : m_pfCySq,
            0
        );

        this->PartialSqForceGauge(
            (const deviceDataGauge*)gauge->GetData(),
            gauge->m_byFieldId,
            (deviceDataGauge*)gaugeforce->GetData(),
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxSqShift : m_pfCxSq,
            1
        );

#pragma endregion

        //gaugeu1->DebugPrintMe();
        checkCudaErrors(cudaDeviceSynchronize());
    }

    void SetBosonOmega(DOUBLE fOmega)
    {
        m_fOmega = fOmega;
        this->UpdatePooledParamters();
    }

protected:

    /**
    * C [partial_{mu}phi partial_{nu}phi + partial_{nu}phi partial_{mu}phi]
    * = C phi^2/a^2  - C (1/4a^2) 2Re[phi(n+mu) to phi(n+nu) - phi(n+mu) to phi(n-nu) - phi(n-mu) to phi(n+nu) + phi(n-mu) to phi(n-nu)]
    * the 2 Re is in fact sum of two terms partial_{mu}phi partial_{nu}phi and partial_{nu}phi partial_{mu}phi which are Hermitian of each other
    *
    * partial _x phi partial _x phi
    * = + (2/a^2)phi^2 - phi(n)[U_{mu} phi(n+mu) + U_{-mu} phi(n-mu)]
    *
    */
    void DFromSource(const CFieldBoson* pSource, INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGaugeFields, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override
    {
        CFieldBosonVN<deviceDataBoson, deviceDataGauge>::DFromSource(pSource, gaugeNum, bosonNum, pGaugeFields, pBoson, eCoeffType, fCoeffReal, fCoeffImg);

        const CFieldGauge* pGauge = this->GetDefaultGauge(gaugeNum, pGaugeFields);

        if (NULL != pGauge && this->VectorN() != pGauge->MatrixN())
        {
            appCrucial(_T("CFieldBosonVU can only play with gauge UN!"));
            return;
        }

        if (NULL == pSource || this->GetFieldType() != pSource->GetFieldType())
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
            return;
        }

        const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* pSourceVN = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(pSource);
        if (NULL == pSourceVN)
        {
            appCrucial(_T("CFieldBosonU1 can only work with CFieldBosonU1!"));
            return;
        }

        //try external gauge field
        if (NULL == pGauge && this->m_byGaugeFieldIds.Num() > 0)
        {
            const CField* externelfield = appGetLattice()->GetFieldById(this->m_byGaugeFieldIds[0]);
            if (NULL != externelfield)
            {
                pGauge = dynamic_cast<const CFieldGauge*>(appGetLattice()->GetFieldById(this->m_byGaugeFieldIds[0]));
                if (pGauge->IsDynamic())
                {
                    appCrucial(_T("CFieldBosonUN: A dynamic field is configured for this UN, but not for the action!\n"));
                    pGauge = NULL;
                }
            }
        }

        Real fRealCoeff = fCoeffReal;
        const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
        if (EOCT_Minus == eCoeffType)
        {
            eCoeffType = EOCT_Real;
            fRealCoeff = F(-1.0);
        }

#pragma region +2 Re[Omega y pt px ]

        //====== +2 Re[Omega y pt px ]
        INT hostPath[2] = { 1, 4 };
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            -m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = 4; hostPath[1] = 1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            -m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = 4; hostPath[1] = -1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = -1; hostPath[1] = 4;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega,
            m_bShiftCenter ? m_pfCyShift : m_pfCy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

#pragma endregion

#pragma region -2 Re[Omega x pt py ]

        //====== -2 Re[Omega x pt py ]
        hostPath[0] = 2; hostPath[1] = 4;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = 4; hostPath[1] = 2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = 4; hostPath[1] = -2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            -m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = -2; hostPath[1] = 4;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            -m_fOmega,
            m_bShiftCenter ? m_pfCxShift : m_pfCx,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

#pragma endregion

#pragma region -2 Re[Omega x y px py ]

        //====== -2 Re[Omega x y px py ]
        hostPath[0] = 1; hostPath[1] = 2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftXYPP : m_pfCxy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = 2; hostPath[1] = 1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftYXPP : m_pfCxy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = 2; hostPath[1] = -1;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            -m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftYXPM : m_pfCxy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

        hostPath[0] = -1; hostPath[1] = 2;
        CFieldBosonVNKernel<deviceDataBoson, deviceDataGauge>::CopyPathBuffer(m_pDevicePath, hostPath, 2);
        this->OneLink(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            -m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxyShiftXYMP : m_pfCxy,
            m_pDevicePath,
            2,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

#pragma endregion

#pragma region +Omega^2 y^2 px^2 + Omega^2 x^2 py^2, square term

        this->DiagnalTerm(
            pSourceVN->m_pDeviceData,
            -F(2.0) * m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxCySqShift : m_pfCxCySq,
            eCoeffType,
            fRealCoeff,
            cCompCoeff);

#pragma endregion

#pragma region +Omega^2 y^2 px^2 + Omega^2 x^2 py^2, shift term

        this->PartialSq(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCySqShift : m_pfCySq,
            0,
            eCoeffType,
            fRealCoeff,
            cCompCoeff
        );

        this->PartialSq(
            pSourceVN->m_pDeviceData,
            NULL == pGauge ? NULL : (const deviceDataGauge*)pGauge->GetData(),
            NULL == pGauge ? 0 : pGauge->m_byFieldId,
            m_fOmega * m_fOmega,
            m_bShiftCenter ? m_pfCxSqShift : m_pfCxSq,
            1,
            eCoeffType,
            fRealCoeff,
            cCompCoeff
        );

#pragma endregion
    }
};

__CLG_REGISTER_HELPER_HEADER(CFieldBosonVNRotationU1)
class CLGAPI CFieldBosonVNRotationU1 : public CFieldBosonVNRotation<CLGComplex, CLGComplex>
{ 
__CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldBosonVNRotationU1)
public: 
    EFieldType GetFieldType() const override { return EFT_BosonComplex; }
    UINT VectorN() const override { return 1; } 
    UINT FloatN() const override { return 2; } 
};

__DEFINE_ROTATION_BOSON(2)
__DEFINE_ROTATION_BOSON(3)
__DEFINE_ROTATION_BOSON(4)
__DEFINE_ROTATION_BOSON(5)
__DEFINE_ROTATION_BOSON(6)
__DEFINE_ROTATION_BOSON(7)
__DEFINE_ROTATION_BOSON(8)

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVNROTATION_H_

//=============================================================================
// END OF FILE
//=============================================================================