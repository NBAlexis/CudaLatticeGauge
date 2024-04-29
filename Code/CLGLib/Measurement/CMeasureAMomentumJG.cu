//=============================================================================
// FILENAME : CMeasureAMomemtumJG.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAMomentumJG)

#pragma region kernles


/**
* calculate momentum, and sum over y and z
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAngularMomentumJG(
    const deviceSU3* __restrict__ pDeviceData,
    Real* pBuffer, 
    Real betaOverN,
    BYTE byFieldId)
{
    intokernalOnlyInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);
    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {

        //======================================================
        //4-chair terms except for the last one
        betaOverN = F(0.125) * betaOverN;

        const Real fX = (sSite4.x - _DC_Centerx);

        //===============
        //+x Omega V412
        const Real fV412 = fX * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

        //===============
        //+x Omega V432
        const Real fV432 = fX * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

        const Real fY = -(sSite4.y - _DC_Centery);

        //===============
        //-y Omega V421
        const Real fV421 = fY * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

        //===============
        //-y Omega V431
        const Real fV431 = fY * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);

        fRes = (fV412 + fV432 + fV421 + fV431) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -fRes);
}

/**
* Copy from _kernelAdd4PlaqutteTermSU3_Test and _kernelAddChairTermSU3_Term5
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAngularMomentumS2(
    const deviceSU3* __restrict__ pDeviceData,
    Real* pBuffer,
    Real betaOverN,
    BYTE byFieldId)
{
    intokernalOnlyInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);

    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {
        const Real betaOverN1over8 = F(0.125) * betaOverN;
        const Real fXYOmega2 = -(sSite4.x - _DC_Centerx) * (sSite4.y - _DC_Centery);

        //===============
        //-Omega^2 xy V132
        const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

        fRes = fV132 * betaOverN1over8;

        Real fXSq = (sSite4.x - _DC_Centerx);
        fXSq = fXSq * fXSq;
        Real fYSq = (sSite4.y - _DC_Centery);
        fYSq = fYSq * fYSq;

        //======================================================
        //4-plaqutte terms
        //Omega^2 x^2 Retr[1 - U_2,3]
        const Real fU23 = fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiN, sSite4, byFieldId);

        //Omega^2 y^2 Retr[1 - U_1,3]
        const Real fU13 = fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiN, sSite4, byFieldId);

        //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
        const Real fU12 = (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiN, sSite4, byFieldId);

        fRes += (fU23 + fU13 + fU12) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes);
}


/**
 * In _deviceChairTerm, the projective plane is already considered
 * So we need only to check the x and y
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAngularMomentumJGProjectivePlane(
    const deviceSU3* __restrict__ pDeviceData,
    Real* pBuffer,
    Real betaOverN,
    BYTE byFieldId)
{
    intokernalOnlyInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);
    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {

        //======================================================
        //4-chair terms except for the last one
        betaOverN = F(0.125) * betaOverN;

        const Real fX = (sSite4.x - _DC_Centerx + F(0.5));

        //===============
        //+x Omega V412
        const Real fV412 = fX * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

        //===============
        //+x Omega V432
        const Real fV432 = fX * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

        const Real fY = -(sSite4.y - _DC_Centery + F(0.5));

        //===============
        //-y Omega V421
        const Real fV421 = fY * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

        //===============
        //-y Omega V431
        const Real fV431 = fY * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);

        fRes = (fV412 + fV432 + fV421 + fV431) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -fRes);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAngularMomentumS2ProjectivePlane(
    const deviceSU3* __restrict__ pDeviceData,
    Real* pBuffer,
    Real betaOverN,
    BYTE byFieldId)
{
    intokernalOnlyInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);

    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {
        const Real betaOverN1over8 = -F(0.125) * betaOverN;
        const Real fXYOmega2 = (sSite4.x - _DC_Centerx + F(0.5)) * (sSite4.y - _DC_Centery + F(0.5));

        //===============
        //-Omega^2 xy V132
        const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

        fRes = fV132 * betaOverN1over8;

        Real fXSq = (sSite4.x - _DC_Centerx + F(0.5));
        fXSq = fXSq * fXSq;
        Real fYSq = (sSite4.y - _DC_Centery + F(0.5));
        fYSq = fYSq * fYSq;

        //======================================================
        //4-plaqutte terms
        //Omega^2 x^2 Retr[1 - U_2,3]
        const Real fU23 = fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiN, sSite4, byFieldId);

        //Omega^2 y^2 Retr[1 - U_1,3]
        const Real fU13 = fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiN, sSite4, byFieldId);

        //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
        const Real fU12 = (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiN, sSite4, byFieldId);

        fRes += (fU23 + fU13 + fU12) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes);
}

/**
 * To calculate spin, there is no site out of lattice
 * we use E_i = -{U_{4i}}_{TA}
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGaugeSpin(
    const deviceSU3* __restrict__ pE, 
    const deviceSU3* __restrict__ piAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    
    if (!site.IsDirichlet())
    {
        //Note iAphys is ia x A_phys
        UINT uiLinkX = _deviceGetLinkIndex(uiSiteIndex, 0);
        UINT uiLinkY = _deviceGetLinkIndex(uiSiteIndex, 1);
        deviceSU3 beforeTrace = pE[uiLinkX].MulC(piAphys[uiLinkY]);
        beforeTrace.Sub(pE[uiLinkY].MulC(piAphys[uiLinkX]));
        CLGComplex cRes = cuCmulf_cr(beforeTrace.Tr(), -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], cRes.x);
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateJGSurf(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGauge,
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (!site.IsDirichlet())
    {
        //const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
        //const BYTE uiDir2 = uiDir * 2;
        const Real fY = static_cast<Real>(sSite4.y - _DC_Centery);
        const Real fX = static_cast<Real>(sSite4.x - _DC_Centerx);

        Real fRes = F(0.0);
        for (BYTE dir = 0; dir < 3; ++dir)
        {
            //p_i A_j = U_i(n) A_j(n+i) U_i^+(n) - U_i^+(n-i) A_j(n-i) U_i(n-i)
            const SSmallInt4 x_p_i_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
            const UINT x_p_i_bi = __idx->_deviceGetBigIndex(x_p_i_site);
            const SIndex& x_p_i_x_idx = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_i_bi * 4 + 0];
            const SIndex& x_p_i_y_idx = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_i_bi * 4 + 1];

            //x p_i A_y : U_i(n) A_y(n+i) U_i^+(n)
            //- y p_i A_x : U_i(n) A_x(n+i) U_i^+(n)
            //U_i(n) [x A_y(n+i) - y A_x(n+i)] U_i^+(n)
            deviceSU3 u(_deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, dir));
            deviceSU3 a(_deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_i_y_idx));
            a.MulReal(fX);
            a.Sub(_deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_i_x_idx).MulRealC(fY));
            a.MulDagger(u);
            u.Mul(a);
            //E_i (x p_i A_y  - y p_i A_x)
            u = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(u);
            fRes += u.ReTr();

            //x p_i A_y : -U_i^+(n-i) A_y(n-i) U_i(n-i)
            //- y p_i A_x : -U_i^+(n-i) A_x(n-i) U_i(n-i)
            /*
            deviceSU3 u(_deviceGetGaugeBCSU3DirOne(pGauge, x_m_i_Gauge, dir));
            a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_i_Gauge, 1);
            a.MulReal(fX);
            a.Sub(_deviceGetGaugeBCSU3DirZero(pAphys, x_m_i_Gauge, 0).MulRealC(fY));
            a.Mul(u);
            u.DaggerMul(a);
            u = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(u);
            fRes -= u.ReTr();
            */
            a = _deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 1);
            a.MulReal(fX);
            a.Sub(_deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 0).MulRealC(fY));
            u = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(a);
            fRes -= u.ReTr();
        }

        //atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -fRes * fBetaOverN * F(0.5));
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes * fBetaOverN);
    }
}

/**
 * Already use m_pDeviceIndexLinkToSIndex
 * So we just concerntrate on the X and Y
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateJGSurfProjectivePlane(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGauge,
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (!site.IsDirichlet())
    {
        //const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
        //const BYTE uiDir2 = uiDir * 2;
        const Real fY = static_cast<Real>(sSite4.y - _DC_Centery + F(0.5));
        const Real fX = static_cast<Real>(sSite4.x - _DC_Centerx + F(0.5));

        Real fRes = F(0.0);
        #pragma unroll
        for (BYTE dir = 0; dir < 3; ++dir)
        {
            //p_i A_j = U_i(n) A_j(n+i) U_i^+(n) - U_i^+(n-i) A_j(n-i) U_i(n-i)
            const SSmallInt4 x_p_i_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
            const UINT x_p_i_bi = __idx->_deviceGetBigIndex(x_p_i_site);
            const SIndex& x_p_i_x_idx = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_i_bi * 4 + 0];
            const SIndex& x_p_i_y_idx = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_i_bi * 4 + 1];

            //x p_i A_y : U_i(n) A_y(n+i) U_i^+(n)
            //- y p_i A_x : U_i(n) A_x(n+i) U_i^+(n)
            //U_i(n) [x A_y(n+i) - y A_x(n+i)] U_i^+(n)
            deviceSU3 u(_deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, dir));
            deviceSU3 a(_deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_i_y_idx));
            a.MulReal(fX);
            a.Sub(_deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_i_x_idx).MulRealC(fY));
            a.MulDagger(u);
            u.Mul(a);
            //E_i (x p_i A_y  - y p_i A_x)
            u = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(u);
            fRes += u.ReTr();

            //x p_i A_y : -U_i^+(n-i) A_y(n-i) U_i(n-i)
            //- y p_i A_x : -U_i^+(n-i) A_x(n-i) U_i(n-i)
            /*
            u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_i_Gauge, dir);
            a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_i_Gauge, 1);
            a.MulReal(fX);
            a.Sub(_deviceGetGaugeBCSU3DirZero(pAphys, x_m_i_Gauge, 0).MulRealC(fY));
            a.Mul(u);
            u.DaggerMul(a);
            u = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(u);
            fRes -= u.ReTr();
            */
            a = _deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 1);
            a.MulReal(fX);
            a.Sub(_deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 0).MulRealC(fY));
            u = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(a);
            fRes -= u.ReTr();
        }

        //atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -fRes * fBetaOverN * F(0.5));
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes * fBetaOverN);
    }
}

/**
 * This is (nabla . E) r x A
 * In z-dir this is (nabla . E) (x Ay - y Ax)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateJGPot(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const UINT uiNablaE = _deviceGetLinkIndex(site.m_uiSiteIndex, 3);
    if (!site.IsDirichlet())
    {
        const Real fY = static_cast<Real>(sSite4.y - _DC_Centery);
        const Real fX = static_cast<Real>(sSite4.x - _DC_Centerx);

        deviceSU3 nablaE(pE[uiNablaE]);
        deviceSU3 a(_deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 1));
        a.MulReal(fX);
        a.Sub(_deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 0).MulRealC(fY));
        nablaE.Mul(a);

        //atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -fRes * fBetaOverN * F(0.5));
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -nablaE.ReTr() * fBetaOverN);
    }
}

/**
 * nablaE = { -sum_i U_{4,i}U_{4,-i} }_{TA}
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateJGPotProjectivePlane(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const UINT uiNablaE = _deviceGetLinkIndex(site.m_uiSiteIndex, 3);
    if (!site.IsDirichlet())
    {
        const Real fY = static_cast<Real>(sSite4.y - _DC_Centery + F(0.5));
        const Real fX = static_cast<Real>(sSite4.x - _DC_Centerx + F(0.5));

        deviceSU3 nablaE(pE[uiNablaE]);
        deviceSU3 a(_deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 1));
        a.MulReal(fX);
        a.Sub(_deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 0).MulRealC(fY));
        nablaE.Mul(a);

        //atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -fRes * fBetaOverN * F(0.5));
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -nablaE.ReTr() * fBetaOverN);
    }
}

/**
 * Use Apure directly
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChenApprox(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pApure,
    const deviceSU3* __restrict__ pAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const Real fmY = -static_cast<Real>(sSite4.y - _DC_Centery);
    const Real fmX = -static_cast<Real>(sSite4.x - _DC_Centerx);

    CLGComplex res = _zeroc;
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            deviceSU3 DxAphys = _deviceDPureMu(pAphys, pApure, sSite4, uiBigIdx, 0, dir, byFieldId);
            DxAphys.MulReal(fmY);
            deviceSU3 DyAphys = _deviceDPureMu(pAphys, pApure, sSite4, uiBigIdx, 1, dir, byFieldId);
            DyAphys.MulReal(fmX);
            DyAphys.Sub(DxAphys);
            res = _cuCaddf(res, _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(DyAphys).Tr());
        }

        res = cuCmulf_cr(res, -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChenApprox2(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pApure,
    const deviceSU3* __restrict__ pAphys,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalOnlyInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    const Real fmY = -static_cast<Real>(sSite4.y - _DC_Centery);
    const Real fmX = -static_cast<Real>(sSite4.x - _DC_Centerx);

    CLGComplex res = _zeroc;
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        for (BYTE dir = 0; dir < uiDir - 1; ++dir)
        {
            deviceSU3 DxAphys = _deviceDPureMu2(pAphys, pApure, sSite4, uiBigIdx, 0, dir, byFieldId);
            DxAphys.MulReal(fmY);
            deviceSU3 DyAphys = _deviceDPureMu2(pAphys, pApure, sSite4, uiBigIdx, 1, dir, byFieldId);
            DyAphys.MulReal(fmX);
            DyAphys.Sub(DxAphys);
            res = _cuCaddf(res, _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir).MulC(DyAphys).Tr());
        }

        res = cuCmulf_cr(res, -fBetaOverN);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
    }
}

/**
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomemtumJGChen(
    const deviceSU3* __restrict__ pE,
    const deviceSU3* __restrict__ pXcrossDpureA,
    Real* pBuffer,
    Real fBetaOverN)
{
    intokernalInt4;
    
    //const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    Real res = F(0.0);
    if (!site.IsDirichlet())
    {
        //only calculate x,y,z
        #pragma unroll
        for (BYTE dir = 0; dir < 3; ++dir)
        {
            deviceSU3 beforeTrace = _deviceGetGaugeBCSU3DirZero(pE, uiBigIdx, dir);
            beforeTrace.Mul(pXcrossDpureA[_deviceGetLinkIndex(uiSiteIndex, dir)]);
            res += beforeTrace.ReTr(); // _cuCaddf(res, beforeTrace.Tr());
        }

        //res = cuCmulf_cr(res, fBetaOverN * F(0.5));
        //res = cuCmulf_cr(res, fBetaOverN);
        //atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], res.x);
        atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], -res * fBetaOverN);
    }
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomentumJGChenDpureA(
    deviceSU3* pXcrossDpureA,
    const deviceSU3* __restrict__ pGauge,
    const deviceSU3* __restrict__ pAphys,
    BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const Real fmY = -static_cast<Real>(sSite4.y - _DC_Centery);
    const Real fmX = -static_cast<Real>(sSite4.x - _DC_Centerx);


    const SSmallInt4 x_p_x_site = _deviceSmallInt4OffsetC(sSite4, 1);
    const SSmallInt4 x_p_y_site = _deviceSmallInt4OffsetC(sSite4, 2);
    const UINT x_p_x_bi4 = __idx->_deviceGetBigIndex(x_p_x_site) * _DC_Dir;
    const UINT x_p_y_bi4 = __idx->_deviceGetBigIndex(x_p_y_site) * _DC_Dir;

    //idir = mu
    for (UINT idir = 0; idir < uiDir - 1; ++idir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pXcrossDpureA[uiLinkIndex] = deviceSU3::makeSU3Zero();

        //should not be idir + 1!
        const SIndex& x_p_x_idir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_x_bi4 + idir];
        const SIndex& x_p_y_idir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_y_bi4 + idir];

        //U_x(n) A_dir(n+x)U_x^+(n)
        deviceSU3 u(_deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, 0));
        deviceSU3 a(_deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_x_idir));
        a.MulDagger(u);
        u.Mul(a);
        u.MulReal(fmY);
        pXcrossDpureA[uiLinkIndex].Add(u);

        //U_x^+(n-x) A_dir(n-x) U_x(n-x)
        //u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_x_Gauge, 0);
        //a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_x_Gauge, idir);
        //a.Mul(u);
        //u.DaggerMul(a);
        //u.MulReal(fmY);
        //pXcrossDpureA[uiLinkIndex].Sub(u);
        a = _deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, idir);
        a.MulReal(fmY - fmX);
        pXcrossDpureA[uiLinkIndex].Sub(a);

        //U_y(n) A_dir(n+y)U_y^+(n)
        u = _deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, 1);
        a = _deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_y_idir);
        a.MulDagger(u);
        u.Mul(a);
        u.MulReal(fmX);
        pXcrossDpureA[uiLinkIndex].Sub(u);

        //U_y^+(n-y) A_dir(n-y) U_y(n-y)
        //u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_y_Gauge, 1);
        //a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_y_Gauge, idir);
        //a.Mul(u);
        //u.DaggerMul(a);
        //u.MulReal(fmX);
        //pXcrossDpureA[uiLinkIndex].Add(u);
    }
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelMomentumJGChenDpureAProjectivePlane(
    deviceSU3* pXcrossDpureA,
    const deviceSU3* __restrict__ pGauge,
    const deviceSU3* __restrict__ pAphys,
    BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const Real fmY = -static_cast<Real>(sSite4.y - _DC_Centery + F(0.5));
    const Real fmX = -static_cast<Real>(sSite4.x - _DC_Centerx + F(0.5));


    const SSmallInt4 x_p_x_site = _deviceSmallInt4OffsetC(sSite4, 1);
    const SSmallInt4 x_p_y_site = _deviceSmallInt4OffsetC(sSite4, 2);
    const UINT x_p_x_bi4 = __idx->_deviceGetBigIndex(x_p_x_site) * _DC_Dir;
    const UINT x_p_y_bi4 = __idx->_deviceGetBigIndex(x_p_y_site) * _DC_Dir;

    //idir = mu
    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pXcrossDpureA[uiLinkIndex] = deviceSU3::makeSU3Zero();

        //Link to SIndex should not have 1!
        const SIndex& x_p_x_idir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_x_bi4 + idir];
        const SIndex& x_p_y_idir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][x_p_y_bi4 + idir];

        //U_x(n) A_dir(n+x)U_x^+(n)
        deviceSU3 u(_deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, 0));
        deviceSU3 a(_deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_x_idir));
        a.MulDagger(u);
        u.Mul(a);
        u.MulReal(fmY);
        pXcrossDpureA[uiLinkIndex].Add(u); //this is -y U_x A_j(n+x) U_x+

        //U_x^+(n-x) A_dir(n-x) U_x(n-x)
        //u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_x_Gauge, 0);
        //a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_x_Gauge, idir);
        //a.Mul(u);
        //u.DaggerMul(a);
        //u.MulReal(fmY);
        //pXcrossDpureA[uiLinkIndex].Sub(u);
        a = _deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, idir);
        a.MulReal(fmY - fmX);
        pXcrossDpureA[uiLinkIndex].Sub(a); //this is (y-x)A_j

        //U_y(n) A_dir(n+y)U_y^+(n)
        u = _deviceGetGaugeBCSU3DirOne(pGauge, uiBigIdx, 1);
        a = _deviceGetGaugeBCSU3DirZeroSIndex(pAphys, x_p_y_idir);
        a.MulDagger(u);
        u.Mul(a);
        u.MulReal(fmX);
        pXcrossDpureA[uiLinkIndex].Sub(u); //this is +x U_y A_j(n+y) U_y+

        //U_y^+(n-y) A_dir(n-y) U_y(n-y)
        //u = _deviceGetGaugeBCSU3DirOne(pGauge, x_m_y_Gauge, 1);
        //a = _deviceGetGaugeBCSU3DirZero(pAphys, x_m_y_Gauge, idir);
        //a.Mul(u);
        //u.DaggerMul(a);
        //u.MulReal(fmX);
        //pXcrossDpureA[uiLinkIndex].Add(u);
    }
}

#pragma endregion

CMeasureAMomentumJG::~CMeasureAMomentumJG()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
    }
    if (NULL != m_pDeviceDataBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBuffer));
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionJG)
    {
        checkCudaErrors(cudaFree(m_pDistributionJG));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionJG)
    {
        free(m_pHostDistributionJG);
    }

    appSafeDelete(m_pE);
    appSafeDelete(m_pDpureA);
}

void CMeasureAMomentumJG::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly));
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("MeasureDist"), iValue);
    m_bMeasureDistribution = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("MeasureSpin"), iValue);
    m_bMeasureSpin = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("MeasureApprox"), iValue);
    m_bMeasureApprox = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("ProjectivePlane"), iValue);
    m_bProjectivePlane = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("NaiveNabla"), iValue);
    m_bNaiveNabla = iValue != 0;

    if (m_bMeasureSpin)
    {
        m_pE = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(GetGaugeFieldIdSingleField())->GetCopy());
        m_pDpureA = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(GetGaugeFieldIdSingleField())->GetCopy());
    }

    if (m_bMeasureDistribution)
    {
        //assuming the center is really at center
        //m_uiMaxR = ((_HC_Lx + 1) / 2) * ((_HC_Lx + 1) / 2)
        //         + ((_HC_Ly + 1) / 2) * ((_HC_Ly + 1) / 2);

        //m_uiEdgeR = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1);
        SetMaxAndEdge(&m_uiMaxR, &m_uiEdgeR, m_bProjectivePlane);

        checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1)));

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistributionJG = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
    }
}

void CMeasureAMomentumJG::OnConfigurationAcceptedSingleField(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

#if !_CLG_DOUBLEFLOAT
    const Real fBetaOverN = static_cast<Real>(CCommonData::m_fBeta / static_cast<DOUBLE>(_HC_SUN));
#else
    const Real fBetaOverN = CCommonData::m_fBeta / static_cast<Real>(_HC_SUN);
#endif

    _ZeroXYPlane(m_pDeviceDataBuffer);

    preparethread;

    if (m_bProjectivePlane)
    {
        _kernelCalculateAngularMomentumJGProjectivePlane << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            m_pDeviceDataBuffer,
            fBetaOverN,
            GetGaugeFieldIdSingleField());
    }
    else
    {
        _kernelCalculateAngularMomentumJG << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            m_pDeviceDataBuffer,
            fBetaOverN,
            GetGaugeFieldIdSingleField());
    }

    _AverageXYPlane(m_pDeviceDataBuffer);

    checkCudaErrors(cudaGetLastError());

    if (m_bMeasureDistribution)
    {
        XYDataToRdistri_R(m_bProjectivePlane, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
            m_uiMaxR, TRUE, GetGaugeFieldIdSingleField());

        checkCudaErrors(cudaGetLastError());

        //extract res
        checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        FillDataWithR_R(
            m_lstJG, &m_lstJGInner, m_lstJGAll, m_lstR, 
            m_pHostDistributionJG, m_pHostDistributionR, 
            m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR,
            F(1.0),
            TRUE
            );

        checkCudaErrors(cudaGetLastError());
    }

    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));

    if (m_bShowResult)
    {
        appDetailed(_T(" === Angular Momentum JG of site y=%d ======\n"), _HC_Centery);
    }
    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            m_lstRes.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
        }
    }

    for (UINT i = 1; i < _HC_Lx; ++i)
    {
        if (m_bShowResult)
        {
            appDetailed(_T("(%d,%d)=%1.6f  "), 
                i, 
                _HC_Centery, 
                m_pHostDataBuffer[i * _HC_Ly + _HC_Centery]);
        }
    }
    if (m_bShowResult)
    {
        appDetailed(_T("\n"));
    }

#pragma region measure s2

    _ZeroXYPlane(m_pDeviceDataBuffer);

    if (m_bProjectivePlane)
    {
        _kernelCalculateAngularMomentumS2ProjectivePlane << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            m_pDeviceDataBuffer,
            fBetaOverN,
            GetGaugeFieldIdSingleField());
    }
    else
    {
        _kernelCalculateAngularMomentumS2 << <block, threads >> > (
            pGaugeSU3->m_pDeviceData,
            m_pDeviceDataBuffer,
            fBetaOverN,
            GetGaugeFieldIdSingleField());
    }

    _AverageXYPlane(m_pDeviceDataBuffer);
    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaGetLastError());

    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            m_lstResJGS2.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
        }
    }

    if (m_bMeasureDistribution)
    {
        XYDataToRdistri_R(
            m_bProjectivePlane, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
            m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());

        checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        FillDataWithR_R(
            m_lstJGS2, &m_lstJGS2Inner, m_lstJGS2All, m_lstR,
            m_pHostDistributionJG, m_pHostDistributionR,
            m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
        );
    }

#pragma endregion


    if (m_bMeasureSpin)
    {
        _ZeroXYPlane(m_pDeviceDataBuffer);

        //projective plane has been considered
        pGaugeSU3->CalculateE_Using_U(m_pE);
        const CFieldGaugeSU3* pESU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pE);
        const CFieldGaugeSU3* pAphysSU3 = dynamic_cast<const CFieldGaugeSU3*>(appGetLattice()->m_pAphys);
        CFieldGaugeSU3* pDpureA = dynamic_cast<CFieldGaugeSU3*>(m_pDpureA);
        if (NULL == pAphysSU3)
        {
            appCrucial(_T("CMeasureAMomentumJG: A phys not calculated\n"));
        }
        else
        {
#pragma region Spin

            _kernelCalculateGaugeSpin << <block, threads >> > (
                pESU3->m_pDeviceData, 
                pAphysSU3->m_pDeviceData, 
                m_pDeviceDataBuffer,
                fBetaOverN);

            _AverageXYPlane(m_pDeviceDataBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaGetLastError());

            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGS.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_bProjectivePlane, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                    m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());
                
                checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGS, &m_lstJGSInner, m_lstJGSAll, m_lstR,
                    m_pHostDistributionJG, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
                );
            }

#pragma endregion

#pragma region Surf

            _ZeroXYPlane(m_pDeviceDataBuffer);
            if (m_bProjectivePlane)
            {
                _kernelCalculateJGSurfProjectivePlane << <block, threads >> > (
                    GetGaugeFieldIdSingleField(),
                    pGaugeSU3->m_pDeviceData,
                    pESU3->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    m_pDeviceDataBuffer,
                    fBetaOverN);
            }
            else
            {
                _kernelCalculateJGSurf << <block, threads >> > (
                    GetGaugeFieldIdSingleField(),
                    pGaugeSU3->m_pDeviceData,
                    pESU3->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    m_pDeviceDataBuffer,
                    fBetaOverN);
            }

            _AverageXYPlane(m_pDeviceDataBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaGetLastError());

            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGSurf.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_bProjectivePlane, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                    m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());

                checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGSurf, &m_lstJGSurfInner, m_lstJGSurfAll, m_lstR,
                    m_pHostDistributionJG, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
                );
            }

#pragma endregion

#pragma region Chen

            _ZeroXYPlane(m_pDeviceDataBuffer);

            if (m_bProjectivePlane)
            {
                _kernelMomentumJGChenDpureAProjectivePlane << <block, threads >> > (
                    pDpureA->m_pDeviceData,
                    pGaugeSU3->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    pAphysSU3->m_byFieldId
                    );
            }
            else
            {
                _kernelMomentumJGChenDpureA << <block, threads >> > (
                    pDpureA->m_pDeviceData,
                    pGaugeSU3->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    pAphysSU3->m_byFieldId
                    );
            }

            _kernelMomemtumJGChen << <block, threads >> > (
                pESU3->m_pDeviceData,
                pDpureA->m_pDeviceData,
                m_pDeviceDataBuffer,
                fBetaOverN);

            _AverageXYPlane(m_pDeviceDataBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaGetLastError());

            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGChen.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_bProjectivePlane, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                    m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());

                //extract res
                checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGChen, &m_lstJGChenInner, m_lstJGChenAll, m_lstR,
                    m_pHostDistributionJG, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
                );

                checkCudaErrors(cudaGetLastError());
            }

#pragma endregion

#pragma region JGPot

            //checkCudaErrors(cudaGetLastError());
            pGaugeSU3->CalculateNablaE_Using_U(m_pE, m_bNaiveNabla);
            //checkCudaErrors(cudaGetLastError());
            _ZeroXYPlane(m_pDeviceDataBuffer);
            if (m_bProjectivePlane)
            {
                _kernelCalculateJGPotProjectivePlane << <block, threads >> > (
                    GetGaugeFieldIdSingleField(),
                    pESU3->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    m_pDeviceDataBuffer,
                    fBetaOverN);
            }
            else
            {
                _kernelCalculateJGPot << <block, threads >> > (
                    GetGaugeFieldIdSingleField(),
                    pESU3->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    m_pDeviceDataBuffer,
                    fBetaOverN);
            }

            _AverageXYPlane(m_pDeviceDataBuffer);
            checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaGetLastError());

            for (UINT i = 1; i < _HC_Ly; ++i)
            {
                for (UINT j = 1; j < _HC_Lx; ++j)
                {
                    m_lstResJGPot.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                }
            }

            if (m_bMeasureDistribution)
            {
                XYDataToRdistri_R(
                    m_bProjectivePlane, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                    m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());

                checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                FillDataWithR_R(
                    m_lstJGPot, &m_lstJGPotInner, m_lstJGPotAll, m_lstR,
                    m_pHostDistributionJG, m_pHostDistributionR,
                    m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
                );
            }

#pragma endregion

            if (m_bMeasureApprox)
            {
                if (m_bProjectivePlane)
                {
                    appGeneral(_T("NOTE!!!: Projective plane does not support Approx"));
                }

                pGaugeSU3->CopyTo(pDpureA);
                pDpureA->TransformToIA();
                //[A, Aphys] = [Apure, Aphys], so we do not need to minus Aphys
                //pDpureA->AxpyMinus(pAphysSU3);

                _ZeroXYPlane(m_pDeviceDataBuffer);

                _kernelMomemtumJGChenApprox << <block, threads >> > (
                    pAphysSU3->m_byFieldId,
                    pESU3->m_pDeviceData,
                    pDpureA->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    m_pDeviceDataBuffer,
                    fBetaOverN
                    );

                _AverageXYPlane(m_pDeviceDataBuffer);
                checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real)* _HC_Lx* _HC_Ly, cudaMemcpyDeviceToHost));
                checkCudaErrors(cudaGetLastError());

                for (UINT i = 1; i < _HC_Ly; ++i)
                {
                    for (UINT j = 1; j < _HC_Lx; ++j)
                    {
                        m_lstResJGChenApprox.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                    }
                }

                if (m_bMeasureDistribution)
                {
                    XYDataToRdistri_R(
                        FALSE, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                        m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());

                    //extract res
                    checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                    FillDataWithR_R(
                        m_lstJGChenApprox, &m_lstJGChenApproxInner, m_lstJGChenApproxAll, m_lstR,
                        m_pHostDistributionJG, m_pHostDistributionR,
                        m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
                    );
                    checkCudaErrors(cudaGetLastError());
                }

                _ZeroXYPlane(m_pDeviceDataBuffer);

                _kernelMomemtumJGChenApprox2 << <block, threads >> > (
                    pAphysSU3->m_byFieldId,
                    pESU3->m_pDeviceData,
                    pDpureA->m_pDeviceData,
                    pAphysSU3->m_pDeviceData,
                    m_pDeviceDataBuffer,
                    fBetaOverN
                    );

                _AverageXYPlane(m_pDeviceDataBuffer);
                checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
                checkCudaErrors(cudaGetLastError());

                for (UINT i = 1; i < _HC_Ly; ++i)
                {
                    for (UINT j = 1; j < _HC_Lx; ++j)
                    {
                        m_lstResJGChenApprox2.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
                    }
                }

                if (m_bMeasureDistribution)
                {
                    XYDataToRdistri_R(
                        FALSE, m_pDeviceDataBuffer, m_pDistributionR, m_pDistributionJG,
                        m_uiMaxR, FALSE, GetGaugeFieldIdSingleField());

                    //extract res
                    checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
                    FillDataWithR_R(
                        m_lstJGChenApprox2, &m_lstJGChenApprox2Inner, m_lstJGChenApprox2All, m_lstR,
                        m_pHostDistributionJG, m_pHostDistributionR,
                        m_uiConfigurationCount, m_uiMaxR, m_uiEdgeR, F(1.0), FALSE
                    );
                    checkCudaErrors(cudaGetLastError());
                }
            }
        }
    }

    ++m_uiConfigurationCount;
}

void CMeasureAMomentumJG::Report()
{
    appPushLogDate(FALSE);

    appGeneral(_T("\n===================================================\n"));
    appGeneral(_T("=========== Angular Momentum JG of sites ==========\n"), _HC_Centerx);
    appGeneral(_T("===================================================\n"));

    ReportDistributionXY_R(m_uiConfigurationCount, m_lstRes);

    appGeneral(_T("===================================================\n"));

    appGeneral(_T("\n===================================================\n"));
    appGeneral(_T("=========== Angular Momentum S2 of sites ==========\n"), _HC_Centerx);
    appGeneral(_T("===================================================\n"));

    ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGS2);

    if (m_bMeasureSpin)
    {
        appGeneral(_T("\n===================================================\n"));
        appGeneral(_T("=========== Angular Momentum JGS of sites ==========\n"), _HC_Centerx);
        appGeneral(_T("===================================================\n"));

        ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGS);

        appGeneral(_T("\n===================================================\n"));
        appGeneral(_T("=========== Angular Momentum JGSurf of sites ==========\n"), _HC_Centerx);
        appGeneral(_T("===================================================\n"));

        ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGSurf);

        appGeneral(_T("\n===================================================\n"));
        appGeneral(_T("=========== Angular Momentum JG Chen of sites ==========\n"), _HC_Centerx);
        appGeneral(_T("===================================================\n"));

        ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGChen);

        appGeneral(_T("\n===================================================\n"));
        appGeneral(_T("=========== Angular Momentum JGPot of sites ==========\n"), _HC_Centerx);
        appGeneral(_T("===================================================\n"));

        ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGPot);

        if (m_bMeasureApprox)
        {
            appGeneral(_T("\n========================================================\n"));
            appGeneral(_T("=========== Angular Momentum JG Chen Approx of sites ==========\n"), _HC_Centerx);
            appGeneral(_T("========================================================\n"));

            ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGChenApprox);

            appGeneral(_T("\n========================================================\n"));
            appGeneral(_T("=========== Angular Momentum JG Chen Approx 2 of sites ==========\n"), _HC_Centerx);
            appGeneral(_T("========================================================\n"));

            ReportDistributionXY_R(m_uiConfigurationCount, m_lstResJGChenApprox2);
        }
    }

    appGeneral(_T("===================================================\n"));
    appGeneral(_T("===================================================\n"));

    appPopLogDate();
}

void CMeasureAMomentumJG::Reset()
{
    CMeasure::Reset();

    m_lstRes.RemoveAll();
    m_lstResJGS.RemoveAll();
    m_lstResJGChen.RemoveAll();
    m_lstResJGChenApprox.RemoveAll();
    m_lstResJGChenApprox2.RemoveAll();
    m_lstResJGSurf.RemoveAll();
    m_lstResJGPot.RemoveAll();
    m_lstResJGS2.RemoveAll();

    m_lstR.RemoveAll();
    m_lstJG.RemoveAll();
    m_lstJGS.RemoveAll();
    m_lstJGChen.RemoveAll();
    m_lstJGChenApprox.RemoveAll();
    m_lstJGChenApprox2.RemoveAll();
    m_lstJGSurf.RemoveAll();
    m_lstJGPot.RemoveAll();
    m_lstJGS2.RemoveAll();

    m_lstJGAll.RemoveAll();
    m_lstJGInner.RemoveAll();
    m_lstJGSAll.RemoveAll();
    m_lstJGSInner.RemoveAll();
    m_lstJGChenAll.RemoveAll();
    m_lstJGChenInner.RemoveAll();
    m_lstJGChenApproxAll.RemoveAll();
    m_lstJGChenApproxInner.RemoveAll();
    m_lstJGChenApprox2All.RemoveAll();
    m_lstJGChenApprox2Inner.RemoveAll();
    m_lstJGSurfAll.RemoveAll();
    m_lstJGSurfInner.RemoveAll();
    m_lstJGPotAll.RemoveAll();
    m_lstJGPotInner.RemoveAll();
    m_lstJGS2All.RemoveAll();
    m_lstJGS2Inner.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================