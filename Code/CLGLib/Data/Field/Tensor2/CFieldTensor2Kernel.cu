//=============================================================================
// FILENAME : CFieldTensor2Kernel.cu
//
// DESCRIPTION:
// The kernels for tensor2 (plaquette) fields
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldTensor2Kernel.h"

__BEGIN_NAMESPACE

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialTensor2(T* pDevicePtr, UINT uiSiteCount, UINT uiPlaqutteCount, EFieldInitialType eInitialType)
{
    intokernal;
    //the layout is component-major, data[plaqutteIndex * siteCount + siteIndex]
    //the seed table has only Dir streams per site, draw the plaquettes from the same stream sequentially
    const UINT uiSeedIndex = _deviceGetLinkIndex(uiSiteIndex, 0);
    for (UINT i = 0; i < uiPlaqutteCount; ++i)
    {
        switch (eInitialType)
        {
        case EFIT_Zero:
            {
                _Zero(pDevicePtr[i * uiSiteCount + uiSiteIndex]);
            }
            break;
        case EFIT_Identity:
            {
                _Id(pDevicePtr[i * uiSiteCount + uiSiteIndex]);
            }
            break;
        case EFIT_RandomGaussian:
        case EFIT_RandomGenerator:
            {
                pDevicePtr[i * uiSiteCount + uiSiteIndex] = _makeGaussian<T>(uiSeedIndex);
            }
            break;
        case EFIT_Random:
            {
                pDevicePtr[i * uiSiteCount + uiSiteIndex] = _makeRandom<T>(uiSeedIndex);
            }
            break;
        default:
            {
                printf("Tensor2 field cannot be initialized with this type! %d\n", eInitialType);
            }
            break;
        }
    }
}

template<typename T>
void CFieldTensor2Kernel<T>::Initial(T* pointer, UINT uiSiteCount, UINT uiPlaqutteCount, EFieldInitialType eInitialType)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelInitialTensor2<T>, block, threads, pointer, uiSiteCount, uiPlaqutteCount, eInitialType);
}

template class CFieldTensor2Kernel<Real>;
template class CFieldTensor2Kernel<CLGComplex>;
template class CFieldTensor2Kernel<deviceSU2>;
template class CFieldTensor2Kernel<deviceSU3>;
template class CFieldTensor2Kernel<deviceZN<3>>;

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
