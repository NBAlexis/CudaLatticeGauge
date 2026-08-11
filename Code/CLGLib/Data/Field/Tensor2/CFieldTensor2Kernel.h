//=============================================================================
// FILENAME : CFieldTensor2Kernel.h
//
// DESCRIPTION:
// The kernels for tensor2 (plaquette) fields.
// The random number seed table has only Dir streams per site, so the initial
// of a tensor2 field must be done site by site, drawing the 6 plaquettes
// of one site from the same stream sequentially.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================

#ifndef _CFIELDTENSOR2KERNEL_H_
#define _CFIELDTENSOR2KERNEL_H_

__BEGIN_NAMESPACE

template<typename T>
#if _CLG_WIN
class __DLL_EXPORT CFieldTensor2Kernel
#else
class CFieldTensor2Kernel
#endif
{
public:
    static void Initial(T* pointer, UINT uiSiteCount, UINT uiPlaqutteCount, EFieldInitialType eInitialType);
};

#if !_CLG_WIN
extern template class CFieldTensor2Kernel<Real>;
extern template class CFieldTensor2Kernel<CLGComplex>;
extern template class CFieldTensor2Kernel<deviceSU2>;
extern template class CFieldTensor2Kernel<deviceSU3>;
extern template class CFieldTensor2Kernel<deviceZN<3>>;
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDTENSOR2KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================
