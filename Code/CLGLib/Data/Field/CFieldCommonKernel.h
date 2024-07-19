//=============================================================================
// FILENAME : CFieldCommonKernel.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [07/20/2024 nbale]
//=============================================================================
#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"

#ifndef _CFIELDCOMMON_KERNEL_H_
#define _CFIELDCOMMON_KERNEL_H_

__BEGIN_NAMESPACE

//to avoid explict instantiation of every fucntion, we just explict instantiation a class
template<typename T>
class __DLL_EXPORT CCommonKernel
{
public:
    static void AllocateBuffer(T*& pointer, UINT count);
    static void FreeBuffer(T*& pointer);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDCOMMON_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================