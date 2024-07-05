//=============================================================================
// FILENAME : CBase.h
// 
// DESCRIPTION:
// This is a simple RTTI implementation, can be use for serialize or factory etc
// Now, only for factory... Not necessary to implement a full version
// 
// Since c++ 11, the virtual is replaced with override
//
// REVISION:
//  [12/14/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CBASE_H_
#define _CBASE_H_

#define __CLGDECLARE_CLASS( class_name) \
public: \
struct __DLL_EXPORT class##class_name : public ClassList \
{ \
    class##class_name() \
    { \
        m_sClassName = _T(#class_name); \
        m_dwObjectSize = (UINT)(sizeof(class class_name)); \
        m_pfnCreateObject = class_name::_StaticCreateObject; \
        GClassGather.AddClass(this); \
    } \
}; \
static class##class_name m_StaticClass; \
static CClass * StaticClass(void) { return &m_StaticClass; } \
static CBase* _StaticCreateObject() { return new class_name(); } \
const CClass* GetClass() const override \
{ return &class_name::m_StaticClass; } \
protected:


//=================================================================
//========= To declare field, we force the copy function ==========

#define __CLGDECLARE_FIELD(class_name) \
__CLGDECLARE_CLASS(class_name) \
public: \
void CopyTo(CField* U) const override; \
CField* GetCopy() const override \
{ \
    class_name* ret = new class_name(); \
    CopyTo(ret); \
    return ret; \
} \
protected:

#define __CLGDECLARE_FIELDWITHOUTCOPYTO(class_name) \
__CLGDECLARE_CLASS(class_name) \
public: \
CField* GetCopy() const override \
{ \
    class_name* ret = new class_name(); \
    CopyTo(ret); \
    return ret; \
} \
protected:

//In Ubuntu, we are using static-link, see the problem at
//https://stackoverflow.com/questions/4767925/how-to-force-gcc-to-link-unreferenced-static-c-objects-from-a-library
//Tested that, the helper must be outside of the class


#define __CLGIMPLEMENT_CLASS(class_name) \
    class_name##helper::class_name##helper() {} \
    __DLL_EXPORT class_name::class##class_name class_name::m_StaticClass;


#define __CLG_REGISTER_HELPER_HEADER(class_name) \
struct CLGAPI class_name##helper \
{ \
    class_name##helper(); \
}; \
static class_name##helper s_##class_name##helper; 


__BEGIN_NAMESPACE

enum { kMaxClassNameLength = 127 };
class CLGAPI CClass
{
protected:

    const TCHAR* m_sClassName;
    UINT  m_dwObjectSize;
    class CBase* (*m_pfnCreateObject)();

public:

    UINT GetSize() const { return m_dwObjectSize; }
    class CBase* Create() const { return m_pfnCreateObject(); }
    const TCHAR* GetName() const { return m_sClassName; }
};

typedef TSimpleDoubleLinkedList<CClass> ClassList;


class CLGAPI CBase
{
public:
    CBase() { }
    virtual ~CBase() { }
    virtual const CClass* GetClass() const { return NULL; }
    virtual CCString GetInfos(const CCString& tab) const
    {
        return tab + _T("Name : ") + GetClass()->GetName() + _T("\n");
    }
};

class CLGAPI CClassGather
{
public:
    void AddClass(ClassList *static_class)
    {
        if (FindClass(static_class->GetName()))
            return;
        static_class->Link(m_pClasses);
    }

    CClass* FindClass(const TCHAR* class_name) const
    {
        ClassList * CurrentList = m_pClasses;
        for (; CurrentList; CurrentList = CurrentList->m_pNext)
        {
            if (appStrcmp(class_name, CurrentList->GetName()) == 0)
                return dynamic_cast<CClass*>(CurrentList);
        }
        return NULL;
    }

    void TraceAllClass() const
    {
        ClassList * CurrentList = m_pClasses;
        for (; CurrentList; CurrentList = CurrentList->m_pNext)
        {
            appGeneral(_T("==we have class: %s\n"), CurrentList->GetName());
        }
    }

    ClassList* m_pClasses;
};

extern CLGAPI CClassGather GClassGather;

inline CBase* appCreate(const CCString& name)
{
    CClass * pClass = GClassGather.FindClass(name);
    return (NULL != pClass) ? pClass->Create() : NULL;
}

__END_NAMESPACE

#endif//#ifndef _CBASE_H_

//=============================================================================
// END OF FILE
//=============================================================================