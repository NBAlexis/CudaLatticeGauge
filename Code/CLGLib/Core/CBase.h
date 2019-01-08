//=============================================================================
// FILENAME : CBase.h
// 
// DESCRIPTION:
// This is a simple RTTI implementation, can be use for serialize or factory etc
// Now, only for factory... Not necessary to implement a full version
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
virtual const CClass* GetClass() const \
{ return &class_name::m_StaticClass; } \
protected:


//=================================================================
//========= To declare field, we force the copy function ==========

#define __CLGDECLARE_FIELD(class_name) \
__CLGDECLARE_CLASS(class_name) \
public: \
virtual void CopyTo(CField* U) const; \
virtual CField* GetCopy() const \
{ \
    class_name* ret = new class_name(); \
    CopyTo(ret); \
    return ret; \
} \
virtual CField* GetZero() const \
{ \
    class_name* ret = new class_name(); \
    ret->InitialField(EFIT_Zero); \
    return ret; \
} \
protected:


#define __CLGIMPLEMENT_CLASS(class_name) \
    __DLL_EXPORT class_name::class##class_name class_name::m_StaticClass;


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
    class CBase* Create() { return m_pfnCreateObject(); }
    const TCHAR* GetName() const { return m_sClassName; }
};

typedef TSimpleDoubleLinkedList<CClass> ClassList;


class CLGAPI CBase
{
public:
    CBase() { ; }
    virtual ~CBase() { ; }
    virtual const CClass* GetClass() const { return NULL; }
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
            printf(_T("==we have class: %s"), CurrentList->GetName());
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