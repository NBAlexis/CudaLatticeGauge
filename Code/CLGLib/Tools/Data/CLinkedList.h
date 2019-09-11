//=============================================================================
// FILENAME : CLinkedList.h
// 
// DESCRIPTION:
//  Very useful tools to gather things
//
// REVISION:
//  [3/12/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CLINKEDLIST_H_
#define _CLINKEDLIST_H_

__BEGIN_NAMESPACE

template< class T > class TSimpleDoubleLinkedList : public T
{
public:
    TSimpleDoubleLinkedList* m_pNext;
    TSimpleDoubleLinkedList** m_ppPrevLink;
    void Unlink()
    {
        if( m_pNext )
            m_pNext->m_ppPrevLink = m_ppPrevLink;
        *m_ppPrevLink = m_pNext;
    }
    void Link( TSimpleDoubleLinkedList*& Before )
    {
        if( Before )
            Before->m_ppPrevLink = &m_pNext;
        m_pNext         = Before;
        m_ppPrevLink    = &Before;
        Before  = this;
    }
};

/*----------------------------------------------------------------------------
SimpleList.
//        Sample:
//            1) TSimpleSingleList<MyStruct*>* MyList = NULL;            //create the list head
//            2) MyList = new TSimpleSingleList<MyStruct*>(MyList);    //add a new item, and make it as head
----------------------------------------------------------------------------*/
template <class ElementType> class TSimpleList : public ElementType
{
public:

    ElementType m_Element;
    TSimpleList<ElementType>* m_pNext;

    // Constructor.

    TSimpleList(ElementType InElement,TSimpleList<ElementType>* pInNext = NULL)
    {
        m_Element = InElement;
        m_pNext = pInNext;
    }
};


__END_NAMESPACE

#endif//#ifndef _CLINKEDLIST_H_

//=============================================================================
// END OF FILE
//=============================================================================