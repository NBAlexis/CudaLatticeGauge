//=============================================================================
// FILENAME : TArray.h
// 
// DESCRIPTION:
// This is a dynamic size array, with data in block of memory
//
// REVISION:
//  [3/13/2018 nbale]
//=============================================================================

#ifndef _TARRAY_H_
#define _TARRAY_H_

__BEGIN_NAMESPACE
///////////////////////////////////////////////////////////////////////////
// TArray<TYPE, ARG_TYPE>
//	1)The implementation is similar to MFC's CArray,
//	  See help of CArray in MFC for details
//	2)some STL interface is added, e.g: clear, empty, size, etc
//	3)memory issues
//		* ReserveSpace	: memory is reserved
//		* SetSize		: memory will be allocated or deallocated according to size and growby
//		* RemoveAll		: memory is freed
//		* FreeExtra		: array is shrinked

template<class TYPE, class ARG_TYPE = const TYPE&>
class __DLL_EXPORT TArray
{
public:

    //----------------------------------------------------------------------------------------
    // Construction
    //----------------------------------------------------------------------------------------

    TArray(void)
     {
         m_pData = NULL;
         m_nSize = m_nMaxSize = m_nGrowBy = 0;
     }

    TArray( const TArray& Other ) 
    {
        m_pData = NULL;
        m_nSize = m_nMaxSize = m_nGrowBy = 0;
        *this = Other;
    }

    ~TArray()
    {
        if (m_pData != NULL)
        {
            for( INT i = 0; i < m_nSize; ++i )
                (m_pData + i)->~TYPE();
            delete[] (BYTE*)(m_pData);
        }
    }

    //----------------------------------------------------------------------------------------
    // Accessing
    //----------------------------------------------------------------------------------------

    inline INT GetSize() const
    {
        return m_nSize;
    }

    inline INT GetCount() const
    {
        return m_nSize;
    }

    inline UBOOL IsEmpty() const
    {
        return m_nSize == 0;
    }

    inline INT GetUpperBound() const
    {
        return m_nSize - 1;
    }

    inline INT Num() const
    { 
        return m_nSize; 
    }

    // Accessing elements
    inline const TYPE& GetAt(INT nIndex) const
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        return m_pData[nIndex]; 
    }

    inline TYPE& GetAt(INT nIndex)
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        return m_pData[nIndex]; 
    }

    inline const TYPE& ElementAt(INT nIndex) const
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        return m_pData[nIndex]; 
    }

    inline TYPE& ElementAt(INT nIndex)
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        return m_pData[nIndex]; 
    }

    // overloaded operator helpers
    inline const TYPE& operator[](INT nIndex) const
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        return m_pData[nIndex]; 
    }

    inline TYPE& operator[](INT nIndex)
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        return m_pData[nIndex]; 
    }

    // Direct Access to the element data (may return NULL)
    inline const TYPE* GetData() const
    {
        return (const TYPE*)(m_pData);
    }

    inline TYPE* GetData()
    {
        return (TYPE*)(m_pData);
    }

    inline INT GetCapacity() const
    { 
        return m_nMaxSize; 
    }

    inline INT FindItemIndex(ARG_TYPE Item)
    {
        INT ret = INDEX_NONE;
        for(INT i = 0; i < m_nSize; ++i)
        {
            if(m_pData[i] == Item)
            {
                return i;
            }
        }
        return ret;
    }

    //----------------------------------------------------------------------------------------
    // Operations
    //----------------------------------------------------------------------------------------

    inline void SetAt(INT nIndex, ARG_TYPE newElement)
    {
        assert(nIndex >= 0 && nIndex < m_nSize);
        m_pData[nIndex] = newElement; 
    }

    inline void SetSize(INT nNewSize, INT nGrowBy = -1)
    {
        assert(nNewSize >= 0);

        if (nGrowBy >= 0)
            m_nGrowBy = nGrowBy;  // set new size

        if (nNewSize == 0)
        {
            // shrink to nothing
            if (m_pData != NULL)
            {
                for( INT i = 0; i < m_nSize; ++i )
                    (m_pData + i)->~TYPE();
                delete[] (BYTE*)(m_pData);
                m_pData = NULL;
            }
            m_nSize = m_nMaxSize = 0;
        }
        else if (m_pData == NULL)
        {
            // create buffer big enough to hold number of requested elements or
            // _nGrowBy elements, whichever is larger.
            INT nAllocSize = appMax(nNewSize, m_nGrowBy);
            m_pData = (TYPE*)(new BYTE[(SIZE_T)(nAllocSize) * sizeof(TYPE)]);
            memset( (void*)(m_pData), 0, (SIZE_T)(nAllocSize) * sizeof(TYPE) );
            for( INT i = 0; i < nNewSize; ++i )
            {
                //placement new
                new( (void*)(m_pData + i) ) TYPE;
            }
            m_nSize = nNewSize;
            m_nMaxSize = nAllocSize;
        }
        else if (nNewSize <= m_nMaxSize)
        {
            // it fits
            if (nNewSize > m_nSize)
            {
                // initialize the new elements
                memset((void*)(m_pData + m_nSize), 0, (SIZE_T)(nNewSize-m_nSize) * sizeof(TYPE));
                for( INT i = 0; i < nNewSize - m_nSize; ++i )
                {
                    new ( (void*)(m_pData + m_nSize + i) ) TYPE;
                }
            }
            else if (m_nSize > nNewSize)
            {
                // destroy the old elements
                for( INT i = 0; i < m_nSize - nNewSize; ++i )
                    (m_pData + nNewSize + i)->~TYPE();
            }
            m_nSize = nNewSize;
        }
        else
        {
            // otherwise, grow array
            nGrowBy = m_nGrowBy;
            if (nGrowBy == 0)
            {
                // heuristically determine growth when nGrowBy == 0
                //  (this avoids heap fragmentation in many situations)
                nGrowBy = m_nSize / 8;
                nGrowBy = (nGrowBy < 4) ? 4 : ((nGrowBy > 1024) ? 1024 : nGrowBy);
            }
            INT nNewMax;
            if (nNewSize < m_nMaxSize + nGrowBy)
                nNewMax = m_nMaxSize + nGrowBy;  // granularity
            else
                nNewMax = nNewSize;  // no slush

            assert(nNewMax >= m_nMaxSize);  // no wrap around

            TYPE* pNewData = (TYPE*)(new BYTE[(SIZE_T)(nNewMax) * sizeof(TYPE)]);

            // copy new data from old
            memcpy(pNewData, m_pData, (SIZE_T)(m_nSize) * sizeof(TYPE));

            // construct remaining elements
            assert(nNewSize > m_nSize);
            memset((void*)(pNewData + m_nSize), 0, (SIZE_T)(nNewSize-m_nSize) * sizeof(TYPE));
            for( INT i = 0; i < nNewSize-m_nSize; ++i )
            {
                new ( (void*)(pNewData + m_nSize + i) ) TYPE;
            }

            // get rid of old stuff (note: no destructors called)
            delete[] (BYTE*)(m_pData);
            m_pData = pNewData;
            m_nSize = nNewSize;
            m_nMaxSize = nNewMax;
        }
    }

    // Clean up
    inline void FreeExtra()
    {
        if (m_nSize != m_nMaxSize)
        {
            // shrink to desired size
            TYPE* pNewData = NULL;
            if (m_nSize != 0)
            {
                pNewData = (TYPE*)(new BYTE[m_nSize * sizeof(TYPE)]);
                // copy new data from old
                memcpy(pNewData, m_pData, m_nSize * sizeof(TYPE));
            }

            // get rid of old stuff (note: no destructors called)
            delete[] (BYTE*)(m_pData);
            m_pData = pNewData;
            m_nMaxSize = m_nSize;
        }
    }

    inline void Shrink(){ FreeExtra(); }

    inline void RemoveAll(UBOOL bFreeMemory = TRUE)
    {
        if(bFreeMemory)
            SetSize(0, -1);
        else
        {
            if (m_pData != NULL)
            {
                for( INT i = 0; i < m_nSize; ++i )
                    (m_pData + i)->~TYPE();
            }
            m_nSize = 0;
        }
    }

    inline void Reset() { RemoveAll(FALSE); }

    // Potentially growing the array
    inline void SetAtGrow(INT nIndex, ARG_TYPE newElement)
    {
        assert(nIndex >= 0);
        if (nIndex >= m_nSize)
            SetSize(nIndex + 1, -1);
        m_pData[nIndex] = newElement;
    }

    inline INT AddSize(INT nSizeToAdd, UBOOL bZeroMemory = FALSE)
    {
        INT nOldSize = GetSize();
        INT nNewSize = nSizeToAdd + nOldSize; 
        SetSize(nNewSize); 
        if(bZeroMemory)
        {
            for ( INT i = nOldSize; i < nNewSize; ++i )
            {
                memset(&(m_pData[i]), 0, sizeof(TYPE));
            }
        }
        return nOldSize;
    }

    inline INT AddZeroed(INT n){ return AddSize(n, TRUE); }

    inline INT AddItem(ARG_TYPE newElement)
    {
        INT nIndex = m_nSize;
        SetAtGrow(nIndex, newElement);
        return nIndex; 
    }

    //Only return first found
    inline INT RemoveItem( ARG_TYPE Item )
    {
        INT ret = INDEX_NONE;
        for(INT i = 0; i < m_nSize; ++i)
        {
            if(m_pData[i] == Item)
            {
                RemoveAt(i);
                return i;
            }
        }
        return ret;
    }

    inline INT Append(const TArray& src)
    {
        //assert(this != &src);   // cannot append to itself (but why?)
        INT nOldSize = m_nSize;
        SetSize(m_nSize + src.m_nSize);
        appCopyElements<TYPE>(m_pData + nOldSize, src.m_pData, src.m_nSize);
        return nOldSize;
    }

    inline void Copy(const TArray& src)
    {
        //assert(this != &src);   // cannot append to itself (but why?)
        if(this != &src)
        {
            SetSize(src.m_nSize);
            appCopyElements<TYPE>(m_pData, src.m_pData, src.m_nSize);
        }
    }

    // Operations that move elements around
    inline void InsertAt(INT nIndex, ARG_TYPE newElement, INT nCount = 1)
    {
        assert(nIndex >= 0);    // will expand to meet need
        assert(nCount > 0);     // zero or negative size not allowed

        if (nIndex >= m_nSize)
        {
            // adding after the end of the array
            SetSize(nIndex + nCount, -1);   // grow so nIndex is valid
        }
        else
        {
            // inserting in the middle of the array
            INT nOldSize = m_nSize;
            SetSize(m_nSize + nCount);  // grow it to new size
            // destroy initial data before copying over it
            for( INT i = 0; i < nCount; i++ )
                (m_pData + nOldSize + i)->~TYPE();
            // shift old data up to fill gap
            memmove(m_pData + nIndex + nCount, m_pData + nIndex,
                (nOldSize-nIndex) * sizeof(TYPE));

            // re-init slots we copied from
            memset( (void*)(m_pData + nIndex), 0, (SIZE_T)(nCount) * sizeof(TYPE));
            for( INT i = 0; i < nCount; ++i )
            {
                new ( (void*)(m_pData + nIndex + i ) ) TYPE;
            }
        }

        // insert new value in the gap
        assert(nIndex + nCount <= m_nSize);
        while (nCount--)
            m_pData[nIndex++] = newElement;
    }

    inline void RemoveAt(INT nIndex, INT nCount = 1)
    {
        assert(nIndex >= 0);
        assert(nCount >= 0);
        assert(nIndex + nCount <= m_nSize);

        // just remove a range
        INT nMoveCount = m_nSize - (nIndex + nCount);
        for( INT i = 0; i < nCount; ++i )
            (m_pData + nIndex + i)->~TYPE();
        if (nMoveCount)
            memmove( m_pData + nIndex, m_pData + nIndex + nCount,
            (SIZE_T)(nMoveCount) * sizeof(TYPE));
        m_nSize -= nCount;
    }

    inline void InsertAt(INT nStartIndex, TArray* pNewArray)
    {
        assert(pNewArray != NULL);
        assert(nStartIndex >= 0);

        if (pNewArray->GetSize() > 0)
        {
            InsertAt(nStartIndex, pNewArray->GetAt(0), pNewArray->GetSize());
            for (INT i = 0; i < pNewArray->GetSize(); ++i)
                SetAt(nStartIndex + i, pNewArray->GetAt(i));
        }
    }

    inline void ReserveSpace(INT nCount)
    {
        INT grow_by = m_nGrowBy;
        SetSize(m_nSize, nCount);
        m_nGrowBy = grow_by;
    }

    inline INT PushBack(ARG_TYPE newElement)
    { 
        SetAtGrow(m_nSize, newElement); 
        return m_nSize; 
    }

    inline TYPE Pop(void)
    {
        assert(m_nSize > 0);
        assert(m_nMaxSize >= m_nSize);
        TYPE Result = ((TYPE*)m_pData)[m_nSize - 1];
        RemoveAt( m_nSize - 1 );
        return Result;
    }

    //----------------------------------------------------------------------------------------
    // Other Operators
    //----------------------------------------------------------------------------------------

    inline TArray& operator=( const TArray& Other ) { Copy(Other); return *this; }

protected:

    TYPE*   m_pData;   // the actual array of data
    INT     m_nSize;     // # of elements (upperBound - 1)
    INT     m_nMaxSize;  // max allocated
    INT     m_nGrowBy;   // grow amount

};

__END_NAMESPACE

#endif //#ifndef _TARRAY_H_
//=============================================================================
// END OF FILE
//=============================================================================
