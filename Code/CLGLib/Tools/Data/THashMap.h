//=============================================================================
// FILENAME : THashMap.h
// 
// DESCRIPTION:
// CCMemStack improved hash map
// REVISION:
//  [3/13/2018 nbale]
//=============================================================================
#pragma once

#ifndef _THASHMAP_H_
#define _THASHMAP_H_

__BEGIN_NAMESPACE

enum { kTMapDefaultHashTableSize = 257 }; //17 be better?

/************************************************************************/
/* MemStack                                                                                                                                               */
/************************************************************************/
class CLGAPI CCMapMemStack : public CMemStack
{
    enum {kMapMemStackChunkSize = (4<<10)};
public:
    CCMapMemStack() { Init(kMapMemStackChunkSize); }
    ~CCMapMemStack() { Reset(); }
};

/************************************************************************/
/* Compare Functions and Hash Function                                                                                             */
/************************************************************************/
template<class ARG_KEY>
inline UINT TMapHashKey(ARG_KEY key)
{
    // default identity hash - works for most primitive values
    return (UINT)((DWORD)(key) >> 4 );
}
template<>
inline UINT TMapHashKey<const CCString&>(const CCString& key)
{
    const TCHAR* name = key;
    DWORD hash = 0;
    while (*name)
    {
        const TCHAR Ch = (TCHAR)appToUpper(*name++);
        const BYTE  B = (BYTE)Ch;
        hash = ((hash << 3) + B) ^ (hash);
#ifdef UNICODE
        B = Ch >> 8;
        hash = ((hash << 3) + B) ^ (hash);
#endif
    }
    return hash & (kTMapDefaultHashTableSize - 1);
}

template<class TYPE, class ARG_TYPE>
inline UBOOL TMapCompareElements(const TYPE* pElement1, const ARG_TYPE* pElement2)
{
    return *pElement1 == *pElement2;
}

/************************************************************************/
/* Tmap                                                                                                                                                        */
/************************************************************************/
template<class KEY, class VALUE, class ARG_KEY = const KEY&, class ARG_VALUE = const VALUE&>
class __DLL_EXPORT THashMap
{
public:
    class TPair
    {
    public:
        const KEY    m_Key;
        VALUE    m_Value;
    protected:
        TPair( ARG_KEY keyval ) : m_Key( keyval )    {}
    };
protected:
    class TAssoc : public TPair
    {
    public:
        friend class THashMap<KEY, ARG_KEY, VALUE, ARG_VALUE>;
        TAssoc* m_pNext;
        UINT     m_nHashValue;  // needed for efficient iteration
    public:
        TAssoc( ARG_KEY key ) : TPair( key ) {}
    };

public:

    typedef UINT(*_pfHash)(ARG_KEY key);
    typedef UBOOL(*_pfCompare)(ARG_KEY key1, ARG_KEY key2);

    //----------------------------------------------------------------------------------------------------------
    // Construction and Destruction
    THashMap(_pfHash pfHash = NULL, _pfCompare pfCompare = NULL, INT HashTableSize = kTMapDefaultHashTableSize)
        : m_pHashTable(NULL)
        , m_nHashTableSize(HashTableSize)
        , m_nCount(0)
        , m_pHashFunction(pfHash)
        , m_pCompareFunction(pfCompare)
    {
        ResetHashTable(m_nHashTableSize);
    }
    THashMap(const THashMap& other, _pfHash pfHash = NULL, _pfCompare pfCompare = NULL, INT HashTableSize = kTMapDefaultHashTableSize)
        : m_pHashTable(NULL)
        , m_nHashTableSize(HashTableSize)
        , m_nCount(0)
        , m_pHashFunction(pfHash)
        , m_pCompareFunction(pfCompare)
    {
        ResetHashTable(m_nHashTableSize);

        TArray<KEY> keys = other.GetAllKeys();
        for (INT i = 0; i < keys.Num(); ++i)
        {
            SetAt(keys[i], other.GetAt(keys[i]));
        }
    }
    ~THashMap()
    {
        RemoveAll();
        assert(0 == m_nCount);
    }
    void RemoveAll()
    {
        if (m_pHashTable != NULL)
        {
            // destroy elements (values and keys)
            for (UINT nHash = 0; nHash < m_nHashTableSize; ++nHash)
            {
                TAssoc* pAssoc;
                for (pAssoc = m_pHashTable[nHash]; pAssoc != NULL; pAssoc = pAssoc->m_pNext)
                    pAssoc->TAssoc::~TAssoc();
            }
        }
        // free hash table
        delete[] m_pHashTable;
        m_pHashTable = NULL;
        m_nCount = 0;
        m_MemStack.Reset();
    }
    void ResetHashTable(UINT nHashSize, UBOOL bAllocNow = TRUE)
    {
        if (m_nCount > 0)
            RemoveAll();
        if (m_pHashTable != NULL)
        {
            // free hash table
            delete[] m_pHashTable;
            m_pHashTable = NULL;
        }
        if (bAllocNow)
        {
            m_pHashTable = new TAssoc*[nHashSize];
            memset(m_pHashTable, 0, sizeof(TAssoc*) * nHashSize);
        }
        m_nHashTableSize = nHashSize;
    }

    //-----------------------------------------------------------------------------------------------------
    // Attributes
    // number of elements
    INT GetCount() const { return m_nCount; }
    INT GetSize() const  { return m_nCount; }
    UBOOL IsEmpty() const { return m_nCount == 0; }
    UINT GetHashTableSize() const { return m_nHashTableSize; }

    //------------------------------------------------------------------------------------------------------
    // Lookup
    UBOOL Exist(ARG_KEY key) const
    {
        UINT nHashBucket, nHashValue;
        TAssoc* pAssoc = GetAssocAt(key, nHashBucket, nHashValue);
        if (pAssoc == NULL) return FALSE;  // not in map
        return TRUE;
    }
    UBOOL Lookup(ARG_KEY key, VALUE& rValue) const
    {
        UINT nHashBucket, nHashValue;
        TAssoc* pAssoc = GetAssocAt(key, nHashBucket, nHashValue);
        if (pAssoc == NULL) return FALSE;  // not in map
        rValue = pAssoc->m_Value;
        return TRUE;
    }
    const TAssoc*PLookup(ARG_KEY key) const
    {
        UINT nHashBucket, nHashValue;
        TAssoc* pAssoc = GetAssocAt(key, nHashBucket, nHashValue);
        return pAssoc;
    }
    TPair *PLookup(ARG_KEY key)
    {
        UINT nHashBucket, nHashValue;
        TAssoc* pAssoc = GetAssocAt(key, nHashBucket, nHashValue);
        return pAssoc;
    }

    //------------------------------------------------------------------------------------------------------------
    // Operations
    // Lookup and add if not there
    VALUE& operator[](ARG_KEY key)
    {
        UINT nHashBucket, nHashValue;
        TAssoc* pAssoc = GetAssocAt(key, nHashBucket, nHashValue);
        if (NULL == pAssoc)
        {
            if (m_pHashTable == NULL)
                ResetHashTable(m_nHashTableSize);

            // it doesn't exist, add a new Association
            pAssoc = NewAssoc(key);
            pAssoc->m_nHashValue = nHashValue;
            //'pAssoc->value' is a constructed object, nothing more

            // put into hash table
            pAssoc->m_pNext = m_pHashTable[nHashBucket];
            m_pHashTable[nHashBucket] = pAssoc;
        }
        return pAssoc->m_Value;  // return new reference
    }

    inline void Copy(const THashMap& src)
    {
        if (this != &src)
        {
            RemoveAll();
            m_nHashTableSize = kTMapDefaultHashTableSize;
            ResetHashTable(m_nHashTableSize);

            TArray<KEY> keys = src.GetAllKeys();
            for (INT i = 0; i < keys.Num(); ++i)
            {
                SetAt(keys[i], src.GetAt(keys[i]));
            }
        }
    }
    inline THashMap& operator=(const THashMap& other)
    {
        Copy(other);
        return *this;
    }
    // removing existing (key, ?) pair
    UBOOL RemoveKey(ARG_KEY key)
    {
        if (m_pHashTable == NULL) return FALSE;  // nothing in the table

        UINT nHashValue;
        TAssoc** ppAssocPrev;
        nHashValue = TMapHashKey<ARG_KEY>(key);
        ppAssocPrev = &m_pHashTable[nHashValue % m_nHashTableSize];
        TAssoc* pAssoc;
        for (pAssoc = *ppAssocPrev; pAssoc != NULL; pAssoc = pAssoc->m_pNext)
        {
            if ( (nHashValue == pAssoc->m_nHashValue) && TMapCompareElements(&pAssoc->m_Key, &key) )
            {
                // remove it
                *ppAssocPrev = pAssoc->m_pNext;  // remove from list
                FreeAssoc(pAssoc);
                return TRUE;
            }
            ppAssocPrev = &pAssoc->m_pNext;
        }
        return FALSE;  // not found
    }

    // add a new (key, value) pair
    void SetAt(ARG_KEY key, ARG_VALUE newValue) { (*this)[key] = newValue; }
    // get a value
    VALUE GetAt(ARG_KEY key) const 
    { 
        UINT nHashBucket, nHashValue;
        TAssoc* pAssoc = GetAssocAt(key, nHashBucket, nHashValue);
        if (NULL == pAssoc)
        {
            return VALUE();
        }
        return pAssoc->m_Value;
    }

    //----------------------------------------
    //Interface for loop
    TArray<KEY> GetAllKeys() const
    {
        TArray<KEY> ret;
        if (m_pHashTable != NULL)
        {
            for (UINT nHash = 0; nHash < m_nHashTableSize; ++nHash)
            {
                TAssoc* pAssoc;
                for (pAssoc = m_pHashTable[nHash]; pAssoc != NULL; pAssoc = pAssoc->m_pNext)
                {
                    ret.AddItem(pAssoc->m_Key);
                }
            }
        }
        return ret;
    }
    //-----------------------------------------------------------------------------------------------------------------------------------
    //some interface like sort
    //     virtual const TPair *PGetFirstAssoc() const  { return NULL; }
    //     virtual TPair *PGetFirstAssoc()  { return NULL; }
    //     virtual const TPair *PGetNextAssoc(const TPair *pAssocRet) const  { return NULL; }
    //     virtual TPair *PGetNextAssoc(const TPair *pAssocRet)  { return NULL; }
    //     virtual const TPair *PGetLastAssoc() const  { return NULL; }
    //     virtual TPair *PGetLastAssoc()  { return NULL; }
    //Currently Not Need

protected:

    TAssoc** m_pHashTable;
    UINT    m_nHashTableSize;
    INT        m_nCount;
    CCMapMemStack m_MemStack;

public:

    _pfHash m_pHashFunction;
    _pfCompare m_pCompareFunction;

protected:

    //------------------------------------------------------------------------------------------------------------
    //Implement Helpers
    TAssoc* NewAssoc(ARG_KEY key)
    {
        BYTE* p = m_MemStack.PushBytes((INT)sizeof(TAssoc));
        memset(p, 0, sizeof(TAssoc));
        TAssoc* pAssoc = (TAssoc*)p;
        //pAssoc->TAssoc::TAssoc(key);
        //calling constructor is non-standard, it only works in MSVC
        //change to replacement new
        pAssoc = new (pAssoc) TAssoc(key);
        ++m_nCount;
        return pAssoc;
    }
    void FreeAssoc(TAssoc*pAssoc)
    {
        pAssoc->TAssoc::~TAssoc();
        --m_nCount;
        assert(m_nCount >= 0);  // make sure we don't underflow
        // if no more elements, cleanup completely
        if (m_nCount == 0) RemoveAll();
    }
    TAssoc* GetAssocAt(ARG_KEY key, UINT& outHashBucket, UINT& outHashValue) const
    {
        outHashValue = TMapHashKey<ARG_KEY>(key);
        outHashBucket = outHashValue % m_nHashTableSize;
        if (NULL == m_pHashTable) return NULL;
        // see if it exists
        TAssoc* pAssoc;
        for (pAssoc = m_pHashTable[outHashBucket];  pAssoc != NULL;  pAssoc = pAssoc->m_pNext)
            if (pAssoc->m_nHashValue == outHashValue && TMapCompareElements(&pAssoc->m_Key, &key))
                return pAssoc;
        return NULL;
    }
};



__END_NAMESPACE

#endif //#ifndef _THASHMAP_H_

//=============================================================================
// END OF FILE
//=============================================================================
