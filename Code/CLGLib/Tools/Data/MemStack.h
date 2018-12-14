//=============================================================================
// FILENAME : MemStack.h
// 
// DESCRIPTION:
// This is a pool, maybe use boost::pool instead?
//
// REVISION:
//  [3/13/2018 nbale]
//=============================================================================
#ifndef _MEMSTACK_H_
#define _MEMSTACK_H_

__BEGIN_NAMESPACE

enum { kDefaultMemoryAlignment = sizeof(void*)};

//============================================================
//	CCMemStack: a simple memory pool, generally use with CMemMark
//		1) basically, you should derive a class, and specify 
//		   the chunk size (e.g 64k)
//			class CMyMemStack : public CCMemStack
//			{
//			public:
//				CMyMemStack() { Init(YourMemChunkSize, AlignmentSize); }
//				~CMyMemStack() { Reset(); }
//			};
//		2) and then, use PushBytes to allocate some memory, the 
//		   memory is aligned with the size that your specified,
//		   if it is -1, use the value specified while calling Init()
//		3) generally, you use CMemMark to free the allocated
//		   chunks
//============================================================
class CLGAPI CMemStack
{
protected:
    // Types.
    struct CTaggedMemory
    {
        CTaggedMemory* m_pNext;
        INT m_iDataSize;
        BYTE m_iData[1];
    };

public:


    // Get bytes.

    //inline function
    BYTE* PushBytes( INT AllocSize, INT Align = -1)
    {
        if(Align < 0) 
            Align = m_iDefaultAlignment;

        // Debug checks.
        assert(AllocSize>=0);
        assert((Align&(Align-1))==0);
        assert(m_pTop<=m_pEnd);

        // Try to get memory from the current chunk.
        BYTE* pResult = (BYTE *)(((PTRINT)(m_pTop)+(Align-1))&~(Align-1));
        m_pTop = pResult + AllocSize;

        // Make sure we didn't overflow.
        if( m_pTop > m_pEnd )
        {
            // We'd pass the end of the current chunk, so allocate a new one.
            AllocateNewChunk( AllocSize + Align );
            pResult = (BYTE *)(((PTRINT)(m_pTop)+(Align-1))&~(Align-1));
            m_pTop    = pResult + AllocSize;
        }
        return pResult;
    }

    // Main functions.
    void Init( INT DefaultChunkSize, INT DefaultAlignment=kDefaultMemoryAlignment)
    {
        m_iDefaultChunkSize = DefaultChunkSize;
        m_iDefaultAlignment = DefaultAlignment;
        m_pTopChunk = NULL;
        m_pEnd = NULL;
        m_pTop = NULL;
    }

    void Reset() { FreeChunks(NULL); }

    INT GetByteCount()
    {
        SIZE_T Count = 0;
        for (CTaggedMemory* pChunk = m_pTopChunk; pChunk; pChunk = pChunk->m_pNext)
        {
            if (pChunk != m_pTopChunk)
                Count += pChunk->m_iDataSize;
            else
                Count += m_pTop - pChunk->m_iData;
        }
        return (INT)(Count);
    }

    // Friends.
    friend class CMemMark;

protected:
    // Variables.
    BYTE*			m_pTop;				// Top of current chunk (Top<=End).
    BYTE*			m_pEnd;				// End of current chunk.
    INT				m_iDefaultChunkSize;	// Maximum chunk size to allocate.
    INT				m_iDefaultAlignment;  // Memory alignment for each allocation
    CTaggedMemory*	m_pTopChunk;			// Only chunks 0..ActiveChunks-1 are valid.

    // Functions.
    BYTE* AllocateNewChunk( INT MinSize )
    {
        CTaggedMemory* pChunk = NULL;
        INT DataSize = appMax(MinSize, m_iDefaultChunkSize - (INT)sizeof(CTaggedMemory));
        pChunk = (CTaggedMemory*)(new BYTE[DataSize + sizeof(CTaggedMemory)]);
        pChunk->m_iDataSize = DataSize;
        pChunk->m_pNext = m_pTopChunk;
        m_pTopChunk = pChunk;
        m_pTop = pChunk->m_iData;
        m_pEnd = m_pTop + pChunk->m_iDataSize;
        return m_pTop;
    }

    void FreeChunks( CTaggedMemory* NewTopChunk )
    {
        while (m_pTopChunk != NewTopChunk)
        {
            CTaggedMemory* RemoveChunk = m_pTopChunk;
            m_pTopChunk = m_pTopChunk->m_pNext;
            delete RemoveChunk;
        }
        m_pTop = NULL;
        m_pEnd = NULL;
        if (m_pTopChunk)
        {
            m_pTop = m_pTopChunk->m_iData;
            m_pEnd = m_pTop + m_pTopChunk->m_iDataSize;
        }
    }
};

//============================================================
//	CMemMark: a simple memory pool, generally use with FMemStack
//			  it pops the memory allocated after FMemMark is created
//		1) FMyMemStack MyStack
//		2) FMemMark MyMark(MyStack)
//		3) MyStack.PushBytes(..)
//		4) MyMark.Pop()
//============================================================
class CLGAPI CMemMark
{
public:
    // Constructors.
    CMemMark(CMemStack& InMem )
    {
        m_pMem          = &InMem;
        m_pTop          = m_pMem->m_pTop;
        m_pSavedChunk   = m_pMem->m_pTopChunk;
    }

    // CMemMark interface.
    //inline
    void Pop()
    {
        // Unlock any new chunks that were allocated.
        if( m_pSavedChunk != m_pMem->m_pTopChunk )
            m_pMem->FreeChunks( m_pSavedChunk );

        // Restore the memory stack's state.
        m_pMem->m_pTop = m_pTop;
    }

protected:
    // Implementation variables.
    //Should not use this construction
    CMemMark(){}

    CMemStack* m_pMem;
    BYTE* m_pTop;
    CMemStack::CTaggedMemory* m_pSavedChunk;
};

__END_NAMESPACE

#endif//#ifndef _MEMSTACK_H_

//=============================================================================
// END OF FILE
//=============================================================================