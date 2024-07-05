//=============================================================================
// FILENAME : VectorsN.h
// 
// DESCRIPTION:
// This is helper for calculate all kinds of vectors
//
//
// REVISION:
//  [07/01/2024 nbale]
//=============================================================================

#ifndef _VECTORSN_H_
#define _VECTORSN_H_

__BEGIN_NAMESPACE

template<INT N, INT NofE>
struct deviceSUNVector
{
public:
    __device__ deviceSUNVector()
    {

    }

    __device__ deviceSUNVector(const deviceSUNVector<N, NofE>& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = other.m_ve[i];
        }
    }

    __device__ __inline__ void DebugPrint() const
    {
        printf("{");
        for (INT i = 0; i < N; ++i)
        {
            printf("%2.3f %s %2.3f i\n",
                m_ve[i].x,
                m_ve[i].y < 0 ? "" : "+",
                m_ve[i].y
            );
            if (i < N - 1)
            {
                printf(", ");
            }
        }
        printf("}\n");
    }

    __device__ __inline__ static deviceSUNVector<N, NofE> makeRandomGaussian(UINT fatIndex)
    {
        deviceSUNVector<N, NofE> ret;
        for (INT i = 0; i < N; ++i)
        {
            ret.m_ve[i] = _deviceRandomGaussC(fatIndex);
        }
        return ret;
    }

    __device__ __inline__ static deviceSUNVector<N, NofE> makeRandomZ4(UINT fatIndex)
    {
        deviceSUNVector<N, NofE> ret;
        for (INT i = 0; i < N; ++i)
        {
            ret.m_ve[i] = _deviceRandomZ4(fatIndex);
        }
        return ret;
    }

    __device__ __inline__ static deviceSUNVector<N, NofE> makeZeroSUNVector()
    {
        deviceSUNVector<N, NofE> ret;
        for (INT i = 0; i < N; ++i)
        {
            ret.m_ve[i] = _make_cuComplex(F(0.0), F(0.0));
        }
        return ret;
    }

    __device__ __inline__ static deviceSUNVector<N, NofE> makeOneSUNVector()
    {
        deviceSUNVector<N, NofE> ret;
        for (INT i = 0; i < N; ++i)
        {
            ret.m_ve[i] = _make_cuComplex(F(1.0), F(0.0));
        }
        return ret;
    }

    __device__ __inline__ static deviceSUNVector<N, NofE> makeOneSUNVectorColor(BYTE byColor)
    {
        deviceSUNVector<N, NofE> ret;
        for (INT i = 0; i < N; ++i)
        {
            ret.m_ve[i] = _make_cuComplex(i == byColor ? F(1.0) : F(0.0), F(0.0));
        }
        return ret;
    }

    /**
    * This is sum _3 (v^* v)
    */
    __device__ __inline__ CLGComplex ConjugateDotC(const deviceSUNVector<N, NofE>& other) const
    {
        CLGComplex ret = _cuCmulf(_cuConjf(m_ve[0]), other.m_ve[0]);
        for (INT i = 1; i < N; ++i)
        {
            ret = _cuCaddf(ret, _cuCmulf(_cuConjf(m_ve[i]), other.m_ve[i]));
        }
        
        return ret;
    }

    __device__ __inline__ void Sub(const deviceSUNVector<N, NofE>& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCsubf(m_ve[i], other.m_ve[i]);
        }
    }

    __device__ __inline__ void SubComp(const CLGComplex& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCsubf(m_ve[i], other);
        }
    }

    __device__ __inline__ void SubReal(Real other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = cuCsubf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void Add(const deviceSUNVector<N, NofE>& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCaddf(m_ve[i], other.m_ve[i]);
        }
    }

    __device__ __inline__ void AddComp(const CLGComplex& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCaddf(m_ve[i], other);
        }
    }

    __device__ __inline__ void AddReal(Real other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = cuCaddf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void MulReal(Real other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = cuCmulf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void MulComp(const CLGComplex& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCmulf(m_ve[i], other);
        }
    }

    __device__ __inline__ void Mul(const deviceSUNVector<N, NofE>& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCmulf(m_ve[i], other.m_ve[i]);
        }
    }

    __device__ __inline__ deviceSUNVector<N, NofE> MulC(const deviceSUNVector<N, NofE>& other) const
    {
        deviceSUNVector<N, NofE> ret(*this);
        ret.Mul(other);
        return ret;
    }

    __device__ __inline__ void DaggerMul(const deviceSUNVector<N, NofE>& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCmulf(_cuConjf(m_ve[i]), other.m_ve[i]);
        }
    }

    __device__ __inline__ deviceSUNVector<N, NofE> DaggerMulC(const deviceSUNVector<N, NofE>& other) const
    {
        deviceSUNVector<N, NofE> ret(*this);
        ret.DaggerMul(other);
        return ret;
    }

    __device__ __inline__ void MulDagger(const deviceSUNVector<N, NofE>& other)
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i] = _cuCmulf(m_ve[i], _cuConjf(other.m_ve[i]));
        }
    }

    __device__ __inline__ deviceSUNVector<N, NofE> MulDaggerC(const deviceSUNVector<N, NofE>& other) const
    {
        deviceSUNVector<N, NofE> ret(*this);
        ret.MulDagger(other);
        return ret;
    }

    /**
    * v = i^k v
    */
    __device__ __inline__ void MulZ4(BYTE byZ4)
    {
        switch (byZ4)
        {
        case 1:
        {
            for (INT i = 0; i < N; ++i)
            {
                m_ve[i] = _make_cuComplex(-m_ve[i].y, m_ve[i].x);
            }
        }
        break;
        case 2:
        {
            for (INT i = 0; i < N; ++i)
            {
                m_ve[i] = _make_cuComplex(-m_ve[i].x, -m_ve[i].y);
            }
        }
        break;
        case 3:
        {
            for (INT i = 0; i < N; ++i)
            {
                m_ve[i] = _make_cuComplex(m_ve[i].y, -m_ve[i].x);
            }
        }
        break;
        }
    }

    __device__ __inline__ void Opposite()
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i].x = -m_ve[i].x;
            m_ve[i].y = -m_ve[i].y;
        }
    }

    __device__ __inline__ void Conjugate()
    {
        for (INT i = 0; i < N; ++i)
        {
            m_ve[i].y = -m_ve[i].y;
        }
    }

    __device__ __inline__ deviceSUNVector<N, NofE> ConjugateC() const
    {
        deviceSUNVector<N, NofE> ret(*this);
        ret.Conjugate();
        return ret;
    }

    __device__ __inline__ deviceSUNVector<N, NofE> SubRealC(Real other) const { deviceSUNVector<N, NofE> ret(*this); ret.SubReal(other); return ret; }
    __device__ __inline__ deviceSUNVector<N, NofE> SubCompC(const CLGComplex& other) const { deviceSUNVector<N, NofE> ret(*this); ret.SubComp(other); return ret; }
    __device__ __inline__ deviceSUNVector<N, NofE> SubC(const deviceSUNVector<N, NofE>& other) const { deviceSUNVector<N, NofE> ret(*this); ret.Sub(other); return ret; }

    __device__ __inline__ deviceSUNVector<N, NofE> AddRealC(const Real& other) const { deviceSUNVector<N, NofE> ret(*this); ret.AddReal(other); return ret; }
    __device__ __inline__ deviceSUNVector<N, NofE> AddCompC(const CLGComplex& other) const { deviceSUNVector<N, NofE> ret(*this); ret.AddComp(other); return ret; }
    __device__ __inline__ deviceSUNVector<N, NofE> AddC(const deviceSUNVector<N, NofE>& other) const { deviceSUNVector<N, NofE> ret(*this); ret.Add(other); return ret; }

    __device__ __inline__ deviceSUNVector<N, NofE> MulRealC(Real other) const { deviceSUNVector<N, NofE> ret(*this); ret.MulReal(other); return ret; }
    __device__ __inline__ deviceSUNVector<N, NofE> MulCompC(const CLGComplex& other) const { deviceSUNVector<N, NofE> ret(*this); ret.MulComp(other); return ret; }
    __device__ __inline__ deviceSUNVector<N, NofE> MulZ4C(BYTE z4) const { deviceSUNVector<N, NofE> ret(*this); ret.MulZ4(z4); return ret; }

    __device__ __inline__ CLGComplex Sum() const
    {
        if (1 == N)
        {
            return m_ve[0];
        }

        CLGComplex res = _cuCaddf(m_ve[0], m_ve[1]);
        for (INT i = 2; i < N; ++i)
        {
            res = _cuCaddf(res, m_ve[i]);
        }
        return res;
    }

    CLGComplex m_ve[NofE];
};

typedef deviceSUNVector<4, 4> deviceSU4Vector;
typedef deviceSUNVector<5, 8> deviceSU5Vector;
typedef deviceSUNVector<6, 8> deviceSU6Vector;
typedef deviceSUNVector<7, 8> deviceSU7Vector;
typedef deviceSUNVector<8, 8> deviceSU8Vector;

__END_NAMESPACE

#endif //#ifndef _VECTORSN_H_

//=============================================================================
// END OF FILE
//=============================================================================
