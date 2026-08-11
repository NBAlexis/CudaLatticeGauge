//=============================================================================
// FILENAME : VectorsN2.h
// 
// DESCRIPTION:
// This is helper for calculate all kinds of vectors
//
//
// REVISION:
//  [mm/dd/yy]
//  [05/06/2025 nbale]
//=============================================================================

#ifndef _VECTORSN2_H_
#define _VECTORSN2_H_

__BEGIN_NAMESPACE

template<INT M, INT N, INT NofE, INT mn = M * N>
struct deviceVectorN2
{
public:
    __device__ deviceVectorN2()
    {

    }

    __device__ deviceVectorN2(const deviceVectorN2<M, N, NofE>& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = other.m_ve[i];
        }
    }

    __device__ __inline__ void DebugPrint(const char* header = NULL) const
    {
        if (NULL == header)
        {
            printf("{");
        }
        else
        {
            printf("%s={", header);
        }
        for (INT i = 0; i < M; ++i)
        {
            printf("{");
            for (INT j = 0; j < N; ++j)
            {
                printf("%2.3f %s %2.3f I",
                    m_ve[i * N + j].x,
                    m_ve[i * N + j].y < 0 ? "" : "+",
                    m_ve[i * N + j].y
                );
                if (j < N - 1)
                {
                    printf(", ");
                }
                else
                {
                    printf("}");
                }
            }
            if (i < M - 1)
            {
                printf(",\n");
            }
            else
            {
                printf("}");
            }
        }
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeRandom(UINT fatIndex)
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            ret.m_ve[i] = _deviceRandomC(fatIndex);
        }
        return ret;
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeRandomGaussian(UINT fatIndex)
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            ret.m_ve[i] = _deviceRandomGaussC(fatIndex);
        }
        return ret;
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeRandomZ4(UINT fatIndex)
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            ret.m_ve[i] = _deviceRandomZ4(fatIndex);
        }
        return ret;
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeZeroVectorN2()
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            ret.m_ve[i] = _make_cuComplex(F(0.0), F(0.0));
        }
        return ret;
    }

    __device__ __inline__ void Zero()
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _make_cuComplex(F(0.0), F(0.0));
        }
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeOneSUNVector()
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            ret.m_ve[i] = _make_cuComplex(F(1.0), F(0.0));
        }
        return ret;
    }

    __device__ __inline__ void Id()
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _make_cuComplex(F(1.0), F(0.0));
        }
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeOneVectorColor(BYTE byColor)
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < M; ++i)
        {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
            for (INT j = 0; j < N; ++j)
            {
                ret.m_ve[i * N + j] = _make_cuComplex(j == byColor ? F(1.0) : F(0.0), F(0.0));
            }
        }
        return ret;
    }

    __device__ __inline__ static deviceVectorN2<M, N, NofE> makeOneVectorColor2(BYTE byColor2)
    {
        deviceVectorN2<M, N, NofE> ret;
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < M; ++i)
        {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
            for (INT j = 0; j < N; ++j)
            {
                ret.m_ve[i * N + j] = _make_cuComplex(i == byColor2 ? F(1.0) : F(0.0), F(0.0));
            }
        }
        return ret;
    }

    /**
    * This is sum _3 (v^* v)
    */
    __device__ __inline__ CLGComplex ConjugateDotC(const deviceVectorN2<M, N, NofE>& other) const
    {
        CLGComplex ret = _cuCmulf(_cuConjf(m_ve[0]), other.m_ve[0]);
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 1; i < mn; ++i)
        {
            ret = _cuCaddf(ret, _cuCmulf(_cuConjf(m_ve[i]), other.m_ve[i]));
        }

        return ret;
    }

    __device__ __inline__ void Sub(const deviceVectorN2<M, N, NofE>& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCsubf(m_ve[i], other.m_ve[i]);
        }
    }

    __device__ __inline__ void SubComp(const CLGComplex& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCsubf(m_ve[i], other);
        }
    }

    __device__ __inline__ void SubReal(Real other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = cuCsubf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void Add(const deviceVectorN2<M, N, NofE>& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCaddf(m_ve[i], other.m_ve[i]);
        }
    }

    __device__ __inline__ void AddComp(const CLGComplex& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCaddf(m_ve[i], other);
        }
    }

    __device__ __inline__ void AddReal(Real other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = cuCaddf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void MulReal(Real other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = cuCmulf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void MulComp(const CLGComplex& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCmulf(m_ve[i], other);
        }
    }

    __device__ __inline__ void DivReal(Real other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = cuCdivf_cr(m_ve[i], other);
        }
    }

    __device__ __inline__ void DivComp(const CLGComplex& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCdivf(m_ve[i], other);
        }
    }

    __device__ __inline__ void Mul(const deviceVectorN2<M, N, NofE>& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCmulf(m_ve[i], other.m_ve[i]);
        }
    }

    __device__ __inline__ deviceVectorN2<M, N, NofE> MulC(const deviceVectorN2<M, N, NofE>& other) const
    {
        deviceVectorN2<M, N, NofE> ret(*this);
        ret.Mul(other);
        return ret;
    }

    __device__ __inline__ void DaggerMul(const deviceVectorN2<M, N, NofE>& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCmulf(_cuConjf(m_ve[i]), other.m_ve[i]);
        }
    }

    __device__ __inline__ deviceVectorN2<M, N, NofE> DaggerMulC(const deviceVectorN2<M, N, NofE>& other) const
    {
        deviceVectorN2<M, N, NofE> ret(*this);
        ret.DaggerMul(other);
        return ret;
    }

    __device__ __inline__ void MulDagger(const deviceVectorN2<M, N, NofE>& other)
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i] = _cuCmulf(m_ve[i], _cuConjf(other.m_ve[i]));
        }
    }

    __device__ __inline__ deviceVectorN2<M, N, NofE> MulDaggerC(const deviceVectorN2<M, N, NofE>& other) const
    {
        deviceVectorN2<M, N, NofE> ret(*this);
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
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
            for (INT i = 0; i < mn; ++i)
            {
                m_ve[i] = _make_cuComplex(-m_ve[i].y, m_ve[i].x);
            }
        }
        break;
        case 2:
        {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
            for (INT i = 0; i < mn; ++i)
            {
                m_ve[i] = _make_cuComplex(-m_ve[i].x, -m_ve[i].y);
            }
        }
        break;
        case 3:
        {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
            for (INT i = 0; i < mn; ++i)
            {
                m_ve[i] = _make_cuComplex(m_ve[i].y, -m_ve[i].x);
            }
        }
        break;
        }
    }

    __device__ __inline__ void Opposite()
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i].x = -m_ve[i].x;
            m_ve[i].y = -m_ve[i].y;
        }
    }

    __device__ __inline__ void Conjugate()
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
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

    __device__ __inline__ void Re()
    {
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 0; i < mn; ++i)
        {
            m_ve[i].y = F(0.0);
        }
    }

    __device__ __inline__ deviceVectorN2<M, N, NofE> SubRealC(Real other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.SubReal(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> SubCompC(const CLGComplex& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.SubComp(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> SubC(const deviceVectorN2<M, N, NofE>& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.Sub(other); return ret; }

    __device__ __inline__ deviceVectorN2<M, N, NofE> AddRealC(const Real& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.AddReal(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> AddCompC(const CLGComplex& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.AddComp(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> AddC(const deviceVectorN2<M, N, NofE>& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.Add(other); return ret; }

    __device__ __inline__ deviceVectorN2<M, N, NofE> MulRealC(Real other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.MulReal(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> MulCompC(const CLGComplex& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.MulComp(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> DivRealC(Real other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.DivReal(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> DivCompC(const CLGComplex& other) const { deviceVectorN2<M, N, NofE> ret(*this); ret.DivComp(other); return ret; }
    __device__ __inline__ deviceVectorN2<M, N, NofE> MulZ4C(BYTE z4) const { deviceVectorN2<M, N, NofE> ret(*this); ret.MulZ4(z4); return ret; }

    __device__ __inline__ CLGComplex Sum() const
    {
        if (1 == mn)
        {
            return m_ve[0];
        }

        CLGComplex res = _cuCaddf(m_ve[0], m_ve[1]);
#if defined(__CUDA_ARCH__)
#pragma unroll
#endif
        for (INT i = 2; i < mn; ++i)
        {
            res = _cuCaddf(res, m_ve[i]);
        }
        return res;
    }

    __device__ __inline__ void Norm()
    {
        Real len = ConjugateDotC(*this).x;
        if (len > _CLG_FLT_MIN_)
        {
            len = __rcp(_sqrt(len));
            MulReal(len);
        }
    }

    __device__ __inline__ Real Abs() const
    {
        Real len = ConjugateDotC(*this).x;
        if (len > _CLG_FLT_MIN_)
        {
            return _sqrt(len);
        }
        return F(0.0);
    }

    CLGComplex m_ve[NofE];
};

__END_NAMESPACE

#endif //#ifndef _VECTORSN2_H_

//=============================================================================
// END OF FILE
//=============================================================================
