//=============================================================================
// FILENAME : Vectors.h
// 
// DESCRIPTION:
// This is helper for calculate all kinds of vectors
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================

#ifndef _VECTORS_H_
#define _VECTORS_H_

#if _CLG_DOUBLEFLOAT
#define __SU3VECTOR_ALIGN 48
#define __WILSONVECTOR_ALIGN 256
#else
#define __SU3VECTOR_ALIGN 24
#define __WILSONVECTOR_ALIGN 128
#endif

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    //NOTE!!! DO NOT ALIGN ME!!! THIS WILL BREAK UNION!!!
    struct deviceSU3Vector
    {
    public:
        __device__ deviceSU3Vector()
        {

        }

        __device__ deviceSU3Vector(const deviceSU3Vector& other)
        {
            m_ve[0] = other.m_ve[0];
            m_ve[1] = other.m_ve[1];
            m_ve[2] = other.m_ve[2];
        }

        __device__ __inline__ void DebugPrint() const
        {
            printf("%2.3f %s %2.3f i, %2.3f %s %2.3f i, %2.3f %s %2.3f i\n",
                m_ve[0].x,
                m_ve[0].y < 0 ? "" : "+",
                m_ve[0].y,

                m_ve[1].x,
                m_ve[1].y < 0 ? "" : "+",
                m_ve[1].y,

                m_ve[2].x,
                m_ve[2].y < 0 ? "" : "+",
                m_ve[2].y
            );
        }

        __device__ __inline__ static deviceSU3Vector makeRandomGaussian(UINT fatIndex)
        {
            deviceSU3Vector ret;
            ret.m_ve[0] = _deviceRandomGaussC(fatIndex);
            ret.m_ve[1] = _deviceRandomGaussC(fatIndex);
            ret.m_ve[2] = _deviceRandomGaussC(fatIndex);
            return ret;
        }

        __device__ __inline__ static deviceSU3Vector makeZeroSU3Vector()
        {
            deviceSU3Vector ret;
            ret.m_ve[0] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_ve[1] = _make_cuComplex(F(0.0), F(0.0));
            ret.m_ve[2] = _make_cuComplex(F(0.0), F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU3Vector makeOneSU3Vector()
        {
            deviceSU3Vector ret;
            ret.m_ve[0] = _make_cuComplex(F(1.0), F(0.0));
            ret.m_ve[1] = _make_cuComplex(F(1.0), F(0.0));
            ret.m_ve[2] = _make_cuComplex(F(1.0), F(0.0));
            return ret;
        }

        __device__ __inline__ static deviceSU3Vector makeOneSU3VectorColor(BYTE byColor)
        {
            deviceSU3Vector ret;
            ret.m_ve[0] = _make_cuComplex(0 == byColor ? F(1.0) : F(0.0), F(0.0));
            ret.m_ve[1] = _make_cuComplex(1 == byColor ? F(1.0) : F(0.0), F(0.0));
            ret.m_ve[2] = _make_cuComplex(2 == byColor ? F(1.0) : F(0.0), F(0.0));
            return ret;
        }

        /**
        * This is sum _3 (v^* v)
        */
        __device__ __inline__ _Complex ConjugateDotC(const deviceSU3Vector& other) const
        {
            _Complex ret = _cuCmulf(_cuConjf(m_ve[0]), other.m_ve[0]);
            ret = _cuCaddf(ret, _cuCmulf(_cuConjf(m_ve[1]), other.m_ve[1]));
            ret = _cuCaddf(ret, _cuCmulf(_cuConjf(m_ve[2]), other.m_ve[2]));
            return ret;
        }

        __device__ __inline__ void Sub(const deviceSU3Vector& other)
        {
            m_ve[0] = _cuCsubf(m_ve[0], other.m_ve[0]);
            m_ve[1] = _cuCsubf(m_ve[1], other.m_ve[1]);
            m_ve[2] = _cuCsubf(m_ve[2], other.m_ve[2]);
        }

        __device__ __inline__ void SubComp(const _Complex& other)
        {
            m_ve[0] = _cuCsubf(m_ve[0], other);
            m_ve[1] = _cuCsubf(m_ve[1], other);
            m_ve[2] = _cuCsubf(m_ve[2], other);
        }

        __device__ __inline__ void SubReal(Real other)
        {
            m_ve[0] = cuCsubf_cr(m_ve[0], other);
            m_ve[1] = cuCsubf_cr(m_ve[1], other);
            m_ve[2] = cuCsubf_cr(m_ve[2], other);
        }

        __device__ __inline__ void Add(const deviceSU3Vector& other)
        {
            m_ve[0] = _cuCaddf(m_ve[0], other.m_ve[0]);
            m_ve[1] = _cuCaddf(m_ve[1], other.m_ve[1]);
            m_ve[2] = _cuCaddf(m_ve[2], other.m_ve[2]);
        }

        __device__ __inline__ void AddComp(const _Complex& other)
        {
            m_ve[0] = _cuCaddf(m_ve[0], other);
            m_ve[1] = _cuCaddf(m_ve[1], other);
            m_ve[2] = _cuCaddf(m_ve[2], other);
        }

        __device__ __inline__ void AddReal(Real other)
        {
            m_ve[0] = cuCaddf_cr(m_ve[0], other);
            m_ve[1] = cuCaddf_cr(m_ve[1], other);
            m_ve[2] = cuCaddf_cr(m_ve[2], other);
        }

        __device__ __inline__ void MulReal(Real other)
        {
            m_ve[0] = cuCmulf_cr(m_ve[0], other);
            m_ve[1] = cuCmulf_cr(m_ve[1], other);
            m_ve[2] = cuCmulf_cr(m_ve[2], other);
        }

        __device__ __inline__ void MulComp(const _Complex& other)
        {
            m_ve[0] = _cuCmulf(m_ve[0], other);
            m_ve[1] = _cuCmulf(m_ve[1], other);
            m_ve[2] = _cuCmulf(m_ve[2], other);
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
                m_ve[0] = _make_cuComplex(-m_ve[0].y, m_ve[0].x);
                m_ve[1] = _make_cuComplex(-m_ve[1].y, m_ve[1].x);
                m_ve[2] = _make_cuComplex(-m_ve[2].y, m_ve[2].x);
            }
            break;
            case 2:
            {
                m_ve[0] = _make_cuComplex(-m_ve[0].x, -m_ve[0].y);
                m_ve[1] = _make_cuComplex(-m_ve[1].x, -m_ve[1].y);
                m_ve[2] = _make_cuComplex(-m_ve[2].x, -m_ve[2].y);
            }
            break;
            case 3:
            {
                m_ve[0] = _make_cuComplex(m_ve[0].y, -m_ve[0].x);
                m_ve[1] = _make_cuComplex(m_ve[1].y, -m_ve[1].x);
                m_ve[2] = _make_cuComplex(m_ve[2].y, -m_ve[2].x);
            }
            break;
            }
        }

        __device__ __inline__ void Opposite()
        {
            m_ve[0].x = -m_ve[0].x;
            m_ve[0].y = -m_ve[0].y;
            m_ve[1].x = -m_ve[1].x;
            m_ve[1].y = -m_ve[1].y;
            m_ve[2].x = -m_ve[2].x;
            m_ve[2].y = -m_ve[2].y;
        }

        __device__ __inline__ deviceSU3Vector SubRealC(Real other) const { deviceSU3Vector ret(*this); ret.SubReal(other); return ret; }
        __device__ __inline__ deviceSU3Vector SubCompC(const _Complex& other) const { deviceSU3Vector ret(*this); ret.SubComp(other); return ret; }
        __device__ __inline__ deviceSU3Vector SubC(const deviceSU3Vector& other) const { deviceSU3Vector ret(*this); ret.Sub(other); return ret; }

        __device__ __inline__ deviceSU3Vector AddRealC(const Real& other) const { deviceSU3Vector ret(*this); ret.AddReal(other); return ret; }
        __device__ __inline__ deviceSU3Vector AddCompC(const _Complex& other) const { deviceSU3Vector ret(*this); ret.AddComp(other); return ret; }
        __device__ __inline__ deviceSU3Vector AddC(const deviceSU3Vector& other) const { deviceSU3Vector ret(*this); ret.Add(other); return ret; }

        __device__ __inline__ deviceSU3Vector MulRealC(Real other) const { deviceSU3Vector ret(*this); ret.MulReal(other); return ret; }
        __device__ __inline__ deviceSU3Vector MulCompC(const _Complex& other) const { deviceSU3Vector ret(*this); ret.MulComp(other); return ret; }
        __device__ __inline__ deviceSU3Vector MulZ4C(BYTE z4) const { deviceSU3Vector ret(*this); ret.MulZ4(z4); return ret; }

        _Complex m_ve[3];
    };

    struct alignas(__WILSONVECTOR_ALIGN) deviceWilsonVectorSU3
    {
    public:
        __device__ deviceWilsonVectorSU3() { ; }

        __device__ deviceWilsonVectorSU3(const deviceWilsonVectorSU3& other)
        {
            m_d[0] = deviceSU3Vector(other.m_d[0]);
            m_d[1] = deviceSU3Vector(other.m_d[1]);
            m_d[2] = deviceSU3Vector(other.m_d[2]);
            m_d[3] = deviceSU3Vector(other.m_d[3]);
        }

        __device__ __inline__ void DebugPrint() const
        {
            printf("=d0: %2.2f %2.2fi, %2.2f %2.2fi, %2.2f %2.2fi\n d1: %2.2f %2.2fi, %2.2f %2.2fi, %2.2f %2.2fi\n d2: %2.2f %2.2fi, %2.2f %2.2fi, %2.2f %2.2fi\n d3: %2.2f %2.2fi, %2.2f %2.2fi, %2.2f %2.2fi\n",

                m_d[0].m_ve[0].x,
                m_d[0].m_ve[0].y,

                m_d[0].m_ve[1].x,
                m_d[0].m_ve[1].y,

                m_d[0].m_ve[2].x,
                m_d[0].m_ve[2].y,

                m_d[1].m_ve[0].x,
                m_d[1].m_ve[0].y,

                m_d[1].m_ve[1].x,
                m_d[1].m_ve[1].y,

                m_d[1].m_ve[2].x,
                m_d[1].m_ve[2].y,

                m_d[2].m_ve[0].x,
                m_d[2].m_ve[0].y,

                m_d[2].m_ve[1].x,
                m_d[2].m_ve[1].y,

                m_d[2].m_ve[2].x,
                m_d[2].m_ve[2].y,

                m_d[3].m_ve[0].x,
                m_d[3].m_ve[0].y,

                m_d[3].m_ve[1].x,
                m_d[3].m_ve[1].y,

                m_d[3].m_ve[2].x,
                m_d[3].m_ve[2].y
            );
        }

        __device__ __inline__ void Opposite()
        {
            m_d[0].Opposite();
            m_d[1].Opposite();
            m_d[2].Opposite();
            m_d[3].Opposite();
        }

        __device__ __inline__ static deviceWilsonVectorSU3 makeRandomGaussian(UINT fatIndex)
        {
            deviceWilsonVectorSU3 ret;
            ret.m_d[0] = deviceSU3Vector::makeRandomGaussian(fatIndex);
            ret.m_d[1] = deviceSU3Vector::makeRandomGaussian(fatIndex);
            ret.m_d[2] = deviceSU3Vector::makeRandomGaussian(fatIndex);
            ret.m_d[3] = deviceSU3Vector::makeRandomGaussian(fatIndex);
            return ret;
        }

        __device__ __inline__ static deviceWilsonVectorSU3 makeZeroWilsonVectorSU3()
        {
            deviceWilsonVectorSU3 ret;
            ret.m_d[0] = deviceSU3Vector::makeZeroSU3Vector();
            ret.m_d[1] = deviceSU3Vector::makeZeroSU3Vector();
            ret.m_d[2] = deviceSU3Vector::makeZeroSU3Vector();
            ret.m_d[3] = deviceSU3Vector::makeZeroSU3Vector();
            return ret;
        }

        __device__ __inline__ static deviceWilsonVectorSU3 makeOneWilsonVectorSU3()
        {
            deviceWilsonVectorSU3 ret;
            ret.m_d[0] = deviceSU3Vector::makeOneSU3Vector();
            ret.m_d[1] = deviceSU3Vector::makeOneSU3Vector();
            ret.m_d[2] = deviceSU3Vector::makeOneSU3Vector();
            ret.m_d[3] = deviceSU3Vector::makeOneSU3Vector();
            return ret;
        }

        __device__ __inline__ static deviceWilsonVectorSU3 makeOneWilsonVectorSU3Spin(BYTE spinIndex)
        {
            deviceWilsonVectorSU3 ret;
            for (BYTE bySp = 0; bySp < 4; ++bySp)
            {
                if (bySp == spinIndex)
                {
                    ret.m_d[bySp] = deviceSU3Vector::makeOneSU3Vector();
                }
                else
                {
                    ret.m_d[bySp] = deviceSU3Vector::makeZeroSU3Vector();
                }
            }
            return ret;
        }

        __device__ __inline__ static deviceWilsonVectorSU3 makeOneWilsonVectorSU3SpinColor(BYTE spinIndex, BYTE colorIndex)
        {
            deviceWilsonVectorSU3 ret;
            for (BYTE bySp = 0; bySp < 4; ++bySp)
            {
                if (bySp == spinIndex)
                {
                    ret.m_d[bySp] = deviceSU3Vector::makeOneSU3VectorColor(colorIndex);
                }
                else
                {
                    ret.m_d[bySp] = deviceSU3Vector::makeZeroSU3Vector();
                }
            }
            return ret;
        }

        /**
        * This is sum of all elements
        * ret = sum _12 (w^* w)
        */
        __device__ __inline__ _Complex ConjugateDotC(const deviceWilsonVectorSU3& other) const
        {
            _Complex ret = m_d[0].ConjugateDotC(other.m_d[0]);
            ret = _cuCaddf(ret, m_d[1].ConjugateDotC(other.m_d[1]));
            ret = _cuCaddf(ret, m_d[2].ConjugateDotC(other.m_d[2]));
            ret = _cuCaddf(ret, m_d[3].ConjugateDotC(other.m_d[3]));
            return ret;
        }

        __device__ __inline__ void Add(const deviceWilsonVectorSU3& other)
        {
            m_d[0].Add(other.m_d[0]);
            m_d[1].Add(other.m_d[1]);
            m_d[2].Add(other.m_d[2]);
            m_d[3].Add(other.m_d[3]);
        }

        __device__ __inline__ void AddVector(const deviceSU3Vector& other)
        {
            m_d[0].Add(other);
            m_d[1].Add(other);
            m_d[2].Add(other);
            m_d[3].Add(other);
        }

        __device__ __inline__ void AddComp(const _Complex& other)
        {
            m_d[0].AddComp(other);
            m_d[1].AddComp(other);
            m_d[2].AddComp(other);
            m_d[3].AddComp(other);
        }

        __device__ __inline__ void AddReal(Real other)
        {
            m_d[0].AddReal(other);
            m_d[1].AddReal(other);
            m_d[2].AddReal(other);
            m_d[3].AddReal(other);
        }

        __device__ __inline__ void Sub(const deviceWilsonVectorSU3& other)
        {
            m_d[0].Sub(other.m_d[0]);
            m_d[1].Sub(other.m_d[1]);
            m_d[2].Sub(other.m_d[2]);
            m_d[3].Sub(other.m_d[3]);
        }

        __device__ __inline__ void SubVector(const deviceSU3Vector& other)
        {
            m_d[0].Sub(other);
            m_d[1].Sub(other);
            m_d[2].Sub(other);
            m_d[3].Sub(other);
        }

        __device__ __inline__ void SubComp(const _Complex& other)
        {
            m_d[0].SubComp(other);
            m_d[1].SubComp(other);
            m_d[2].SubComp(other);
            m_d[3].SubComp(other);
        }

        __device__ __inline__ void SubReal(Real other)
        {
            m_d[0].SubReal(other);
            m_d[1].SubReal(other);
            m_d[2].SubReal(other);
            m_d[3].SubReal(other);
        }

        __device__ __inline__ void MulReal(Real other)
        {
            m_d[0].MulReal(other);
            m_d[1].MulReal(other);
            m_d[2].MulReal(other);
            m_d[3].MulReal(other);
        }

        __device__ __inline__ void MulComp(const _Complex& other)
        {
            m_d[0].MulComp(other);
            m_d[1].MulComp(other);
            m_d[2].MulComp(other);
            m_d[3].MulComp(other);
        }

        __device__ __inline__ void MulZ4(BYTE z4)
        {
            m_d[0].MulZ4(z4);
            m_d[1].MulZ4(z4);
            m_d[2].MulZ4(z4);
            m_d[3].MulZ4(z4);
        }

        __device__ __inline__ deviceWilsonVectorSU3 AddRealC(Real other) const { deviceWilsonVectorSU3 ret(*this); ret.AddReal(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 AddCompC(const _Complex& other) const { deviceWilsonVectorSU3 ret(*this); ret.AddComp(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 AddVectorC(const deviceSU3Vector& other) const { deviceWilsonVectorSU3 ret(*this); ret.AddVector(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 AddC(const deviceWilsonVectorSU3& other) const { deviceWilsonVectorSU3 ret(*this); ret.Add(other); return ret; }

        __device__ __inline__ deviceWilsonVectorSU3 SubRealC(Real other) const { deviceWilsonVectorSU3 ret(*this); ret.SubReal(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 SubCompC(const _Complex& other) const { deviceWilsonVectorSU3 ret(*this); ret.SubComp(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 SubVectorC(const deviceSU3Vector& other) const { deviceWilsonVectorSU3 ret(*this); ret.SubVector(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 SubC(const deviceWilsonVectorSU3& other) const { deviceWilsonVectorSU3 ret(*this); ret.Sub(other); return ret; }

        __device__ __inline__ deviceWilsonVectorSU3 MulRealC(Real other) const { deviceWilsonVectorSU3 ret(*this); ret.MulReal(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 MulCompC(const _Complex& other) const { deviceWilsonVectorSU3 ret(*this); ret.MulComp(other); return ret; }
        __device__ __inline__ deviceWilsonVectorSU3 MulZ4C(BYTE z4) const { deviceWilsonVectorSU3 ret(*this); ret.MulZ4(z4); return ret; }

        union
        {
            deviceSU3Vector m_d[4];
            _Complex m_me[12];
        };
        //deviceSU3Vector m_d[4];
    };

#if defined(__cplusplus)
}
#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _VECTORS_H_

//=============================================================================
// END OF FILE
//=============================================================================
