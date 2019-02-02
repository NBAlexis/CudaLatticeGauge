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

__BEGIN_NAMESPACE

struct CLGAPI deviceSU3Vector
{
public:
    __device__ deviceSU3Vector() 
    { 
        m_ve[0] = _make_cuComplex(0, 0);
        m_ve[1] = _make_cuComplex(0, 0);
        m_ve[2] = _make_cuComplex(0, 0);
    }

    __device__ deviceSU3Vector(const deviceSU3Vector& other)
    {
        m_ve[0] = other.m_ve[0];
        m_ve[1] = other.m_ve[1];
        m_ve[2] = other.m_ve[2];
    }

    __device__ __inline__ void MakeRandomGaussian(UINT fatIndex)
    {
        m_ve[0] = _deviceRandomGaussC(fatIndex);
        m_ve[1] = _deviceRandomGaussC(fatIndex);
        m_ve[2] = _deviceRandomGaussC(fatIndex);
    }

    __device__ __inline__ void MakeZero()
    {
        m_ve[0] = _make_cuComplex(0, 0);
        m_ve[1] = _make_cuComplex(0, 0);
        m_ve[2] = _make_cuComplex(0, 0);
    }

    __device__ __inline__ _Complex ConjugateDotC(const deviceSU3Vector& other) const
    {
        _Complex ret = _make_cuComplex(0, 0);
        ret = _cuCaddf(ret, _cuCmulf(_cuConjf(m_ve[0]), other.m_ve[0]));
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

    __device__ __inline__ void Sub(const _Complex& other)
    {
        m_ve[0] = _cuCsubf(m_ve[0], other);
        m_ve[1] = _cuCsubf(m_ve[1], other);
        m_ve[2] = _cuCsubf(m_ve[2], other);
    }

    __device__ __inline__ void Sub(Real other)
    {
        m_ve[0] = _cuCsubf(m_ve[0], other);
        m_ve[1] = _cuCsubf(m_ve[1], other);
        m_ve[2] = _cuCsubf(m_ve[2], other);
    }

    __device__ __inline__ void Add(const deviceSU3Vector& other)
    {
        m_ve[0] = _cuCaddf(m_ve[0], other.m_ve[0]);
        m_ve[1] = _cuCaddf(m_ve[1], other.m_ve[1]);
        m_ve[2] = _cuCaddf(m_ve[2], other.m_ve[2]);
    }

    __device__ __inline__ void Add(const _Complex& other)
    {
        m_ve[0] = _cuCaddf(m_ve[0], other);
        m_ve[1] = _cuCaddf(m_ve[1], other);
        m_ve[2] = _cuCaddf(m_ve[2], other);
    }

    __device__ __inline__ void Add(Real other)
    {
        m_ve[0] = _cuCaddf(m_ve[0], other);
        m_ve[1] = _cuCaddf(m_ve[1], other);
        m_ve[2] = _cuCaddf(m_ve[2], other);
    }

    __device__ __inline__ void Mul(Real other)
    {
        m_ve[0] = _cuCmulf(m_ve[0], other);
        m_ve[1] = _cuCmulf(m_ve[1], other);
        m_ve[2] = _cuCmulf(m_ve[2], other);
    }

    __device__ __inline__ void Mul(const _Complex& other)
    {
        m_ve[0] = _cuCmulf(m_ve[0], other);
        m_ve[1] = _cuCmulf(m_ve[1], other);
        m_ve[2] = _cuCmulf(m_ve[2], other);
    }

    __device__ __inline__ void Mul(const cuComplexI& other)
    {
        m_ve[0] = cuCmulI(m_ve[0], other);
        m_ve[1] = cuCmulI(m_ve[1], other);
        m_ve[2] = cuCmulI(m_ve[2], other);
    }

    __device__ __inline__ deviceSU3Vector SubC(Real other) const { deviceSU3Vector ret(*this); ret.Sub(other); return ret; }
    __device__ __inline__ deviceSU3Vector SubC(const _Complex& other) const { deviceSU3Vector ret(*this); ret.Sub(other); return ret; }
    __device__ __inline__ deviceSU3Vector SubC(const deviceSU3Vector& other) const { deviceSU3Vector ret(*this); ret.Sub(other); return ret; }

    __device__ __inline__ deviceSU3Vector AddC(const Real& other) const { deviceSU3Vector ret(*this); ret.Add(other); return ret; }
    __device__ __inline__ deviceSU3Vector AddC(const _Complex& other) const { deviceSU3Vector ret(*this); ret.Add(other); return ret; }
    __device__ __inline__ deviceSU3Vector AddC(const deviceSU3Vector& other) const { deviceSU3Vector ret(*this); ret.Add(other); return ret; }

    __device__ __inline__ deviceSU3Vector MulC(Real other) const { deviceSU3Vector ret(*this); ret.Mul(other); return ret; }
    __device__ __inline__ deviceSU3Vector MulC(const _Complex& other) const { deviceSU3Vector ret(*this); ret.Mul(other); return ret; }
    __device__ __inline__ deviceSU3Vector MulC(const cuComplexI& other) const { deviceSU3Vector ret(*this); ret.Mul(other); return ret; }

    _Complex m_ve[3];
};

struct CLGAPI deviceWilsonVectorSU3
{
public:
    __device__ deviceWilsonVectorSU3() { ; }

    __device__ deviceWilsonVectorSU3(const deviceWilsonVectorSU3& other) 
    { 
        m_d[0] = other.m_d[0];
        m_d[1] = other.m_d[1];
        m_d[2] = other.m_d[2];
        m_d[3] = other.m_d[3];
    }

    __device__ __inline__ void MakeRandomGaussian(UINT fatIndex)
    { 
        m_d[0].MakeRandomGaussian(fatIndex);
        m_d[1].MakeRandomGaussian(fatIndex);
        m_d[2].MakeRandomGaussian(fatIndex);
        m_d[3].MakeRandomGaussian(fatIndex);
    }

    __device__ __inline__ void MakeZero()
    {
        m_d[0].MakeZero();
        m_d[1].MakeZero();
        m_d[2].MakeZero();
        m_d[3].MakeZero();
    }

    __device__ __inline__ _Complex ConjugateDotC(const deviceWilsonVectorSU3& other) const
    {
        _Complex ret = _make_cuComplex(0, 0);
        ret = _cuCaddf(ret, m_d[0].ConjugateDotC(other.m_d[0]));
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

    __device__ __inline__ void Add(const deviceSU3Vector& other)
    {
        m_d[0].Add(other);
        m_d[1].Add(other);
        m_d[2].Add(other);
        m_d[3].Add(other);
    }

    __device__ __inline__ void Add(const _Complex& other)
    {
        m_d[0].Add(other);
        m_d[1].Add(other);
        m_d[2].Add(other);
        m_d[3].Add(other);
    }

    __device__ __inline__ void Add(Real other)
    {
        m_d[0].Add(other);
        m_d[1].Add(other);
        m_d[2].Add(other);
        m_d[3].Add(other);
    }

    __device__ __inline__ void Sub(const deviceWilsonVectorSU3& other)
    {
        m_d[0].Sub(other.m_d[0]);
        m_d[1].Sub(other.m_d[1]);
        m_d[2].Sub(other.m_d[2]);
        m_d[3].Sub(other.m_d[3]);
    }

    __device__ __inline__ void Sub(const deviceSU3Vector& other)
    {
        m_d[0].Sub(other);
        m_d[1].Sub(other);
        m_d[2].Sub(other);
        m_d[3].Sub(other);
    }

    __device__ __inline__ void Sub(const _Complex& other)
    {
        m_d[0].Sub(other);
        m_d[1].Sub(other);
        m_d[2].Sub(other);
        m_d[3].Sub(other);
    }

    __device__ __inline__ void Sub(Real other)
    {
        m_d[0].Sub(other);
        m_d[1].Sub(other);
        m_d[2].Sub(other);
        m_d[3].Sub(other);
    }

    __device__ __inline__ void Mul(Real other)
    {
        m_d[0].Mul(other);
        m_d[1].Mul(other);
        m_d[2].Mul(other);
        m_d[3].Mul(other);
    }

    __device__ __inline__ void Mul(const _Complex& other)
    {
        m_d[0].Mul(other);
        m_d[1].Mul(other);
        m_d[2].Mul(other);
        m_d[3].Mul(other);
    }

    __device__ __inline__ deviceWilsonVectorSU3 AddC(Real other) const { deviceWilsonVectorSU3 ret(*this); ret.Add(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 AddC(const _Complex& other) const { deviceWilsonVectorSU3 ret(*this); ret.Add(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 AddC(const deviceSU3Vector& other) const { deviceWilsonVectorSU3 ret(*this); ret.Add(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 AddC(const deviceWilsonVectorSU3& other) const { deviceWilsonVectorSU3 ret(*this); ret.Add(other); return ret; }

    __device__ __inline__ deviceWilsonVectorSU3 SubC(Real other) const { deviceWilsonVectorSU3 ret(*this); ret.Sub(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 SubC(const _Complex& other) const { deviceWilsonVectorSU3 ret(*this); ret.Sub(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 SubC(const deviceSU3Vector& other) const { deviceWilsonVectorSU3 ret(*this); ret.Sub(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 SubC(const deviceWilsonVectorSU3& other) const { deviceWilsonVectorSU3 ret(*this); ret.Sub(other); return ret; }

    __device__ __inline__ deviceWilsonVectorSU3 MulC(Real other) const { deviceWilsonVectorSU3 ret(*this); ret.Mul(other); return ret; }
    __device__ __inline__ deviceWilsonVectorSU3 MulC(const _Complex& other) const { deviceWilsonVectorSU3 ret(*this); ret.Mul(other); return ret; }

    //union
    //{
    //    deviceSU3Vector m_d[4];
    //    _Complex m_me[12];
    //};
    deviceSU3Vector m_d[4];
};


__END_NAMESPACE

#endif //#ifndef _VECTORS_H_

//=============================================================================
// END OF FILE
//=============================================================================
