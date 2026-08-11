//=============================================================================
// FILENAME : DN.h
//
// DESCRIPTION:
// deviceDN<N>: D_N dihedral group element (order 2N)
// Stored as 2x2 complex matrix in standard irreducible representation.
//   Rotation r^k:    [[e^{i\theta k}, 0          ], [0,            e^{-i\theta k}]]
//   Reflection r^k s: [[0,           e^{i\theta k}], [e^{-i\theta k}, 0          ]]
// where theta = 2\pi / N.
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================

#ifndef _DN_H_
#define _DN_H_

__BEGIN_NAMESPACE

template<INT N, INT EN=2*N>
struct deviceDN
{

public:

    static constexpr INT matrixN = 2;
    static constexpr INT elementN = EN;

    // Performance: default constructor leaves m_me uninitialized.
    // Caller must explicitly initialize via Id(), Zero(), or makeDN*() before use.
    __device__ deviceDN()
    {
    }

    __device__ deviceDN(const deviceDN<N>& other)
    {
        memcpy(m_me, other.m_me, sizeof(CLGComplex) * 4);
    }

    __device__ deviceDN(const CLGComplex* __restrict__ other)
    {
        memcpy(m_me, other, sizeof(CLGComplex) * 4);
    }

    __device__ void DebugPrint(const char* header = NULL) const
    {
        printf("%s%s{{%f%s%fi, %f%s%fi},\n {%f%s%fi, %f%s%fi}};\n",
            NULL == header ? "" : header, NULL == header ? "" : "=",
            m_me[0].x, m_me[0].y < 0 ? "" : "+", m_me[0].y,
            m_me[1].x, m_me[1].y < 0 ? "" : "+", m_me[1].y,
            m_me[2].x, m_me[2].y < 0 ? "" : "+", m_me[2].y,
            m_me[3].x, m_me[3].y < 0 ? "" : "+", m_me[3].y);
    }

#pragma region creation

    __device__ __inline__ static deviceDN<N> makeDNId()
    {
        deviceDN<N> ret;
        ret.m_me[0] = _make_cuComplex(F(1.0), F(0.0));
        ret.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
        ret.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
        ret.m_me[3] = _make_cuComplex(F(1.0), F(0.0));
        return ret;
    }

    __device__ __inline__ static deviceDN<N> makeAsK(UINT k)
    {
        deviceDN<N> ret;
        ret.SetFromIndex(k);
        return ret;
    }

    __device__ __inline__ static deviceDN<N> makeDNRandom(UINT fatIndex)
    {
        INT idx = static_cast<INT>(_deviceRandomF(fatIndex) * (2 * N)) % (2 * N);
        deviceDN<N> ret;
        ret.SetFromIndex(idx);
        return ret;
    }

    __device__ __inline__ static deviceDN<N> makeDNZero()
    {
        deviceDN<N> ret;
        ret.m_me[0] = _make_cuComplex(F(0.0), F(0.0));
        ret.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
        ret.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
        ret.m_me[3] = _make_cuComplex(F(0.0), F(0.0));
        return ret;
    }

    __device__ __inline__ void Id()
    {
        m_me[0] = _make_cuComplex(F(1.0), F(0.0));
        m_me[1] = _make_cuComplex(F(0.0), F(0.0));
        m_me[2] = _make_cuComplex(F(0.0), F(0.0));
        m_me[3] = _make_cuComplex(F(1.0), F(0.0));
    }

    __device__ __inline__ void Zero()
    {
        m_me[0] = _make_cuComplex(F(0.0), F(0.0));
        m_me[1] = _make_cuComplex(F(0.0), F(0.0));
        m_me[2] = _make_cuComplex(F(0.0), F(0.0));
        m_me[3] = _make_cuComplex(F(0.0), F(0.0));
    }

    __device__ __inline__ void SetFromIndex(UINT idx)
    {
        UINT k = ((idx % (2 * N)) + (2 * N)) % (2 * N);
        UINT parity = k / N;
        k = k % N;
        Real theta = F(2.0) * PI * k / N;
        CLGComplex eip = _make_cuComplex(_cos(theta), _sin(theta));
        CLGComplex eim = _make_cuComplex(_cos(theta), -_sin(theta));
        if (0 == parity)
        {
            m_me[0] = eip;
            m_me[1] = _make_cuComplex(F(0.0), F(0.0));
            m_me[2] = _make_cuComplex(F(0.0), F(0.0));
            m_me[3] = eim;
        }
        else
        {
            m_me[0] = _make_cuComplex(F(0.0), F(0.0));
            m_me[1] = eip;
            m_me[2] = eim;
            m_me[3] = _make_cuComplex(F(0.0), F(0.0));
        }
    }

#pragma endregion

#pragma region operators

    __device__ __inline__ void Add(const deviceDN<N>& right)
    {
        m_me[0] = _cuCaddf(m_me[0], right.m_me[0]);
        m_me[1] = _cuCaddf(m_me[1], right.m_me[1]);
        m_me[2] = _cuCaddf(m_me[2], right.m_me[2]);
        m_me[3] = _cuCaddf(m_me[3], right.m_me[3]);
    }

    __device__ __inline__ void Sub(const deviceDN<N>& right)
    {
        m_me[0] = _cuCsubf(m_me[0], right.m_me[0]);
        m_me[1] = _cuCsubf(m_me[1], right.m_me[1]);
        m_me[2] = _cuCsubf(m_me[2], right.m_me[2]);
        m_me[3] = _cuCsubf(m_me[3], right.m_me[3]);
    }

    __device__ __inline__ void Mul(const deviceDN<N>& right)
    {
        CLGComplex t0 = _cuCaddf(_cuCmulf(m_me[0], right.m_me[0]), _cuCmulf(m_me[1], right.m_me[2]));
        CLGComplex t1 = _cuCaddf(_cuCmulf(m_me[0], right.m_me[1]), _cuCmulf(m_me[1], right.m_me[3]));
        CLGComplex t2 = _cuCaddf(_cuCmulf(m_me[2], right.m_me[0]), _cuCmulf(m_me[3], right.m_me[2]));
        CLGComplex t3 = _cuCaddf(_cuCmulf(m_me[2], right.m_me[1]), _cuCmulf(m_me[3], right.m_me[3]));
        m_me[0] = t0;
        m_me[1] = t1;
        m_me[2] = t2;
        m_me[3] = t3;
    }

    __device__ __inline__ void MulComp(const CLGComplex& right)
    {
        m_me[0] = _cuCmulf(m_me[0], right);
        m_me[1] = _cuCmulf(m_me[1], right);
        m_me[2] = _cuCmulf(m_me[2], right);
        m_me[3] = _cuCmulf(m_me[3], right);
    }

    __device__ __inline__ void MulReal(Real right)
    {
        m_me[0] = cuCmulf_cr(m_me[0], right);
        m_me[1] = cuCmulf_cr(m_me[1], right);
        m_me[2] = cuCmulf_cr(m_me[2], right);
        m_me[3] = cuCmulf_cr(m_me[3], right);
    }

    __device__ __inline__ void Dagger()
    {
        CLGComplex t1 = m_me[1];
        m_me[1] = _cuConjf(m_me[2]);
        m_me[2] = _cuConjf(t1);
        m_me[0] = _cuConjf(m_me[0]);
        m_me[3] = _cuConjf(m_me[3]);
    }

    __device__ __inline__ void DaggerMul(const deviceDN<N>& right)
    {
        Dagger();
        Mul(right);
    }

    __device__ __inline__ void MulDagger(const deviceDN<N>& right)
    {
        deviceDN<N> rDag(right);
        rDag.Dagger();
        Mul(rDag);
    }

    __device__ __inline__ void Opposite()
    {
        m_me[0].x = -m_me[0].x; m_me[0].y = -m_me[0].y;
        m_me[1].x = -m_me[1].x; m_me[1].y = -m_me[1].y;
        m_me[2].x = -m_me[2].x; m_me[2].y = -m_me[2].y;
        m_me[3].x = -m_me[3].x; m_me[3].y = -m_me[3].y;
    }

    __device__ __inline__ deviceDN<N> AddC(const deviceDN<N>& right) const { deviceDN<N> ret(*this); ret.Add(right); return ret; }
    __device__ __inline__ deviceDN<N> SubC(const deviceDN<N>& right) const { deviceDN<N> ret(*this); ret.Sub(right); return ret; }
    __device__ __inline__ deviceDN<N> MulC(const deviceDN<N>& right) const { deviceDN<N> ret(*this); ret.Mul(right); return ret; }
    __device__ __inline__ deviceDN<N> MulCompC(const CLGComplex& right) const { deviceDN<N> ret(*this); ret.MulComp(right); return ret; }
    __device__ __inline__ deviceDN<N> MulRealC(const Real& right) const { deviceDN<N> ret(*this); ret.MulReal(right); return ret; }
    __device__ __inline__ deviceDN<N> DaggerC() const { deviceDN<N> ret(*this); ret.Dagger(); return ret; }
    __device__ __inline__ deviceDN<N> OppositeC() const { deviceDN<N> ret(*this); ret.Opposite(); return ret; }
    __device__ __inline__ deviceDN<N> DaggerMulC(const deviceDN<N>& right) const { deviceDN<N> ret(*this); ret.DaggerMul(right); return ret; }
    __device__ __inline__ deviceDN<N> MulDaggerC(const deviceDN<N>& right) const { deviceDN<N> ret(*this); ret.MulDagger(right); return ret; }

    __device__ __inline__ void AddDagger(const deviceDN<N>& right)
    {
        deviceDN<N> rDag(right);
        rDag.Dagger();
        Add(rDag);
    }

    __device__ __inline__ void SubDagger(const deviceDN<N>& right)
    {
        deviceDN<N> rDag(right);
        rDag.Dagger();
        Sub(rDag);
    }

    __device__ __inline__ void AddReal(Real right)
    {
        m_me[0].x += right;
        m_me[3].x += right;
    }

    __device__ __inline__ void AddComp(const CLGComplex& right)
    {
        m_me[0] = _cuCaddf(m_me[0], right);
        m_me[3] = _cuCaddf(m_me[3], right);
    }

    __device__ __inline__ void SubReal(Real right)
    {
        m_me[0].x -= right;
        m_me[3].x -= right;
    }

    __device__ __inline__ void SubComp(const CLGComplex& right)
    {
        m_me[0] = _cuCsubf(m_me[0], right);
        m_me[3] = _cuCsubf(m_me[3], right);
    }

    __device__ __inline__ void DivComp(const CLGComplex& right)
    {
        m_me[0] = _cuCdivf(m_me[0], right);
        m_me[1] = _cuCdivf(m_me[1], right);
        m_me[2] = _cuCdivf(m_me[2], right);
        m_me[3] = _cuCdivf(m_me[3], right);
    }

    __device__ __inline__ void DivReal(const Real& right)
    {
        m_me[0] = cuCdivf_cr(m_me[0], right);
        m_me[1] = cuCdivf_cr(m_me[1], right);
        m_me[2] = cuCdivf_cr(m_me[2], right);
        m_me[3] = cuCdivf_cr(m_me[3], right);
    }

    __device__ __inline__ deviceDN<N> AddCompC(const CLGComplex& right) const { deviceDN<N> ret(*this); ret.AddComp(right); return ret; }
    __device__ __inline__ deviceDN<N> AddRealC(const Real& right) const { deviceDN<N> ret(*this); ret.AddReal(right); return ret; }
    __device__ __inline__ deviceDN<N> SubCompC(const CLGComplex& right) const { deviceDN<N> ret(*this); ret.SubComp(right); return ret; }
    __device__ __inline__ deviceDN<N> SubRealC(const Real& right) const { deviceDN<N> ret(*this); ret.SubReal(right); return ret; }
    __device__ __inline__ deviceDN<N> DivCompC(const CLGComplex& right) const { deviceDN<N> ret(*this); ret.DivComp(right); return ret; }
    __device__ __inline__ deviceDN<N> DivRealC(const Real& right) const { deviceDN<N> ret(*this); ret.DivReal(right); return ret; }

    __device__ __inline__ void MulOnMe(const deviceDN<N>& left)
    {
        CLGComplex t0 = _cuCaddf(_cuCmulf(left.m_me[0], m_me[0]), _cuCmulf(left.m_me[1], m_me[2]));
        CLGComplex t1 = _cuCaddf(_cuCmulf(left.m_me[0], m_me[1]), _cuCmulf(left.m_me[1], m_me[3]));
        CLGComplex t2 = _cuCaddf(_cuCmulf(left.m_me[2], m_me[0]), _cuCmulf(left.m_me[3], m_me[2]));
        CLGComplex t3 = _cuCaddf(_cuCmulf(left.m_me[2], m_me[1]), _cuCmulf(left.m_me[3], m_me[3]));
        m_me[0] = t0;
        m_me[1] = t1;
        m_me[2] = t2;
        m_me[3] = t3;
    }

    __device__ __inline__ void MulOnMeDN(const deviceDN<N>& left)
    {
        deviceDN<N> lDag(left);
        lDag.Dagger();
        lDag.MulOnMe(*this);
        *this = lDag;
    }

    __device__ __inline__ void MulOnMeND(const deviceDN<N>& left)
    {
        Dagger();
        MulOnMe(left);
    }

#pragma endregion

#pragma region useful functions

    __device__ __inline__ CLGComplex Tr() const
    {
        return _cuCaddf(m_me[0], m_me[3]);
    }

    __device__ __inline__ Real ReTr() const
    {
        return m_me[0].x + m_me[3].x;
    }

    __device__ __inline__ Real ImTr() const
    {
        return m_me[0].y + m_me[3].y;
    }

    __device__ __inline__ void Re()
    {
        m_me[0].y = F(0.0);
        m_me[1].y = F(0.0);
        m_me[2].y = F(0.0);
        m_me[3].y = F(0.0);
    }

    __device__ __inline__ deviceDN<N> ReC() const { deviceDN<N> ret(*this); ret.Re(); return ret; }

    __device__ __inline__ void Im()
    {
        m_me[0].x = m_me[0].y; m_me[0].y = F(0.0);
        m_me[1].x = m_me[1].y; m_me[1].y = F(0.0);
        m_me[2].x = m_me[2].y; m_me[2].y = F(0.0);
        m_me[3].x = m_me[3].y; m_me[3].y = F(0.0);
    }

    __device__ __inline__ deviceDN<N> ImC() const { deviceDN<N> ret(*this); ret.Im(); return ret; }

    __device__ __inline__ void AddId()
    {
        m_me[0].x += F(1.0);
        m_me[3].x += F(1.0);
    }

    /**
     * Project to nearest D_N element.
     * D_N elements in 2D irrep:
     *   Rotation r^k:    diag(e^{i\theta k}, e^{-i\theta k}), k=0..N-1
     *   Reflection r^k s: anti-diag(e^{i\theta k}, e^{-i\theta k}), k=0..N-1
     * where theta = 2\pi / N.
     *
     * Finds closest element by maximizing Re[Tr(V^dag * this)] over all V in D_N.
     */
    __device__ __inline__ void Proj(INT ite = 0)
    {
        UN_USE(ite);
        Real theta = F(2.0) * PI / N;
        Real bestVal = F(-1e30);
        INT bestK = 0;
        INT bestParity = 0;

        for (INT k = 0; k < N; ++k)
        {
            Real thk = theta * k;
            Real cosk = _cos(thk);
            Real sink = _sin(thk);

            // Rotation (parity=0): overlap with diag(e^{i\theta k}, e^{-i\theta k})
            // Re[Tr(rot^dag * this)] = Re[e^{-i\theta k}*U_00 + e^{i\theta k}*U_11]
            Real rotVal = cosk * m_me[0].x + sink * m_me[0].y
                        + cosk * m_me[3].x - sink * m_me[3].y;

            if (rotVal > bestVal)
            {
                bestVal = rotVal;
                bestK = k;
                bestParity = 0;
            }

            // Reflection (parity=1): overlap with anti-diag(e^{i\theta k}, e^{-i\theta k})
            // Re[Tr(ref^dag * this)] = Re[e^{-i\theta k}*U_10 + e^{i\theta k}*U_01]
            Real refVal = cosk * m_me[2].x + sink * m_me[2].y
                        + cosk * m_me[1].x - sink * m_me[1].y;

            if (refVal > bestVal)
            {
                bestVal = refVal;
                bestK = k;
                bestParity = 1;
            }
        }

        Real thk = theta * bestK;
        CLGComplex eip = _make_cuComplex(_cos(thk), _sin(thk));
        CLGComplex eim = _make_cuComplex(_cos(thk), -_sin(thk));
        if (0 == bestParity)
        {
            m_me[0] = eip;
            m_me[1] = _make_cuComplex(F(0.0), F(0.0));
            m_me[2] = _make_cuComplex(F(0.0), F(0.0));
            m_me[3] = eim;
        }
        else
        {
            m_me[0] = _make_cuComplex(F(0.0), F(0.0));
            m_me[1] = eip;
            m_me[2] = eim;
            m_me[3] = _make_cuComplex(F(0.0), F(0.0));
        }
    }

    __device__ __inline__ void Norm()
    {
        Proj();
    }

    __device__ __inline__ deviceDN<N> Exp(const CLGComplex& a, BYTE uiPrecision = N + 1) const
    {
        deviceDN<N> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            const Real exp_factor = __rcp(uiPrecision - i);
            CLGComplex alpha = cuCmulf_cr(a, exp_factor);
            deviceDN<N> aUoN(*this);
            aUoN.MulComp(alpha);
            if (0 == i) { tmp = aUoN; }
            else { tmp.Mul(aUoN); }
            tmp.AddId();
        }
        return tmp;
    }

    __device__ __inline__ deviceDN<N> ExpReal(Real a, BYTE uiPrecision = N + 1) const
    {
        deviceDN<N> tmp;
        for (BYTE i = 0; i < uiPrecision; ++i)
        {
            deviceDN<N> aUoN(*this);
            aUoN.MulReal(a * __rcp(uiPrecision - i));
            if (0 == i) { tmp = aUoN; }
            else { tmp.Mul(aUoN); }
            tmp.AddId();
        }
        tmp.Norm();
        return tmp;
    }

    __device__ __inline__ void CalculateEigenValues(CLGComplex& c1, CLGComplex& c2, deviceDN<N>& vectors) const
    {
        CLGComplex tr = _cuCaddf(m_me[0], m_me[3]);
        CLGComplex det = _cuCsubf(_cuCmulf(m_me[0], m_me[3]), _cuCmulf(m_me[1], m_me[2]));

        CLGComplex disc = _cuCsubf(_cuCmulf(tr, tr), _cuCmulf(_make_cuComplex(F(4.0), F(0.0)), det));
        CLGComplex sqrt_disc = __cuCsqrtf(disc);

        c1 = _cuCmulf(_cuCaddf(tr, sqrt_disc), _make_cuComplex(F(0.5), F(0.0)));
        c2 = _cuCmulf(_cuCsubf(tr, sqrt_disc), _make_cuComplex(F(0.5), F(0.0)));

        vectors.m_me[0] = m_me[1];
        vectors.m_me[2] = _cuCsubf(c1, m_me[0]);

        Real norm1_sq = __cuCabsSqf(vectors.m_me[0]) + __cuCabsSqf(vectors.m_me[2]);
        if (norm1_sq < _CLG_FLT_MIN)
        {
            vectors.m_me[0] = _cuCsubf(c1, m_me[3]);
            vectors.m_me[2] = m_me[2];
            norm1_sq = __cuCabsSqf(vectors.m_me[0]) + __cuCabsSqf(vectors.m_me[2]);
            if (norm1_sq < _CLG_FLT_MIN)
            {
                vectors.m_me[0] = _make_cuComplex(F(1.0), F(0.0));
                vectors.m_me[2] = _make_cuComplex(F(0.0), F(0.0));
                norm1_sq = F(1.0);
            }
        }

        vectors.m_me[1] = m_me[1];
        vectors.m_me[3] = _cuCsubf(c2, m_me[0]);

        Real norm2_sq = __cuCabsSqf(vectors.m_me[1]) + __cuCabsSqf(vectors.m_me[3]);
        if (norm2_sq < _CLG_FLT_MIN)
        {
            vectors.m_me[1] = _cuCsubf(c2, m_me[3]);
            vectors.m_me[3] = m_me[2];
            norm2_sq = __cuCabsSqf(vectors.m_me[1]) + __cuCabsSqf(vectors.m_me[3]);
            if (norm2_sq < _CLG_FLT_MIN)
            {
                vectors.m_me[1] = _make_cuComplex(F(0.0), F(0.0));
                vectors.m_me[3] = _make_cuComplex(F(1.0), F(0.0));
                norm2_sq = F(1.0);
            }
        }

        Real inv_norm1 = __rcp(_sqrt(norm1_sq));
        Real inv_norm2 = __rcp(_sqrt(norm2_sq));

        vectors.m_me[0] = cuCmulf_cr(vectors.m_me[0], inv_norm1);
        vectors.m_me[2] = cuCmulf_cr(vectors.m_me[2], inv_norm1);
        vectors.m_me[1] = cuCmulf_cr(vectors.m_me[1], inv_norm2);
        vectors.m_me[3] = cuCmulf_cr(vectors.m_me[3], inv_norm2);
    }

    __device__ __inline__ deviceDN<N> Log() const
    {
        CLGComplex c1, c2;
        deviceDN<N> v;
        CalculateEigenValues(c1, c2, v);

        deviceDN<N> ret;
        ret.m_me[0] = __cuClogf(c1);
        ret.m_me[1] = _zeroc;
        ret.m_me[2] = _zeroc;
        ret.m_me[3] = __cuClogf(c2);

        ret.MulDagger(v);
        ret = v.MulC(ret);
        return ret;
    }

    __device__ __inline__ deviceDN<N> StrictExp() const
    {
        CLGComplex c1, c2;
        deviceDN<N> v;
        CalculateEigenValues(c1, c2, v);

        deviceDN<N> ret;
        ret.m_me[0] = __cuCexpf(c1);
        ret.m_me[1] = _zeroc;
        ret.m_me[2] = _zeroc;
        ret.m_me[3] = __cuCexpf(c2);

        ret.MulDagger(v);
        ret = v.MulC(ret);
        return ret;
    }

    __device__ __inline__ void Ta()
    {
        deviceDN<N> uDag(*this);
        uDag.Dagger();
        Sub(uDag);
        MulReal(F(0.5));
        Real halfTrIm = F(0.5) * (m_me[0].y + m_me[3].y);
        m_me[0].y -= halfTrIm;
        m_me[3].y -= halfTrIm;
    }

    __device__ __inline__ void Th()
    {
        deviceDN<N> uDag(*this);
        uDag.Dagger();
        Add(uDag);
        MulReal(F(0.5));
        Real halfTrRe = F(0.5) * (m_me[0].x + m_me[3].x);
        m_me[0].x -= halfTrRe;
        m_me[3].x -= halfTrRe;
        m_me[0].y = F(0.0);
        m_me[1].y = F(0.0);
        m_me[2].y = F(0.0);
        m_me[3].y = F(0.0);
    }

    __device__ __inline__ void Re2()
    {
        Re();
        MulReal(F(2.0));
    }

    __device__ __inline__ void iIm2()
    {
        Im();
        MulReal(F(2.0));
    }

    __device__ __inline__ static Real TrIm(const deviceDN<N>& a, const deviceDN<N>& b)
    {
        return a.m_me[0].x * b.m_me[0].y - a.m_me[0].y * b.m_me[0].x
             + a.m_me[1].x * b.m_me[1].y - a.m_me[1].y * b.m_me[1].x
             + a.m_me[2].x * b.m_me[2].y - a.m_me[2].y * b.m_me[2].x
             + a.m_me[3].x * b.m_me[3].y - a.m_me[3].y * b.m_me[3].x;
    }

    __device__ __inline__ CLGComplex Det() const
    {
        return _cuCsubf(_cuCmulf(m_me[0], m_me[3]), _cuCmulf(m_me[1], m_me[2]));
    }

#pragma endregion

#pragma region global device helpers

    CLGComplex m_me[4];
};


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _DN_H_

//=============================================================================
// END OF FILE
//=============================================================================
