> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 08c. KS/HISQ 质量约定

## KS/HISQ rational-approximation mass convention

For staggered / HISQ fermions, CLGLib supports two equivalent YAML conventions. Let

```text
D = D0 + 2am
```

where `D0` is often called `dslash` in the literature, and `D` is sometimes called `M`.

### Convention 1: put the mass in `Mass`

Use rational approximations of

```text
f(D^\dagger D)
```

and set

```yaml
Mass: <2am>
```

### Convention 2: absorb the mass into rational coefficients

Use rational approximations of

```text
f(D0^\dagger D0)
```

and set

```yaml
Mass: 0
```

In this convention, the physical mass is encoded in the `MC` / `MD` rational-approximation coefficients. The HISQ block in `RotationImprovedFermion.yaml` uses this convention: HISQ fermion fields set `Mass: 0.0`, while the comments' `2am` values enter coefficients such as `(x + 4m^2)^p`.

Do not mix the two conventions. If `Mass` is already `2am`, the rational coefficients should correspond to `D^\dagger D`; if the rational coefficients already contain mass terms such as `(x + 4m^2)^p` or mass-preconditioned ratios, set `Mass: 0` to avoid adding the mass twice.

---


---

[< 返回 08. 费米子场](08-fermion-fields.md) | [< 返回 yaml-reference 目录](home.md)
