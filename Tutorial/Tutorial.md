# CLGLib 入门教程

> 本教程基于 `Code/Applications/CLGExample` 这一最简示例，手把手教你用 CLGLib 完成四种典型模拟任务，以及如何在已保存的组态（configuration）上进行物理量测量。
>
> 阅读前请确保你已能成功编译 CLGLib（参见项目根目录 `README.md` 的编译说明）。
>
> **作者**：Kimi-k2.6

---

## 目录

| 章节 | 内容 | 文件 |
|------|------|------|
| 1 | [CLGLib 程序的基本骨架](01-skeleton.md) | `01-skeleton.md` |
| 2 | [纯 SU(3) 规范场的 HMC 模拟](02-pure-su3.md) | `02-pure-su3.md` |
| 3 | [纯 Z₂ 规范场的 Heatbath 模拟](03-pure-z2.md) | `03-pure-z2.md` |
| 4 | [SU(3) + Wilson Dirac 费米子的 HMC 模拟](04-su3-wilson.md) | `04-su3-wilson.md` |
| 5 | [SU(3) + Staggered (KS) 费米子的 HMC 模拟](05-su3-ks.md) | `05-su3-ks.md` |
| 6 | [对已保存组态进行测量](06-measurements.md) | `06-measurements.md` |
| 7 | [如何把你的项目加入 CMake 编译](07-cmake.md) | `07-cmake.md` |
| 8 | [常见问题速查](08-faq.md) | `08-faq.md` |

---

## 快速导航

**第一次使用？** 从 [1. 基本骨架](01-skeleton.md) 开始，了解一个最小 CLGLib 程序的构成。

**只想跑纯规范场？** 看 [2. 纯 SU(3)](02-pure-su3.md) 或 [3. 纯 Z₂](03-pure-z2.md)。

**需要加费米子？** [4. Wilson Dirac](04-su3-wilson.md) 或 [5. Staggered/KS](05-su3-ks.md)。

**已有组态，要做测量？** 直接跳 [6. 测量](06-measurements.md)。

---

*本教程由 Kimi 2.6 自动生成，如有疑问请参考 `wiki/home.md` 或项目根目录 `README.md`。*
