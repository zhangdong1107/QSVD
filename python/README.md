# Quaternion SVD in Python

This NumPy implementation uses the same complex representation as the MATLAB package. It accepts dense matrices and returns one singular value per quaternion direction, in descending order. Tested with Python 3.12 and NumPy 2.3.5.

Extract the ZIP, open the directory containing `qsvd.py`, and install NumPy with `python -m pip install numpy`. Run `python -m unittest -v` there to check the installation.

```python
import numpy as np
from qsvd import quaternion_to_complex, qsvd, qsvds

rng = np.random.default_rng(8)
a0, a1, a2, a3 = rng.normal(size=(4, 16, 12))
ac = quaternion_to_complex(a0 + 1j*a1, a2 + 1j*a3)

u, s, v = qsvd(ac)
np.testing.assert_allclose((u * np.tile(s, 2)) @ v.conj().T, ac, atol=1e-11)
sigma_max = qsvds(ac, 1, values_only=True)[0]
uk, sk, vk = qsvds(ac, 3)
approx = (uk * np.tile(sk, 2)) @ vk.conj().T
```

For k quaternion triplets the returned complex factors have 2k columns. The first k columns are representatives and the last k are their quaternion partners. The middle factor is represented by a one-dimensional array of k real singular values. Use `v.conj().T` for the adjoint. Full QSVD returns economy factors, including for zero input.

Full QSVD computes a complex economy SVD and recovers paired quaternion bases. Python QSVDS uses restarted block Krylov iteration and thin projected QSVD with NumPy only. When k is close to the full dimension, a two-sided subspace iteration fits the smaller search space. For k equal to the full quaternion dimension, QSVDS calls full QSVD. For smaller k, the routine checks both singular equations and the orthogonality of the paired factors before returning. The `values_only=True` partial path performs the same convergence checks internally.

The optional arguments `tol=1e-10`, `maxiter=1000`, and `seed=0` control partial iteration. A failed convergence check raises `numpy.linalg.LinAlgError`. The software page reports timings for the MATLAB implementations. Python timings depend on the NumPy build, matrix, and hardware.

[Software page and MATLAB package](https://zhangdong1107.github.io/qsvd/) · [Accompanying paper](https://doi.org/10.1016/j.aml.2026.110047)
