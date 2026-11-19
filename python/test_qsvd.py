"""Run with python -m unittest -v from this directory."""
import unittest
from unittest.mock import patch

import numpy as np
from qsvd import qsvd, qsvds, quaternion_to_complex


def rotation(n, rng):
    # Real orthogonal mixing and diagonal unit quaternions are unitary.
    mix, _ = np.linalg.qr(rng.normal(size=(n, n)))
    phases = rng.normal(size=(4, n))
    phases /= np.linalg.norm(phases, axis=0)
    return quaternion_to_complex((phases[0]+1j*phases[1])[:, None]*mix,
                                 (phases[2]+1j*phases[3])[:, None]*mix)


class TestQSVD(unittest.TestCase):
    def check_factors(self, a, result, expected):
        u, s, v = result
        k = len(s)
        np.testing.assert_allclose(s, expected[:k], atol=1e-12*expected[0], rtol=1e-10)
        for factor in (u, v):
            mid = factor.shape[0]//2
            partner = np.vstack((-factor[mid:, :k].conj(), factor[:mid, :k].conj()))
            np.testing.assert_allclose(factor[:, k:], partner, atol=1e-12)
            np.testing.assert_allclose(factor.conj().T @ factor, np.eye(2*k), atol=1e-10)
        twice = np.tile(s, 2)
        scale = expected[0] or 1
        np.testing.assert_allclose((a/scale) @ v, u*(twice/scale), atol=2e-10)
        np.testing.assert_allclose((a/scale).conj().T @ u, v*(twice/scale), atol=2e-10)
        if k == min(a.shape)//2:
            np.testing.assert_allclose((u*(twice/scale)) @ v.conj().T, a/scale, atol=1e-11)

    def test_known_spectra(self):
        rng = np.random.default_rng(41)
        for m, n in [(1, 4), (12, 9), (9, 12), (13, 13)]:
            p = min(m, n)
            for s in [np.linspace(1, .85, p), np.ones(p), np.zeros(p),
                      np.r_[1., np.zeros(p-1)], np.geomspace(1, 1e-13, p)]:
                d = np.zeros((m, n))
                np.fill_diagonal(d, s)
                a = rotation(m, rng) @ quaternion_to_complex(d, 0*d) @ rotation(n, rng).conj().T
                with self.subTest(m=m, n=n, s=s):
                    self.check_factors(a, qsvd(a), s)
                    for k in sorted({1, min(3,p), max(1,p-1), p}):
                        self.check_factors(a, qsvds(a, k), s)
                    np.testing.assert_allclose(qsvd(a, values_only=True), s, atol=1e-14)

    def test_partial_sizes_and_scales(self):
        rng = np.random.default_rng(421)
        parts = rng.normal(size=(4, 90, 72))
        a = quaternion_to_complex(parts[0]+1j*parts[1], parts[2]+1j*parts[3])
        expected = np.linalg.svd(a, compute_uv=False)[::2]
        dense = np.linalg.svd
        def thin_only(x, *args, **kwargs):
            self.assertLess(min(x.shape), min(a.shape))
            return dense(x, *args, **kwargs)
        with patch("numpy.linalg.svd", side_effect=thin_only):
            for scale in (1e-180, 1e180):
                self.check_factors(a*scale, qsvds(a*scale, 3), expected*scale)
            np.testing.assert_allclose(qsvds(a, 1, values_only=True), expected[:1])
            with self.assertRaises(np.linalg.LinAlgError):
                qsvds(a, 3, maxiter=1)

    def test_invalid_input(self):
        for k in (0, -1, 3, True, 1.5):
            with self.assertRaises((ValueError, TypeError)):
                qsvds(np.eye(4), k)
        for x in (np.ones((3, 4)), np.full((4, 4), np.nan), np.ones((4, 4))):
            with self.assertRaises(ValueError):
                qsvd(x)

    def test_independent_real_embedding(self):
        a, b, c, d = np.random.default_rng(107).normal(size=(4, 21, 16))
        real = np.block([[a, -b, -c, -d], [b, a, -d, c],
                         [c, d, a, -b], [d, -c, b, a]])
        expected = np.linalg.svd(real, compute_uv=False)[::4]
        ac = quaternion_to_complex(a+1j*b, c+1j*d)
        self.check_factors(ac, qsvd(ac), expected)
        self.check_factors(ac, qsvds(ac, 5), expected)


if __name__ == "__main__":
    unittest.main()
