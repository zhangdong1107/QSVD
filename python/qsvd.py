"""Quaternion SVD on dense complex representations. Requires NumPy."""
import operator

import numpy as np


def quaternion_to_complex(b1, b2):
    """Encode b1 + b2*j as [[b1,b2],[-conj(b2),conj(b1)]]."""
    b1, b2 = (np.asarray(b, dtype=np.complex128) for b in (b1, b2))
    if b1.ndim != 2 or b1.shape != b2.shape:
        raise ValueError("Expected two equally sized matrices")
    return np.block([[b1, b2], [-b2.conj(), b1.conj()]])


def _partner(x):
    top, bottom = np.split(x.conj(), 2, axis=0)
    return np.concatenate((-bottom, top), axis=0)


def _pair(x):
    return np.column_stack((x, _partner(x)))


def _scaled(ac):
    a = np.asarray(ac, dtype=np.complex128)
    if a.ndim != 2 or min(a.shape) == 0 or any(n % 2 for n in a.shape):
        raise ValueError("Expected a nonempty 2m by 2n complex representation")
    if not np.isfinite(a).all():
        raise ValueError("Input contains nonfinite entries")
    scale = max(abs(a.real).max(), abs(a.imag).max())
    if scale:
        a = a.real / scale + 1j * (a.imag / scale)
    m, n = (d // 2 for d in a.shape)
    if not (np.allclose(a[m:, n:], a[:m, :n].conj(), atol=1e-14, rtol=1e-13)
            and np.allclose(a[m:, :n], -a[:m, n:].conj(), atol=1e-14, rtol=1e-13)):
        raise ValueError("Input does not have quaternion block structure")
    return a, scale


def _extend(basis, candidates, count):
    """Select independent quaternion directions with two projection passes."""
    selected = []
    q = _pair(basis)
    for x in candidates.T:
        y = x.copy()
        for _ in range(2):
            y -= q @ (q.conj().T @ y)
        norm = np.linalg.norm(y)
        if norm > 1e-12 * np.linalg.norm(x):
            y /= norm
            selected.append(y)
            q = np.column_stack((q, y, _partner(y)))
            if len(selected) == count:
                break
    return np.column_stack(selected) if selected else basis[:, :0]


def _dense(a):
    """Complex economy SVD with joint left/right rotation in each cluster."""
    left, values, vh = np.linalg.svd(a, full_matrices=False)
    right = vh.conj().T
    u, v = left[:, :0], right[:, :0]
    s = values[::2].copy()
    # Choose representative vectors before adjoining their exact partners.
    # In a repeated positive subspace U and V need the same rotation.
    band = 64 * np.finfo(float).eps * values[0]
    start = 0
    while start < len(s):
        stop = start + 1
        while stop < len(s) and s[start] - s[stop] <= band:
            stop += 1
        cols = slice(2 * start, 2 * stop)
        count = stop - start
        vg = right[:, cols]
        new_v = _extend(v, vg, count)
        if s[start] <= band:
            new_u = _extend(u, left[:, cols], count)
        else:
            # Include the component removed by orthogonalization against
            # earlier pairs (relevant for almost coincident singular values).
            mapped = left @ (right.conj().T @ new_v)
            new_u = _extend(u, mapped, count)
        if new_u.shape[1] != count or new_v.shape[1] != count:
            raise np.linalg.LinAlgError("Quaternion pair recovery failed")
        u, v = np.column_stack((u, new_u)), np.column_stack((v, new_v))
        start = stop
    return _pair(u), s, _pair(v)


def _zero(a, k):
    def identity_columns(size):
        return _pair(np.eye(size, k, dtype=complex))
    return identity_columns(a.shape[0]), np.zeros(k), identity_columns(a.shape[1])


def qsvd(ac, *, values_only=False):
    """All quaternion singular values or economy factors U, s, V."""
    a, scale = _scaled(ac)
    p = min(a.shape) // 2
    if values_only:
        return np.linalg.svd(a, compute_uv=False)[::2] * scale if scale else np.zeros(p)
    u, s, v = _dense(a) if scale else _zero(a, p)
    return u, s * scale, v


def qsvds(ac, k=1, *, values_only=False, tol=1e-10, maxiter=1000, seed=0):
    """Leading k quaternion triplets via restarted paired block Krylov iteration.

    Returns U, s, V, or just s. Raises LinAlgError on nonconvergence.
    For k < min(m,n), only thin projected matrices undergo dense SVD.
    """
    a, scale = _scaled(ac)
    p = min(a.shape) // 2
    if isinstance(k, (bool, np.bool_)) or not 1 <= operator.index(k) <= p:
        raise ValueError("k must be an integer between 1 and min(m,n)")
    if not np.isfinite(tol) or tol <= 0 or operator.index(maxiter) < 1:
        raise ValueError("Expected positive tol and maxiter")
    if k == p:
        return qsvd(ac, values_only=values_only)
    if not scale:
        result = _zero(a, k)
        return result[1] if values_only else result
    rng = np.random.default_rng(seed)
    n = a.shape[1]
    width = min(p - 1, k + max(4, min(k, 20)))
    capacity = min(p - 1, max(60, 12*k + 20))
    def random_block(cols):
        return rng.normal(size=(n, cols)) + 1j * rng.normal(size=(n, cols))
    empty = np.empty((n, 0), dtype=complex)
    w = _extend(empty, random_block(width), width)
    fresh = w
    ah = a.conj().T
    for _ in range(maxiter):
        r = w.shape[1]
        u, s, rot = _dense(_pair(a @ w))
        choose = np.r_[0:k, r:r+k]
        uk, vk = u[:, choose], _pair(w) @ rot[:, choose]
        s = s[:k]
        both = np.tile(s, 2)
        res = np.hypot(np.linalg.norm(a @ vk - uk * both, axis=0),
                       np.linalg.norm(ah @ uk - vk * both, axis=0))
        orth = max(np.linalg.norm(f.conj().T @ f - np.eye(2*k)) for f in (uk, vk))
        if res.max() <= tol * s[0] and orth <= max(tol, 1e-12)*10:
            return s * scale if values_only else (uk, s * scale, vk)
        if r == capacity:
            w = _pair(w) @ rot[:, :width]
            fresh = w
            if width == capacity:
                # With no room to expand, perform a two-sided power step.
                w = _extend(empty, ah @ (a @ w), width)
                if w.shape[1] < width:
                    w = np.column_stack((w, _extend(w, random_block(width), width-w.shape[1])))
                fresh = w
                continue
        count = min(width, capacity-w.shape[1])
        fresh = _extend(w, ah @ (a @ fresh), count)
        if not fresh.shape[1]:
            fresh = _extend(w, random_block(count), count)
        w = np.column_stack((w, fresh))
    raise np.linalg.LinAlgError(f"QSVDS did not converge in {maxiter} iterations")
