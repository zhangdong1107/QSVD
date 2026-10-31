# QSVD / QSVDS: Quaternion Singular Value Decomposition in MATLAB

[Algorithm project page](https://zhangdong1107.github.io/qsvd/) · [Download the MATLAB package](https://zhangdong1107.github.io/downloads/qsvd-aml-matlab.zip?v=20261031)

MATLAB code accompanying **A fast universal algorithm for quaternion singular value decomposition**, by Dong Zhang, Chuan Jiang, V. I. Vasil'ev, and Tongsong Jiang. *Applied Mathematics Letters*, 182 (2026), 110047.

Paper: https://doi.org/10.1016/j.aml.2026.110047

This package computes a quaternion singular value decomposition using ordinary complex MATLAB matrices. You do not need to install a quaternion class or a third-party toolbox.

## Requirements and validation

Validated with MATLAB R2025a on Windows. Full-QSVD, compatibility, structured-spectrum, independent real-representation, backend and failure-recovery suites are bundled. Their measured results and the large-matrix benchmark are documented in BENCHMARK.txt. These checks do not prove correctness or speed for every matrix. No quaternion class or external toolbox is required; GNU Octave has not been tested.

## QSVD and QSVDS: selectable backends (31 October 2026)

The original full algorithm remains available as `QSVD`. QSVDS now defaults to MATLAB `svds`, followed by quaternion pair recovery, thin QSVD refinement, zero-mode basis completion and checks of both singular equations and orthogonality. This choice follows the owner's request to use `svds` if it is faster after benchmarking. The software no longer claims that the default QSVDS is a built-in-free solver.

```matlab
[U,E,V] = QSVD(AC,m,n);                    % Original full decomposition.
[Uk,Ek,Vk,info] = QSVDS(AC,m,n,k);          % Default: SVDS + structured recovery.
opts = struct('backend','krylov');
[Uk,Ek,Vk,info] = QSVDS(AC,m,n,k,opts);     % Pure restarted block Krylov.
s = QSVDS(AC,m,n,k);                       % k descending quaternion values.
sigma_max = QSVDS(AC,m,n,1);               % Largest value.
```

- `backend='svds'` (default): request 2k complex triplets, recover k independent quaternion directions, then use the original QSVD on a thin projection. Recovery/refinement cost is part of the backend. The internal SVDS tolerance is `max(eps,tol/10)`. It does not simply reinterpret unstructured complex factors as quaternion factors.
- `backend='krylov'`: self-contained paired, thick-restarted block Krylov iteration, with twice-reorthogonalized bases and thin QSVD Ritz extraction. No `svds` or `eigs` call occurs in this backend. It retains multiple search blocks, unlike the previous power iteration. When k is near the full dimension and no expansion fits, paired two-sided iteration is used.
- `backend='subspace'`: retain the previous block power/subspace solver for reproducibility; flat spectra may fail to converge within 1000 iterations.

If the default SVDS backend fails or its recovered factors do not pass its residual check, the wrapper attempts the pure Krylov backend and records `info.fallback_reason`. It does not fall back to a full-input decomposition. `info.used_builtin_svds` reports actual use, including a failed native attempt. For an explicit k=min(m,n) request, full QSVD is used and recorded as `info.full_decomposition=true`. Zero inputs use explicit compact structured factors.

The outputs are complex representations with sizes `2*m` by `2*k`, `2*k` by `2*k`, and `2*n` by `2*k`. k counts quaternion singular values; true quaternion multiplicities are retained. `Uk*Ek*Vk'` is a rank-k approximation. Inspect `info.flag`, `info.relative_residuals` and `info.orthogonality` before using results. A failed tolerance check warns and returns flag=1; it is not claimed as successful convergence.

```matlab
opts = struct('backend','krylov','tol',1e-10,'maxit',1000, ...
              'oversampling',8,'maxspace',80,'seed',0);
[Uk,Ek,Vk,info] = QSVDS(AC,m,n,k,opts);
s = diag(Ek(1:k,1:k));
```

Options: `tol` (default 1e-10), `maxit` (1000), `seed` (0), `cluster_tol` (1e-12), `oversampling` (`min(20,max(4,k))`), and `maxspace` (`max(60,12*k+20)`). The latter two control the pure backend; dimensions are capped below min(m,n) for partial solves. maxspace is at least the retained block size. Krylov `iterations` counts block iterations; SVDS's internal iteration count is not exposed, so it is reported as NaN. The caller's RNG state is restored.

The compatibility helpers `qsvd_values(AC,m,n,k)` and `qsvd_complex(AC,m,n,[],k)` dispatch through the same QSVDS wrapper. The fourth argument of `qsvd_complex` remains the original clustering tolerance. Without k, its original output conventions are retained. New uppercase QSVD/QSVDS interfaces return values when called with one output.

See ALGORITHM.md for details and BENCHMARK.txt for the measured speed/accuracy comparison. Speed depends on the spectrum and machine. In particular, a converged partial solve is not guaranteed to beat full QSVD on a flat spectrum.

## Start here: run the example

1. Extract `QSVD-AML-MATLAB.zip`.
2. In MATLAB, set **Current Folder** to the extracted `QSVD-AML-MATLAB` folder, where the `.m` files are located.
3. Enter the following commands in the Command Window:

```matlab
demo_qsvd
test_qsvd
test_qsvd_partial
test_qsvds
test_qsvds_real_reference
test_qsvds_backends
test_qsvds_real_reference('krylov')
test_qsvds_fallback
```

The demo prints the full and largest singular values and reconstruction/orthogonality errors. The tests report full-QSVD and top-K checks. Small numerical differences between computers are normal.

If you prefer to work in another folder, first add this package to your path:

```matlab
addpath('C:/your/path/QSVD-AML-MATLAB'); % Replace with your extracted folder.
demo_qsvd
```

## Files

| File | Purpose |
| --- | --- |
| `QSVD.m` | Full decomposition interface; original core. |
| `QSVDS.m` | Validated public interface and backend selection. |
| `qsvds_native.m` | SVDS, paired recovery and thin refinement. |
| `qsvds_krylov.m` | Self-contained restarted paired block Krylov. |
| `qsvds_subspace.m` | Previous subspace solver for comparison. |
| `qsvds_repair_left_null.m` | Paired left-zero-mode completion for rank-deficient thin factors. |
| `test_qsvds_backends.m` | Backend, rank deficiency, options and flat-spectrum checks. |
| `test_qsvds_fallback.m` | Injected native failure and invalid-factor recovery checks. |
| `benchmark_qsvds_optimized.m` | Reproducible speed/accuracy benchmark, including factor recovery cost. |
| `BENCHMARK.txt` | Measured backend speed and precision results. |
| `qsvd_complex.m` | Original full core; fifth-argument calls dispatch to QSVDS. |
| `qsvd_values.m` | Compatibility values helper, delegated to QSVDS. |
| `qsvd_validate_inputs.m` | Shared input and representation validation. |
| `quaternion_to_complex.m` | Construct the complex representation from two complex blocks. |
| `demo_qsvd.m` | Complete, reproducible example that needs no external data. |
| `test_qsvd.m` | Eight numerical checks, including rectangular and zero matrices. |
| `test_qsvd_partial.m` | Top-K values, optimal approximation error, factor shapes, repetition, zero input, scales, and invalid inputs. |
| `test_qsvds.m` | Rotated spectra, convergence status, RNG and forbidden-built-in checks. |
| `test_qsvds_real_reference.m` | Independent real representation comparisons and an analytic rank-one value. |
| `ALGORITHM.md` | Iteration, residual criterion, scope and references. |
| `README.md` | This guide. |
| `CITATION.bib` | BibTeX entry for the accompanying paper. |

## Use your own matrix

Represent an `m`-by-`n` quaternion matrix as

```text
A = A0 + A1*i + A2*j + A3*k
```

Here `A0`, `A1`, `A2`, and `A3` are real MATLAB matrices of the same size. The symbols `i`, `j`, and `k` above are quaternion units; MATLAB's numeric `1i` below is the ordinary complex unit. Build two complex blocks:

```matlab
% Replace these four real matrices with your data.
A0 = [1 2; 0 3; 2 -1];
A1 = [1 0; 0 -1; 0 2];
A2 = [0 0; 1 2; 0 1];
A3 = [0 1; 0 0; -1 0];

B1 = A0 + 1i*A1;
B2 = A2 + 1i*A3;
[m,n] = size(B1);

AC = quaternion_to_complex(B1,B2);
rng(0,'twister');                  % Reproducible repeated-value processing.
[UC,EC,VC] = qsvd_complex(AC,m,n); % Default clustering tolerance: 1e-12.

p = min(m,n);
singular_values = diag(EC(1:p,1:p));
reconstruction_error = norm(AC-UC*EC*VC','fro') / max(norm(AC,'fro'),eps);
disp(singular_values);
disp(reconstruction_error);
```

If you already have complex blocks `B1` and `B2`, start with `[m,n] = size(B1)`. The exact representation expected by the routine is

```matlab
AC = [B1, B2; -conj(B2), conj(B1)];
```

Do not supply a different quaternion-to-complex convention without converting it first.

## Function reference

```matlab
[UC,EC,VC] = qsvd_complex(AC,m,n)
[UC,EC,VC] = qsvd_complex(AC,m,n,tol)
s = qsvd_values(AC,m,n,k)
s = qsvd_complex(AC,m,n,tol,k)
[UC,EC,VC] = qsvd_complex(AC,m,n,tol,k)
```

### Inputs

- `AC`: a finite, dense, double-precision complex representation matrix of size `2*m` by `2*n`, assembled as above.
- `m`, `n`: positive integer dimensions of the **original quaternion matrix**, not the doubled complex matrix.
- `tol`: optional positive singular-value clustering tolerance, default `1e-12`. Clustering is performed on singular values of `AC / norm(AC,'fro')`. Start with the default; increasing it can merge distinct nearby singular values and affect accuracy.

Inputs are validated. Empty arrays, NaN/Inf values, invalid dimensions/counts/tolerance, and nonconforming complex representations are rejected.

### Outputs and reconstruction

The outputs are **complex representations** of quaternion factors, not MATLAB quaternion objects. Always reconstruct with

```matlab
AC_reconstructed = UC * EC * VC';
```

`VC'` is the complex conjugate transpose. Do not replace it with `VC.'`.

For a nonzero input without K, let `p = min(m,n)`:

| Output | Size |
| --- | --- |
| `UC` | `2*m` by `2*p` |
| `EC` | `2*p` by `2*p` |
| `VC` | `2*n` by `2*p` |

This is the compact decomposition. Columns of `UC` and `VC` are orthonormal. The quaternion singular values are the first `p` diagonal entries of `EC`; the second diagonal block repeats them.

**Exact-zero special case:** the original implementation returns full identity factors: `UC` is `2*m` by `2*m`, `EC` is `2*m` by `2*n`, and `VC` is `2*n` by `2*n`. The same reconstruction expression still applies. Determine factor widths from their actual sizes when writing reusable code.

### Recover quaternion component matrices

For either output shape, extract the left factor's four real components as follows:

```matlab
pu = size(UC,2)/2;
U_B1 = UC(1:m,1:pu);
U_B2 = UC(1:m,pu+1:2*pu);
U0 = real(U_B1); U1 = imag(U_B1);
U2 = real(U_B2); U3 = imag(U_B2);
% Quaternion factor: U = U0 + U1*i + U2*j + U3*k.
```

For the right factor, use `VC`, `n`, and `pv = size(VC,2)/2` in the same pattern. Use the complex representation for multiplication unless you have separately implemented quaternion arithmetic.

## Reproducibility and checks

Repeated singular values allow multiple valid singular-vector bases. The routine can call `randn` internally when processing those subspaces. Set `rng(seed,'twister')` before a call when you need reproducible runs in the same environment. The demo and tests set the random seed and therefore change MATLAB's global random-number state.

Compare reconstruction errors, orthogonality, and singular values, rather than expecting singular vectors to match entry by entry across machines. `test_qsvd` checks random square/tall/wide matrices, repeated values, a rank-deficient matrix, tall/wide zero matrices, and a scalar quaternion. It also verifies output sizes and the quaternion block structure. These checks cover the included examples; they are not a proof for every possible input.

## Troubleshooting

- **Unrecognized function or variable:** change Current Folder to this package or run `addpath` as above. Check `which qsvd_complex -all` to rule out another copy on your path.
- **Matrix dimensions must agree:** `B1` and `B2` must have the same size. Derive `[m,n] = size(B1)` before constructing `AC`.
- **Unexpected duplicate singular values:** the complex representation duplicates each quaternion singular value. Use `diag(EC(1:p,1:p))` to obtain one copy.
- **Unexpected output dimensions for a zero matrix:** see the exact-zero special case above.
- **Large reconstruction error:** check the block signs/conjugation, `m` and `n`, and the use of `VC'`. Begin with the default tolerance and run `demo_qsvd`.

## Names and compatibility

The full-QSVD numerical core is preserved. The October 2026 update includes QSVD/QSVDS interfaces, selectable native and pure iterative backends, structured recovery and input validation. This distribution replaces `CC.m` with `quaternion_to_complex.m`, and `svd_quat_fix_complex_xin.m` with `qsvd_complex.m`. Update existing scripts to call these English names. No old-name aliases are included.

The QSVDS tests audit SVD call sites by making instrumented copies of the pure Krylov helper and full QSVD core, replacing only function names and the SVD call target. The instrumented and original outputs are compared. This avoids MATLAB's cached binding of built-ins defeating path-based shadowing. On the audited partial problem with complex minimum dimension 50, the largest intercepted dense SVD had dimension 12; a full-input SVD would have raised an error.

## Citation and contact

Please cite the accompanying paper when using this code in academic work; a BibTeX entry is included in `CITATION.bib`.

- Dong Zhang: zhangdong@qfnu.edu.cn
- Alternative email: dz_zhangdong@sina.com
