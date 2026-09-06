# QSVD-AML-MATLAB

MATLAB code accompanying **A fast universal algorithm for quaternion singular value decomposition**, by Dong Zhang, Chuan Jiang, V. I. Vasil'ev, and Tongsong Jiang. *Applied Mathematics Letters*, 182 (2026), 110047.

Paper: https://doi.org/10.1016/j.aml.2026.110047

This package computes a quaternion singular value decomposition using ordinary complex MATLAB matrices. You do not need to install a quaternion class or a third-party toolbox.

## Requirements and validation

Validated with MATLAB R2025a (25.1.0.2943329) on Windows. The package uses built-in MATLAB functions only. The demo and all eight bundled checks passed; the demo relative reconstruction error was approximately `1.14e-15`. Other MATLAB releases and GNU Octave have not been tested for this distribution.

## Start here: run the example

1. Extract `QSVD-AML-MATLAB.zip`.
2. In MATLAB, set **Current Folder** to the extracted `QSVD-AML-MATLAB` folder, where the `.m` files are located.
3. Enter the following commands in the Command Window:

```matlab
demo_qsvd
test_qsvd
```

The demo prints the singular values and reconstruction/orthogonality errors, then prints `Demo passed.` The checks finish with `All 8 QSVD checks passed.` Small numerical differences between computers are normal.

If you prefer to work in another folder, first add this package to your path:

```matlab
addpath('C:/your/path/QSVD-AML-MATLAB'); % Replace with your extracted folder.
demo_qsvd
```

## Files

| File | Purpose |
| --- | --- |
| `qsvd_complex.m` | Main QSVD routine. |
| `quaternion_to_complex.m` | Construct the complex representation from two complex blocks. |
| `demo_qsvd.m` | Complete, reproducible example that needs no external data. |
| `test_qsvd.m` | Eight numerical checks, including rectangular and zero matrices. |
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
```

### Inputs

- `AC`: a finite, dense, double-precision complex representation matrix of size `2*m` by `2*n`, assembled as above.
- `m`, `n`: positive integer dimensions of the **original quaternion matrix**, not the doubled complex matrix.
- `tol`: optional positive singular-value clustering tolerance, default `1e-12`. Clustering is performed on singular values of `AC / norm(AC,'fro')`. Start with the default; increasing it can merge distinct nearby singular values and affect accuracy.

The routine expects correctly shaped inputs. Empty arrays, NaN/Inf values, and nonconforming complex matrices are outside the documented input contract.

### Outputs and reconstruction

The outputs are **complex representations** of quaternion factors, not MATLAB quaternion objects. Always reconstruct with

```matlab
AC_reconstructed = UC * EC * VC';
```

`VC'` is the complex conjugate transpose. Do not replace it with `VC.'`.

For a nonzero input, let `p = min(m,n)`:

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

The supplied numerical operations are preserved. This distribution replaces `CC.m` with `quaternion_to_complex.m`, and `svd_quat_fix_complex_xin.m` with `qsvd_complex.m`. Update existing scripts to call these English names. The original header's dimensions have been corrected to describe the implementation's compact and exact-zero branches. No old-name aliases are included.

## Citation and contact

Please cite the accompanying paper when using this code in academic work; a BibTeX entry is included in `CITATION.bib`.

- Dong Zhang: zhangdong@qfnu.edu.cn
- Alternative email: dz_zhangdong@sina.com
