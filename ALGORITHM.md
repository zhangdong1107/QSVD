# QSVD and selectable QSVDS backends

The original full complex-representation algorithm remains in qsvd_complex.m and QSVD.m. QSVDS is a checked partial-decomposition interface. Its default backend now uses MATLAB svds, as authorized by the owner after performance comparison. The optional Krylov backend is self-contained and never calls svds/eigs. This software extension is not a claim that the original QSVD paper already presented these partial algorithms.

## Common representation and acceptance checks

Input is A=[B1 B2; -conj(B2) conj(B1)], of size 2m by 2n. Let K_d(x)=[-conj(x_bottom);conj(x_top)]. K_d^2=-I and A K_n(x)=K_m(Ax). A quaternion basis is represented by [W,K(W)]. k counts quaternion singular values, with genuine multiplicities preserved.

All backends scale by alpha=norm(A,'fro') and restore singular-value scale on return. Both singular equations are checked on the original scaled operator. For each representative triplet, require

sqrt(norm(A_s*v_j-s_j*u_j)^2+norm(A_s'*u_j-s_j*v_j)^2)/max(s_1,eps) <= tol.

The public interface additionally checks finiteness and orthogonality of both full paired factors. info.flag=0 requires residual convergence and Frobenius orthogonality error <= max(1e-10,10*tol). Failure is warned and reported; approximate output is not presented as successful convergence.

Thin complex SVD nullspaces need not be K-invariant. For numerically zero modes, qsvds_repair_left_null completes representative left vectors orthogonally to all positive modes and their K companions. The actual A'*U residual is checked afterwards, so completion alone is not assumed to prove a null vector is valid.

## Default backend: svds with structured recovery

1. Request 2k complex triplets from MATLAB svds, internally using max(eps,tol/10) and the specified iteration budget.
2. Sort the right vectors by singular value. Twice reorthogonalize candidates against both previously accepted representatives and K companions. Skip numerically dependent candidates and retain k independent quaternion directions.
3. Form the thin matrix Y=A_s W and its K companions. Apply the original QSVD to this 2m by 2k projection, not to the full original input. Lift its right factors through [W,K(W)].
4. Complete zero left modes when required, then test both singular equations and orthogonality. The wrapper can attempt the pure Krylov backend if the native backend raises an error or fails its residual check. The final public flag includes orthogonality as well.

This recovery/refinement cost is included in QSVDS timings. A raw complex svds timing does not include it and is labeled separately. info.used_builtin_svds=true records a native attempt even if the final result came from Krylov fallback. info.method identifies the actual solver; info.requested_backend identifies the requested one. MATLAB may make its own internal algorithm decisions; the wrapper makes no claim that its native backend is built-in-free.

## Pure backend: thick-restarted paired block Krylov

With k<p=min(m,n), retain ell=min(p-1,k+oversampling) quaternion directions. The search capacity is min(p-1,max(ell,maxspace)). Construct representatives with twice-reorthogonalized Gram-Schmidt against all stored representatives and their K companions.

Each expansion applies the operator A_s'*(A_s*W_last) without forming A_s'*A_s. Keep the older search blocks rather than discarding them after each power sweep. Cache A_s times the representatives. At regular checkpoints, apply the original QSVD to the tall thin projected matrix, extract the leading Ritz triplets and test their actual residuals. When capacity is reached, restart with the leading ell Ritz directions.

The full SVD inside the original core is limited to the thin projection: its minimum complex dimension is at most 2*capacity < 2*p. No full-input SVD is used for a partial pure solve. Twice reorthogonalization and K completion maintain the quaternion structure. Rank breakdown can be completed with new random directions. If capacity=ell (k near p), no expansion fits below p; the implementation uses the paired two-sided subspace iteration so that the retained search space continues changing instead of restarting to an unchanged basis.

Default settings: oversampling=min(20,max(4,k)), maxspace=max(60,12*k+20), maxit=1000, tol=1e-10. Increasing capacity costs more memory and reduced decomposition work. Iteration counts measure block steps and are not directly comparable to native svds internal iterations. Slow/clustered spectra can still be expensive; this backend has no universal speedup guarantee.

## Compatibility and limits

- backend='subspace' retains the previous two-sided power method for comparison. It may fail on flat spectra.
- k=p explicitly requests the full result; QSVD is used and full_decomposition=true. This is not a hidden fallback after partial failure.
- Exact zero input has explicit compact structured factors. Finite, dense, double complex representations are supported; sparse/operator input and smallest/interior values are not implemented.
- Every public call saves and restores the caller RNG. seed controls internal randomness.
- Random initialization and small triplet residuals are not a deterministic proof that no larger singular value was missed in every possible matrix. Independent reference spectra, multiple seeds, orthogonality and residual tests are used in the bundled validation.
- Instrumented pure-backend tests replace only function names and SVD call targets to audit projection dimensions; they also compare instrumented and ordinary outputs. This avoids relying on path shadowing for MATLAB's cached built-in SVD binding.

## Sources and attribution

Original full QSVD: https://doi.org/10.1016/j.aml.2026.110047.

MATLAB svds documentation: https://www.mathworks.com/help/matlab/ref/svds.html.

Background on restarted block Krylov/Lanczos: J. Baglama and L. Reichel, Restarted block Lanczos bidiagonalization methods, https://digitalcommons.uri.edu/math_facpubs/79/; pure implementation here is not a verbatim implementation of their published algorithm. General quaternion randomized subspace context: https://arxiv.org/abs/2011.03368. These references motivate the framework, not an unproved all-input accuracy or speed claim for this code.
