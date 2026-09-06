function [UC, EC, VC] = qsvd_complex(AC, m, n, tol)
%QSVD_COMPLEX Quaternion SVD through its complex representation.
%   [UC, EC, VC] = QSVD_COMPLEX(AC, M, N) uses TOL = 1e-12.
%   [UC, EC, VC] = QSVD_COMPLEX(AC, M, N, TOL) sets the relative
%   singular-value clustering tolerance (applied after scaling AC).
%
%   AC = [B1, B2; -conj(B2), conj(B1)] is a 2*M-by-2*N complex matrix.
%   Construct AC with QUATERNION_TO_COMPLEX(B1, B2).
%   M and N are the ORIGINAL quaternion matrix dimensions.
%
%   For nonzero input, P = min(M,N):
%       UC is 2*M-by-2*P; EC is 2*P-by-2*P; VC is 2*N-by-2*P.
%   The exact-zero branch returns full factors:
%       UC is 2*M-by-2*M; EC is 2*M-by-2*N; VC is 2*N-by-2*N.
%   In both cases, AC is approximately UC * EC * VC'.
%   The prime denotes the complex conjugate transpose.
%
%   The first P diagonal entries of EC are the quaternion singular values.
%   Repeated singular values may use random bases; use RNG for reproducibility.
%   See README.md and DEMO_QSVD for examples and interpretation.
%
%   Core numerical operations preserved from the author-supplied implementation.
%   Paper: https://doi.org/10.1016/j.aml.2026.110047

    if nargin < 4 || isempty(tol)
        tol = 1e-12;
    end

    % Scale to improve numerical stability
    alpha = norm(AC,'fro');
    if alpha == 0
        % Zero matrix: return trivial factors
        UC = eye(2*m);
        VC = eye(2*n);
        EC = zeros(2*m,2*n);
        return;
    end
    ACs = AC / alpha;

    % Economy SVD of complex representation
    [U,S,V] = svd(ACs,'econ');
    s = diag(S);

    % Working copies
    U2 = U;
   % size(U2)
    V2 = V;
    S2 = S * alpha;

    % Anti-linear "K" operators: K(x) = J * conj(x)
    Jm = [zeros(m), -eye(m); eye(m), zeros(m)];   % 2m-by-2m
    Jn = [zeros(n), -eye(n); eye(n), zeros(n)];   % 2n-by-2n
    Km = @(x) Jm * conj(x);
    Kn = @(x) Jn * conj(x);

    % Zero threshold (deciding the 0-singular-value cluster)
    smax = max(s);
    zero_thr = 100 * eps(max(1,smax)) * max(1,smax);

    N = length(s);
    used = false(N,1);

    % ---- helper: polar projection to nearest unitary ----
    function T = polar_unitary(M)
        % M is k-by-k
        [Um,~,Vm] = svd(M);
        T = Um*Vm';
    end

    % ---- helper: build K-structured orthonormal basis inside span(G) ----
    function Q = build_structured_basis_pair(G, K)
        % Build a K-structured orthonormal basis Q inside span(G)
        % using pairs [w, proj_span(G) K(w)].
        k = size(G,2);
        Q = zeros(size(G));
        qcol = 0;

        % Projector onto span(G)
        % (G is assumed orthonormal in our usage: columns from U/V)
        % But we use G*(G'*x) anyway.
        while qcol < k
            c = randn(k,1) + 1i*randn(k,1);
            w = G*c;

            % orthogonalize against current Q
            if qcol > 0
                w = w - Q(:,1:qcol)*(Q(:,1:qcol)'*w);
            end
            nw = norm(w);
            if nw < 1e-14
                continue;
            end
            w = w/nw;

            % companion via anti-linear map, projected back
            wK = G*(G' * K(w));

            B = [w, wK];
            if qcol > 0
                B = B - Q(:,1:qcol)*(Q(:,1:qcol)'*B);
            end

            [Bq,~] = qr(B,0);  % complex QR, columns orthonormal

            take = min(2, k-qcol);
            Q(:,qcol+1:qcol+take) = Bq(:,1:take);
            qcol = qcol + take;
        end
    end

    % ---- Main loop over singular value clusters ----
    % IMPORTANT update:
    %   - Because singular values of AC come in (at least) pairs,
    %     the "true" problematic case for quaternion QSVD extraction is when
    %     a quaternion singular value repeats, i.e., AC has >=4 equal singular values.
    %   - Thus, we ONLY post-process clusters with k >= 4.
    for p = 1:N
        if used(p), continue; end

        idx = find(abs(s - s(p)) <= tol*max(1,abs(s(p))));
        used(idx) = true;
        k = numel(idx);

        % Only handle true repetition (>=4 in complex representation)
        if k < 4
            continue;
        end

        Ug = U(:,idx);
        Vg = V(:,idx);

        is_zero_cluster = (max(abs(s(idx))) <= zero_thr);

        if ~is_zero_cluster
            % nonzero cluster: use the SAME unitary rotation for U and V
            Qu = build_structured_basis_pair(Ug, Km);
            T  = polar_unitary(Ug' * Qu);
            U2(:,idx) = Ug * T;
            V2(:,idx) = Vg * T;
        else
            % zero cluster (>=4 zeros): left and right nullspaces are independent
            Qu = build_structured_basis_pair(Ug, Km);
            Tu = polar_unitary(Ug' * Qu);
            U2(:,idx) = Ug * Tu;

            Qv = build_structured_basis_pair(Vg, Kn);
            Tv = polar_unitary(Vg' * Qv);
            V2(:,idx) = Vg * Tv;
        end
    end

    % ---- Extract quaternion factors from (possibly) post-processed U2,V2 ----
    % Select one vector per K-pair (odd indices 1,3,5,...)
    p = min(m,n);
    u = U2(:,1:2:2*p);
    v = V2(:,1:2:2*p);
    e = S2(1:2:2*p,1:2:2*p);

    % Form complex blocks for Cayley--Dickson recombination
    u1 = u(1:m,:);
    u2 = -conj(u(m+1:2*m,:));
    v1 = v(1:n,:);
    v2 = -conj(v(n+1:2*n,:));

    % Assemble Chi(Uq), Chi(Vq), Chi(E)
    function X = assemble_complex_blocks(B1,B2)
        X = [B1, B2;
             -conj(B2), conj(B1)];
    end

    UC = assemble_complex_blocks(u1,u2);
    VC = assemble_complex_blocks(v1,v2);

    % Quaternion diagonal E corresponds to Chi(E) = [E 0; 0 E]
    % Here e is p-by-p (real nonnegative on diagonal, numerically),
    % so we place it and zeros.
    EC = assemble_complex_blocks(e, zeros(p,p));
end
