function [UC, EC, VC, info] = qsvds_subspace(AC,m,n,k,opts)
%QSVDS Iterative largest-K quaternion singular triplets, without SVDS/EIGS.
%   [UC,EC,VC,INFO] = qsvds_subspace(AC,M,N,K,OPTS) uses structure-preserving
%   two-sided block subspace iteration and the original QSVD on a small core.
%   S = QSVDS(...) returns a K-by-1 vector of largest singular values.
%   K defaults to 1 and counts quaternion singular values, not complex pairs.
%
%   OPTS fields: tol (residual tolerance, 1e-10), maxit (1000),
%   oversampling (min(20,max(4,K))), seed (0), cluster_tol (1e-12).
%   With K < P=min(M,N), the working quaternion dimension is
%   L=min(P-1,K+oversampling). Neither the input nor a P-by-P core is sent
%   to full SVD. SVD is used only inside the original QSVD on L-by-L cores.
%   K=P explicitly requests the full decomposition and delegates to QSVD.
%
%   UC: 2*M-by-2*K, EC: 2*K-by-2*K, VC: 2*N-by-2*K.
%   INFO.flag=0 means the two singular equations met the residual tolerance;
%   flag=1 means the iteration limit was reached (a warning is also issued).
%   There is NO fallback to a full input SVD after nonconvergence.
%   Random initialization targets the dominant subspace; residual convergence
%   alone is not a deterministic certificate that no larger value was missed.
%   Closely spaced spectra may converge slowly; no speedup is guaranteed.

if nargin < 4 || isempty(k), k = 1; end
if nargin < 5 || isempty(opts), opts = struct(); end
p = qsvd_validate_inputs(AC,m,n,k);
validateattributes(opts,{'struct'},{'scalar'},mfilename,'opts');
defaults = struct('tol',1e-10,'maxit',1000,'oversampling',min(20,max(4,k)), ...
                  'seed',0,'cluster_tol',1e-12);
names = fieldnames(opts);
for j = 1:numel(names)
    if ~isfield(defaults,names{j})
        error('qsvds:UnknownOption','Unknown QSVDS option: %s',names{j});
    end
    defaults.(names{j}) = opts.(names{j});
end
opts = defaults;
validateattributes(opts.tol,{'numeric'},{'real','scalar','finite','positive','<',1});
validateattributes(opts.cluster_tol,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(opts.maxit,{'numeric'},{'real','scalar','finite','integer','positive'});
validateattributes(opts.oversampling,{'numeric'},{'real','scalar','finite','integer','nonnegative'});
validateattributes(opts.seed,{'numeric'},{'real','scalar','finite','integer','nonnegative','<=',2^32-1});
state = rng; restore_rng = onCleanup(@() rng(state)); %#ok<NASGU>
rng(opts.seed,'twister');
info = struct('method','structured-subspace','converged',true,'flag',0, ...
    'iterations',0,'relative_residuals',zeros(k,1),'subspace_dimension',0, ...
    'largest_projected_svd',0,'full_decomposition',false,'fallback_reason','');
alpha = norm(AC,'fro');
if alpha == 0
    UC = quaternion_to_complex(eye(m,k),zeros(m,k));
    VC = quaternion_to_complex(eye(n,k),zeros(n,k));
    EC = zeros(2*k);
    info.method = 'zero';
elseif k == p
    [UC,EC,VC] = QSVD(AC,m,n,opts.cluster_tol);
    info.method = 'full-requested';
    info.full_decomposition = true;
    info.subspace_dimension = p;
    info.largest_projected_svd = p;
    A = AC/alpha; D = EC(1:k,1:k)/alpha;
    Rleft = A*VC(:,1:k)-UC(:,1:k)*D;
    Rright = A'*UC(:,1:k)-VC(:,1:k)*D;
    info.relative_residuals = sqrt(sum(abs(Rleft).^2,1)+sum(abs(Rright).^2,1)).'/max(D(1,1),eps);
    info.converged = all(info.relative_residuals <= opts.tol);
    info.flag = double(~info.converged);
else
    ell = min(p-1,k+opts.oversampling);
    info.subspace_dimension = ell;
    info.largest_projected_svd = ell;
    A = AC/alpha;
    Q = paired_orth(randn(2*n,ell)+1i*randn(2*n,ell),n,ell);
    info.converged = false;
    for iteration = 1:opts.maxit
        % Orthonormalize after EACH multiplication; never form A'*A.
        P = paired_orth(A*Q(:,1:ell),m,ell);
        Q = paired_orth(A'*P(:,1:ell),n,ell);
        B = P'*A*Q;
        % Remove only roundoff-level drift from the quaternion block pattern.
        B = quaternion_to_complex(B(1:ell,1:ell),B(1:ell,ell+1:2*ell));
        % Reuse the author's full algorithm only on the reduced core.
        [Ub,Eb,Vb] = qsvd_complex(B,ell,ell,opts.cluster_tol);
        columns = [1:k,ell+(1:k)];
        UC = P*Ub(:,columns);
        VC = Q*Vb(:,columns);
        s = real(diag(Eb(1:k,1:k)));
        D = diag(s);
        Rleft = A*VC(:,1:k)-UC(:,1:k)*D;
        Rright = A'*UC(:,1:k)-VC(:,1:k)*D;
        residuals = sqrt(sum(abs(Rleft).^2,1)+sum(abs(Rright).^2,1)).';
        residuals = residuals/max(s(1),eps);
        info.iterations = iteration;
        info.relative_residuals = residuals;
        if all(residuals <= opts.tol)
            info.converged = true;
            break;
        end
    end
    EC = quaternion_to_complex(diag(s*alpha),zeros(k));
    if ~info.converged
        info.flag = 1;
    end
end
if ~info.converged
    warning('qsvds:NoConvergence', ...
        'QSVDS stopped after %d iterations; max relative residual %.3e exceeds %.3e. Inspect INFO.', ...
        info.iterations,max(info.relative_residuals),opts.tol);
end
if nargout <= 1, UC = diag(EC(1:k,1:k)); end
end

function Q = paired_orth(X,d,r)
% K(x)=[-conj(x_bottom);conj(x_top)]. Paired, twice-reorthogonalized
% Gram-Schmidt keeps span(Q) K-invariant, including repeated-value spaces.
W = complex(zeros(2*d,r));
for j = 1:r
    accepted = false;
    for attempt = 1:20
        if attempt == 1, w = X(:,j); else, w = randn(2*d,1)+1i*randn(2*d,1); end
        if j > 1
            previous = W(:,1:j-1);
            G = [previous,[-conj(previous(d+1:2*d,:));conj(previous(1:d,:))]];
            for pass = 1:2, w = w-G*(G'*w); end
        end
        nw = norm(w);
        if nw > 1e-14
            W(:,j) = w/nw;
            accepted = true;
            break;
        end
    end
    if ~accepted, error('qsvds:BasisBreakdown','Unable to complete a paired orthonormal basis.'); end
end
Q = [W,[-conj(W(d+1:2*d,:));conj(W(1:d,:))]];
end
