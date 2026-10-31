function [U,E,V,info] = QSVDS(AC,m,n,k,opts)
%QSVDS Largest K quaternion singular triplets with checked structured factors.
%   [U,E,V,INFO]=QSVDS(AC,M,N,K,OPTS), S=QSVDS(...).
%   AC is [B1,B2;-conj(B2),conj(B1)], M,N are quaternion dimensions.
%   Default backend 'svds' uses MATLAB SVDS, then recovers quaternion pairs
%   and refines the factors on a thin projection using the original QSVD.
%   Set opts.backend='krylov' for a self-contained thick-restarted block
%   Krylov solver that calls neither SVDS nor EIGS. 'subspace' retains the
%   previous two-sided power/subspace iteration for comparison.
%   Options: tol=1e-10, maxit=1000, seed=0, cluster_tol=1e-12,
%   oversampling=min(20,max(4,K)), maxspace=max(60,12*K+20).
%   maxspace/oversampling apply to Krylov (oversampling also to subspace).
%   INFO.used_builtin_svds records actual use, including failed attempts.
%   INFO.flag=0 requires both singular equations AND orthogonality to pass.
%   No partial backend falls back to a full-input SVD. An explicit K=P
%   request delegates to full QSVD. QSVDS preserves the caller's RNG state.
if nargin<4 || isempty(k), k=1; end
if nargin<5 || isempty(opts), opts=struct(); end
p=qsvd_validate_inputs(AC,m,n,k);
validateattributes(opts,{'struct'},{'scalar'},mfilename,'opts');
d=struct('backend','svds','tol',1e-10,'maxit',1000,'seed',0,'cluster_tol',1e-12, ...
 'oversampling',min(20,max(4,k)),'maxspace',max(60,12*k+20));
for name=fieldnames(opts).'
    if ~isfield(d,name{1}), error('qsvds:UnknownOption','Unknown option: %s',name{1}); end
    d.(name{1})=opts.(name{1});
end
opts=d;
opts.backend=validatestring(opts.backend,{'svds','krylov','subspace'},mfilename,'backend');
validateattributes(opts.tol,{'numeric'},{'real','scalar','finite','positive','<',1});
validateattributes(opts.maxit,{'numeric'},{'real','scalar','finite','integer','positive'});
validateattributes(opts.seed,{'numeric'},{'real','scalar','finite','integer','nonnegative','<=',2^32-1});
validateattributes(opts.cluster_tol,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(opts.oversampling,{'numeric'},{'real','scalar','finite','integer','nonnegative'});
validateattributes(opts.maxspace,{'numeric'},{'real','scalar','finite','integer','positive'});
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
used_native=false; reason=''; requested=opts.backend;
basic=rmfield(opts,{'backend','maxspace'});
if ~any(AC(:)) || k==p
    [U,E,V,info]=qsvds_subspace(AC,m,n,k,basic);
elseif strcmp(requested,'subspace')
    [U,E,V,info]=qsvds_subspace(AC,m,n,k,basic);
elseif strcmp(requested,'krylov')
    [U,E,V,info]=qsvds_krylov(AC,m,n,k,rmfield(opts,'backend'));
else
    used_native=true;
    native_opts=rmfield(opts,{'backend','maxspace','oversampling'});
    try
        [U,E,V,info]=qsvds_native(AC,m,n,k,native_opts);
        if ~info.converged, reason='SVDS/recovered factors did not meet the requested tolerance.'; end
    catch err
        reason=['SVDS backend failed: ',err.identifier];
    end
    if ~isempty(reason)
        [U,E,V,info]=qsvds_krylov(AC,m,n,k,rmfield(opts,'backend'));
    end
end
% Recheck the public result rather than trusting the backend status alone.
alpha=norm(AC,'fro');
if alpha>0
    U=qsvds_repair_left_null(U,real(diag(E(1:k,1:k))),m);
    V=qsvds_repair_left_null(V,real(diag(E(1:k,1:k))),n);
    A=AC/alpha; D=E(1:k,1:k)/alpha;
    R1=A*V(:,1:k)-U(:,1:k)*D; R2=A'*U(:,1:k)-V(:,1:k)*D';
    info.relative_residuals=sqrt(sum(abs(R1).^2,1)+sum(abs(R2).^2,1)).'/max(real(D(1,1)),eps);
end
info.orthogonality=max(norm(U'*U-eye(2*k),'fro'),norm(V'*V-eye(2*k),'fro'));
info.converged=info.converged && all(isfinite(info.relative_residuals)) ...
    && all(info.relative_residuals<=opts.tol) && info.orthogonality<=max(1e-10,10*opts.tol);
info.flag=double(~info.converged);
info.requested_backend=requested; info.used_builtin_svds=used_native; info.fallback_reason=reason;
if info.flag~=0
    warning('qsvds:NoConvergence','QSVDS factors did not pass tolerance; residual %.3e, orthogonality %.3e.',max(info.relative_residuals),info.orthogonality);
end
if nargout<=1, U=diag(E(1:k,1:k)); end
end
