function [U,E,V,info]=qsvds_native(AC,m,n,k,opts)
% SVDS backend plus quaternion pair recovery and reduced QSVD refinement.
if nargin<5, opts=struct(); end
p=qsvd_validate_inputs(AC,m,n,k);
d=struct('tol',1e-10,'maxit',1000,'seed',0,'cluster_tol',1e-12);
for name=fieldnames(opts).', d.(name{1})=opts.(name{1}); end
opts=d;
if ~any(AC(:)) || k==p
    [U,E,V,info]=qsvds_subspace(AC,m,n,k); return;
end
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(opts.seed,'twister'); alpha=norm(AC,'fro'); A=AC/alpha;
[~,D,Z,flag]=svds(A,2*k,'largest',struct('tol',max(eps,opts.tol/10),'maxit',opts.maxit));
[~,order]=sort(real(diag(D)),'descend'); Z=Z(:,order);
W=complex(zeros(2*n,0));
for j=1:size(Z,2)
    x=Z(:,j); Q=[W,K(W,n)];
    for pass=1:2, x=x-Q*(Q'*x); end
    nx=norm(x);
    if nx>1e-7
        W(:,end+1)=x/nx; %#ok<AGROW>
        if size(W,2)==k, break; end
    end
end
if size(W,2)<k, error('qsvds:PairRecovery','SVDS did not supply k independent quaternion directions.'); end
Y=A*W; YC=[Y,K(Y,m)];
[U,E,Vsmall]=qsvd_complex(YC,m,k,opts.cluster_tol);
V=[W,K(W,n)]*Vsmall; s=real(diag(E(1:k,1:k)));
U=qsvds_repair_left_null(U,s,m);
R1=A*V(:,1:k)-U(:,1:k)*diag(s); R2=A'*U(:,1:k)-V(:,1:k)*diag(s);
res=sqrt(sum(abs(R1).^2,1)+sum(abs(R2).^2,1)).'/max(s(1),eps);
orth=max(norm(U'*U-eye(2*k),'fro'),norm(V'*V-eye(2*k),'fro'));
E=quaternion_to_complex(diag(s*alpha),zeros(k));
info=struct('method','matlab-svds','converged',flag==0 && all(res<=opts.tol) && orth<=max(1e-10,10*opts.tol), ...
 'flag',double(flag~=0 || any(res>opts.tol) || orth>max(1e-10,10*opts.tol)),'iterations',NaN,'relative_residuals',res, ...
 'subspace_dimension',k,'largest_projected_svd',k,'full_decomposition',false,'fallback_reason','', ...
 'native_flag',flag);
% The public wrapper may recover using Krylov before emitting any warning.
if nargout<=1, U=diag(E(1:k,1:k)); end
end
function y=K(x,n)
y=[-conj(x(n+1:2*n,:));conj(x(1:n,:))];
end
