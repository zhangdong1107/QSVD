function [U,E,V,info] = qsvds_krylov(AC,m,n,k,opts)
% Thick-restarted, paired block Krylov iteration; no SVDS/EIGS.
if nargin < 5, opts=struct(); end
p=qsvd_validate_inputs(AC,m,n,k);
defaults=struct('tol',1e-10,'maxit',1000,'oversampling',min(20,max(4,k)), ...
 'seed',0,'cluster_tol',1e-12,'maxspace',max(60,12*k+20));
for name=fieldnames(opts).', defaults.(name{1})=opts.(name{1}); end
opts=defaults;
if ~any(AC(:)) || k==p
    basic=rmfield(opts,'maxspace');
    [U,E,V,info]=qsvds_subspace(AC,m,n,k,basic); return;
end
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(opts.seed,'twister'); alpha=norm(AC,'fro'); A=AC/alpha;
ell=min(p-1,k+opts.oversampling); capacity=min(p-1,max(ell,opts.maxspace));
if capacity==ell
    % When k is almost full, no Krylov expansion fits below p. A paired
    % two-sided iteration still changes the retained space each sweep.
    [U,E,V,info]=qsvds_subspace(AC,m,n,k,rmfield(opts,'maxspace'));
    info.method='structured-krylov-compact'; return;
end
W=append_pairs([],randn(2*n,ell)+1i*randn(2*n,ell),n,ell,0);
Y=A*W; last=1:size(W,2); restarts=0;
info=struct('method','structured-krylov','converged',false,'flag',1, ...
 'iterations',0,'relative_residuals',Inf(k,1),'subspace_dimension',0, ...
 'largest_projected_svd',0,'full_decomposition',false,'fallback_reason','');
for iteration=1:opts.maxit
    r=size(W,2);
    % Stable singular solve of a tall, thin projection, not a Gram eigensolve.
    if iteration==1 || mod(iteration,2)==0 || r==capacity
        YC=quaternion_to_complex(Y(1:m,:),-conj(Y(m+1:2*m,:)));
        [Ub,Eb,Vb]=qsvd_complex(YC,m,r,opts.cluster_tol);
        cols=[1:k,r+(1:k)]; Q=[W,K(W,n)];
        U=Ub(:,cols); V=Q*Vb(:,cols); s=real(diag(Eb(1:k,1:k)));
        U=qsvds_repair_left_null(U,s,m);
        R1=A*V(:,1:k)-U(:,1:k)*diag(s);
        R2=A'*U(:,1:k)-V(:,1:k)*diag(s);
        residual=sqrt(sum(abs(R1).^2,1)+sum(abs(R2).^2,1)).'/max(s(1),eps);
        info.iterations=iteration; info.relative_residuals=residual;
        info.subspace_dimension=r;
        info.largest_projected_svd=max(info.largest_projected_svd,r);
        if all(residual<=opts.tol), info.flag=0; info.converged=true; break; end
    end
    if r==capacity
        W=Q*Vb(:,1:ell); Y=A*W; last=1:ell; restarts=restarts+1;
    else
        Z=A'*Y(:,last);
        fresh=append_pairs(W,Z,n,min(ell,capacity-r),1e-12);
        if isempty(fresh)
            fresh=append_pairs(W,randn(2*n,ell)+1i*randn(2*n,ell),n,min(ell,capacity-r),0);
        end
        W=[W,fresh]; Y=[Y,A*fresh]; last=r+(1:size(fresh,2));
    end
end
% An odd iteration limit still needs a final extraction from the latest basis.
if info.iterations~=opts.maxit && ~info.converged
    r=size(W,2); YC=[Y,K(Y,m)];
    [Ub,Eb,Vb]=qsvd_complex(YC,m,r,opts.cluster_tol); cols=[1:k,r+(1:k)];
    U=Ub(:,cols); V=[W,K(W,n)]*Vb(:,cols); s=real(diag(Eb(1:k,1:k)));
    U=qsvds_repair_left_null(U,s,m);
    R1=A*V(:,1:k)-U(:,1:k)*diag(s); R2=A'*U(:,1:k)-V(:,1:k)*diag(s);
    info.relative_residuals=sqrt(sum(abs(R1).^2,1)+sum(abs(R2).^2,1)).'/max(s(1),eps);
    info.converged=all(info.relative_residuals<=opts.tol); info.flag=double(~info.converged);
    info.iterations=opts.maxit; info.largest_projected_svd=max(info.largest_projected_svd,r);
end
info.restarts=restarts; E=quaternion_to_complex(diag(s*alpha),zeros(k));
if ~info.converged, warning('qsvds:NoConvergence','Krylov QSVDS did not converge: %.3e.',max(info.relative_residuals)); end
if nargout<=1, U=diag(E(1:k,1:k)); end
end

function fresh=append_pairs(W,X,n,count,drop)
fresh=complex(zeros(2*n,0)); Q=[W,K(W,n)];
for j=1:size(X,2)
    x=X(:,j); scale=norm(x);
    for pass=1:2, x=x-Q*(Q'*x); end
    nx=norm(x);
    if nx>max(realmin,drop*scale)
        x=x/nx; fresh(:,end+1)=x; Q=[Q,x,K(x,n)]; %#ok<AGROW>
        if size(fresh,2)==count, return; end
    end
end
end

function y=K(x,n)
if isempty(x), y=complex(zeros(2*n,0)); return; end
y=[-conj(x(n+1:2*n,:));conj(x(1:n,:))];
end
