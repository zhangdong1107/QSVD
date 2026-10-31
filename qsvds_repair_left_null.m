function U=qsvds_repair_left_null(U,s,m)
% A thin complex SVD need not return a K-invariant left nullspace.
% Complete zero-mode representatives in pairs against all nonzero modes.
% Actual A'*U residuals are subsequently checked by the calling solver.
k=numel(s); zero=find(s<=100*eps*max(s(1),realmin));
if isempty(zero), return; end
positive=find(s>100*eps*max(s(1),realmin));
W=U(:,1:k); G=W(:,positive); Q=[G,K(G,m)];
for j=zero(:).'
    accepted=false;
    for attempt=1:20
        x=randn(2*m,1)+1i*randn(2*m,1);
        for pass=1:2, x=x-Q*(Q'*x); end
        nx=norm(x);
        if nx>1e-12
            x=x/nx; W(:,j)=x; Q=[Q,x,K(x,m)]; accepted=true; break; %#ok<AGROW>
        end
    end
    if ~accepted, error('qsvds:NullBasis','Unable to complete left zero modes.'); end
end
U=[W,K(W,m)];
end
function y=K(x,m)
y=[-conj(x(m+1:2*m,:));conj(x(1:m,:))];
end
