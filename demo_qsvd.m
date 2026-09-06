function demo_qsvd
%DEMO_QSVD Run a small reproducible example; no input files are required.
%   In MATLAB, change Current Folder to this package and enter: demo_qsvd

rng(0, 'twister');
m = 3;
n = 2;
B1 = [1+1i, 2; 0, 3-1i; 2, -1+2i];
B2 = [0, 1i; 1, 2; -1i, 1];
AC = quaternion_to_complex(B1, B2);
[UC, EC, VC] = qsvd_complex(AC, m, n, 1e-12);

% Use conjugate transpose (apostrophe), NOT nonconjugate transpose (dot-prime).
relative_error = norm(AC - UC*EC*VC', 'fro') / norm(AC, 'fro');
u_error = norm(UC'*UC - eye(size(UC,2)), 'fro');
v_error = norm(VC'*VC - eye(size(VC,2)), 'fro');
p = min(m,n);
singular_values = diag(EC(1:p,1:p));

fprintf('QSVD example: %d-by-%d quaternion matrix\n', m, n);
fprintf('Relative reconstruction error: %.3e\n', relative_error);
fprintf('Left orthogonality error:       %.3e\n', u_error);
fprintf('Right orthogonality error:      %.3e\n', v_error);
fprintf('Quaternion singular values:\n');
disp(singular_values);
assert(max([relative_error,u_error,v_error]) < 1e-10, ...
    'Demo validation failed. Check the input representation and MATLAB path.');
fprintf('Demo passed.\n');
end
