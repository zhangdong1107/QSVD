function AC = quaternion_to_complex(B1, B2)
%QUATERNION_TO_COMPLEX Build the complex representation of a quaternion matrix.
%   AC = QUATERNION_TO_COMPLEX(B1, B2) requires same-size complex blocks.
%   For A = A0 + A1*i + A2*j + A3*k, set B1 = A0 + 1i*A1 and
%   B2 = A2 + 1i*A3. The result has twice as many rows and columns.
%   No quaternion class or external toolbox is required.

AC = [B1, B2; -conj(B2), conj(B1)];
end
