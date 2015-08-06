%% BIMOND_extract_pp.m
% M. Henry Linder (mhlinder@gmail.com)
% 
% This M-file formats the output from the BIMOND algorithms as a
% ppform piecewise-polynomial for use with functions from the Curve
% Fitting Toolbox.

%% BIMOND_extract_pp
function pp2d = BIMOND_extract_pp(x, y, p, px, py, pxy);
nx = length(x);
ny = length(y);

% For explanation of how to structure inputs to `ppmak`, see
% http://cn.mathworks.com/matlabcentral/answers/77100-reconstruct-multivariate-spline-from-csapi
coefs = zeros(1, nx-1, 4, ny-1, 4);
for i = 1:(nx-1)
    for j = 1:(ny-1)
        alpha = bicubic_polynomial(abs(x(i+1) - x(i)), ...
                                   abs(y(j+1) - y(j)), ...
                                   p(i:i+1, j:j+1), ...
                                   px(i:i+1, j:j+1), ...
                                   py(i:i+1, j:j+1), ...
                                   pxy(i:i+1, j:j+1));
        for k = 4:-1:1
            ix = (4-k)*4 + (1:4);
            coefs(1, i, :, j, k) = flipud(alpha(ix));
        end
    end
end

pp2d = ppmak({x, y}, reshape(coefs, [(nx-1)*4 (ny-1)*4]));
end % function extract_pp

    
%% bicubic_polynomial
% This function calculates the coefficients for a bicubic
% polynomial from the function and derivative values.
function alpha = bicubic_polynomial(h, k, p, px, py, pxy);
% Inputs p, px, py, pxy are assumed to be 2x2, with rows indexing x
% and columns indexing y, i.e.,
% p = [f(0,0) f(0,k);
%      f(h,0) f(h,k)]

%% Calculate polynomial coefficients
% For alpha = [a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33]
% and x = [f(0,0) f(h,0) f(0,k) f(h,k) fx(0,0) fx(h,0) fx(0,k) fx(h,k)
%          fy(0,0) fy(h,0) fy(0,k) fy(h,k) fxy(0,0) fxy(h,0) fxy(0,k) fxy(h,k)]
% we have A*alpha = x, so we can solve by alpha = (A^-1)*x
x = [p(:); px(:); py(:); pxy(:)];
A = [1 0 0   0     0 0   0     0       0   0     0       0         0     0       0       0;           ...
     1 h h^2 h^3   0 0   0     0       0   0     0       0         0     0       0       0;           ...
     1 0 0   0     k 0   0     0       k^2 0     0       0         k^3   0       0       0;           ...
     1 h h^2 h^3   k h*k h^2*k h^3*k   k^2 h*k^2 h^2*k^2 h^3*k^2   k^3   h*k^3   h^2*k^3 h^3*k^3;     ...
     0 1 0   0     0 0   0     0       0   0     0       0         0     0       0       0;           ...
     0 1 2*h 3*h^2 0 0   0     0       0   0     0       0         0     0       0       0;           ...
     0 1 0   0     0 k   0     0       0   k^2   0       0         0     k^3     0       0;           ...
     0 1 2*h 3*h^2 0 k   2*h*k 3*h^2*k 0   k^2   2*h*k^2 3*h^2*k^2 0     k^3     2*h*k^3 3*h^2*k^3;   ...
     0 0 0   0     1 0   0     0       0   0     0       0         0     0       0       0;           ...
     0 0 0   0     1 h   h^2   h^3     0   0     0       0         0     0       0       0;           ...
     0 0 0   0     1 0   0     0       2*k 0     0       0         3*k^2 0       0       0;           ...
     0 0 0   0     1 h   h^2   h^3     2*k 2*h*k 2*h^2*k 2*h^3*k   3*k^2 3*h*k^2 3*h^2*k^2 3*h^3*k^3; ...
     0 0 0   0     0 1   0     0       0   0     0       0         0     0       0       0;           ...
     0 0 0   0     0 1   2*h   3*h^2   0   0     0       0         0     0       0       0;           ...
     0 0 0   0     0 1   0     0       0   2*k   0       0         0     3*k^2   0       0;           ...
     0 0 0   0     0 1   2*h   3*h^2   0   2*k   4*h*k   6*h^2*k   0     3*k^2   6*h*k^2 9*h^2*k^2];

alpha = inv(A) * x;
end % function bicubic_polynomial
