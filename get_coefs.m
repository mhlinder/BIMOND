%% get_coefs.m
% M. Henry Linder (mhlinder@gmail.com)
%
% This M-file extracts the coefficients of a bicubic polynomial
% from its pp-form. Each entry in the cell `coefs` is a vector
% alpha = [a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33]
% for the cubic polynomial
% f(x, y) = \sum_{i=0}^3 \sum_{j=0}^3 a_ij x^i y^j.

%% get_coefs
function coefs = get_coefs(nx, ny, incoefs);

incoefs = reshape(incoefs, [1 nx-1 4 ny-1 4]);
for i = 1:(nx-1)
    for j = 1:(ny-1)
        coefs{i, j} = nan(16, 1);
        for k = 4:-1:1
            ix = (4-k)*4 + (1:4);
            coefs{i, j}(ix) = flipud(squeeze(incoefs(1, i, :, j, k)));
        end
    end
end
end
