function alphas = extract_pp(x, y, p, px, py, pxy);
% Inputs are NOT verified

nx = length(x);
ny = length(y);

breaks = cell(1,2);
breaks{1} = x;
breaks{2} = y;

alphas = cell(nx-1, ny-1);
for i = 1:(nx-1)
    for j = 1:(ny-1)
        alpha = bicubic_polynomial(abs(x(i+1) - x(i)), ...
                                   abs(y(j+1) - y(j)), ...
                                   p(i:i+1, j:j+1), ...
                                   px(i:i+1, j:j+1), ...
                                   py(i:i+1, j:j+1), ...
                                   pxy(i:i+1, j:j+1));
        alphas{i, j} = alpha;
    end
end

% For explanation of how to structure inputs to `ppmak`, see
% http://cn.mathworks.com/matlabcentral/answers/77100-reconstruct-multivariate-spline-from-csapi
coefs = zeros(1, nx-1, 4, ny-1, 4);
for i = 1:(nx-1)
    for j = 1:(ny-1)
        alpha = alphas{i, j};
        for k = 4:-1:1
            ix = (4-k)*4 + (1:4);
            coefs(1, i, :, j, k) = flipud(alpha(ix));
        end
    end
end

end
