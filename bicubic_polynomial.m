function coefs = bicubic_polynomial(h, k, p, px, py, pxy)

if ~isscalar(h) | ~isscalar(k)
    error('Surface grid width and height must both be scalars.');
end

if ~(all(size(p) == size(px)) ...
     & all(size(p) == size(py)) ...
     & all(size(p) == size(pxy)))
    error('Input derivative matrices must have identical dimensions.');
end
