%% BIMOND3.m
% M. Henry Linder (mhlinder@gmail.com)
% 
% This M-file implements the BIMOND3 algorithm presented in "An
% Algorithm For Monotone Piecewise Bicubic Interpolation" (Carlson and
% Fritsch, 1989).
% 
% Recreating the results of Carlson and Fritsch requires the following
% functions be available in your path. These files can be obtained as
% part of the `slatec` MATLAB package, available at
% http://www.mathworks.com/matlabcentral/fileexchange/14535-slatec
% * pchci.m
% * pchcs.m
% * pchic.m
% * pchst.m
% * pchsw.m
% * r1mach.m


%% BIMOND3
function pp2d = BIMOND3(x, y, p);

%% Verify inputs
nx = length(x);
ny = length(y);

% Check that input is properly 2-dimensional
if ~(isvector(x) & isvector(y)) | ~(nx >= 4 & ny >= 4)
    error(['Input data must be 2-dimensional with at least 4 points ' ...
           'along each input axis (i.e., `length(x) >= 4 & length(y) ' ...
           '>= 4`']);
end

% Check that dimensions match. Note that x indexes rows and y
% indexes columns
if ~all([nx ny] == size(p))
    error('Input value dimensions do not match function value dimensions.');
end


%% Step 1: verify monotonicity, calculate sx and sy, and other
% variable setup

% Verify monotonicity in both directions
if ~(all(all(diff(p, 1, 1) > 0)) | all(all(diff(p, 1, 1) < 0))) | ...
        ~(all(all(diff(p, 1, 2) > 0)) | all(all(diff(p, 1, 2) < 0)))
    error(['Input function values are not monotone along both axes. ' ...
           'BIMOND3 requires that the interpolated function be ' ...
           'monotone in both directions.'])
end
sx = sign(p(2, 1) - p(1, 1));
sy = sign(p(1, 2) - p(1, 1));

h = [x(2:end) - x(1:(end-1)), nan];
k = [y(2:end) - y(1:(end-1)), nan];

del1 = nan(nx, ny);
del2 = nan(nx, ny);
for i = 1:nx
    for j = 1:(ny-1)
        del1(i, j) = (p(i, j+1) - p(i, j)) / k(j);
    end
    del1(i, end) = nan;
end
for j = 1:ny
    for i = 1:(nx-1)
        del2(i, j) = (p(i+1, j) - p(i, j)) / h(i);
    end
    del2(end, j) = nan;
end

a = nan(size(del1));
b = nan(size(del2));
for i = 1:nx
    for j = 1:(ny-1)
        a(i, j) = 3 * k(j) * sy * del1(i, j);
    end
end
for i = 1:(nx-1)
    for j = 1:ny
        b(i, j) = 3 * h(i) * sx * del2(i, j);
    end
end


%% Step 2: Initialize partial derivatives px, py

% Intitialize inputs to PCHIC. The value of `SWITCH` ensures
% numerically identical results to those in the paper. Specifically,
% this "smooths" the interpolating function, preventing extrema from
% having partial derivatives equal to zero.
IC = [0 0];
VC = [0 0];
SWITCH = 1;
INCFD = 1;
IERR = 0;

% Initialize partial derivative matrices
px = nan(size(p));
py = nan(size(p));

% Set partials with respect to x-axis
for i = 1:length(y)
    [ic,vc,switchml,n,outx,f,d,incfd,wk,nwk,ierr] = pchic(IC, VC, SWITCH, ...
                                                      nx, x, p(:, i), ...
                                                      zeros(1, nx), INCFD, ...
                                                      zeros(1, 2*(nx-1)), ...
                                                      2*(nx-1), IERR);
    px(:, i) = d;
end

% Set partials with respect to y-axis
for i = 1:length(x)
    [ic,vc,switchml,n,outy,f,d,incfd,wk,nwk,ierr] = pchic(IC, VC, SWITCH, ...
                                                      ny, y, p(i, :), ...
                                                      zeros(1, ny), INCFD, ...
                                                      zeros(1, 2*(ny-1)), ...
                                                      2*(ny-1), IERR);
    py(i, :) = d;
end


%% Step 3: Adjust partial derivatives as needed
% This is omitted as unnecessary (see Carlson and Fritsch).


%% Step 4: Adjust partial derivatives as needed
% All sweeps are done with the original values of px, py
pxx = px;
pyy = py;

% Step 4a: Adjust x partial derivatives
D = nan(nx, ny);
U = nan(nx, ny);
for i = 1:nx
    for j = 1:(ny-1)
        D(i, j) = (3 * k(j) * sy * del1(i, j) ...
                   - k(j) * max([sy*pyy(i, j+1), sy*pyy(i, j)])) / h(i);
        if i > 1
            U(i, j) = (3 * k(j) * sy * del1(i, j) ...
                       - k(j) * max([sy*pyy(i, j+1), sy*pyy(i, j)])) / h(i-1);
        end
    end
end

for i = 1:nx
    px(i, :) = sx * sweep_bimond3(abs(px(i, :)), D(i, :), U(i, :));
end

% Step 4b: Adjust y partial derivatives
L = nan(nx, ny);
R = nan(nx, ny);
for j = 1:ny
    for i = 1:(nx-1)
        L(i, j) = (3 * h(i) * sx * del2(i, j) ...
                   - h(i) * max([sx*pxx(i+1, j), sx*pxx(i, j)])) / k(j);
        if j > 1
            R(i, j) = (3 * h(i) * sx * del2(i, j) ...
                       - h(i) * max([sx*pxx(i+1, j), sx*pxx(i, j)])) / k(j-1);
        end
    end
end
for j = 1:ny
    py(:, j) = sy * sweep_bimond3(abs(py(:, j)), L(:, j), R(:, j));
end


%% Step 5: Compute values of crossed derivatives pxy
% Three point difference formulae from
% http://www.sitmo.com/article/numerical-differentiation/
%
% For first derivative, one equation is:
% d/dy px(x, y) = (px(x, y + h) - px(x, y - h)) / (2 * h)
% d/dy px(x, y) = (-px(x, y + 2h), + 4 * px(x, y + h)
%                   - 3 * px(x,y)) / (2 * h)

partialxy = nan(size(p));
for i = 1:nx
    for j = 1:ny
        if j > 1 && j < ny
            partialxy(i, j) = (px(i, j+1) - px(i, j-1)) / 2;
        elseif j == 1
            partialxy(i, j) = (-1*px(i, j+2) + 4*px(i, j+1) - 3*px(i, j)) ...
                / 2;
        else % j == ny
            partialxy(i, j) = (-1*px(i, j-2) + 4*px(i, j-1) - 3*px(i, j)) ...
                / -2;
        end
    end
end

partialyx = nan(size(p));
for j = 1:ny
    for i = 1:nx
        if i > 1 && i < nx
            partialyx(i, j) = (py(i+1, j) - py(i-1, j)) / 2;
        elseif i == 1
            partialyx(i, j) = (-1*py(i+2, j) + 4*py(i+1, j) - 3*py(i, j)) ...
                / 2;
        else % i == nx
            partialyx(i, j) = (-1*py(i-2, j) + 4*py(i-1, j) - 3*py(i, j)) ...
                / -2;
        end
    end
end

pxy = (partialxy + partialyx) / 2;

% 2.3
pyplus = nan(size(p));
pyminus = nan(size(p));
hlag = [nan h(1:(end-1))];
for i = 1:nx
    pyplus(i, :) = sy * py(i, :) / h(i);
    pyminus(i, :) = sy * py(i, :) / hlag(i);
end
pxplus = nan(size(p));
pxminus = nan(size(p));
klag = [nan k(1:(end-1))];
for j = 1:ny
    pxplus(:, j) = sx * px(:, j) / k(j);
    pxminus(:, j) = sx * px(:, j) / klag(j);
end

del1prime = nan(size(p));
del2prime = nan(size(p));
for i = 1:nx
    for j = 1:(ny-1)
        del1prime(i, j) = (px(i, j+1) - px(i, j)) / k(j);
    end
    del1prime(i, end) = nan;
end
for j = 1:ny
    for i = 1:(nx-1)
        del2prime(i, j) = (py(i+1, j) - py(i, j)) / h(i);
    end
    del2prime(end, j) = nan;
end

A = nan(size(p));
B = nan(size(p));
C = nan(size(p));
D = nan(size(p));
for i = 1:nx
    for j = 1:ny
        A(i, j) = sy * (del1prime(i, j) - 3 * del1(i, j) / ...
                        hlag(i));
        B(i, j) = sy * (del1prime(i, j) + 3 * del1(i, j) / h(i));
        C(i, j) = sx * (del2prime(i, j) - 3 * del2(i, j) / ...
                        klag(j));
        D(i, j) = sx * (del2prime(i, j) + 3 * del2(i, j) / k(j));
    end
end

boundsly1 = nan(size(p)); % left y bounds
boundsry1 = nan(size(p)); % right y bounds
boundsly2 = nan(size(p));
boundsry2 = nan(size(p));
boundsly3 = nan(size(p));
boundsry3 = nan(size(p));
for i = 1:nx
    for j = 1:ny
        boundsly1(i, j) = -3 * pyplus(i, j);
        boundsly2(i, j) = 3 * (A(i, j) + pyminus(i, j));
        boundsry1(i, j) = 3 * (B(i, j) - pyplus(i, j));
        boundsry3(i, j) = 3 * pyminus(i, j);
        if j > 1
            boundsly3(i, j) = 3 * (A(i, j-1) + pyminus(i, j));
            boundsry2(i, j) = 3 * (B(i, j-1) - pyplus(i, j));
        end
    end
end

boundslx1 = nan(size(p));
boundsrx1 = nan(size(p));
boundslx2 = nan(size(p));
boundsrx2 = nan(size(p));
boundslx3 = nan(size(p));
boundsrx3 = nan(size(p));
for j = 1:ny
    for i = 1:nx
        boundslx1(i, j) = -3 * pxplus(i, j);
        boundslx2(i, j) = 3 * (C(i, j) + pxminus(i, j));
        boundsrx1(i, j) = 3 * (D(i, j) - pxplus(i, j));
        boundsrx3(i, j) = 3 * pxminus(i, j);
        if i > 1
            boundslx3(i, j) = 3 * (C(i-1, j) + pxminus(i, j));
            boundsrx2(i ,j) =3 * (D(i-1, j) - pxplus(i, j));
        end
    end
end
boundsly = max( max(boundsly1, boundsly2), ...
                boundsly3);
boundsry = min( min(boundsry1, boundsry2), ...
                boundsry3);
boundslx = max( max(boundslx1, boundslx2), ...
                boundslx3);
boundsrx = min( min(boundsrx1, boundsrx2), ...
                boundsrx3);

for i = 1:nx
    for j = 1:ny
        if sy * pxy(i, j) < boundsly(i, j)
            pxy(i, j) = boundsly(i, j);
        end
        if sy * pxy(i, j) > boundsry(i, j)
            pxy(i, j) = boundsry(i, j);
        end
        if sx * pxy(i, j) < boundslx(i, j)
            pxy(i, j) = boundslx(i, j);
        end
        if sx * pxy(i, j) > boundsrx(i, j)
            pxy(i, j) = boundsrx(i, j);
        end
    end
end

pp2d = extract_pp(x, y, p, px, py, pxy);
end % function BIMOND3
    

%% sweep_bimond3
function d = sweep_bimond3(d, lhs, rhs)
% This function does NOT verify inputs
n = length(d);

% upsweep
for k = 1:(n-1)
    if d(k+1) - d(k) > rhs(k)
        d(k+1) = d(k) + rhs(k);
    end
end

% downsweep
for k = flip(1:(n-1))
    if d(k+1) - d(k) < -1*lhs(k)
        d(k) = d(k+1) + lhs(k);
    end
end
end % function sweep_bimond3
