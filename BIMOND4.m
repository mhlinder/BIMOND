%% BIMOND4.m
% M. Henry Linder (mhlinder@gmail.com)
% 
% This M-file implements the BIMOND4 algorithm presented in "A
% Bivariate Interpolation Algorithm For Data That Are Monotone In One
% Variable" (Carlson and Fritsch, 1991).
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


%% BIMOND4
function [px, py, pxy] = BIMOND4(x, y, p);

%% Verify inputs
nx = length(x);
ny = length(y);

% Check that input is properly 2-dimensional
if ~(nx >= 4 & ny >= 4)
    error(['Input data must be 2-dimensional with at least 4 points ' ...
           'along each input axis (i.e., `length(x) >= 4 & length(y) ' ...
           '>= 4`']);
end

% Check that dimensions match. Note that x indexes rows and y
% indexes columns
if ~all([nx ny] == size(p))
    error('Input value dimensions do not match function value dimensions.');
end


%% Step 1: Verify monotonicity, compute sx, and other variable setup

% Verify monotonicity in at least one direction
if ~(all(all(diff(p, 1, 1) > 0)) | all(all(diff(p, 1, 1) < 0)))
    % If monotone along y axis, flip axes so monotone in x direction
    if all(all(diff(p, 1, 2) > 0)) | all(all(diff(p, 1, 2) < 0))
        flipped = 1;
        input_x = x;
        x = y;
        y = input_x;
        p = p';
    else
        error(['Input function values are not monotone along either ' ...
               'axis. BIMOND4 requires that the interpolated function ' ...
               'be monotone along at least one axis.']);
    end
else
    flipped = 0;
end
sx = sign(p(2, 1) - p(1, 1));

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
        a(i, j) = 3 * k(j) * del1(i, j);
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


%% Step 3: Adjust y partial derivatives as needed
pxx = px;
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
    py(:, j) = sweep_bimond4(py(:, j), L(:, j), R(:, j));
end


%% Step 4: Compute values of crossed derivates pxy
% Three point difference formulae from
% http://www.sitmo.com/article/numerical-differentiation/
%
% For first derivative, these equations are:
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

% Adjust to satisfy inequalities
pxplus = nan(size(p));
pxminus = nan(size(p));
klag = [nan k(1:(end-1))];
for j = 1:ny
    pxplus(:, j) = sx * px(:, j) / k(j);
    pxminus(:, j) = sx * px(:, j) / klag(j);
end

del2prime = nan(size(p));
for j = 1:ny
    for i = 1:(nx-1)
        del2prime(i, j) = (py(i+1, j) - py(i, j)) / h(i);
    end
    del2prime(end, j) = nan;
end

C = nan(size(p));
D = nan(size(p));
for i = 1:nx
    for j = 1:ny
        C(i, j) = sx * (del2prime(i, j) - 3 * del2(i, j) / ...
                        klag(j));
        D(i, j) = sx * (del2prime(i, j) + 3 * del2(i, j) / k(j));
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
        boundsrx2(i, j) = 3 * (D(i, j) - pxplus(i, j));
        boundsrx1(i, j) = 3 * pxminus(i, j);
        if i > 1
            boundslx3(i, j) = 3 * (C(i-1, j) + pxminus(i, j));
            boundsrx3(i ,j) = 3 * (D(i-1, j) - pxplus(i, j));
        end
    end
end
boundslx = max( max(boundslx1, boundslx2), ...
                boundslx3);
boundsrx = min( min(boundsrx1, boundsrx2), ...
                boundsrx3);

for i = 1:nx
    for j = 1:ny
        if sx * pxy(i, j) < boundslx(i, j)
            pxy(i, j) = boundslx(i, j);
        end
        if sx * pxy(i, j) > boundsrx(i, j)
            pxy(i, j) = boundsrx(i, j);
        end
    end
end


%% If input data is monotone along y-axis, revert to original axes
if flipped
    px_old = px;
    px = py;
    py = px_old;
    pxy = pxy';
end
end % function BIMOND4


%% sweep_bimond4
function d = sweep_bimond4(d, lhs, rhs)
% This function does NOT verify inputs
n = length(d);

for i = 2:(n-1)
    if sign(d(i)) ~= sign(d(i+1))
        sd = sign(d(i));
        sd1 = sign(d(i+1));
        if sd == 0
            if sd1 == -1
                sd = 1;
            else % sd1 == 1
                sd = -1;
            end
        end
        if sd == 1
            if d(i+1) - d(i) < -1*lhs(i)
                rho = -1*lhs(i) / (d(i+1) - d(i));
                d(i) = rho*d(i);
                d(i+1) = rho*d(i+1);
            end
        else % sd == -1
            if d(i+1) - d(i) > rhs(i)
                rho = rhs(i) / (d(i+1) - d(i));
                d(i) = rho*d(i);
                d(i+1) = rho*d(i+1);
            end
        end
    end
end

% Find strings of same sign
signs = sign(d);
sign_locs = find(signs(1:(end-1))' ~= signs(2:end)');
ix = [1; reshape(sign_locs, [length(sign_locs), 1]); n];
ixpairs = [ix(1:(end-1)) ix(2:end)];
ixpairs(:, 1) = ixpairs(:, 1) + 1;
ixpairs(1, 1) = ixpairs(1, 1) - 1;
% Rows of ixpairs give [start end] for strings

for i = 1:size(ixpairs, 1)
    if signs(ixpairs(i, 1) == 1)
        % upsweep
        for j = ixpairs(i, 1):(ixpairs(i, 2)-1)
            if d(j+1) - d(j) >  rhs(j)
                d(j+1) = d(j) + rhs(j);
            end
        end
        % downsweep
        for j = flip(ixpairs(i, 1):(ixpairs(i, 2)-1))
            if d(j+1) - d(j) < -1*lhs(j)
                d(j) = d(j+1) + lhs(j);
            end
        end
    else
        % upsweep
        for j = ixpairs(i, 1):(ixpairs(i, 2)-1)
            if d(j+1) - d(j) < -1*lhs(j)
                d(j+1) = d(j) - lhs(j);
            end
        end
        % downsweep
        for j = flip(ixpairs(i, 1):(ixpairs(i, 2)-1))
            if d(j+1) - d(j) > rhs(j)
                d(j) = d(j+1) - rhs(j);
            end
        end
    end
end
end % function sweep_bimond4
    
