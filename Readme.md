# BIMOND: Piecewise Bicubic Interpolation

These files provide MATLAB implementations of the BIMOND3 and BIMOND4
algorithms presented by Carlson and Fritsch in [1] and [2],
respectively. The BIMOND algorithms are a 2D extension of PCHIP
(Piecewise Cubic Hermite Interpolating Polynomial), preserving
monotonicity in the underlying data.

This package depends on the MATLAB Curve Fitting
Toolbox. Additionally, these functions assume that the following files
are available in the user's MATLAB path:

  * `pchci.m`
  * `pchcs.m`
  * `pchic.m`
  * `pchst.m`
  * `pchsw.m`
  * `r1mach.m`
  
These files are included in this repository but can also be obtained
as part of the
[slatec](http://www.mathworks.com/matlabcentral/fileexchange/14535-slatec)
MATLAB package.

## Input

Both `BIMOND3` and `BIMOND4` take as input vectors `x` and `y`, and a
matrix of function values `p` where the rows index `x` and the columns
index`y`.

The functions output a ppform spline that can be evaluated with
functions from the MATLAB Curve Fitting Toolbox.

## BIMOND3

The BIMOND3 algorithm presented in [1] provides two-dimensional
interpolation of data that is monotone in both
variables. `BIMOND3_example.m` reproduces the numeric example given in
Appendix B.

## BIMOND4

The BIMOND4 algorithm presented in [2] provides two-dimensional
interpolation of data that is monotone in only one
variable. `BIMOND4_example.m` reproduces the numeric example given in
Section 5.

## References

[1] Carlson, R. E., & Fritsch, F. N. (1989). An algorithm for monotone
piecewise bicubic interpolation. *SIAM Journal on Numerical Analysis*,
26(1), 230-238.

[2] Carlson, R. E., & Fritsch, F. N. (1991). A bivariate interpolation
algorithm for data that are monotone in one variable. *SIAM Journal on
Scientific and Statistical Computing*, 12(4), 859-866.

