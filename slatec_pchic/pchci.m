function [n,h,slope,d,incfd]=pchci(n,h,slope,d,incfd);
%***BEGIN PROLOGUE  PCHCI
%***SUBSIDIARY
%***PURPOSE  Set interior derivatives for PCHIC
%***LIBRARY   SLATEC (PCHIP)
%***TYPE      SINGLE PRECISION (PCHCI-S, DPCHCI-D)
%***AUTHOR  Fritsch, F. N., (LLNL)
%***DESCRIPTION
%
%          PCHCI:  PCHIC Initial Derivative Setter.
%
%    Called by PCHIC to set derivatives needed to determine a monotone
%    piecewise cubic Hermite interpolant to the data.
%
%    Default boundary conditions are provided which are compatible
%    with monotonicity.  If the data are only piecewise monotonic, the
%    interpolant will have an extremum at each point where monotonicity
%    switches direction.
%
%    To facilitate two-dimensional applications, includes an increment
%    between successive values of the D-array.
%
%    The resulting piecewise cubic Hermite function should be identical
%    (within roundoff error) to that produced by PCHIM.
%
% ----------------------------------------------------------------------
%
%  Calling sequence:
%
%        PARAMETER  (INCFD = ...)
%        INTEGER  N
%        REAL  H(N), SLOPE(N), D(INCFD,N)
%
%        CALL  PCHCI (N, H, SLOPE, D, INCFD)
%
%   Parameters:
%
%     N -- (input) number of data points.
%           If N=2, simply does linear interpolation.
%
%     H -- (input) real array of interval lengths.
%     SLOPE -- (input) real array of data slopes.
%           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
%                  H(I) =  X(I+1)-X(I),
%              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
%
%     D -- (output) real array of derivative values at the data points.
%           If the data are monotonic, these values will determine a
%           a monotone cubic Hermite function.
%           The value corresponding to X(I) is stored in
%                D(1+(I-1)*INCFD),  I=1(1)N.
%           No other entries in D are changed.
%
%     INCFD -- (input) increment between successive values in D.
%           This argument is provided primarily for 2-D applications.
%
%    -------
%    WARNING:  This routine does no validity-checking of arguments.
%    -------
%
%  Fortran intrinsics used:  ABS, MAX, MIN.
%
%***SEE ALSO  PCHIC
%***ROUTINES CALLED  PCHST
%***REVISION HISTORY  (YYMMDD)
%   820218  DATE WRITTEN
%   820601  Modified end conditions to be continuous functions of
%           data when monotonicity switches in next interval.
%   820602  1. Modified formulas so end conditions are less prone
%             to over/underflow problems.
%           2. Minor modification to HSUM calculation.
%   820805  Converted to SLATEC library version.
%   890411  Added SAVE statements (Vers. 3.2).
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890831  Modified array declarations.  (WRB)
%   890831  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900328  Added TYPE section.  (WRB)
%   910408  Updated AUTHOR section in prologue.  (WRB)
%   930503  Improved purpose.  (FNF)
%***end PROLOGUE  PCHCI
%
%  Programming notes:
%     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
%        either argument is zero, +1 if they are of the same sign, and
%        -1 if they are of opposite sign.
%**end
%
%  DECLARE ARGUMENTS.
%
persistent del1 del2 dmax dmin drat1 drat2 firstCall hsum hsumt3 i nless1 three w1 w2 zero ; if isempty(firstCall),firstCall=1;end; 

h_shape=size(h);h=reshape(h,1,[]);
slope_shape=size(slope);slope=reshape(slope,1,[]);
d_shape=size(d);d=reshape([d(:).',zeros(1,ceil(numel(d)./prod([incfd])).*prod([incfd])-numel(d))],incfd,[]);
%
%  DECLARE LOCAL VARIABLES.
%
if isempty(i), i=0; end;
if isempty(nless1), nless1=0; end;
if isempty(del1), del1=0; end;
if isempty(del2), del2=0; end;
if isempty(dmax), dmax=0; end;
if isempty(dmin), dmin=0; end;
if isempty(drat1), drat1=0; end;
if isempty(drat2), drat2=0; end;
if isempty(hsum), hsum=0; end;
if isempty(hsumt3), hsumt3=0; end;
if isempty(three), three=0; end;
if isempty(w1), w1=0; end;
if isempty(w2), w2=0; end;
if isempty(zero), zero=0; end;
%
%  INITIALIZE.
%
if firstCall,   zero=[0.];  end;
if firstCall,  three=[3.];  end;
firstCall=0;
%***FIRST EXECUTABLE STATEMENT  PCHCI
nless1 = fix(n - 1);
del1 = slope(1);
%
%  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
%
if( nless1>1 )
%
%  NORMAL CASE  (N .GE. 3).
%
del2 = slope(2);
%
%  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
%     SHAPE-PRESERVING.
%
hsum = h(1) + h(2);
w1 =(h(1)+hsum)./hsum;
w2 = -h(1)./hsum;
d(1,1) = w1.*del1 + w2.*del2;
if( pchst(d(1,1),del1)<=zero )
d(1,1) = zero;
elseif( pchst(del1,del2)<zero ) ;
%        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
dmax = three.*del1;
if( abs(d(1,1))>abs(dmax) )
d(1,1) = dmax;
end;
end;
%
%  LOOP THROUGH INTERIOR POINTS.
%
for i = 2 : nless1;
if( i~=2 )
%
hsum = h(i-1) + h(i);
del1 = del2;
del2 = slope(i);
end;
%
%        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
%
d(1,i) = zero;
if( pchst(del1,del2)>zero )
%
%        use BRODLIE MODIFICATION OF BUTLAND FORMULA.
%
hsumt3 = hsum + hsum + hsum;
w1 =(hsum+h(i-1))./hsumt3;
w2 =(hsum+h(i))./hsumt3;
dmax = max(abs(del1),abs(del2));
dmin = min(abs(del1),abs(del2));
drat1 = del1./dmax;
drat2 = del2./dmax;
d(1,i) = dmin./(w1.*drat1+w2.*drat2);
end;
%
end; i = fix(nless1+1);
%
%  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
%     SHAPE-PRESERVING.
%
w1 = -h(n-1)./hsum;
w2 =(h(n-1)+hsum)./hsum;
d(1,n) = w1.*del1 + w2.*del2;
if( pchst(d(1,n),del2)<=zero )
d(1,n) = zero;
elseif( pchst(del1,del2)<zero ) ;
%        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
dmax = three.*del2;
if( abs(d(1,n))>abs(dmax) )
d(1,n) = dmax;
end;
end;
else;
d(1,1) = del1;
d(1,n) = del1;
end;
%
%  NORMAL RETURN.
%
%------------- LAST LINE OF PCHCI FOLLOWS ------------------------------
h_shape=zeros(h_shape);h_shape(:)=h(1:numel(h_shape));h=h_shape;
slope_shape=zeros(slope_shape);slope_shape(:)=slope(1:numel(slope_shape));slope=slope_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
end
%DECK PCHCM
