function [ic,vc,switchml,n,x,f,d,incfd,wk,nwk,ierr]=pchic(ic,vc,switchml,n,x,f,d,incfd,wk,nwk,ierr);
%***BEGIN PROLOGUE  PCHIC
%***PURPOSE  Set derivatives needed to determine a piecewise monotone
%            piecewise cubic Hermite interpolant to given data.
%            User control is available over boundary conditions and/or
%            treatment of points where monotonicity switches direction.
%***LIBRARY   SLATEC (PCHIP)
%***CATEGORY  E1A
%***TYPE      SINGLE PRECISION (PCHIC-S, DPCHIC-D)
%***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
%             PCHIP, PIECEWISE CUBIC INTERPOLATION,
%             SHAPE-PRESERVING INTERPOLATION
%***AUTHOR  Fritsch, F. N., (LLNL)
%             Lawrence Livermore National Laboratory
%             P.O. Box 808  (L-316)
%             Livermore, CA  94550
%             FTS 532-4275, (510) 422-4275
%***DESCRIPTION
%
%         PCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
%
%     Sets derivatives needed to determine a piecewise monotone piece-
%     wise cubic interpolant to the data given in X and F satisfying the
%     boundary conditions specified by IC and VC.
%
%     The treatment of points where monotonicity switches direction is
%     controlled by argument SWITCH.
%
%     To facilitate two-dimensional applications, includes an increment
%     between successive values of the F- and D-arrays.
%
%     The resulting piecewise cubic Hermite function may be evaluated
%     by PCHFE or PCHFD.
%
% ----------------------------------------------------------------------
%
%  Calling sequence:
%
%        PARAMETER  (INCFD = ...)
%        INTEGER  IC(2), N, NWK, IERR
%        REAL  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
%
%        CALL  PCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
%
%   Parameters:
%
%     IC -- (input) integer array of length 2 specifying desired
%           boundary conditions:
%           IC(1) = IBEG, desired condition at beginning of data.
%           IC(2) = IEND, desired condition at end of data.
%
%           IBEG = 0  for the default boundary condition (the same as
%                     used by PCHIM).
%           If IBEG.NE.0, then its sign indicates whether the boundary
%                     derivative is to be adjusted, if necessary, to be
%                     compatible with monotonicity:
%              IBEG.GT.0  if no adjustment is to be performed.
%              IBEG.LT.0  if the derivative is to be adjusted for
%                     monotonicity.
%
%           Allowable values for the magnitude of IBEG are:
%           IBEG = 1  if first derivative at X(1) is given in VC(1).
%           IBEG = 2  if second derivative at X(1) is given in VC(1).
%           IBEG = 3  to use the 3-point difference formula for D(1).
%                     (Reverts to the default b.c. if N.LT.3 .)
%           IBEG = 4  to use the 4-point difference formula for D(1).
%                     (Reverts to the default b.c. if N.LT.4 .)
%           IBEG = 5  to set D(1) so that the second derivative is con-
%              tinuous at X(2). (Reverts to the default b.c. if N.LT.4.)
%              This option is somewhat analogous to the 'not a knot'
%              boundary condition provided by PCHSP.
%
%          NOTES (IBEG):
%           1. An error return is taken if ABS(IBEG).GT.5 .
%           2. Only in case  IBEG.LE.0  is it guaranteed that the
%              interpolant will be monotonic in the first interval.
%              If the returned value of D(1) lies between zero and
%              3*SLOPE(1), the interpolant will be monotonic.  This
%              is **NOT** checked if IBEG.GT.0 .
%           3. If IBEG.LT.0 and D(1) had to be changed to achieve mono-
%              tonicity, a warning error is returned.
%
%           IEND may take on the same values as IBEG, but applied to
%           derivative at X(N).  In case IEND = 1 or 2, the value is
%           given in VC(2).
%
%          NOTES (IEND):
%           1. An error return is taken if ABS(IEND).GT.5 .
%           2. Only in case  IEND.LE.0  is it guaranteed that the
%              interpolant will be monotonic in the last interval.
%              If the returned value of D(1+(N-1)*INCFD) lies between
%              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
%              This is **NOT** checked if IEND.GT.0 .
%           3. If IEND.LT.0 and D(1+(N-1)*INCFD) had to be changed to
%              achieve monotonicity, a warning error is returned.
%
%     VC -- (input) real array of length 2 specifying desired boundary
%           values, as indicated above.
%           VC(1) need be set only if IC(1) = 1 or 2 .
%           VC(2) need be set only if IC(2) = 1 or 2 .
%
%     SWITCH -- (input) indicates desired treatment of points where
%           direction of monotonicity switches:
%           Set SWITCH to zero if interpolant is required to be mono-
%           tonic in each interval, regardless of monotonicity of data.
%             NOTES:
%              1. This will cause D to be set to zero at all switch
%                 points, thus forcing extrema there.
%              2. The result of using this option with the default boun-
%                 dary conditions will be identical to using PCHIM, but
%                 will generally cost more compute time.
%                 This option is provided only to facilitate comparison
%                 of different switch and/or boundary conditions.
%           Set SWITCH nonzero to use a formula based on the 3-point
%              difference formula in the vicinity of switch points.
%           If SWITCH is positive, the interpolant on each interval
%              containing an extremum is controlled to not deviate from
%              the data by more than SWITCH*DFLOC, where DFLOC is the
%              maximum of the change of F on this interval and its two
%              immediate neighbors.
%           If SWITCH is negative, no such control is to be imposed.
%
%     N -- (input) number of data points.  (Error return if N.LT.2 .)
%
%     X -- (input) real array of independent variable values.  The
%           elements of X must be strictly increasing:
%                X(I-1) .LT. X(I),  I = 2(1)N.
%           (Error return if not.)
%
%     F -- (input) real array of dependent variable values to be inter-
%           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
%
%     D -- (output) real array of derivative values at the data points.
%           These values will determine a monotone cubic Hermite func-
%           tion on each subinterval on which the data are monotonic,
%           except possibly adjacent to switches in monotonicity.
%           The value corresponding to X(I) is stored in
%                D(1+(I-1)*INCFD),  I=1(1)N.
%           No other entries in D are changed.
%
%     INCFD -- (input) increment between successive values in F and D.
%           This argument is provided primarily for 2-D applications.
%           (Error return if  INCFD.LT.1 .)
%
%     WK -- (scratch) real array of working storage.  The user may wish
%           to know that the returned values are:
%              WK(I)     = H(I)     = X(I+1) - X(I) ;
%              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
%           for  I = 1(1)N-1.
%
%     NWK -- (input) length of work array.
%           (Error return if  NWK.LT.2*(N-1) .)
%
%     IERR -- (output) error flag.
%           Normal return:
%              IERR = 0  (no errors).
%           Warning errors:
%              IERR = 1  if IBEG.LT.0 and D(1) had to be adjusted for
%                        monotonicity.
%              IERR = 2  if IEND.LT.0 and D(1+(N-1)*INCFD) had to be
%                        adjusted for monotonicity.
%              IERR = 3  if both of the above are truemlv.
%           'Recoverable' errors:
%              IERR = -1  if N.LT.2 .
%              IERR = -2  if INCFD.LT.1 .
%              IERR = -3  if the X-array is not strictly increasing.
%              IERR = -4  if ABS(IBEG).GT.5 .
%              IERR = -5  if ABS(IEND).GT.5 .
%              IERR = -6  if both of the above are truemlv.
%              IERR = -7  if NWK.LT.2*(N-1) .
%             (The D-array has not been changed in any of these cases.)
%               NOTE:  The above errors are checked in the order listed,
%                   and following arguments have **NOT** been validated.
%
%***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
%                 Package, Report UCRL-87285, Lawrence Livermore Nation-
%                 al Laboratory, July 1982.  [Poster presented at the
%                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
%               2. F. N. Fritsch and J. Butland, A method for construc-
%                 ting local monotone piecewise cubic interpolants, SIAM
%                 Journal on Scientific and Statistical Computing 5, 2
%                 (June 1984), pp. 300-304.
%               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
%                 cubic interpolation, SIAM Journal on Numerical Ana-
%                 lysis 17, 2 (April 1980), pp. 238-246.
%***ROUTINES CALLED  PCHCE, PCHCI, PCHCS, XERMSG
%***REVISION HISTORY  (YYMMDD)
%   820218  DATE WRITTEN
%   820804  Converted to SLATEC library version.
%   870813  Updated Reference 2.
%   890411  Added SAVE statements (Vers. 3.2).
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890703  Corrected category record.  (WRB)
%   890831  Modified array declarations.  (WRB)
%   890831  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
%   920429  Revised format and order of references.  (WRB,FNF)
%***end PROLOGUE  PCHIC
%  Programming notes:
%
%     To produce a doubleprecision version, simply:
%        a. Change PCHIC to DPCHIC wherever it occurs,
%        b. Change PCHCE to DPCHCE wherever it occurs,
%        c. Change PCHCI to DPCHCI wherever it occurs,
%        d. Change PCHCS to DPCHCS wherever it occurs,
%        e. Change the real declarations to doubleprecision, and
%        f. Change the constant  ZERO  to doubleprecision.
%
%  DECLARE ARGUMENTS.
%
persistent firstCall i ibeg iend nless1 zero ; if isempty(firstCall),firstCall=1;end; 

x_shape=size(x);x=reshape(x,1,[]);
f_shape=size(f);f=reshape([f(:).',zeros(1,ceil(numel(f)./prod([incfd])).*prod([incfd])-numel(f))],incfd,[]);
d_shape=size(d);d=reshape([d(:).',zeros(1,ceil(numel(d)./prod([incfd])).*prod([incfd])-numel(d))],incfd,[]);
%
%  DECLARE LOCAL VARIABLES.
%
if isempty(i), i=0; end;
if isempty(ibeg), ibeg=0; end;
if isempty(iend), iend=0; end;
if isempty(nless1), nless1=0; end;
if isempty(zero), zero=0; end;
if firstCall,   zero=[0.];  end;
firstCall=0;
%
%  VALIDITY-CHECK ARGUMENTS.
%
%***FIRST EXECUTABLE STATEMENT  PCHIC
if( n<2 )
%
%  ERROR RETURNS.
%
%     N.LT.2 RETURN.
ierr = -1;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','NUMBER OF DATA POINTS LESS THAN TWO',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
elseif( incfd<1 ) ;
%
%     INCFD.LT.1 RETURN.
ierr = -2;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','INCREMENT LESS THAN ONE',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
else;
for i = 2 : n;
if( x(i)<=x(i-1) )
%
%     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','X-ARRAY NOT STRICTLY INCREASING',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
end;
end; i = fix(n+1);
%
ibeg = fix(ic(1));
iend = fix(ic(2));
ierr = 0;
if( abs(ibeg)>5 )
ierr = fix(ierr - 1);
end;
if( abs(iend)>5 )
ierr = fix(ierr - 2);
end;
if( ierr<0 )
%
%     IC OUT OF RANGE RETURN.
ierr = fix(ierr - 3);
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','IC OUT OF RANGE',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
else;
%
%  function DEFINITION IS OK -- GO ON.
%
nless1 = fix(n - 1);
if( nwk<2.*nless1 )
%
%     NWK .LT. 2*(N-1)  RETURN.
ierr = -7;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','WORK ARRAY TOO SMALL',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
else;
%
%  SET UP H AND SLOPE ARRAYS.
%
for i = 1 : nless1;
wk(i) = x(i+1) - x(i);
wk(nless1+i) =(f(1,i+1)-f(1,i))./wk(i);
end; i = fix(nless1+1);
%
%  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
%
if( nless1>1 )
%
%  NORMAL CASE  (N .GE. 3) .
%
%
%  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
%
%     --------------------------------------
[n,dumvar2,dumvar3,d,incfd]=pchci(n,wk(sub2ind(size(wk),max(1,1)):end),wk(sub2ind(size(wk),max(n,1)):end),d,incfd);   dumvar2i=find((wk(sub2ind(size(wk),max(1,1)):end))~=(dumvar2));dumvar3i=find((wk(sub2ind(size(wk),max(n,1)):end))~=(dumvar3));   wk(1-1+dumvar2i)=dumvar2(dumvar2i); wk(n-1+dumvar3i)=dumvar3(dumvar3i); 
%     --------------------------------------
%
%  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
%
if( switchml~=zero )
%     ----------------------------------------------------
[switchml,n,dumvar3,dumvar4,d,incfd,ierr]=pchcs(switchml,n,wk(sub2ind(size(wk),max(1,1)):end),wk(sub2ind(size(wk),max(n,1)):end),d,incfd,ierr);   dumvar3i=find((wk(sub2ind(size(wk),max(1,1)):end))~=(dumvar3));dumvar4i=find((wk(sub2ind(size(wk),max(n,1)):end))~=(dumvar4));   wk(1-1+dumvar3i)=dumvar3(dumvar3i); wk(n-1+dumvar4i)=dumvar4(dumvar4i); 
%     ----------------------------------------------------
if( ierr~=0 )
%
%     ERROR RETURN FROM PCHCS.
ierr = -8;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','ERROR RETURN FROM PCHCS',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
end;
end;
else;
d(1,1) = wk(2);
d(1,n) = wk(2);
end;
%
%  SET END CONDITIONS.
%
if((ibeg~=0) ||(iend~=0) )
%     -------------------------------------------------------
[ic,vc,n,x,dumvar5,dumvar6,d,incfd,ierr]=pchce(ic,vc,n,x,wk(sub2ind(size(wk),max(1,1)):end),wk(sub2ind(size(wk),max(n,1)):end),d,incfd,ierr);   dumvar5i=find((wk(sub2ind(size(wk),max(1,1)):end))~=(dumvar5));dumvar6i=find((wk(sub2ind(size(wk),max(n,1)):end))~=(dumvar6));   wk(1-1+dumvar5i)=dumvar5(dumvar5i); wk(n-1+dumvar6i)=dumvar6(dumvar6i); 
%     -------------------------------------------------------
if( ierr<0 )
%
%     ERROR RETURN FROM PCHCE.
%   *** THIS CASE SHOULD NEVER OCCUR ***
ierr = -9;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHIC','ERROR RETURN FROM PCHCE',ierr,1);
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
end;
end;
%
%  NORMAL RETURN.
%
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
end;
end;
end;
%------------- LAST LINE OF PCHIC FOLLOWS ------------------------------
x_shape=zeros(x_shape);x_shape(:)=x(1:numel(x_shape));x=x_shape;
f_shape=zeros(f_shape);f_shape(:)=f(1:numel(f_shape));f=f_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
end
%DECK PCHID
