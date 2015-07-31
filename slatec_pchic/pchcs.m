function [switchml,n,h,slope,d,incfd,ierr]=pchcs(switchml,n,h,slope,d,incfd,ierr);
%***BEGIN PROLOGUE  PCHCS
%***SUBSIDIARY
%***PURPOSE  Adjusts derivative values for PCHIC
%***LIBRARY   SLATEC (PCHIP)
%***TYPE      SINGLE PRECISION (PCHCS-S, DPCHCS-D)
%***AUTHOR  Fritsch, F. N., (LLNL)
%***DESCRIPTION
%
%         PCHCS:  PCHIC Monotonicity Switch Derivative Setter.
%
%     Called by  PCHIC  to adjust the values of D in the vicinity of a
%     switch in direction of monotonicity, to produce a more 'visually
%     pleasing' curve than that given by  PCHIM .
%
% ----------------------------------------------------------------------
%
%  Calling sequence:
%
%        PARAMETER  (INCFD = ...)
%        INTEGER  N, IERR
%        REAL  SWITCH, H(N), SLOPE(N), D(INCFD,N)
%
%        CALL  PCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
%
%   Parameters:
%
%     SWITCH -- (input) indicates the amount of control desired over
%           local excursions from data.
%
%     N -- (input) number of data points.  (assumes N.GT.2 .)
%
%     H -- (input) real array of interval lengths.
%     SLOPE -- (input) real array of data slopes.
%           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
%                  H(I) =  X(I+1)-X(I),
%              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
%
%     D -- (input) real array of derivative values at the data points,
%           as determined by PCHCI.
%          (output) derivatives in the vicinity of switches in direction
%           of monotonicity may be adjusted to produce a more 'visually
%           pleasing' curve.
%           The value corresponding to X(I) is stored in
%                D(1+(I-1)*INCFD),  I=1(1)N.
%           No other entries in D are changed.
%
%     INCFD -- (input) increment between successive values in D.
%           This argument is provided primarily for 2-D applications.
%
%     IERR -- (output) error flag.  should be zero.
%           If negative, trouble in PCHSW.  (should never happen.)
%
%    -------
%    WARNING:  This routine does no validity-checking of arguments.
%    -------
%
%  Fortran intrinsics used:  ABS, MAX, MIN.
%
%***SEE ALSO  PCHIC
%***ROUTINES CALLED  PCHST, PCHSW
%***REVISION HISTORY  (YYMMDD)
%   820218  DATE WRITTEN
%   820617  Redesigned to (1) fix  problem with lack of continuity
%           approaching a flat-topped peak (2) be cleaner and
%           easier to verify.
%           Eliminated subroutines PCHSA and PCHSX in the process.
%   820622  1. Limited fact to not exceed one, so computed D is a
%             convex combination of PCHCI value and PCHSD value.
%           2. Changed fudge from 1 to 4 (based on experiments).
%   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
%   820805  Converted to SLATEC library version.
%   870813  Minor cosmetic changes.
%   890411  Added SAVE statements (Vers. 3.2).
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890831  Modified array declarations.  (WRB)
%   890831  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900328  Added TYPE section.  (WRB)
%   910408  Updated AUTHOR section in prologue.  (WRB)
%   930503  Improved purpose.  (FNF)
%***end PROLOGUE  PCHCS
%
%  Programming notes:
%     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
%        either argument is zero, +1 if they are of the same sign, and
%        -1 if they are of opposite sign.
%**end
%
%  DECLARE ARGUMENTS.
%
persistent del dext dfloc dfmx fact firstCall fudge h1 h2 i indx k nless1 one s1 s2 slmax wtave zero ; if isempty(firstCall),firstCall=1;end; 

h_shape=size(h);h=reshape(h,1,[]);
slope_shape=size(slope);slope=reshape(slope,1,[]);
d_shape=size(d);d=reshape([d(:).',zeros(1,ceil(numel(d)./prod([incfd])).*prod([incfd])-numel(d))],incfd,[]);
%
%  DECLARE LOCAL VARIABLES.
%
if isempty(i), i=0; end;
if isempty(indx), indx=0; end;
if isempty(k), k=0; end;
if isempty(nless1), nless1=0; end;
if isempty(del), del=zeros(1,3); end;
if isempty(dext), dext=0; end;
if isempty(dfloc), dfloc=0; end;
if isempty(dfmx), dfmx=0; end;
if isempty(fact), fact=0; end;
if isempty(fudge), fudge=0; end;
if isempty(one), one=0; end;
if isempty(slmax), slmax=0; end;
if isempty(wtave), wtave=zeros(1,2); end;
if isempty(zero), zero=0; end;
%
%  DEFINE INLINE FUNCTION FOR WEIGHTED AVERAGE OF SLOPES.
%
% pchsd= @(s1,s2,h1,h2) (h2/(h1+h2))*s1 +(h1/(h1+h2))*s2;real :: pchsd ;
if isempty(s1), s1=0; end;
if isempty(s2), s2=0; end;
if isempty(h1), h1=0; end;
if isempty(h2), h2=0; end;
%
%  INITIALIZE.
%
if firstCall,   zero=[0.];  end;
if firstCall,  one=[1.];  end;
if firstCall,   fudge=[4.];  end;
firstCall=0;
pchsd= @(s1,s2,h1,h2) (h2./(h1+h2)).*s1 +(h1./(h1+h2)).*s2;
%***FIRST EXECUTABLE STATEMENT  PCHCS
ierr = 0;
nless1 = fix(n - 1);
%
%  LOOP OVER SEGMENTS.
%
for i = 2 : nless1;
if( pchst(slope(i-1),slope(i))<0 )
%             --------------------------
%
%
%....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
%
%           DO NOT CHANGE D IF 'UP-DOWN-UP'.
if( i>2 )
if( pchst(slope(i-2),slope(i))>zero )
continue;
end;
%                   --------------------------
end;
if( i<nless1 )
if( pchst(slope(i+1),slope(i-1))>zero )
continue;
end;
%                   ----------------------------
end;
%
%   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
%
pchsd= @(s1,s2,h1,h2) (h2./(h1+h2)).*s1 +(h1./(h1+h2)).*s2;
dext = pchsd(slope(i-1),slope(i),h(i-1),h(i));
%
%   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
%
if( pchst(dext,slope(i-1))<0 )
%                -----------------------
%
%              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
%                        EXTREMUM IS IN (X(I-1),X(I)).
k = fix(i - 1);
%              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
wtave(2) = dext;
pchsd= @(s1,s2,h1,h2) (h2./(h1+h2)).*s1 +(h1./(h1+h2)).*s2;
if( k>1 )
wtave(1) = pchsd(slope(k-1),slope(k),h(k-1),h(k));
end;
elseif( pchst(dext,slope(i-1))==0 ) ;
continue;
else;
%
%              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
%                        EXTREMUM IS IN (X(I),X(I+1)).
k = fix(i);
%              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
wtave(1) = dext;
pchsd= @(s1,s2,h1,h2) (h2./(h1+h2)).*s1 +(h1./(h1+h2)).*s2;
if( k<nless1 )
wtave(2) = pchsd(slope(k),slope(k+1),h(k),h(k+1));
end;
end;
elseif( pchst(slope(i-1),slope(i))==0 ) ;
%
%
%....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
%                     CHECK FOR FLAT-TOPPED PEAK .......................
%
if( i==nless1 )
continue;
end;
if( pchst(slope(i-1),slope(i+1))>=zero )
continue;
end;
%                -----------------------------
%
%           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
k = fix(i);
%           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
pchsd= @(s1,s2,h1,h2) (h2./(h1+h2)).*s1 +(h1./(h1+h2)).*s2;
wtave(1) = pchsd(slope(k-1),slope(k),h(k-1),h(k));
pchsd= @(s1,s2,h1,h2) (h2./(h1+h2)).*s1 +(h1./(h1+h2)).*s2;
wtave(2) = pchsd(slope(k),slope(k+1),h(k),h(k+1));
else;
continue;
end;
%
%
%....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
%        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
%           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
%                    IF K.GT.1
%           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
%                    IF K.LT.N-1
%
slmax = abs(slope(k));
if( k>1 )
slmax = max(slmax,abs(slope(k-1)));
end;
if( k<nless1 )
slmax = max(slmax,abs(slope(k+1)));
end;
%
if( k>1 )
del(1) = slope(k-1)./slmax;
end;
del(2) = slope(k)./slmax;
if( k<nless1 )
del(3) = slope(k+1)./slmax;
end;
%
if((k>1) &&(k<nless1) )
%           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
fact = fudge.*abs(del(3).*(del(1)-del(2)).*(wtave(2)./slmax));
d(1,k) = d(1,k) + min(fact,one).*(wtave(1)-d(1,k));
fact = fudge.*abs(del(1).*(del(3)-del(2)).*(wtave(1)./slmax));
d(1,k+1) = d(1,k+1) + min(fact,one).*(wtave(2)-d(1,k+1));
else;
%           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
%                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
fact = fudge.*abs(del(2));
d(1,i) = min(fact,one).*wtave(i-k+1);
%              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
%                        I-K+1 = 2 IF K=I-1(=1).
end;
%
%
%....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
%
if( switchml>zero )
%
dfloc = h(k).*abs(slope(k));
if( k>1 )
dfloc = max(dfloc,h(k-1).*abs(slope(k-1)));
end;
if( k<nless1 )
dfloc = max(dfloc,h(k+1).*abs(slope(k+1)));
end;
dfmx = switchml.*dfloc;
indx = fix(i - k + 1);
%        INDX = 1 IF K=I, 2 IF K=I-1.
%        ---------------------------------------------------------------
[dfmx,indx,d(1,k),d(1,k+1),h(k),slope(k),ierr]=pchsw(dfmx,indx,d(1,k),d(1,k+1),h(k),slope(k),ierr);
%        ---------------------------------------------------------------
if( ierr~=0 )
h_shape=zeros(h_shape);h_shape(:)=h(1:numel(h_shape));h=h_shape;
slope_shape=zeros(slope_shape);slope_shape(:)=slope(1:numel(slope_shape));slope=slope_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
return;
end;
end;
%
%....... end OF SEGMENT LOOP.
%
end; i = fix(nless1+1);
%
%------------- LAST LINE OF PCHCS FOLLOWS ------------------------------
h_shape=zeros(h_shape);h_shape(:)=h(1:numel(h_shape));h=h_shape;
slope_shape=zeros(slope_shape);slope_shape(:)=slope(1:numel(slope_shape));slope=slope_shape;
d_shape=zeros(d_shape);d_shape(:)=d(1:numel(d_shape));d=d_shape;
end
%DECK PCHDF
