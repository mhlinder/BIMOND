function [dfmax,iextrm,d1,d2,h,slope,ierr]=pchsw(dfmax,iextrm,d1,d2,h,slope,ierr);
%***BEGIN PROLOGUE  PCHSW
%***SUBSIDIARY
%***PURPOSE  Limits excursion from data for PCHCS
%***LIBRARY   SLATEC (PCHIP)
%***TYPE      SINGLE PRECISION (PCHSW-S, DPCHSW-D)
%***AUTHOR  Fritsch, F. N., (LLNL)
%***DESCRIPTION
%
%         PCHSW:  PCHCS Switch Excursion Limiter.
%
%     Called by  PCHCS  to adjust D1 and D2 if necessary to insure that
%     the extremum on this interval is not further than DFMAX from the
%     extreme data value.
%
% ----------------------------------------------------------------------
%
%  Calling sequence:
%
%        INTEGER  IEXTRM, IERR
%        REAL  DFMAX, D1, D2, H, SLOPE
%
%        CALL  PCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
%
%   Parameters:
%
%     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
%           the cubic determined by derivative values D1,D2.  (assumes
%           DFMAX.GT.0.)
%
%     IEXTRM -- (input) index of the extreme data value.  (assumes
%           IEXTRM = 1 or 2 .  Any value .NE.1 is treated as 2.)
%
%     D1,D2 -- (input) derivative values at the ends of the interval.
%           (Assumes D1*D2 .LE. 0.)
%          (output) may be modified if necessary to meet the restriction
%           imposed by DFMAX.
%
%     H -- (input) interval length.  (Assumes  H.GT.0.)
%
%     SLOPE -- (input) data slope on the interval.
%
%     IERR -- (output) error flag.  should be zero.
%           If IERR=-1, assumption on D1 and D2 is not satisfied.
%           If IERR=-2, quadratic equation locating extremum has
%                       negative discriminant (should never occur).
%
%    -------
%    WARNING:  This routine does no validity-checking of arguments.
%    -------
%
%  Fortran intrinsics used:  ABS, SIGN, SQRT.
%
%***SEE ALSO  PCHCS
%***ROUTINES CALLED  R1MACH, XERMSG
%***REVISION HISTORY  (YYMMDD)
%   820218  DATE WRITTEN
%   820805  Converted to SLATEC library version.
%   870707  Replaced DATA statement for SMALL with a use of R1MACH.
%   890411  1. Added SAVE statements (Vers. 3.2).
%           2. Added REAL R1MACH for consistency with D.P. version.
%   890411  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
%   900328  Added TYPE section.  (WRB)
%   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
%   920526  Eliminated possible divide by zero problem.  (FNF)
%   930503  Improved purpose.  (FNF)
%***end PROLOGUE  PCHSW
%
%**end
%
%  DECLARE ARGUMENTS.
%
%
%  DECLARE LOCAL VARIABLES.
%
persistent cp fact firstCall hphi lambda nu one phi radcal rho sigma small that third three two zero ; if isempty(firstCall),firstCall=1;end; 

if isempty(cp), cp=0; end;
if isempty(fact), fact=0; end;
if isempty(hphi), hphi=0; end;
if isempty(lambda), lambda=0; end;
if isempty(nu), nu=0; end;
if isempty(one), one=0; end;
if isempty(phi), phi=0; end;
if isempty(radcal), radcal=0; end;
if isempty(rho), rho=0; end;
if isempty(sigma), sigma=0; end;
if isempty(small), small=0; end;
if isempty(that), that=0; end;
if isempty(third), third=0; end;
if isempty(three), three=0; end;
if isempty(two), two=0; end;
if isempty(zero), zero=0; end;
%
if firstCall,   zero=[0.];  end;
if firstCall,  one=[1.];  end;
if firstCall,  two=[2.];  end;
if firstCall,  three=[3.];  end;
if firstCall,  fact=[100.];  end;
%        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
if firstCall,   third=[0.33333];  end;
firstCall=0;
%
%  NOTATION AND GENERAL REMARKS.
%
%     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
%     LAMBDA IS THE RATIO OF D2 TO D1.
%     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
%     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
%           WHERE  THAT = (XHAT - X1)/H .
%        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
%     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
%
%      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
%***FIRST EXECUTABLE STATEMENT  PCHSW
small = fact.*r1mach(4);
%
%  DO MAIN CALCULATION.
%
while (1);
if( d1~=zero )
%
rho = slope./d1;
lambda = -d2./d1;
if( d2==zero )
%
%           SPECIAL CASE -- D2.EQ.ZERO .
%
%             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
if( rho>=third )
%
%  NORMAL RETURN.
%
ierr = 0;
return;
else;
cp = two - three.*rho;
nu = one - two.*rho;
that = one./(three.*nu);
end;
elseif( lambda<=zero ) ;
break;
else;
%
%           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
%
nu = one - lambda - two.*rho;
sigma = one - rho;
cp = nu + sigma;
if( abs(nu)>small )
radcal =(nu-(two.*rho+one)).*nu + sigma.^2;
if( radcal<zero )
%
%     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
ierr = -2;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHSW','NEGATIVE RADICAL',ierr,1);
return;
else;
that =(cp-sqrt(radcal))./(three.*nu);
end;
else;
that = one./(two.*sigma);
end;
end;
phi = that.*((nu.*that-cp).*that+one);
%
%          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
if( iextrm~=1 )
phi = phi - rho;
end;
%
%          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
hphi = h.*abs(phi);
if( hphi.*abs(d1)>dfmax )
%           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
d1 = (abs(dfmax./hphi).*sign(d1));
d2 = -lambda.*d1;
end;
ierr = 0;
return;
%
%        SPECIAL CASE -- D1.EQ.ZERO .
%
%          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
elseif( d2~=zero ) ;
%
rho = slope./d2;
%          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO .GE. 1/3 .
if( rho<third )
that =(two.*(three.*rho-one))./(three.*(two.*rho-one));
phi = that.^2.*((three.*rho-one)./three);
%
%          CONVERT TO DISTANCE FROM F2 IF IEXTRM.NE.1 .
if( iextrm~=1 )
phi = phi - rho;
end;
%
%          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
hphi = h.*abs(phi);
%           AT THIS POINT, HPHI.GT.0, SO DIVIDE IS OK.
if( hphi.*abs(d2)>dfmax )
d2 = (abs(dfmax./hphi).*sign(d2));
end;
end;
ierr = 0;
return;
end;
break;
end;
%
%  ERROR RETURNS.
%
%     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
ierr = -1;
[dumvar1,dumvar2,dumvar3,ierr]=xermsg('SLATEC','PCHSW','D1 AND/OR D2 INVALID',ierr,1);
return;
%------------- LAST LINE OF PCHSW FOLLOWS ------------------------------
end %subroutine pchsw
%DECK PCOEF
