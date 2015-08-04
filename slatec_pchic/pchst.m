function [pchstresult,arg1,arg2]=pchst(arg1,arg2);

pchstresult=[];
persistent firstCall one zero ; if isempty(firstCall),firstCall=1;end; ;
%***BEGIN PROLOGUE  PCHST
%***SUBSIDIARY
%***PURPOSE  PCHIP Sign-Testing Routine
%***LIBRARY   SLATEC (PCHIP)
%***TYPE      SINGLE PRECISION (PCHST-S, DPCHST-D)
%***AUTHOR  Fritsch, F. N., (LLNL)
%***DESCRIPTION
%
%         PCHST:  PCHIP Sign-Testing Routine.
%
%     Returns:
%        -1. if ARG1 and ARG2 are of opposite sign.
%         0. if either argument is zero.
%        +1. if ARG1 and ARG2 are of the same sign.
%
%     The object is to do this without multiplying ARG1*ARG2, to avoid
%     possible over/underflow problems.
%
%  Fortran intrinsics used:  SIGN.
%
%***SEE ALSO  PCHCE, PCHCI, PCHCS, PCHIM
%***ROUTINES CALLED  (NONE)
%***REVISION HISTORY  (YYMMDD)
%   811103  DATE WRITTEN
%   820805  Converted to SLATEC library version.
%   870813  Minor cosmetic changes.
%   890411  Added SAVE statements (Vers. 3.2).
%   890411  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900328  Added TYPE section.  (WRB)
%   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
%   930503  Improved purpose.  (FNF)
%***end PROLOGUE  PCHST
%
%**end
%
%  DECLARE ARGUMENTS.
%
%
%  DECLARE LOCAL VARIABLES.
%
if isempty(one), one=0; end;
if isempty(zero), zero=0; end;
if firstCall,   zero=[0.];  end;
if firstCall,  one=[1.];  end;
firstCall=0;
%
%  PERFORM THE TEST.
%
%***FIRST EXECUTABLE STATEMENT  PCHST
pchstresult = (abs(one).*sign(arg1)).*(abs(one).*sign(arg2));
if((arg1==zero) ||(arg2==zero) )
pchstresult = zero;
end;
%
%------------- LAST LINE OF PCHST FOLLOWS ------------------------------
csnil=dbstack(1); csnil=csnil(1).name(1)~='@';
if csnil&&~isempty(inputname(2)), assignin('caller','FUntemp',arg2); evalin('caller',[inputname(2),'=FUntemp;']); end
if csnil&&~isempty(inputname(1)), assignin('caller','FUntemp',arg1); evalin('caller',[inputname(1),'=FUntemp;']); end
end
%DECK PCHSW
