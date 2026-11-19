function [a, b, k] = agm(a_, b_, d, N)
% AGM - Compute the arithmetic-geometric mean of two numbers. 
%
% arithmetic_geometric_mean(a_, b_)
% 
% [a, b] = agm(1,1/sqrt(2))
%
%
% http://www.ma.huji.ac.il/~tamarz/agm-talk.pdf
% a > b > 0
% >> [a, b, k] = agm(1, 1/sqrt(2))
% a = 0.8472130847939792
% b = 0.8472130847939791
% k = 3
% 
% >> M=1/1.202,x=agm(1,M),x=.5*pi/ellipke(1-M^2),K=pi/(2*agm(1,M)),K=ellipke(1-M^2)
% M = 0.8319467554076539
% x = 0.9140412453718602
% x = 0.9140412453718602
% K = 1.718517993305486
% K = 1.718517993305486
%
% >> M=1/1.202,x=agm(1,sqrt(1-M*M)),x=.5*pi/ellipke(M^2),K=pi/(2*agm(1,sqrt(1-M*M))),K=ellipke(M^2)
% M = 0.8319467554076539
% x = 0.7610700109767787
% x = 0.7610700109767787
% K = 2.063931444071607
% K = 2.063931444071607
%
%
% The above shows that the OCTAVE complete elliptic integral function 
% ELLIPKE is related to these as:
%
% ellipke(1-M^2)  ==  pi/(2*agm(1,M))
% ellipke(M^2)    ==  pi/(2*agm(1,sqrt(1-M*M)))
%
% M=1/sqrt(2), 2*agm(1,M)*ellipke(1-M^2) << cheats to find pi because ellipke is written including pi.
% ------------------------------------------------------------------
% Also relate to Daniels
% STEP 1.43 of PART 1 of Fig. 5.12
% XL = FH/FB;
% K = 1/XL; 
% KK   = daniels_fig5d12_part5_ellipk(K);                %  K(xₗ⁻¹)
%% >> K = .5,KK=daniels_fig5d12_part5_ellipk(K),KK=ellipke(K*K),KK=pi/(2*agm(1,sqrt(1-K*K)))
% K = 0.5
% KK = 1.685750354812596 % Daniels
% KK = 1.685750354812596 % Octave 5.1
% KK = 1.685750354812596 % agm()
% 
% >> K = 20/26,KK=daniels_fig5d12_part5_ellipk(K),KK=ellipke(K*K),KK=pi/(2*agm(1,sqrt(1-K*K)))
% K = 0.7692307692307693
% KK = 1.940716710940836 % Daniels
% KK = 1.940716710940837 % Octave 5.1
% KK = 1.940716710940837 % agm() 
% ------------------------------------------------------------------
% 
%
% >> 
% >> k=1/1.202,[x0,y0,count]=agm(1, k)
% x0 = 0.9140412453718602
% y0 = 0.9140412453718602
% count = 3
% >> [x1,y1,count]=agm(1, sqrt(1 - k*k))
% x1 = 0.7610700109767787
% y1 = 0.7610700109767786
% count = 3
% >> qofk = exp(1)^(-pi*x1/x0)
% qofk = 0.07310744898944521
% 
%
% Stop at 1 iteration:
% >> [x,y]=agm(1, 1/sqrt(2), [], 1)
% x = 0.8535533905932737
% y = 0.8408964152537145
% Stop at 2 iterations:
% >> [x,y]=agm(1, 1/sqrt(2), [], 2)
% x = 0.8472249029234942
% y = 0.8472012667468914
% Stop at 3 iterations:
% >> [x,y]=agm(1, 1/sqrt(2), [], 3)
% x = 0.8472130848351929
% y = 0.8472130847527654
% Stop at 4 iterations:
% >> [x,y]=agm(1, 1/sqrt(2), [], 4)
% x = 0.8472130847939792
% y = 0.8472130847939791
% Only 4 iterations are needed for default convergence value:
% >> [x,y,k]=agm(1, 1/sqrt(2), [], 5)
% x = 0.8472130847939792
% y = 0.8472130847939791
% k = 4
% Only 4 iterations are needed for default convergence value:
% >> [x,y,k]=agm(1, 1/sqrt(2), [], 5)
% x = 0.8472130847939792
% y = 0.8472130847939791
% k = 4
% A "too low" convergence value will cause a loop break due to maximum count.
% >> [x,y,k]=agm(1, 1/sqrt(2), eps-.0000000000000000000000000000001)
% warning: Iteration limit of 100 reached. Check convergence.
% warning: called from
    % agm at line 140 column 5
% x = 0.8472130847939792
% y = 0.8472130847939791
% k = 100
% y = 0.8472130847939791
% k = 100
% But in this case, eps is a big enough convergence value to stop looping 
% before maximum count:
% >> [x,y,k]=agm(1, 1/sqrt(2), eps)
% x = 0.8472130847939792
% y = 0.8472130847939791
% k = 4
% In 4 iterations it has converged as much as possible for double precision math.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% http://rosettacode.org/wiki/Arithmetic-geometric_mean#MATLAB_.2F_Octave
%%while (1)
%%    a_ = a;
%%    a = (a_ + b)/2;
%%    b = sqrt(a_*b);
%%    if (abs(a_ - a) < a*eps) break; end; % This one simply checks if a is not changing much anymore.
%%end
%
% See also: dimo_ex4d2_p179_00


% do not swap
if (a_ <= b_) |  (b_ <= 0)
    error('Inputs a and b of agm(a, b) must be greater than zero and a > b.')
end

if ~(exist('d')==1)
    d = [];
end
if isempty(d)
    d = 4*eps; % default 1-ratio difference level, meaning if the two are are essentially converged, quit iterating.
end
if ~(exist('N')==1)
    N = [];
end
if isempty(N)
    N = 100;   % default loop exit level
end

k = 0;
while(k < N)   % quit after N iterations (k goes from 0:N-1 and the loop does not perform k = N)
   a = (a_ + b_)/2; % arithmetic mean
   b = sqrt(a_*b_); % geometic mean
   k = k + 1; % Computation done, inc by 1.
   if abs(1 - a/b) <= d  % This one checks if a and b are about the same. Not sure what happens if a and b are huge.
                         % The "rosetta code" version I started with uses a different method: it asks if "a" is no longer changing...
      break
   end
   a_ = a;
   b_ = b;
end
if k==N
    warning(sprintf('Iteration limit of %d reached. Check convergence.', N))
end

 % elfun18 version 1.2.0.4 (18.8 MB) by milan batista 
 % elfun18   milan batista 