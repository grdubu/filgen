function [a, b, k] = agm(a_, b_, d, N)
% AGM - Compute the arithmetic-geometric mean of two numbers. 
%
% Syntax:  [a, b, k] = agm(a_, b_, d, N)
% 
%  a_ := an integer greater than 1 and greater than b_.
%  b_ := an integer greater than 1 and less than a_.
%  a  := converged arithmetic mean.
%  b  := converged geometric mean; effectively equal to "a" when 
%        convergent. (It is 'within the ratiometric "error" value 
%        "d" of "a".)
%  d  := ratiometric allowed distance between "a" and "b". 
%        Default is:    d = 4*eps     
%        Satisfactory convergence condition:   (abs(1 - a/b)) <= d 
%                                               TRUE is logical 1
%        Note |a/b| is a ratio that should approach unity. 
%  N  := number of iterations allowed before search terminates.
%        Default is:    N = 100
%        A warning is produced if the limit is reached.
%  k  := number of iterations performed during search.
% 
% Conditions: 
%              a_ > 0
%              b_ > 0
%              a_ > b_
% 
% 
% Example 1: 
% >> [a, b, k] = agm(1, 1/sqrt(2))
% a = 0.8472130847939792
% b = 0.8472130847939791
% k = 3
% 
% 
% Example 2: Demonstrate relationship of agm() to built-in function 
%            ellipke() using test values.
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
% The above shows that the OCTAVE complete elliptic integral function 
% ELLIPKE is related to these as:
%
% ellipke(1-M^2)  ==  pi/(2*agm(1,M))
% ellipke(M^2)    ==  pi/(2*agm(1,sqrt(1-M*M)))
% 
% 
% Example 3: Set a lower N iteration allowance while not adjusting d.
% >> [x,y,k]=agm(1, 1/sqrt(2), [], 5)
% x = 0.8472130847939792
% y = 0.8472130847939791
% k =                4
% 
% 
% See also:  1. http://www.ma.huji.ac.il/~tamarz/agm-talk.pdf
%            2. http://rosettacode.org/wiki/Arithmetic-geometric_mean#MATLAB_.2F_Octave
%            3. p.170 of Dimopoulos, Hercules G., Analog Electronic Filters - 
%               Theory, Design and Synthesis, © Springer 2012, New York, 
%               ISBN 978-94-007-2189-0, DOI 10.1007/978-94-007-2190-6


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
   if abs(1 - a/b) <= d  % This one checks if a and b are about the same. What if a and b are huge?
      break
   end
   a_ = a;
   b_ = b;
end
if k==N
    warning(sprintf('Iteration limit of %d reached. Check convergence.', N))
end

% http://rosettacode.org/wiki/Arithmetic-geometric_mean#MATLAB_.2F_Octave
% uses a different compare method:
%%while (1)
%%    a_ = a;
%%    a = (a_ + b)/2;
%%    b = sqrt(a_*b);
%%    if (abs(a_ - a) < a*eps) break; end; % This one simply checks if a is not changing much anymore.
%%end
