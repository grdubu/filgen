function [g, R, g_alt] = chebyshevLCnorm(n, r, ripdB, suppress_warn_r)

% chebyshevLCnorm - Compute the normalized "g" (LC) values for a Chebyshev filter. 
%
% This function will compute the so-called g-values for a Chebyshev passband 
% prototype analog (passive) lowpass filter. The load (termination) of Port 1 is 
% always 1 Ω and the corner frequency is always 1 rad/s. These conditions make 
% this a normalized prototype. The g-values are LC (inductor and capacitor) 
% values of the analog filter in ladder form. 
% 
% The termination value of Port 2 does not need to be equal to that of Port 1. 
% (Under certain conditions they cannot be equal.) This feature is not always 
% present in popularly published algorithms. This is to say the design of 
% transducer loss in the filter can be greater than the minimum.
% 
% The recursion algorithm in the m-coding is interpreted from: 
% M. G. Ellis, Sr., /Electronic Filter Analysis and Synthesis/. Boston, MA: 
% Artech House, 1994, pp.99-102.
% who in turn references: 
% J. Porter, "Chebyshev Filters with Arbitrary Source and Load Resistances", 
% /RF Design Magazine/, August 1989, pp.63-68.
% 
% Syntax:
%         [g, R, g_alt] = chebyshevLCnorm(n, r, ripdB, suppress_warn_r)
% 
% n      := filter order
% r      := ratio of termination loads (resistance)
% ripdB  := filter ripple in dB
% g      := normalized L or C values (they alternate in a ladder)
% R      := normalized Resistance values (2 elements, reciprocals)
% g_alt  := alternate normalized L or C values when loads are mismatched 
% 
% Both single-ended or doubly-terminated design can be produced.
%
%   • r = 0 will produce a single-ended design where Port 2 is terminated with 
%     a voltage source. 
%   • r = inf will produce a single-ended design where Port 2 is terminated 
%     with a current source. 
% 
% Any other value of r (from the prior two cases) will produce a doubly-terminated 
% design for the r-ratio with one exception. The exception occurs when n is even 
% and the load ratio falls within the disallowed range of
%                           
%               tanh(beta/4)^2  <  r  <  coth(beta/4)^2,
% where
%               ep = sqrt( 10^(0.1*ripdB) - 1 ),
%               beta = 2*asinh(1/ep).
%
% If n is even and a disallowed value of r is entered, the function will 
% produce a warning, internally change r to tanh(beta/4)^2, and produce the  
% g-values accordingly. In this case R will be returned as 
%
%               R = [tanh(beta/4)^2, coth(beta/4)^2].
%
% Thus error-checking is left to the caller. That is, if the value of R(1) does 
% not equal the inputted value of r, then this indicates an invalid request 
% was made and the returned R values are the closest  valid values.
% 
% This feature is provided so the user can avoid determining the boundary values. 
% For example, the user enters n as even and also r = 1 while knowing this is 
% invalid but with the intent of getting the boundary values. The returned value 
% of R can then be used to recall the function and obtain the g-values for this 
% boundary case. In fact, these boundary values have already been provided on 
% the first call but not without a warning. To suppress the "first call" warning, 
% the parameter suppress_warn_r can be set to 1 or 'y'. Thus we can get the 
% boundary design without first needing to know the boundary values. Again, we do
% this by entering a known disallowed value for r. r=1 is guaranteed to be "wrong"  
% for an even order filter.
%
% The g_alt parameter is only different from g when the tranducer loss is 
% greater than the minimum for the Chebyshev passband. This difference occurs when 
% the filter is doubly terminated and
%
%     r =/= 1  &  n%2==1, or 
%     r < tanh(beta/4)^2  &  n%2==0, or  
%     r > coth(beta/4)^2  &  n%2==0. 
%
% "n%2" is simply the modulus-2 even-odd order detector. This added feature of
% providing g_alt is not part of the Ellis-Porter algorithm.
% 
% Some examples may be useful and are now performed.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Singly-terminated example:
% >> g = chebyshevLCnorm(3, 0, 1)
% g =                              
%        1.011796320945272         
%        1.333241652417787         
%        1.508847543107054         
% 
% From the g-values we can build two circuits:
% 
%     0 ├───R───┤  1 Ω                     0 ├───R───┤  1 Ω                
%       │       │                            │       │
%     1 ├───C───┤  1.011796320945272 F     1 │       L  1.011796320945272 H
%     2 │       L  1.333241652417787 H     2 ├───C───┤  1.333241652417787 F
%     3 ├───C───┤  1.508847543107054 F     3 │       L  1.508847543107054 H
%       │       │                            │       │
%       │ ┌───┐ │                            │ ┌───┐ │
%     4 ├─┤-V+├─┤  plot i(r0)/i(v4)        4 ├─┤-I+├─┤  plot v(r0)/v(i4)
%         └───┘                                └───┘  
% 
% Note that although we only asked for r=0, we actually get both because the 
% circuits are duals. 
% 
% See also: lpf_cheby_n3_1radps_r1dB_r0_00_00.asc (LTspice)
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Doubly-terminated example:
% >> g = chebyshevLCnorm(5, 1, 0.043648)
% g =
%       0.9732090774393125
%        1.372275955363987
%          1.8031709459949
%        1.372275955363987
%       0.9732090774393128
% 
% 
%     0 ├───R───┤  1 Ω                     0 ├───R───┤  1 Ω                
%       │       │                            │       │
%     1 ├───C───┤  .9732090774393125 F     1 │       L  .9732090774393125 H
%     2 │       L  1.372275955363987 H     2 ├───C───┤  1.372275955363987 F
%     3 ├───C───┤  1.803170945994900 F     3 │       L  1.803170945994900 H
%     4 │       L  1.372275955363987 H     4 ├───C───┤  1.372275955363987 F
%     5 ├───C───┤  .9732090774393128 F     5 │       L  .9732090774393128 H
%       │       │                            │       │
%     6 ├───R───┤  1 Ω                     6 ├───R───┤  1 Ω  
% 
% See also: lpf_cheby_n5_1radps_r20dB_rL1_00_00.asc (LTspice)
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Doubly-terminated example:
% >> [g, R] = chebyshevLCnorm(4, 1, 0.043648, 1) % get R; suppress disallowed load warning.
% g =
%       0.9332324828624304
%        1.292330756694625
%        1.579515172807973
%        0.763553944597803
% R =
%       0.8181819199604095       1.222222070182617
% 
%     0 ├───R───┤  1 Ω                     0 ├───R───┤  1 Ω                
%       │       │                            │       │
%     1 ├───C───┤  .9332324828624304 F     1 │       L  .9332324828624304 H
%     2 │       L  1.292330756694625 H     2 ├───C───┤  1.292330756694625 F
%     3 ├───C───┤  1.579515172807973 F     3 │       L  1.579515172807973 H
%     4 │       L  0.763553944597803 H     4 ├───C───┤  0.763553944597803 F
%       │       │                            │       │
%     5 ├───R───┤  .8181819199604095 Ω     5 ├───R───┤  1.222222070182617 Ω  
% 
% See also: lpf_cheby_n4_1radps_r20dB_rLx_00_00.asc (LTspice)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Doubly-terminated example:
% >> [g, R, g_alt] = chebyshevLCnorm(4, .25, 0.5); g, g_alt % also get alternate
% g =
%        5.029256047800175
%        0.464554325160485
%        5.637699415403508
%        0.2503632542676826
% 
% g_alt =
%        1.00145301707073
%        1.409424853850877
%        1.85821730064194
%        1.257314011950044
% 
%
%     0 ├───R───┤  1 Ω                     0 ├───R───┤  1 Ω                 
%       │       │                            │       │                     
%     1 ├───C───┤  5.029256047800175 F     1 │       L  5.029256047800175 F
%     2 │       L  0.464554325160485 H     2 ├───C───┤  0.464554325160485 H
%     3 ├───C───┤  5.637699415403508 F     3 │       L  5.637699415403508 F
%     4 │       L  0.250363254267683 H     4 ├───C───┤  0.250363254267683 H
%       │       │                            │       │                     
%     5 ├───R───┤  0.25 Ω                  5 ├───R───┤  4 Ω    
%             
%
%     0 ├───R───┤  1 Ω                     0 ├───R───┤  1 Ω                
%       │       │                            │       │
%     1 ├───C───┤  1.00145301707073 F      1 │       L  1.00145301707073 H
%     2 │       L  1.40942485385088 H      2 ├───C───┤  1.40942485385088 F
%     3 ├───C───┤  1.85821730064194 F      3 │       L  1.85821730064194 H
%     4 │       L  1.25731401195004 H      4 ├───C───┤  1.25731401195004 F
%       │       │                            │       │
%     5 ├───R───┤  0.25 Ω                  5 ├───R───┤  4 Ω  
% 
% The alternate realization is possible because the reflection zeros are not 
% on the jω-axis, unlike minimum loss filters for which all the zeros reside 
% on the jω-axis.
% 
% See also: lpf_cheby_n4_1radps_r0d5dB_rL4_00_00.asc (LTspice)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Copyright (c) 2025 Gregory R. White
% I claim it is mine to give to the public domain, which I do.
% If you do use this code, a credit to me is appreciated, but not required.
% 
% This is free and unencumbered software released into the public domain.
% 
% Anyone is free to copy, modify, publish, use, compile, sell, or
% distribute this software, either in source code form or as a compiled
% binary, for any purpose, commercial or non-commercial, and by any
% means.
% 
% In jurisdictions that recognize copyright laws, the author or authors
% of this software dedicate any and all copyright interest in the
% software to the public domain. We make this dedication for the benefit
% of the public at large and to the detriment of our heirs and
% successors. We intend this dedication to be an overt act of
% relinquishment in perpetuity of all present and future rights to this
% software under copyright law.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
% 
% For more information, please refer to <http://unlicense.org/>

if (round(n) ~= n)  |  (n <= 0)  |  (sum(size(n)) ~= 2)
    error('Filter order "n" must be a single positive integer')
end
if (r < 0)  |  (sum(size(n)) ~= 2)
    error('Load (resistor) ratio "r" must be a single positive number')
end
if (ripdB <= 0)  |  (sum(size(ripdB)) ~= 2)
    error('Ripple (dB) must be a single positive number')
end
if ~(exist('suppress_warn_r')==1)
    suppress_warn_r = [];
end
if isempty(suppress_warn_r)
    suppress_warn_r = 0;
end
if lower(suppress_warn_r(1)) == 'y'
    suppress_warn_r = 1;
end

Ap = 10^(ripdB/10);
ep = sqrt(Ap-1);
beta = 2*asinh(1/ep);
gamma = sinh(beta/(2*n));
pdn = pi/n;
k = 1:n;
a = sin((k - 0.5)*pdn);
g = zeros(n,1); 
if r==0  |  isinf(r) % Singly terminated, short (voltage source) or open (current source)
    b = (gamma^2 + sin(0.5*k*pdn).^2).*cos(0.5*k*pdn).^2; % b(end) not used
    g(1) = a(1)/gamma;
    for k0=2:n
        g(k0) = a(k0)*a(k0-1)/(b(k0-1)*g(k0-1));
    end
    if r==0 
        R = [0, inf];
    else
        R = [inf, 0];
    end
else % Doubly-terminated; was: elseif ~(r==0  |  isinf(r))
    rL_even_lo = tanh(beta/4)^2;
    rL_even_hi = coth(beta/4)^2;
    if ~rem(n,2) & (r > rL_even_lo & r < rL_even_hi ) % If true, then disallowed r was inputted.
        % Reassign r to rL_even_lo or rL_even_hi; either is fine.
        r = rL_even_lo;
        % r = rL_even_hi; 
        R = [rL_even_lo, rL_even_hi]; 
        if ~suppress_warn_r
            warning('r not valid between %.16f and %.16f; reverting to those limits.',rL_even_lo,rL_even_hi)
        end
    else
        R = [r, 1/r]; 
    end
    At = 4*r/(r + 1)^2;
    if ~rem(n,2)
        At = At*Ap; % tiny residual errors can make At slightly greater than 1.
    end
    if At >= 1
        At = 1; % clamped; "At" won't be needed in this case because we can assign d=0.
        d = 0;
    else
        d = sinh(asinh(sqrt(1 - At)/ep)/n); % note At>1 is intolerable.
    end
    b = gamma^2 + d^2 - 2*gamma*d*cos(k*pdn) + sin(k*pdn).^2; % b(end) not used
    g(1) = 2*a(1)/(gamma-d);
    for k0=2:n
        g(k0) = 4*a(k0)*a(k0-1)/(b(k0-1)*g(k0-1));
    end
end
g_alt = g; 
if nargout >= 3 & ~(r==0 | isinf(r)) & ~(rem(n,2) & r==1) & ~(~rem(n,2) & (r == rL_even_lo | r == rL_even_hi))  
    % g_alt: Alternate reflection zeros only meaningful when there is more than minimum loss.
    if r < 1
        g_alt(1:2:end) = r*g_alt(1:2:end); 
        g_alt(2:2:end) = (1/r)*g_alt(2:2:end); 
    elseif r > 1
        g_alt(1:2:end) = (1/r)*g_alt(1:2:end); 
        g_alt(2:2:end) = r*g_alt(2:2:end); 
    end
    g_alt = flipud(g_alt);
end

