function g = butterworthLCnorm(n,r)

%BUTTERWORTHLCNORM  Compute the normalized "g" (LC) values for a Butterworth lowpass filter.
% <---   To do: Option for a "not 3 dB" corner frequency.   --->
% Syntax: 
%                 g = butterworthlcnorm(n, r)
%                 g = butterworthlcnorm(n, 1/r)
% 
%
% These are the "g" values for a Butterworth LC lowpass prototype ladder 
% type filter with a 3.01 dB corner frequency of w = 1 (rad/s). 
%
% If doubly terminated but not equal (r =/= 1), then entering r or 
% 1/r produce alternate realizations. (See the examples.) 
%
% The recursion algorithm in the m-coding is interpreted from: 
% M. G. Ellis, Sr., /Electronic Filter Analysis and Synthesis/. Boston, MA: 
% Artech House, 1994, pp.99-102.
% who in turn references: 
% J. Porter, "Chebyshev Filters with Arbitrary Source and Load Resistances", 
% /RF Design Magazine/, August 1989, pp.63-68.
% 
% Fig. 7.5 on p.102 of Ellis has an error in the Butterworth equation for b_k.
% It should be b_k = cos(.5*k*pi/n)^2, not b_k = cos(k*pi/n)^2.
% 
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Some examples now follow: 

% >> g = butterworthLCnorm(3,0),g = butterworthLCnorm(3,inf)
% g =
%       0.4999999999999999
%        1.333333333333333
%                      1.5
% g =
%       0.4999999999999999
%        1.333333333333333
%                      1.5
% 
%   0 ├───R───┤  1                   0 ├───R───┤  1               
%     │       │                        │       │
%   1 ├───C───┤  1.011796320945272   1 │       L  1.011796320945272
%   2 │       L  1.333241652417787   2 ├───C───┤  1.333241652417787
%   3 ├───C───┤  1.508847543107054   3 │       L  1.508847543107054
%     │       │                        │       │
%     │ ┌───┐ │  r = 0                 │ ┌───┐ │  r = inf
%   4 ├─┤-V+├─┤  plot i(r0)/i(v4)    4 ├─┤-I+├─┤  plot v(r0)/v(i4)
%       └───┘                            └───┘      
% 
%  Example LTspice: lpf_butter_n3_1radps_rL0_00_00.asc
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% >> g = butterworthLCnorm(4,0),g = butterworthLCnorm(4,inf)
% g =
%       0.3826834323650898
%        1.082392200292394
%        1.577161014949475
%        1.530733729460359
% g =
%       0.3826834323650898
%        1.082392200292394
%        1.577161014949475
%        1.530733729460359
% 
%   g                                g
%     ├───R───┤  1.0                   ├───R───┤  1.0
%   1 ├───C───┤  0.382683432365090   1 │       L  0.382683432365090
%   2 │       L  1.082392200292394   2 ├───C───┤  1.082392200292394
%   3 ├───C───┤  1.577161014949475   3 │       L  1.577161014949475
%   4 │       L  1.530733729460359   4 ├───C───┤  1.530733729460359
%     │       │                        │       │
%     │ ┌───┐ │  r = 0                 │ ┌───┐ │  r = inf
%   5 ├─┤-V+├─┤  plot i(r0)/i(v4)    5 ├─┤-I+├─┤  plot v(r0)/v(i4)
%       └───┘                            └───┘      
%
%  Example LTspice: lpf_butter_n4_1radps_rL1_00_00.asc 
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% >> g = butterworthLCnorm(4,1)
% g =
%        .7653668647301796
%        1.847759065022574
%        1.847759065022573
%        .7653668647301798
% 
%  g                                g
%    ├───R───┤  1.0                   ├───R───┤  1.0
%  1 ├───C───┤  .7653668647301796   1 │       L  .7653668647301796
%  2 │       L  1.847759065022574   2 ├───C───┤  1.847759065022574
%  3 ├───C───┤  1.847759065022573   3 │       L  1.847759065022573
%  4 │       L  .7653668647301798   4 ├───C───┤  .7653668647301798
%    ├───R───┤  1.0                   ├───R───┤  1.0     
% 
%  Example LTspice: lpf_butter_n4_1radps_rL1_00_00.asc
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% >> g = butterworthLCnorm(4,.5)
% g =
%        3.186846750345509
%        0.882623595140887
%        2.452375708643579
%        0.2174540699936987     
%  g                                g
%    ├───R───┤  1.0                   ├───R───┤  1.0
%  1 ├───C───┤  3.186846750345509   1 │       L  3.186846750345509
%  2 │       L  0.882623595140887   2 ├───C───┤  0.882623595140887
%  3 ├───C───┤  2.452375708643579   3 │       L  2.452375708643579
%  4 │       L  0.217454069993699   4 ├───C───┤  0.217454069993699
%    ├───R───┤  0.5                   ├───R───┤  2.0
%       
% and 
% 
% >> g = butterworthLCnorm(4,2)
% g =
%        0.4349081399873973
%        1.226187854321789
%        1.765247190281774
%        1.593423375172754
%  g                                g
%    ├───R───┤  1.0                   ├───R───┤  1.0
%  1 ├───C───┤  0.434908139987397   1 │       L  0.434908139987397
%  2 │       L  1.226187854321789   2 ├───C───┤  1.226187854321789
%  3 ├───C───┤  1.765247190281774   3 │       L  1.765247190281774
%  4 │       L  1.593423375172754   4 ├───C───┤  1.593423375172754
%    ├───R───┤  0.5                   ├───R───┤  2.0
% 
%  Example LTspice: lpf_butter_n4_1radps_rL2_00_00.asc
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Doing the same for an odd order filter, we have
% 
% >> g = butterworthlcnorm(5,1/2)
% g =
%           3.13311812799952
%          0.923711519119942
%           3.05095872940666
%          0.495521963437276
%          0.685660109978754
%  g                                g
%    ├───R───┤  1.0                   ├───R───┤  1.0
%  1 ├───C───┤  3.13311812799952    1 │       L  3.13311812799952
%  2 │       L  0.923711519119942   2 ├───C───┤  0.923711519119942
%  3 ├───C───┤  3.05095872940666    3 │       L  3.05095872940666
%  4 │       L  0.495521963437276   4 ├───C───┤  0.495521963437276
%  5 ├───C───┤  0.685660109978754   5 │       L  0.685660109978754
%    ├───R───┤  0.5                   ├───R───┤  2.0
%       
% and 
%  
%  >> g = butterworthlcnorm(5,2)
% g =
%          0.342830054989377
%          0.991043926874551
%           1.52547936470333
%           1.84742303823988
%           1.56655906399976
%  g                                g
%    ├───R───┤  1.0                   ├───R───┤  1.0
%  1 ├───C───┤  0.342830054989377   1 │       L  0.342830054989377
%  2 │       L  0.991043926874551   2 ├───C───┤  0.991043926874551
%  3 ├───C───┤  1.52547936470333    3 │       L  1.52547936470333
%  4 │       L  1.84742303823988    4 ├───C───┤  1.84742303823988
%  5 ├───C───┤  1.56655906399976    5 │       L  1.56655906399976
%    ├───R───┤  2.0                   ├───R───┤  0.5
%     
%  Example LTspice: lpf_butter_n5_1radps_rL2_00_00.asc 
%      
% Within the respective groupings, the above networks have identical 
% transmission characteristics. The diagrams also give an idea how the 
% "g" values are to be used in assigning them to real networks.  
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
    error('Termination ratio (load resistor) "r" must be a single positive number')
end
pdn = pi/n;
k = 1:n;
a = sin((k - 0.5)*pdn);
g = zeros(n,1); 
if r==0  |  isinf(r) % Singly terminated, short (voltage source) or open (current source) 
    b = cos(0.5*k*pdn).^2; % b(end) not used; Error in Ellis. DEN is 2n, not n.
    g(1) = a(1);
    for k0=2:n
        g(k0) = a(k0)*a(k0-1)/(b(k0-1)*g(k0-1));
    end
else % Doubly-terminated
    At = 4*r/(r + 1)^2;
    if At >= 1
        At = 1; % clamped; "At" won't be needed in this case because we can assign d=0.
        d = 0;
    else
        d = (1 - At)^(0.5/n); % note At>1 is intolerable.
    end
    b = 1 + d^2 - 2*d*cos(k*pdn); % b(end) not used
    g(1) = 2*a(1)/(1 - d);
    for k0=2:n
        g(k0) = 4*a(k0)*a(k0-1)/(b(k0-1)*g(k0-1));
    end
end
if (r>1) & ~isinf(r)
    % Provides the dual for the unequal resistance. Amounts to z-scale and flipping the network.
    g(1:2:end) = (1/r)*g(1:2:end); 
    g(2:2:end) = r*g(2:2:end); 
    g = flipud(g);
end
