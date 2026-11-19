function [z,p,k] = ellipap1(n, rp, rs)
%ELLIPAP1 Elliptic analog lowpass filter prototype
% [Z,P,K] = ELLIPAP1(N, Rp, Rs) returns the zeros, poles, and 
% gain of an N-th order normalized prototype elliptic analog 
% lowpass filter with Rp decibels of ripple in the passband and 
% a stopband Rs decibels down.
% 
% This program is a faster and more accurate replacement for
% the standard Octave (or matlab) program ELLIPAP.
% 
% Authors: H. J. Orchard and A. N. Willson, Jr., 9-25-95
% Copyright (c) 1995
% 
% Reference (m-code is given):
% [1] H. J. Orchard and A. N. Willson, Jr., "Elliptic Functions
% for Filter Design," IEEE Transactions on Circuits and
% Systems--I: Fundamental Theory and Applications, April 1997.
if n == 1, % special case; for n==1, reduces to Chebyshev type 1
    z = [];
    p = -sqrt(1/(10^(rp/10) - 1));
    k = -p;
    return
end 
dbn = log(10)/20;
no  = rem(n,2); 
n3  = (n - no)/2; 
apn = dbn*rp; 
asn = dbn*rs; 
e(1) = sqrt(2*exp(apn)*sinh(apn));
g(1) = e(1)/sqrt(exp(2*asn) - 1);
v = g(1);
m2 = 1;
while v > 1.0e-150
    v = (v/(1 + sqrt(1 - v^2)))^2; 
    m2 = m2 + 1;
    g(m2) = v;
end
for index = 0:10
    m1 = m2 + index;
    ek(m1) = 4*(g(m2)/4)^((2^index)/n);
    if ek(m1) < 1.0e-14
        break
    end
end
for en = m1:-1:2
    ek(en-1) = 2*sqrt(ek(en))/(1 + ek(en));
end
%compute poles and zeros
for en = 2:m2
    a = (1 + g(en))*e(en - 1)/2;
    e(en) = a + sqrt(a^2 + g(en));
end
u2 = log((1 + sqrt(1 + e(m2)^2))/e(m2))/n;
z = [];
p = [];
for index = 1:n3
    u1 = (2*index - 1)*pi/(2*n);
    c = -i/cos((-u1 + u2*i));
    d = 1/cos(u1);
    for en = m1:-1:2
        c = (c - ek(en)/c)/(1 + ek(en));
        d = (d + ek(en)/d)/(1 + ek(en));
    end
    af(index) = 1/c;
    df(index) = d/ek(1);
    p = [conj(af(index)); af(index); p]; 
    z = [-df(index)*i; df(index)*i; z];
end
if no == 1
    a = 1/sinh(u2);
    for en = m1:-1:2
        a = (a - ek(en)/a)/(1 + ek(en));
    end
    p = [p; -1/a];
end
%gain
k = real(prod(-p)/prod(-z));
if ~rem(n,2) % n is even order so patch the gain
    % k = k/sqrt((1 + epsilon^2)); % This original line is an error.
    k = k*10^(-rp/20);             % See addendum note.
end

% - - - - - - - - - - - - - - - - - - - - - - - -
% Addendum - Note regarding the variable epsilon and k for n even:
%if ~rem(n,2) % n is even order so patch gain
%    % k = k/sqrt((1 + epsilon^2));
%    % ^ This is as shown in the paper. However, epsilon is 
%    %   not otherwise defined and it causes an error.
%    %   I believe it is the standard ϵ² "ripple factor."   
%    %         ϵ² = 10^(passRipple/10)  -  1
%    % So I solved this and I think it is correct, 
%    % as it matches other results:
%    % epssq = 10^(0.1*rp) - 1;
%    % k = k/sqrt(1 + epssq);
%    % k = k/sqrt(10^(0.1*rp));
%    k = k*10^(-rp/20);   
%    % matlab has it as epsilon = sqrt(10^(0.1*rp)-1)
%    % -- which is the same as my guess.
%end
%
% - - - - - - - - - - - - - - - - - - - - - - - - 
