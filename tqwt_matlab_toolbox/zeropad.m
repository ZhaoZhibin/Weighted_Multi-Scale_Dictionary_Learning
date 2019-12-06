function y = zeropad(x,N);
% y = zeropad(x,N)
% Zero pads vector x to length N. (Appends zeros to end of x.)
%  - x: vector
%  - N: scalar
%  - Erorr if length(x) < N

if ~isvector(x)
    error('Error in zeropad: x must be a vector.')
end

if ~isscalar(N)
    error('Error in zeropad: N must be a scalar.')
end

L = length(x);

if L > N
    error('Warning in zeropad: N must be >= length(x).')
elseif L == N;
    y = x;
else
    y = x;
    y(N) = 0;       % peform zero-padding.
end

    