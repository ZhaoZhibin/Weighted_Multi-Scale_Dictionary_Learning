function Y = lps(X,N0)
% Y = lps(X,N0)
% Low-pass scaling
% Notes
% - output Y will be length N0
% - length(X) should be even

% Reference: 'Wavelet Transform with Tunable Q-Factor'
% http://taco.poly.edu/selesi/TQWT/
% Ivan Selesnick
% selesi@poly.edu
% Polytechnic Institute of NYU

% November, 2010


N = length(X);

Y = zeros(1,N0);

% Add 1 to indices because Matlab indexing starts at 1 (not 0)

switch 1
    case N0 <= N
        
        k = 0:N0/2-1;
        Y(k +1) = X(k +1);
        
        Y(N0/2 +1) = X(N/2 + 1);
        
        k = 1:N0/2-1;
        Y(N0-k +1) = X(N-k +1);
        
    case N0 >= N

        k = 0:N/2-1;
        Y(k +1) = X(k +1);
        
        k = N/2:N0/2-1;
        Y(k +1) = 0;
        
        Y(N0/2 +1) = X(N/2 +1);
        
        Y(N0-k +1) = 0;
        
        k = 1:N/2-1;
        Y(N0-k +1) = X(N-k +1); 
        
end
