function C = getConstant(n,k)

% function C = getConstant(n), function C = getConstant(n,k)
% 
% Gets constant for OMP denoising given the dimension of the signal n.
% If w is the gaussian random vector, then this function returns constant C
% such that P( ||w|| <= sqrt(n)*C*sigma ) = k (=0.93 default), given by the
% Rayleigh distribution.
% 
% More details in Mairal, J. Elad, M. Sapiro, G, Sparse Representation for 
% Color Image Restoration, and K. S. Miller, Multidimensional Gaussian 
% Distributions. New York: Wiley, 1964.
% 
% Jeremias Sulam
% CS - Technion
% 9/14

if nargin < 2
    k = 0.93;
end

a = n/2;
x = gammaincinv(k,a);
C = sqrt(x*2/n);

end


