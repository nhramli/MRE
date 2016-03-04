function dyvaldt = NL_rhs( ~ , yval, parvec )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NL_rhs: Evaluates the `nonlinear' terms on the 
% RHS of Eq.(11) consisting of the
% coefficients C_k of the Hermite functions (used by etdrk4all)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: yval   - gridded C_k(z,t) 
%         parvec - contains profiles tau(z), sigma_w(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: dyvaldt - gridded (d/dt) C_k(z,t)          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M = size(yval,2)-1; N = size(yval,1); % M should be odd
dz = 1/N; OM = 3:2:M; EM = 2:2:M;     % odd and even integers excluding 0,1

% first extract the stuff we need
azero = yval(:,1); ak = yval(:,2:M+1); sigma = parvec(:,1);

% M+1 matrix of k*sigma*ak, and skak(:,M+1)=0 :
skak=(sigma*(1:M)).*ak; skak=[skak,zeros(N,1)]; 
sigrep=repmat(sigma,1,M); dakdt=0*ak;

dazerodt = -cFD(skak,1,N,dz);
dakdt(:,1) = - cFD(skak,2,N,dz) -0.5.*sigma.*...
    ([azero(2:N)',azero(N,1)]'-[azero(1),azero(1:N-1,1)']')/dz;
dakdt(:,EM) = -cFD(skak,EM+1,N,dz) - sigrep(:,EM-1).*cFD(ak,EM-1,N,dz);
dakdt(:,OM) = -cFD(skak,OM+1,N,dz) - sigrep(:,OM-1).*cFD(ak,OM-1,N,dz);

dyvaldt = [dazerodt,dakdt];

end