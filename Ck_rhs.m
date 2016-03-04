function [ dyvaldt ] = Ck_rhs( ~ , yval, parvec )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ak_rhs: Evaluates the RHS of Eq.(11) consisting the
% coefficients C_k of the Hermite functions
% Note: brute force (very small time steps) on exponentials! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: yval   - gridded C_k(z,t) 
%         parvec - contains profiles tau(z), sigma_w(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: dyvaldt   - gridded (d/dt) C_k(z,t)          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = size(yval,2)-1;         % K will be odd
N = size(yval,1); dz = 1/N;
OM = 3:2:K; EM = 2:2:K;     % odd and even integers excluding 0,1
OMI=1./OM; EMI=1./EM;       % odd and even integers inverse

% first extract the stuff we need
azero = yval(:,1); ak = yval(:,2:K+1);
sigma = parvec(:,1); tau = parvec(:,2);

% M+1 matrix of k*sigma*ak, and skak(:,M+1)=0 :
skak=(sigma*(1:K)).*ak; skak=[skak,zeros(N,1)]; 
sigrep=repmat(sigma,1,K); dakdt=0*ak;

dazerodt = -cFD(skak,1,N,dz);
dakdt(:,1) = -ak(:,1)./tau - cFD(skak,2,N,dz) ...
   -0.5.*sigma.*([azero(2:N)',azero(N,1)]'-[azero(1),azero(1:N-1,1)']')/dz;

dakdt(:,EM) = -ak(:,EM)./(tau*EMI) -cFD(skak,EM+1,N,dz)...
    - sigrep(:,EM-1).*cFD(ak,EM-1,N,dz);

dakdt(:,OM) = -ak(:,OM)./(tau*OMI) -cFD(skak,OM+1,N,dz)...
    - sigrep(:,OM-1).*cFD(ak,OM-1,N,dz);

dyvaldt = [dazerodt,dakdt];

end
