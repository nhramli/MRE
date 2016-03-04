function ynext = etdrk4all(NLfunc,tstart,ynow,time,nstep,parvec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-stepping code using exponential time-differencing 
% fourth-order Runge Kutta (ETDRK4),
% adapted from Kassam and Trefethen 2005.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = time/nstep;
N = size(ynow,1); M = size(ynow,2);
tau = parvec(:,2);

% Linear function & precomputed various ETDRK4 scalar quantities: 
Ms = 0:M-1; MI = repmat(Ms,N,1); taurep = repmat(tau,1,M);
L = -MI./taurep; E = exp(h*L); E2 = exp(h*L/2);

R = 16; r = exp(1i*pi*((1:R)-0.5)/R); % no. of points for complex means
L = L(:); LR = h*L(:,ones(R,1)) + r(ones(N*M,1),:);  % roots of unity
Q = h*real(mean((exp(LR/2)-1)./LR,2));
f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3,2));
f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3,2));
f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3,2));
f1 = reshape(f1,N,M); f2 = reshape(f2,N,M); f3 = reshape(f3,N,M); 
Q = reshape(Q,N,M);

for n = 1:nstep
    t = tstart + n*h;
    Nv = NLfunc(t,ynow,parvec); a = E2.*ynow + Q.*Nv;
    Na = NLfunc(t+h/2,a,parvec); b = E2.*ynow + Q.*Na;
    Nb = NLfunc(t+h/2,b,parvec); c = E2.*a + Q.*(2*Nb-Nv);
    Nc = NLfunc(t+h,c,parvec);
    ynext = E.*ynow + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;  
    ynow = ynext;
end
