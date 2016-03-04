function [Cknew,z] = FP_main(Nz,M,z0,sigz,zB,epsil,pflag,nstep,tout,tflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB routine to solve the Fokker-Planck equation
% using a Hermite function expansion 
% by HMR and JGE (started Oct 14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
%   Nz   - no. of z-grid points
%   M    - M no. of Hermite functions (M is even, M=K+1)
%   z0   - mean level of initial conc. profile 
%   sigz - standard dev. of initial conc. profile
%   zB   - regularisation param. for profiles                  
%  epsil - param. for neutral profile
%  pflag - flag for profile 
%          (1 - cons. tau, 2 - stable, 3 - neutral)
%  nstep - no. time-steps (dt=tout/nstep)
%   tout - final time for profile
%  tflag - flag for time-stepping (1- rk4, 2 - etdrk4)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS: 
%  Cknew - profiles C_k(z,tout)  k=0,...,M-1=K  
%     z  - z-grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set up grid, initialise profiles

z = (-0.5+(1:Nz)')/Nz;  % staggered half-grids in 0<z<1
Ck = zeros(Nz,M-1);
Czero = (1/(sqrt(2*pi)*sigz)).*exp(-0.5*(z-z0).^2/sigz^2);

% set up profile

switch pflag
    case 1  % constant tau
        tau = 0.1.*ones(Nz,1); sigmaw = 0.5.*(1 + z);   
    case 2  % modified stable profile
        Zm = zB + z.*(1-2*zB);
        sigmaw = 1.3*(1-Zm); tau = 0.1*Zm.^(4/5)./sigmaw; 
    case 3  % modified neutral profile
        Zm = zB + z.*(1-2*zB);
        sigmaw = 1.3*exp(-2.*Zm./epsil); 
        tau = 0.5*Zm./(sigmaw.*(1+(15.*Zm./epsil)));
end

% time-stepping

parvec = [sigmaw,tau]; 
switch tflag
    case 1 
    Cknew = rk4all(@Ck_rhs,0,[Czero,Ck],tout,nstep,parvec);
    case 2
    Cknew = etdrk4all(@NL_rhs,0,[Czero,Ck],tout,nstep,parvec); 
end
      
end


