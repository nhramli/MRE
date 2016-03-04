function ynext = rk4all(func,tstart,yval,time,nstep,parvec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth order Runge-kutta function used to take nsteps
% rk4 steps over an interval time length dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = time/nstep;
hh = 0.5*h;
h6 = h/6;

tval = tstart;
for steps = 1:nstep
    k1 = func(tval, yval, parvec);
    k2 = func(tval+hh, yval+hh*k1, parvec);
    k3 = func(tval+hh, yval+hh*k2, parvec);
    k4 = func(tval+h, yval+k3*h, parvec);
    ynext = yval + h6*(k1 + 2*k2 + 2*k3 + k4); 
    
    tval = tstart + steps*h;
    yval = ynext;
end

