
The MATLAB(R) routines contained within this directory solves vertical 
dimension Fokker-Planck equation using the Hermite function expansion.
This file describes examples of running the main program 'FP_main.m' 
for the three different test cases as described in the paper by 
H. Mohd Ramli and J. G. Esler (2016).

Examples
(i) 	For constant tau profile for time=0.5(h/u*) with Nz=2^10 
        number of z-grid points, using the fourth order Runge-Kutta(rk4)
        [Ck,zgrid] = FP_main(2^10,20,0.5,5e-2,5e-2,0,1,5e3,0.5,1);
	
(ii) 	For modified stable profile for time=1(h/u*), using the 
        exponential time-differencing fourth order Runge-Kutta(etdrk4)
        [Ck,zgrid] = FP_main(2^10,20,0.5,5e-2,5e-2,0,2,5e3,1,2);
	
(iii) 	For modified nuetral profile for time=3(h/u*), with epsilon 
        parameter of 0.8
        [Ck,zgrid] = FP_main(2^10,20,0.5,5e-2,5e-2,0.8,3,5e3,3,2);
	
To plot the particle concentration, use
	plot(Ck(:,1),zgrid); 
	ylabel('Height (z/h)'); 
	xlabel('Particle concentration c(z,t)') 
