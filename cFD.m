
function dakdz = cFD(AKmat,kvec,N,dz)

% Central finite difference including boundary conditions

if mod(kvec,2) == 0         % when k even
  dakdz = 0.5*([AKmat(2:N,kvec)',AKmat(N,kvec)']'...
                -[AKmat(1,kvec)',AKmat(1:N-1,kvec)']')/dz;
else                        % when k odd
  dakdz = 0.5*([AKmat(2:N,kvec)',-AKmat(N,kvec)']'...
                -[-AKmat(1,kvec)',AKmat(1:N-1,kvec)']')/dz;
end

end