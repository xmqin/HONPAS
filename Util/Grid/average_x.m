% Planar Averages of Grid quantities in Octave
%
% You need the octcdf package available from the Octave sourceforge
%
function [ x, fx ] = average_x(file)

nc=netcdf(file,"r");

cell=nc{"cell"}(:)
celvol=det(cell)
%
n1=nc("n1")(:)
n2=nc("n2")(:)
n3=nc("n3")(:)
%
% The array is read in C's row-major order, as it is stored in the netCDF file...
% We have to permute
%
% Note that this does not take into account the spin
%
rhoC = nc{"gridfunc"}(:);                     % Get whole array
rho_spin1=squeeze(rhoC(1,:,:,:));             % Select first spin and flatten spin dimension
rho=permute(rho_spin1,[3,2,1]);               % Permute
%
% Fourier transform
%
rhog=fftn(rho)/numel(rho);
%
% Now, depending on the orientation of the slab, we have to get "x" or "z" averaged
%
% This is for an X-oriented slab
%
fxg(1:n1)=rhog(1:n1,1,1);
fx=n1*ifft(fxg);                  % Proper normalization
%
L1 = cell(1,1)
x=linspace(0.,L1, n1);
%%plot(x,fx)
%
% Check: The integral of fx over the range, times the surface, has to give Q
%        There should be a better way to do the integral
%
Average = rhog(1,1,1)
Average_times_Volume = rhog(1,1,1)*celvol     % This is Q for the charge case
dx = x(2) - x(1);
surf = celvol/L1
sum(fx)*dx*surf

%
end



