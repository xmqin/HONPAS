% Planar Averages of Grid quantities in Octave
%
% You need the octcdf package available from the Octave sourceforge
%
function [ z, fz ] = average_z(file)

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
% This is for a Z-oriented slab
%
fzg(1:n3)=rhog(1,1,1:n3);
fz=n3*ifft(fzg);                  % Proper normalization
%
L3 = cell(3,3)
z=linspace(0.,L3, n3);
%%plot(z,fz)
%
% Check: The integral of fz over the range, times the surface, has to give Q
%        There should be a better way to do the integral
%
Average = rhog(1,1,1)
Average_times_Volume = rhog(1,1,1)*celvol     % This is Q for the charge case
dz = z(2) - z(1);
surf = celvol/L3
sum(fz)*dz*surf

%
end



