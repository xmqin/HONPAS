% Sparse matrix reading in Octave
%
% You need the octcdf package available from the Octave sourceforge
%
function [ dm, norbs, m, nnzs ] = read_dm(file)

nc=netcdf(file,"r");

norbs=nc("norbs")(:)
m=nc("no_s")(:)
nnzs=nc("nnzs")(:)
nspin=nc("nspin")(:)
%
% The array is read in C's row-major order, as it is stored in the netCDF file...
%
% Note that this does not take into account the spin
%
row_pointer = nc{"row_pointer"}(:);                     % Get whole array
column=nc{"column"}(:);
dm1=nc{"dm"}(:);
dmval=squeeze(dm1(1,1,:));             % Select first step and spin and flatten those dimensions

% 
%
% Construct row index (we have a sort of packed row index)
%
cdr = [row_pointer; nnzs];   % we are missing the last element
ri = [];  
for i=1:norbs
  ni = cdr(i+1)-cdr(i);
  ri = [ ri; i*ones(ni,1)];
endfor

dm=sparse(ri,column,dmval,norbs,m);
%
end



