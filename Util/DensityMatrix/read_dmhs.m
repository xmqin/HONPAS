% Sparse matrix reading in Octave
%
% You need the octcdf package available from the Octave sourceforge
%
function [ dmin, dmout, h, s, norbs, m, nnzs ] = read_dmhs(file,step)

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
dmin1=nc{"dm_in"}(:);
dminval=squeeze(dmin1(step,1,:));             % Select first step and spin and flatten those dimensions
sval=nc{"overlap"}(:);
dmout1=nc{"dm_out"}(:);
dmoutval=squeeze(dmout1(step,1,:));             % Select first step and spin and flatten those dimensions
h1=nc{"h"}(:);
hval=squeeze(h1(step,1,:));             % Select first step and spin and flatten those dimensions

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

dmin=sparse(ri,column,dminval,norbs,m);
dmout=sparse(ri,column,dmoutval,norbs,m);
s=sparse(ri,column,sval,norbs,m);
h=sparse(ri,column,hval,norbs,m);
%
end



