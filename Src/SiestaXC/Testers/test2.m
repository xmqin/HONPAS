% test2.m
% plots the xc potential calculated by
% Siesta/Src/SiestaXC/Testsers/test2.f90
% J.M.Soler. Jul.2012

clear all
data_drsll = load('test2vxc.drsll');
%data_vv    = load('test2vxc.vv');
data_vv    = load('test2vxc.vv.sqr_rho');
r   = data_drsll(:,2);
rho = data_drsll(:,3);
vxc_drsll = data_drsll(:,4);
vxc_vv    = data_vv(:,4);

figure(1)
plot(r,vxc_drsll, r,vxc_vv)
legend('DRSLL','VV','Location','SouthEast')
xlabel('r'), ylabel('V_{xc}')
