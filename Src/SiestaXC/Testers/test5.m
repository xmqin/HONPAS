% test5.m
% reads and plots the output of Siesta/Src/SiestaXC/Testers/test5.f90, 
% of the fit of the Vydrov-vanVoorhis vdW functional kernel.
% J.M.Soler, Jul.2012

% Read data
data = load('n1n2phi.table');
r  = data(:,1);
kf = data(:,2);
kg = data(:,3);
n1n2phi_exact = data(:,4);
n1n2phi_interp = data(:,5);
r = unique(r);
kf = unique(kf);
kg = unique(kg);
nr = numel(r);
nkf = numel(kf);
nkg = numel(kg);
n1n2phi_exact = reshape(n1n2phi_exact,nkg,nkf,nr);
n1n2phi_interp = reshape(n1n2phi_interp,nkg,nkf,nr);

% Plot data
[kf,kg] = meshgrid(kf,kg);
for ir = 1:nr
    figure(ir)
    surf(kf,kg,n1n2phi_exact(:,:,ir)*27.2e6), hold on
    surf(kf,kg,n1n2phi_interp(:,:,ir)*27.2e6)
%    surf( kf, kg, (n1n2phi_interp(:,:,ir)-n1n2phi_exact(:,:,ir))*27.2e6 )
%    surf( kf, kg, n1n2phi_interp(:,:,ir)./n1n2phi_exact(:,:,ir) )
    xlabel('k_F'), ylabel('\nabla(n) / n')
    zlabel('n_1 n_2 \Phi (k_1,k_2,r_{12} (\mueV / bohr^6)')
    title(['r_{12} = ',num2str(r(ir))])
    hold off
end
