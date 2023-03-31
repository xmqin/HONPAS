% test6.m
% Reads and plots the output of test6.f90
% J.M.Soler, May.2014

% Set mesh parameters
clear all
smin = 0;
smax = 10;
kfmin = 0.001;
kfmax = 3;
kgmin = 0;
kgmax = 3;
ds = 0.1;
dkf = 0.05;
dkg = 0.05;
rs = [0.01,2.0,10.0,1.e2];

% Read Fx(s) data
fx = load('fx.out');
ns = size(fx,1);
nf = size(fx,2)-1;
smesh = fx(:,1);
fxMesh = fx(:,2:nf+1);

% Read epsxc(kf,kg) data
epsx = load('epsx.out');
epsc = load('epsc.out');
fxc = load('fxc.out');
kfmesh = unique(epsx(:,1));
kgmesh = unique(epsx(:,2));
zmesh = unique(epsx(:,3));
nf = size(epsx,2)-3;
nkf = numel(kfmesh);
nkg = numel(kgmesh);
nz = numel(zmesh);
iz0 = find(zmesh==0);
epsx = reshape(epsx(:,1+3:nf+3),nkf,nkg,nz,nf);
epsc = reshape(epsc(:,1+3:nf+3),nkf,nkg,nz,nf);
fxc  = reshape( fxc(:,1+3:nf+3),nkf,nkg,nz,nf);

% Read functional names
filein = fopen('funcs.out');
funcs = fscanf(filein,'%c');
fclose(filein);
n = numel(funcs)-1;
funcs = funcs(1:n);
l = n/nf;
funcs = reshape(funcs,l,nf)'

% Interpolate data
s = (smin:ds:smax)';
fx = zeros(numel(s),nf);
for jf = 1:nf
    fxDat = spline1D( fxMesh(:,jf), smesh, 0, Inf );
    fx(:,jf) = splint1D( s, fxDat );
end

% Plot all enhancement factors
figure(1)
plot(s,fx)
xlabel('s')
ylabel('Fx')
legend(funcs,'location','northwest')

% Replot Fx(s) for specific functionals
hold on
for jf = 1:nf
    if strcmp(strtrim(funcs(jf,:)),'PW91') | ...
       strcmp(strtrim(funcs(jf,:)),'B88') | ...
       strcmp(strtrim(funcs(jf,:)),'AM05')
        plot(smesh,fxMesh(:,jf),'o')
    end
end
hold off
axis([0,smax,0,4])
grid on

% Plot xc enhancement factors as a function of kf,kg
if true   % on|off switch
kf = (kfmin:dkf:kfmax);
kg = (kgmin:dkg:kgmax);
nkf = numel(kf);
nkg = numel(kg);
nkfg = nkf*nkg;
[kf,kg] = meshgrid(kf,kg);
kfg = [reshape(kf,1,nkfg);reshape(kg,1,nkfg)];
fxcMesh = fxc;
figure(2)
for jf = 1:nf
    fxcDat = spline2D( fxcMesh(:,:,iz0,jf), kfmesh, kgmesh );
    fxc = splint2D( kfg, fxcDat );
    surf(kf,kg,reshape(fxc,nkg,nkf))
    axis([0,kfmax,0,kgmax,0,3])
    xlabel('kf')
    ylabel('kg')
    zlabel('Fxc')
    title(strtrim(funcs(jf,:)))
    input('Press Enter for next functional\n')
end
end

% Plot xc enhancement factors as a function of s and rs
if false
ns = numel(s);
nrs = numel(rs);
rho = 3/4/pi./rs.^3;
kf = (3*pi^2*rho).^(1/3);
fxcMesh = fxc;
fxc = zeros(ns,nrs);
figure(3)
for jf = 1:nf
    fxcDat = spline2D( fxcMesh(:,:,iz0,jf), kfmesh, kgmesh );
    for irs = 1:numel(rs)
        kfg = kf(irs)*[ones(size(s)),2*s];
        fxc(:,irs) = splint2D( kfg, fxcDat );
    end
    plot(s,fxc)
    legend(num2str(rs'))
    axis([0,smax,1,2])
    grid on
    xlabel('s')
    ylabel('Fxc')
    title(strtrim(funcs(jf,:)))
    input('Press Enter for next functional\n')
end
end


    