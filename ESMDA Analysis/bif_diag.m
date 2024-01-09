%create bifurcation diagram
addpath('C:\Users\Julie\Documents\MATH\MCRN\Vortex\Matlab\Matlab')

T = 24.0*3600; %seconds per day 
L = 6369426.73; %Earth radius in meters

%% Give the boundaries on which this potential should be defined
Xmin = -0.06;
Xmax = 0.04; %-2e7*T/L**2,2e7*T/L^2
Ymin = -0.05;
Ymax = 0.05; %-2e7*T/L**2,2e7*T/L^2
Umin = 0;
Umax = 0.8; %15*T/L,40*T/L

%% Get Variables    
params_d = default_params();
params_nd = nondimensionalize(params_d);
coeffs = params2coeffs(params_nd);

tau1 = coeffs(1);
R = coeffs(2); %Dont care about this if Lam(t)
s = coeffs(3);
xi = coeffs(4);
dw = coeffs(5);
zeta = coeffs(6);
dw = coeffs(7); %Not sure why dw is returned twice but its okay they're equal
URR = coeffs(8); %Dont care about this if Lam(t)
% URR = 0.47; %% Ligia's value
tau2 = coeffs(9);
eta = coeffs(10);
dlam = coeffs(11);
Lamdott = coeffs(12); %Dont care about this if Lam(t)
h = coeffs(13);
hdot = coeffs(14);
Psi0 = coeffs(15);
Psi0dot = coeffs(16);
URR = coeffs(17); %Dont care about this if Lam(t)
f0_d = params_d(2);
N_d = params_d(4);
f0 = f0_d*T;
N = N_d*T;
H = params_nd(1);
zT = params_nd(5);
k = params_nd(6);
ell = params_nd(7);
beta = params_nd(8);
UB = params_nd(13);
fN = f0^2/N^2;
oh = 1.0/(2*H);
eps = 8.0/(3*pi);
a1 = k^2 + ell^2 + fN*(oh^2 + 8.0/zT^2);
a3 = -k*beta;
a4 = eps*k*fN*(2.0/zT-oh);
a5 = eps*k*fN*(4/zT^2+1.0/(H*zT));

h =  40/L;
hs = linspace(0,250,10000)/L;
lambs = -1*10^(-3);
%lambs = linspace(-3*10^(-3),3*10^(-3),10);

epsil = 0;
equilibs = zeros(length(hs)*2,1);


for i = 1:length(lambs)
    for j = 1:length(hs)
    h = hs(j);
%% Solve ODE
    %%% for Lam, h constant use f in ode45 command

    URR = UB + lambs(i)*T*zT/2;
    R = -(a3+a4*lambs(i)*T+a5*UB)/a1;

%    x' = -x/122.6276 - r(L)*y+ 1.9638*y*u - 1.7488*h + 70.8437*0
%    y' = -y/122.6276 + r(L)*x- 1.9638*x*u + 240.5361*h*u
%    u' = -(u-p(L))/30.3713 - 91310*h*y
%f = @(t,a)[-a(1)/122.6276 - R*a(2)+1.9638*a(2)*a(3)-1.7488*h;
%    -a(2)/122.6276+R*a(1)-1.9638*a(1)*a(3)+240.5361*h*a(3);
%    -(a(3)-URR)/30.3713-91310*h*a(2)];


    f = @(t,a) [-a(1)/tau1 - R*a(2)+ s*a(2)*a(3) - xi*h + dw*hdot; 
    -a(2)/tau1 + R*a(1)- s*a(1)*a(3) + zeta*h*a(3);
    -(a(3)-URR)/tau2 - eta*h*a(2)-dlam*Lamdott;];

    xt0 = [0,0,-1*T/L]; %starting values around equilibrium point?
    [tm,a] = ode45(f,[0 365.25*35],xt0);

    %save endpoint
    equilibs(2*j-1) = a(size(a,1),3);

    xt0 = [0,0,-20*T/L]; %starting values around equilibrium point?
    [~,a] = ode45(f,[0 365.25*35],xt0);
    equilibs(2*j) = a(size(a,1),3);
    end
end

%% Plot bifurcation diagram
figure
scatter(repelem(hs,2)*L,equilibs*L/T)
%writematrix([repelem(hs,2)',equilibs*L/T], 'heqpts_n1.csv')

pts = readtable("l0p6equil.dat");
scatter(pts.Var1*L, pts.Var2*L/T,[], {'1','2'},'filled',pts.Var4)
xlim([0,250])