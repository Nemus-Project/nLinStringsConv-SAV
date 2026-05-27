%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             Cubic string SAV split form
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear
close all

% Flags
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
SAV = false;
    plotPsi = true;         

initType = 1;               %1=first mode, 2=raised cosine

dampOn = false;             % enable damping
nLinOn = true;              % enable nonlinearity

play = false;               % set if to play at the end
realTimeDraw = false;       % plot in real time
computeEnergy = true;       % computes and plots the energy
    gridOn = true;
plotSpect = true;

small_shift = false;

% Custom Parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
OSfac = 1024;
SR = OSfac*44100;           % sample rate
k = 1 / SR;                 % time step
durSec = 0.06;              % length of simulation in seconds
timeSamples = SR*durSec;     

amplitude = 2e-3;           % initial displacement for initial conditions

% potential shift
if small_shift
shiftV = eps;  
else
    shiftV = 1000;
end

%-- string parameters
rho = 8.05*10^3;            % material density [kg/m^3] 
T0 = 75;                    % tension [N] 
radius = 3.5560e-04;        % radius (0.016 gauge) [m] 
E = 174e9;                  % Young modulus [Pa]
L = 1;                      % length in meters

Area = pi*radius^2;         % area of string section
rA = rho*Area;
I = (pi*radius^4)/ 4;       % moment of Inertia
K = sqrt(E*I/rA);           % Stiffness parameter
c = sqrt(T0/rA); 

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Damping
omega1 = 0*2*pi;
omega2 = 1000*2*pi;
T601 = 15; T602 = 10;
param1 = (-c^2 + sqrt(c^4 + 4*K^2*omega1^2))/(2*K^2);
param2 = (-c^2 + sqrt(c^4 + 4*K^2*omega2^2))/(2*K^2);

sigma0 = (6*log(10)/(param2 - param1))*(param2/T601 - param1/T602); 
sigma1 = (6*log(10)/(param2 - param1))*(-1/T601 + 1/T602);

if ~dampOn sigma0 = 0; sigma1 = 0; end

sigmaL = sigma0;

% grid spacing
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
epsilon = 0.01;          %Deviation from stability condition
h = sqrt((c^2*k^2 + 4*sigma1*k + sqrt((c^2*k^2 + 4*sigma1*k)^2+((16*K^2*k^2))))/2)*(1 + epsilon);

N = 2*floor(0.5*L/h);
h = L/N;
k = (-(h^2*(2*sigma1 + 2*epsilon*sigma1 - (4*K^2*epsilon^2 + 8*K^2*epsilon + 4*K^2 + c^2*h^2 + 4*epsilon^2*sigma1^2 + 8*epsilon*sigma1^2 + 4*sigma1^2)^(1/2)))/(4*K^2*epsilon^3 + 12*K^2*epsilon^2 + 12*K^2*epsilon + 4*K^2 + c^2*epsilon*h^2 + c^2*h^2));
SR = floor(1/k);
timeSamples = floor(SR*durSec);

nPoints = N - 1;

% output position
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
outPos = floor((nPoints + 1)/2);

%initialisations 1
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
Id = speye(nPoints);
Dxm = (diag(ones(1,N)) + diag(-1*ones(1,N-1),-1));
Dxm = sparse(Dxm(:,1:end-1)/h);
Dxp = -1*Dxm.';
Dxx = Dxp*Dxm;
Dxxxx = Dxx.'*Dxx;


fac = k^2*(E*Area - T0)/4;
B = 2*rA*Id + (T0*k^2 + 2*rA*sigma1*k)*Dxx - E*I*k^2*Dxxxx;
M1 = (rA*sigma0*k - rA)*Id - 2*rA*sigma1*k*Dxx;
invM = (1/(rA+rA*sigma0*k))*ones(nPoints,1); %Computing inverse of M for Sherman-Morrison
M = sparse((rA+rA*sigma0*k)*Id);

Delta2 = sparse(zeros(N));

%initial conditions
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
u = zeros(nPoints,1);
uPrev = zeros(nPoints,1);
uNext = zeros(nPoints,1);

u1 = zeros(nPoints,1);
uPrev1 = zeros(nPoints,1);
uNext1 = zeros(nPoints,1);

switch initType
    case 1
        uPrev = (amplitude*sin(pi*(1:nPoints)/N)).';
    case 2
        widthMeter = 0.25;  
        width = floor(widthMeter*nPoints);
        excitPos = outPos;
        initDistr = zeros(nPoints,1);
        for i = 1:nPoints
            if abs(i-excitPos)<=width
                initDistr(i) = 0.5*(1+cos(pi*(i-excitPos)/width));
            else
                initDistr(i) = 0;
            end
        end
        uPrev = initDistr*amplitude;
end
r = Dxm*uPrev;
u = uPrev + 0.5*k^2*(h*c^2*Dxx*uPrev - h*K^2*Dxxxx*uPrev + nLinOn*h*(0.5*(E*Area - T0))*(Dxp*(r.*r.*r)));
u1 = u;
uPrev1 = uPrev;

% initialisations 2
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
uMid = 0.5*(u + uPrev);

r = (Dxm*uMid).^2;
norm = h*(r.'*r);
V = norm*(E*Area - T0)/8;

psiPrev = sqrt(2*(V + shiftV));
psiNext = 0;

Out = zeros(timeSamples+2,1); 

Out(1) = uPrev(outPos);
Out(2) = u(outPos);


Kin = zeros(1,timeSamples);
PotLin = zeros(1,timeSamples);
PotNlin = zeros(1,timeSamples);
En = zeros(1,timeSamples);
Psi = zeros(1,timeSamples+1);
Loss = 0;

Psi(1) = psiPrev;

g = zeros(nPoints,1);

% simulation
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
tic
for n=1:timeSamples
    lin = +B*u + (M1*uPrev);
    if SAV
    %%%% SAV 
        r = (Dxm*u);
        norm = h*(r.*r).'*(r.*r);
        V = norm*(E*Area - T0)/8;
        coeff = 1/sqrt(2*(V + shiftV));
        g = - coeff*h*Dxp*(r.*r.*r)*(E*Area - T0)/2;
        %g = - coeff*Dxp*(r.*r.*r)*(E*Area - T0)/2;
            
        b1 = g.'*uPrev;
 
        b = lin + (0.25*k^2/h)*g*b1 - (k^2/h)*g*psiPrev;
        chi = invM .* ((0.25*k^2/h)*g);

        % b = lin + (0.25*k^2*h)*g*b1 - (k^2)*g*psiPrev;
        % chi = invM .* ((0.25*k^2*h)*g);

        csi = invM .* b; 
        c1 = g.'*csi;
        c2 = (1 + g.'*chi);
        uNext = csi - chi*(c1/c2);
        
        psiNext = psiPrev + 0.5*g.'*(uNext - uPrev);
        %psiNext = psiPrev + 0.5*h*g.'*(uNext - uPrev);

        Psi(n+1) = psiNext;

        if computeEnergy
            Kin(n) = 0.5*rA*h*(uNext-u).'*(uNext-u)/k^2;
            PotLin(n) = - 0.5*T0*h*uNext.'*(Dxx*u) + 0.5*E*I*h*(Dxx*uNext).'*(Dxx*u);
            PotNlin(n) = 0.5*psiNext^2;

            deltaD = 0.5*(uNext - uPrev)/k; deltaM = (u - uPrev)/k;
            Loss = (Loss + 2*rA*h*(sigma1*(deltaD.'*(Dxx*deltaM)) - sigma0*(deltaD.'*deltaD)));
            En(n) = Kin(n) + PotLin(n) + PotNlin(n) - Loss*k;
        end

        psiPrev = psiNext;
    else
    %%%%% Implicit scheme
        if nLinOn
            Delta2(1:(N+1):N^2) = (Dxm*u).^2;
        end
        DxDelta = fac*Dxp*Delta2*Dxm;
        Aexpl = M - DxDelta;
        Bexpl = lin + DxDelta*uPrev;
        uNext = Aexpl\Bexpl;
        
        if computeEnergy
            Kin(n) = 0.5*rA*h*(uNext-u).'*(uNext-u)/k^2;
            PotLin(n) = - 0.5*T0*h*uNext.'*(Dxx*u) + 0.5*E*I*h*uNext.'*(Dxxxx*u);
            if nLinOn
                PotNlin(n) = (h*(E*Area - T0)/8)*((Dxm*uNext).^2).'*((Dxm*u).^2);
            end
            deltaD = 0.5*(uNext - uPrev)/k; deltaM = (u - uPrev)/k;
            Loss = (Loss + 2*rA*h*(sigma1*(deltaD.'*(Dxx*deltaM)) - sigma0*(deltaD.'*deltaD)));
            En(n) = Kin(n) + PotLin(n) + PotNlin(n) - Loss*k;
        end
            
        uMid = 0.5*(u + uNext);
        r = (Dxm*uMid).^2;
        norm = h*(r.'*r);
        V = norm*(E*Area - T0)/8;
        psiNext = sqrt(2*(V + shiftV));
        Psi(n+1) = psiNext;
    end
    
    Out(n+2)= uNext(outPos);
       

    if realTimeDraw
        plot(uNext)
        ylim([-amplitude,amplitude]);
        drawnow;
    end
    
    uPrev = u;
    u = uNext;
end
realTimeFrac = toc/(timeSamples/SR)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Play & Plot
if play 
    soundsc(diff(Out),SR); 
end

figure(1)
plot((1:timeSamples+2)*k, Out)
title('Displacement')

if plotSpect
windowLength = 2^nextpow2(SR/8);
if windowLength < timeSamples
    figure(2)
    spectrogram(Out,hann(windowLength),floor(windowLength/1.1),windowLength,SR,'onesided','yaxis');
    ax=gca;
    ylim(ax, [0,1.5]);
    colormap hot
    colorbar off
end
end

if computeEnergy
    figure(3)
    pow2 = 2^(nextpow2(En(1))-1);
    plot(diff(En)/pow2, 'k.')
    if gridOn
        hold on
        yline(0);
        hold on
        for count = 1:10
            lineHeigh = eps*count;
            yline(lineHeigh);
            hold on
            yline(-lineHeigh);
            hold on
        end
        hold off
    end

    title("Total Energy")
end

%if SAV
figure
plot(Psi)
title("\Psi")
%end