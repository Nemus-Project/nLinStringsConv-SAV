%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%    Geometrically exact string form B SAV non split form
%                    Riccardo Russo
%                 University of Bologna
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear
close all

% Flags
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
SAV = true;
    plotPsi = true;         

initType = 2;               %1=first mode, 2=raised cosine

dampOn = false;             % enable damping

hTypeLong = true;          % selects longitudinal grid spacing

play = false;               % set if to play at the end
realTimeDraw = false;       % plot in real time
computeEnergy = true;       % computes and plots the energy
    gridOn = true;         % plots grid of multiples of eps

plotSpect = true;

small_shift = true;

% Custom Parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

OSfac = 256;
SR = OSfac*44100;           % sample rate
k = 1 / SR;                 % time step
durSec = 0.06;              % length of simulation in seconds
timeSamples = SR*durSec;     

amplitude = 2e-3;           % initial displacement for initial conditions

% potential shift
if small_shift
    shiftV = eps;  
else
    shiftV = 1e3;  
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
%++++++++++++++++++++++++++++++++++`+++++++++++++++++++++++++++++++++++++++%
epsilon = 0.01;          %Deviation from stability condition
h1 = sqrt((c^2*k^2 + 4*sigma1*k + sqrt((c^2*k^2 + 4*sigma1*k)^2+((16*K^2*k^2))))/2)*(1 + epsilon);
h2  = sqrt(E / rho) * k ;

if hTypeLong
    h = h2;
else
    h = h1;
end

N = 4*floor(0.25*L/h);
h = L / N;

if hTypeLong
    k = 0.999*h/sqrt(E / rho);
    SR = floor(1/k);
    timeSamples = floor(SR*durSec);
else
    k = 0.999*(-(h^2*(2*sigma1 + 2*epsilon*sigma1 - (4*K^2*epsilon^2 + 8*K^2*epsilon + 4*K^2 + c^2*h^2 + 4*epsilon^2*sigma1^2 + 8*epsilon*sigma1^2 + 4*sigma1^2)^(1/2)))/(4*K^2*epsilon^3 + 12*K^2*epsilon^2 + 12*K^2*epsilon + 4*K^2 + c^2*epsilon*h^2 + c^2*h^2));
    SR = floor(1/k);
    timeSamples = floor(SR*durSec);
end
nPoints = N - 1;

% output positions
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
outPos = floor((nPoints + 1)/2);
outPosL = floor((outPos + 1)/2);

% initialisations 1
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
Id = speye(nPoints);
Dxm = (diag(ones(1,N)) + diag(-1*ones(1,N-1),-1));
Dxm = sparse(Dxm(:,1:end-1)/h);
Dxp = -1*Dxm.';
Dxx = Dxp*Dxm;
DxxPr = Dxm*Dxp;
Dxxxx = Dxx.'*Dxx;

IId  = speye(2*nPoints) ;
zeroMat = sparse(nPoints,nPoints);

facNlin = h*0.5*(E*Area - T0);
facSAV1 = (0.25*k^2/h);
facSAV2 = (k^2/h);

D2 = [c^2*Dxx, zeroMat; zeroMat, c^2*Dxx];
D4 = [K^2*Dxxxx, zeroMat; zeroMat, zeroMat];

B = [2*rA*Id + (T0*k^2 + 2*rA*sigma1*k)*Dxx - E*I*k^2*Dxxxx, zeroMat;
    zeroMat, 2*rA*Id + T0*k^2*Dxx];
M1 = [(rA*sigma0*k - rA)*Id - 2*rA*sigma1*k*Dxx, zeroMat; 
     zeroMat, (rA*sigmaL*k - rA)*Id];
invM = [(1/(rA+rA*sigma0*k))*ones(nPoints,1);
       (1/(rA+rA*sigmaL*k))*ones(nPoints,1)];

Dxm2 = [Dxm, sparse(N, nPoints);
       sparse(N, nPoints), Dxm];


% initial conditions
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
u = zeros(nPoints,1);
uPrev = zeros(nPoints,1);
uNext = zeros(nPoints,1);

v = zeros(nPoints,1);
vPrev = zeros(nPoints,1);
vNext = zeros(nPoints,1);

x = [u; v];
xPrev = [uPrev; vPrev];

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
xPrev = [uPrev; vPrev];
q = Dxm*uPrev;
p  = Dxm*vPrev;
pot = sqrt( (1 + p).^2 + (q).^2);
nlin = [Dxp*(2*facNlin*(pot - 1).*q./pot); Dxp*(2*facNlin*(pot - 1).*(1+p)./pot)];
x = xPrev + 0.5*k^2*(h*D2*xPrev - h*D4*xPrev + invM.*(nlin));

xPrev = [uPrev; vPrev];
u = x(1:nPoints); v = x(nPoints+1:end);
    
% initialisations 2
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
g1   = zeros(nPoints,1) ;
g2   = zeros(nPoints,1) ;
g = [g1;g2];

nlin = zeros(nPoints*2, 1);

q = Dxm*u;
p  = Dxm*v;
qPrev = Dxm*uPrev;
pPrev = Dxm*vPrev;

V = 0 ;
for n = 1 : N
    pot = sqrt((1 + 0.5*( p(n)+pPrev(n) ) )^2 + (0.5*(q(n) + qPrev(n)))^2)   ;
    V = V + (pot-1)^2  ;
end
V = facNlin*V;
psiPrev = sqrt(2*V + shiftV);

Kin = zeros(1,timeSamples);
PotLin = zeros(1,timeSamples);
PotNlin = zeros(1,timeSamples);
En = zeros(1,timeSamples);
Psi = zeros(1,timeSamples+1);
Loss = 0;

Out = zeros(timeSamples+2,1);
OutL = zeros(timeSamples+2,1);

Out(1) = uPrev(outPos);
Out(2) = u(outPos);
OutL(1) = vPrev(outPosL);
OutL(2) = v(outPosL);

Psi(1) = psiPrev;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Simulation
tic
for n = 1 : timeSamples
    lin = B*x + (M1*xPrev);

    q = Dxm*x(1:nPoints);
    p = Dxm*x(nPoints+1:end);
    
    if SAV
    %%%% SAV 
        pot = sqrt( (1 + p).^2 + (q).^2);
        g1 = 2*facNlin*(pot - 1).*q./pot;
        g2 = 2*facNlin*(pot - 1).*(1+p)./pot;
        Vtemp = facNlin*sum((pot-1).^2);
        V = sqrt(2*Vtemp + shiftV);

        g1 =  -Dxp*g1/V  ;
        g2 =  -Dxp*g2/V  ;
        g(1:nPoints) = g1 ;
        g(nPoints+1:end) = g2;
         
        b1 = g.'*xPrev;
         
        b = lin + facSAV1*g*b1 - facSAV2*g*psiPrev;
        chi = invM .* (facSAV1*g);
    
        csi = invM .* b; 
        c1 = g.'*csi;
        c2 = (1 + g.'*chi);
        xNext = csi - chi*(c1/c2);
        
        psiNext = psiPrev + 0.5*g.'*(xNext - xPrev);
        
        if plotPsi
            Psi(n+1) = psiNext;
        end
    
        psiPrev = psiNext ;
    
        if computeEnergy
            uNext = xNext(1:nPoints); u = x(1:nPoints); uPrev = xPrev(1:nPoints);
            vNext = xNext(nPoints+1:end); v = x(nPoints+1:end); vPrev = xPrev(nPoints+1:end);
    
            Kin(n) = 0.5*rA*h*(xNext-x).'*(xNext-x)/k^2;
            PotLin(n) = - 0.5*T0*h*uNext.'*(Dxx*u) - 0.5*T0*h*vNext.'*(Dxx*v) + 0.5*E*I*h*(Dxx*uNext).'*(Dxx*u);
            PotNlin(n) = 0.5*psiNext^2;
    
            deltaDu = 0.5*(uNext - uPrev)/k; deltaM = (u - uPrev)/k; deltaDv = 0.5*(vNext - vPrev)/k;
            Loss = (Loss + 2*rA*h*(sigma1*(deltaDu.'*(Dxx*deltaM)) - sigma0*(deltaDu.'*deltaDu) - sigmaL*(deltaDv.'*deltaDv)));
            En(n) = Kin(n) + PotLin(n) + PotNlin(n) - Loss*k;
        end
    else
    %%%% StormerVerlet 
        pot = sqrt( (1 + p).^2 + (q).^2);
        nlin = [Dxp*(2*facNlin*(pot - 1).*q./pot); Dxp*(2*facNlin*(pot - 1).*(1+p)./pot)];
        
        xNext = invM.*(lin + (k^2/h)*nlin);
    end

    xPrev = x ;
    x = xNext ;

    if ~SAV
        qNext = Dxm*xNext(1:nPoints);
        pNext = Dxm*xNext(nPoints+1:end);
        
        V = facNlin*sum((sqrt((1 + 0.5*(p + pNext)).^2 + (0.5*(q + qNext)).^2) - 1).^2);
        psiNext = sqrt(2*V + shiftV);
        Psi(n+1) = psiNext;
    end

    Out(n+2)= xNext(outPos);
    OutL(n+2) = xNext(nPoints + outPosL);

    if isnan(Out(n+2))
        Out(n+2);
        Out(end-2)=NaN;
        OutL(end-2)=NaN;
        break;
    end

    if realTimeDraw
        plot(xNext(1:nPoints))
        ylim([-amplitude,amplitude]);
        drawnow;
    end
end
realTimeFrac = toc/(timeSamples/SR)

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Play & Plot

if play 
    soundsc(resample(diff(Out),1,OSfac),44100);
end

figure
plot((1:timeSamples+2)*k, Out)
title('Displacement transv')
figure
plot((1:timeSamples+2)*k, OutL)
title('Displacement long')

if plotSpect
windowLength = 2^nextpow2(SR/8);
if windowLength < timeSamples
    figure
    spectrogram(Out,hann(windowLength),floor(windowLength/1.1),windowLength,SR,'onesided','yaxis');
    ax=gca;
    ylim(ax, [0,4]); 
    colormap hot
end
end

if computeEnergy && SAV
    figure
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

if plotPsi && SAV
    figure
    plot(Psi)
    title("Psi")
end