%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%    Kirchhoff-Carrier string SAV split form
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
nLinOn = true;              % enable nonlinearity

play = false;               % set if to play at the end
realTimeDraw = false;       % plot in real time
computeEnergy = true;       % computes and plots the energy

plotSpect = true;

% Custom Parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
OSfac = 2;
SR = OSfac*44100;           % sample rate
k = 1 / SR;                 % time step
durSec = 0.06;              % length of simulation in seconds
timeSamples = SR*durSec;     

amplitude = 2e-3;           % initial displacement for initial conditions

% potential shift
shiftV = eps;  

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

% initialisations 1
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Offline initialization of matrices
Id = speye(nPoints);
Dxm = (diag(ones(1,N)) + diag(-1*ones(1,N-1),-1));
Dxm = sparse(Dxm(:,1:end-1)/h);
Dxp = -1*Dxm.';
Dxx = Dxp*Dxm;
Dxxxx = Dxx.'*Dxx;

B = 2*rA*Id + (T0*k^2 + 2*rA*sigma1*k)*Dxx - E*I*k^2*Dxxxx;
M1 = (rA*sigma0*k - rA)*Id - 2*rA*sigma1*k*Dxx;
invM = (1/(rA+rA*sigma0*k))*Id;

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
u = uPrev + 0.5*k^2*(h*c^2*Dxx*uPrev - h*K^2*Dxxxx*uPrev + nLinOn*h*(0.5*E*h/(L*rho))*(-uPrev.'*Dxx*uPrev)*Dxx*uPrev);

u1 = u;
uPrev1 = uPrev;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%%Initializations

% initializing potential
uMid = 0.5*(u + uPrev);

V = (E*h^2*Area/8/L)*(-uMid.'*Dxx*uMid)^2;

psiPrev = sqrt(2*V);
psiNext = 0;

Out = zeros(timeSamples + 2,1); 

Out(1) = uPrev(outPos);
Out(2) = u(outPos);

Kin = zeros(1,timeSamples);
PotLin = zeros(1,timeSamples);
PotNlin = zeros(1,timeSamples);
En = zeros(1,timeSamples);
Loss = 0;
Psi = zeros(1,timeSamples+1);

Psi(1) = psiPrev;

g = zeros(nPoints,1);
q = zeros(nPoints,1);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%%Simulation
tic
for n=1:timeSamples
    if SAV
    %%%% SAV 
        lin = +B*u + (M1*uPrev);
        g = - sqrt(E*Area*h^2/L)*Dxx*u;      
    
        b1 = g.'*uPrev;
        b = lin - (k^2/4/h)*g*b1;% - k^2*g*psiPrev/h;
        chi = invM * (k^2/4/h)*g;
    
        csi = invM * b; 
        c1 = g.'*csi;
        c2 = (1 + g.'*chi);
        uNext = csi - chi*(c1/c2);
        
        psiNext = psiPrev + 0.5*g.'*(uNext - uPrev);
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
        %%%%% Explicit scheme
    
        lin = +B*u + (M1*uPrev);
        
        q = (k*sqrt(h*E*Area/L)*0.5)*Dxx*u;
    
        coeff0 = q.'*uPrev;
        coeff1 = lin - q*coeff0;
    
        coeff2 = invM*coeff1;
        coeff3 = invM*q;
        coeff4 = 1+q.'*coeff3;
        coeff5 = q.'*coeff2;
        uNext = coeff2 - (1/coeff4)*coeff3*coeff5;
    
       if ~SAV
            uMid = 0.5*(u + uNext);
            V = (E*h^2*Area/8/L)*(-uMid.'*Dxx*uMid)^2;
            psiNext = sqrt(2*V);
            Psi(n+1) = psiNext;
        end
        
        if computeEnergy
            Kin(n) = 0.5*rA*h*(uNext-u).'*(uNext-u)/k^2;
            PotLin(n) = - 0.5*T0*h*uNext.'*(Dxx*u) + 0.5*E*I*h*uNext.'*(Dxxxx*u);
            deltaD = 0.5*(uNext - uPrev)/k; deltaM = (u - uPrev)/k;
            Loss = (Loss + 2*rA*h*(sigma1*(deltaD.'*(Dxx*deltaM)) - sigma0*(deltaD.'*deltaD)));
            if nLinOn
                PotNlin(n) = (E*Area*h^2/8/L)*(-uNext.'*Dxx*u)^2;
            end
            En(n) = Kin(n) + PotLin(n) + PotNlin(n) - Loss*k;
        end
    end

    Out(n+2)= uNext(outPos);
    
    if realTimeDraw
        plot(uNext)
        drawnow;
    end
    uPrev = u;
    u = uNext; 
end
realTimeFrac = toc

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%%%%% Play & Plot
lineWidth = 1.5;
fontSize = 14;

if play 
    soundsc(Out,SR); 
end

figure(1)
plot((1:timeSamples+2)*k*1000, Out,'LineWidth',lineWidth,'Color','k')
set(gca,'FontSize',fontSize);
xlabel('Time (ms)');
ylabel('u(x_o) (m)')

title('Displacement')

windowLength = 2^nextpow2(SR/8);
if windowLength < timeSamples
    figure(2)
    spectrogram(Out,hann(windowLength),floor(windowLength/1.1),windowLength,SR,'onesided','yaxis');
    colorbar off
    ax=gca;
    set(gca,'FontSize',fontSize);
    ylim(ax, [0,4]); 
    colormap hot
end

if computeEnergy
    figure(3)
    pow2 = 2^(nextpow2(En(1))-1);
    plot((1:timeSamples-1)*k*1000,diff(En)/pow2, 'k.','LineWidth',lineWidth,'Color','k')
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
    set(gca,'FontSize',fontSize);
    xlabel('Time (ms)');
    ylabel('\Delta H')
    xlim([0 20])

    if SAV
        figure(4)
        plot(Psi)
        title("Psi")
    end
end