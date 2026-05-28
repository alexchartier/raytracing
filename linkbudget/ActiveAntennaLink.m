c = 299792458;              % Speed of Light
k = 1.380649e-23;           % Boltzmann Constant
f_MHz = linspace(3,30,41);
lambda = c./(f_MHz*1e6);

P_t=10;                  %Tx Power in Watt
LTx=6;                   %Length of Tx Antenna

ecable=10^(-0.25/10);    %Efficiency of coax
efilter=10^(-0.25/10);   %Efficiency of bandpass filter
ebalun=10^(-0.25/10);    %Efficiency of balun
Gamp_dB=30;              %Gain of amp
Gamp=10^(Gamp_dB/10);
Famp_dB = 1.5;           %Amplifier Noise Figure
Famp=10^(Famp_dB/10);  
Fx_dB = 5;               %Receiver Noise Figure
Fx=10^(Fx_dB/10);  
T0 = 290;

%List of ranges
rangeList = [200e3 500e3 1000e3 2000e3 3000e3];

LRx=5;                      % Rx Dipole Length
radiusRx=0.0219/2;            % Dipole Radius
B=10;                       % 10 Hz bandwidth
RL=1e6;                     % Load resistor in receive mode - Use megaohm for Voc and/or high impedance amplifier
Tp = 273.15+100;            % Physical Temperature

polarization_loss = 3;  %Polarization loss in dB
ionospheric_loss = 3;   %Path Loss in dB
pointing_loss = 3;      %Loss due to not being on peak of beam
P_loss=10^((-polarization_loss-ionospheric_loss-pointing_loss)/10);


% Galactic Noise Temperature - Cane (1979)
% Measured in space
% https://academic.oup.com/mnras/article/189/3/465/1001579
% Model fit to Table 2
p=[     0.254660151205025
        -0.419609666259272
         0.476761966430168
         -1.66747017319054
        -0.576048677208052
           7.3310482980969];
T_gal_Cane=10.^(polyval(p,log10(f_MHz)));
% Cane has a little high noise at lower frequencies that the galactic ITU model
% https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.372-11-201309-S!!PDF-E.pdf
% ITU Noise
% Coefficents from Table 1 in https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7909994
c=52; d=23; T0=290;
Fam = c-d*log10(f_MHz);
T_gal_ITU = 10.^(Fam/10)*T0;

% Use MATLAB Antenna Toolbox to estimate realized gain of Tx Antenna
% LTx meter dipole on 1 m ground
% Assume a resistive N:1 balun

ant=monopole(Height=LTx,GroundPlaneLength=1, GroundPlaneWidth=1);
G_dB=pattern(ant, f_MHz*1e6, 0, 0, CoordinateSystem="rectangular",Type="gain")';
Za=impedance(ant,f_MHz*1e6);

% Lossy N^2:1 transformer w/ Rtransformer ohm equivalent resistance
N=2;
Rtransformer=400;
tp=TwoPort(f_MHz*1e6);
tp=tp.addComponent(tp.transformer(1/N));
tp=tp.addComponent(tp.parallel_r(Rtransformer));
[RL_dB, ML_dB, Pnet_dB, PL_dB] = tp.networkLoss(50, tp.cascade(), Za);
G_t_dB = G_dB+PL_dB'; % Gain + power delivered to antenna

G_t=10.^(G_t_dB/10);

% Small Antenna Analysis from Steve Best
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7909994
dipoleZ=dipoleCylindrical(Length=LRx,Radius=radiusRx, ClosedEnd=1); %Remake dipole without load
mBeCu=metal(Name='BeCu',Conductivity=0.25*59600000,Thickness=radius);
dipoleZ.Conductor=mBeCu;
G_r_dB=pattern(dipoleZ, f_MHz*1e6, 0, 0, CoordinateSystem="rectangular",Type="gain")';
Zdipole=impedance(dipoleZ,f_MHz*1e6);
ecd=efficiency(dipoleZ,f_MHz*1e6);

cnt=1;
figure(10);
hold off
for R=rangeList
Efield=P_loss.*sqrt(30*P_t*G_t)/R;
  
S11=(RL-Zdipole)./(RL+Zdipole);  
ML=(1-abs(S11).^2);              

eta_free = 120*pi;               %Impedance of free space
Pd=1/2*Efield.^2/eta_free;       %Best Paragraph 1 of page 2

Si = Pd.*lambda.^2.*10.^(G_r_dB/10)/(4*pi); %Best Eq 34
Ta=ML.*ecd.*T_gal_Cane+ML.*(1-ecd).*Tp;     %Best Eq 33

Tpc= T0;
Tr = (Fx-1)*T0/(ecable*Gamp*efilter*ebalun)+(1-ecable)*Tpc/(ecable*Gamp*efilter*ebalun) +...
     (Famp-1)*T0/(efilter*ebalun)+(1-efilter)*Tpc/(efilter*ebalun)+(1-ebalun)*Tp/ebalun;

% SNR of Two Active Antenna Designs
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7972627
% Eq 3
kappa=ecd.*ML.*Gamp;
SNR2=kappa.*Si ./ (k.*Gamp.*Ta*B+k.*B.*Tr); 

plot(f_MHz,10*log10(SNR2),'LineWidth',3)
hold on
legendstr{cnt}=sprintf('%g km',R/1e3);
cnt=cnt+1;

end

plot([3 30],[12 12],'g-')
xlabel('Frequency (MHz)')
ylabel('SNR (dB)')
title(sprintf('%g m Active Dipole SNR, %g m Monopole, PTx=%d Watt',LRx,LTx,P_t))
grid on
axis([3 30 0 40])
legend(legendstr)
