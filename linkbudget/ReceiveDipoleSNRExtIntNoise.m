c = 299792458;              % Speed of Light
k = 1.380649e-23;           % Boltzmann Constant
f_MHz = linspace(3,30,41);
lambda = c./(f_MHz*1e6);

L=5;                      % Dipole Length
radius=0.0219/2;            % Dipole Radius
B=10;                       % 10 Hz bandwidth
RL=1e6;                     % Load resistor in receive mode - Use megaohm for Voc and/or high impedance amplifier
Efield = 100e-9;            % V/m - Received field strength
Tp = 273.15+100;            % Physical Temperature

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

% Small Antenna Analysis from Steve Best
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7909994
dipoleZ=dipoleCylindrical(Length=L,Radius=radius, ClosedEnd=1); %Remake dipole without load
mBeCu=metal(Name='BeCu',Conductivity=0.25*59600000,Thickness=radius);
dipoleZ.Conductor=mBeCu;
G_r_dB=pattern(dipoleZ, f_MHz*1e6, 0, 0, CoordinateSystem="rectangular",Type="gain")';
Zdipole=impedance(dipoleZ,f_MHz*1e6);
ecd=efficiency(dipoleZ,f_MHz*1e6);

S11=(RL-Zdipole)./(RL+Zdipole);  
ML=(1-abs(S11).^2);              

eta_free = 120*pi;               %Impedance of free space
Pd=1/2*Efield.^2/eta_free;       %Best Paragraph 1 of page 2

Si = Pd.*lambda.^2.*10.^(G_r_dB/10)/(4*pi); %Best Eq 34
Ta=ML.*ecd.*T_gal_Cane+ML.*(1-ecd).*Tp;     %Best Eq 33

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
Tpc= T0;
Tr = (Fx-1)*T0/(ecable*Gamp*efilter*ebalun)+(1-ecable)*Tpc/(ecable*Gamp*efilter*ebalun) +...
     (Famp-1)*T0/(efilter*ebalun)+(1-efilter)*Tpc/(efilter*ebalun)+(1-ebalun)*Tp/ebalun;

SNR = ML.*Si ./ (k*B*(Ta+Tr));    %Best Eq 35

figure(2)
hold off
plot(f_MHz,10*log10(SNR),'LineWidth',3)

% SNR of Two Active Antenna Designs
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7972627
% Eq 3
kappa=ecd.*ML.*Gamp;
SNR2=kappa.*Si ./ (k.*Gamp.*Ta*B+k.*B.*Tr); 

hold on
plot(f_MHz,10*log10(SNR2),'LineWidth',3)
grid on
xlim([0 30]);
xlabel('Frequency (MHz)')
ylabel('SNR (dB)')
legend('Steve Best','Konovalenko, et al')