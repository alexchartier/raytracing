c = 299792458;              % Speed of Light
k = 1.380649e-23;           % Boltzmann Constant
f_MHz = linspace(3,30,41);
lambda = c./(f_MHz*1e6);

L=5;                        % Dipole Length
radius=0.0219/2;            % Dipole Radius
B=10;                       % 10 Hz bandwidth
RL=1e6;                     %Load resistor in receive mode - Use megaohm for Voc and/or high impedance amplifier
Efield = 100e-9;            %V/m - Received field strength

%Note that dipole configuration was verifed to have a nearly identical
%response to what FEKO provides in receiving mode

%Dipole Unloaded with RL
dipoleZ=dipoleCylindrical(Length=L,Radius=radius, ClosedEnd=1); 
G_r_dB=pattern(dipoleZ, f_MHz*1e6, 0, 0, CoordinateSystem="rectangular",Type="gain")';
Zdipole=impedance(dipoleZ,f_MHz*1e6);

%Receive dipole loaded with RL
dipoleZ=dipoleCylindrical(Length=L,Radius=radius, ClosedEnd=1);
le=lumpedElement(Impedance=RL);
dipoleZ.Load=le;

%Apply plane wave excitation to dipole
dirP = [0 1 0]';            %Direction vector
polP = [0 0 Efield]';       %Polarization vector
pZ = planeWaveExcitation(Element=dipoleZ,Direction=dirP, Polarization=polP);

%Determine voltage on dipole across load
I = feedCurrent(pZ,f_MHz*1e6);
Pr=1/2*abs(I).^2*RL;
VL = I*RL;
% PrV=1/2*abs(Voc).^2/RL;  %same as Pr of course - sanity check

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


% SNR of High Z system
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4458269&tag=1
r=0.5;   cm=0.9;   Tmin=290;
cphase=-acos( (cm^2+r^2-1) / (2*cm*r));
c=cm*exp(1i*cphase);
d=2*k*Tmin/(sqrt(1-imag(c)^2)+real(c));
Tn=abs(Zdipole)*d/(2*real(Zdipole)*k)*(1+abs(c)*cos(cphase+angle(Zdipole)));
Tn=0; %Internal noise set to zero for comparison

SNR = (abs(VL).^2/2)  ./ (4*k*B*(T_gal_Cane+Tn).*real(Zdipole)); % Eq 7
figure(1)
hold off
plot(f_MHz,10*log10(SNR),'LineWidth',3)

% Alternate SNR measurement using
% Small Antenna Analysis from Steve Best
% https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7909994

% Note that this analysis ignores system noise
% In another MATLAB this analysis is performed and in that case
% shows that the system noise temperature dominates
% This is a contradition to other paper that I need to resolve

%In any case ML cancels out in the SNR equation and is not needed here, so
%RL does not matter in this case
S11=(RL-Zdipole)./(RL+Zdipole);  %Very close to 1 since RL is huge
ML=(1-abs(S11).^2);              %Very close to 0 since S11 is near 1

eta_free = 120*pi;               %Impedance of free space
Pd=1/2*Efield.^2/eta_free;       %Best Paragraph 1 of page 2

Si = Pd.*lambda.^2.*10.^(G_r_dB/10)/(4*pi); %Best Eq 34
%Note that Si=ML.*Pr;

SNR = ML.*Si ./ (k*ML.*T_gal_Cane*B);       %Best Eq 35

hold on
plot(f_MHz,10*log10(SNR),'LineWidth',3)
grid on
xlim([0 30]);
xlabel('Frequency (MHz)')
ylabel('SNR (dB)')
legend('Warnick / Jensen Analysis','Steve Best Analysis')


%Analysis from:
%https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9928059
%Much higher SNR!! - lets just ignore this one for now
SNR=lambda.^2*abs(Efield).^2./(2*eta_free*k*B*T_gal_Cane); %Eq 18
