f=linspace(2e6,75e6,19);
ant=dipoleCylindrical;
ant.Length=12;
ant.Radius=0.04;

%D=pattern(ant, f, 0, 0, CoordinateSystem="rectangular",Type="directivity")
G_dB=pattern(ant, f, 0, 0, CoordinateSystem="rectangular",Type="gain")';
G_real_dB=pattern(ant, f, 0, 0, CoordinateSystem="rectangular",Type="realizedgain")';

%I cannot figure out how to get realized gain with a different feed
%impedance automatically in matlab, so I use a mismatch loss instead
RL_dB=returnLoss(ant,f,50);
RL=10.^(-RL_dB/20);
ML=10*log10(1-RL.^2);

%Validation of method - compare G+ML to G_realized from matlab
figure(2);
plot(f/1e6,G_dB+ML,'b',LineWidth=3)
hold on
plot(f/1e6,G_real_dB,'r:',LineWidth=3)

%We can only really make 1:n^2 baluns, so lets assume those impedances
figure(3)
clf
for n=1:7
    RL_dB=returnLoss(ant,f,50*n^2);
    RL=10.^(-RL_dB/20);
    ML=10*log10(1-RL.^2);    
    plot(f/1e6,G_dB+ML,LineWidth=3)
    hold on
    legendstr{n}=sprintf('1:%00d',n^2);
end

grid on
xlabel('Frequency (MHz)')
ylabel('Realized Gain (dBi)')
legend(legendstr)
