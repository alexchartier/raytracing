% Three 12m dipoles in Matlab at 12 GHz along cartiesian axes
% https://apps.dtic.mil/sti/tr/pdf/ADA539999.pdf

f=12e6;
antZ=dipoleCylindrical; antZ.Length=12; antZ.Radius=0.04;
antX=dipoleCylindrical; antX.Length=12; antX.Radius=0.04;
antX.Tilt=90; antX.TiltAxis=[0 1 0];
antY=dipoleCylindrical; antY.Length=12; antY.Radius=0.04;
antY.Tilt=270;antY.TiltAxis=[1 0 0];

% simple Y perfectly circular plane wave
% dirP = [0 1 0]';               %direction vector
% polP = [1*exp(1j*pi/2) 0 1]';  %Polarization vector

% Just do random directions
dirP = 2*rand(3,1)-1;         %Random direction vector
dirP = dirP/norm(dirP);       %Normalize vector

% Find an orthgonal polarization to dirP
% We don't really care about the tilt, so we just project X and Y axis
% and correct the cases in the dirP direction
sign = -1;                                %Change sign for RHCP vs LHCP
ang  = pi/2;                             %pi/2 is circular, 0 is linear 
% Note that alogithm need elliptical pol in general or it may fail

pX   = [1 0 0]';              %Start with polarization vector along X
pY   = [0 1 0]';              %  and along Y
if (dirP == pX)      % Fails if in pX direction
    pX   = [0 0 1]'; pY   = [0 1 0]';              
elseif (dirP == pY)  % Fails if in pY direction
    pX   = [1 0 0]'; pY   = [0 0 1]';              
end
projX = pX-dot(dirP,pX)/dot(pX,pX)*dirP;  %Project X onto plane in dirP
projY = pY-dot(dirP,pY)/dot(pY,pY)*dirP;  %Project Y onto plane in dirP
polP  = projX + exp(sign*1j*ang)*projY;   %circular pol in dirP

%Determine plane wave excitation in direction
pX = planeWaveExcitation(Element=antX,Direction=dirP, Polarization=polP);
pY = planeWaveExcitation(Element=antY,Direction=dirP, Polarization=polP);
pZ = planeWaveExcitation(Element=antZ,Direction=dirP, Polarization=polP);

%Determine voltage on dipole
IX = feedCurrent(pX,f); IY = feedCurrent(pY,f); IZ = feedCurrent(pZ,f);
ZX = impedance(antX,f); ZY = impedance(antY,f); ZZ = impedance(antZ,f); %Assumes no matching network, does not matter for now
VX = ZX*IX; VY = ZY*IY; VZ = ZZ*IZ;
V=[VX;VY;VZ]; 

%Add in simple noise for 1 trial
SNR = 1000;   %Noise in dB        
noiseV_r=2*rand(3,1)-1;         %Random direction noise vector - uniform white noise should use gaussian for multiple trials
noiseV_i=2*rand(3,1)-1;         %Random direction noise vector
noiseV=noiseV_r+1i*noiseV_i;    %Noise phaser
noiseV=10^(-SNR/20)*noiseV/norm(noiseV);
V=V/norm(V)+noiseV;


vDOA=imag(cross(V,conj(V)));
vDOA=vDOA/norm(vDOA);
% Assume perfect knowledge for now to resolve 180 deg ambiguity
if norm(dirP+vDOA) < 1
    vDOA=-vDOA;
end


fprintf('DOA Vector Estimate: [%g, %g, %g]\n',vDOA(1),vDOA(2),vDOA(3))
fprintf('Vector Direction:    [%g, %g, %g]\n',dirP(1),dirP(2),dirP(3))

thetaDir = 180/pi*atan2(dirP(2),dirP(1));
phiDir = 180/pi*acos(dirP(3)/norm(dirP));
thetaDOA = 180/pi*atan2(vDOA(2),vDOA(1));
phiDOA = 180/pi*acos(vDOA(3)/norm(vDOA));

fprintf('DOA Spherical Estimate: theta=%g, phi=%g\n', thetaDOA,phiDOA)
fprintf('Spherical Direction:    theta=%g, phi=%g\n', thetaDir,phiDir)



% Perform analysis in time for just two samples, requires many in practice
% Vm=abs(V);
% Vphi=angle(V);
% 
% a1=Vm.*cos(Vphi);
% a2=Vm.*cos(Vphi+pi/4); %Assume next step in time is a pi/4
% 
% aDOA=cross(a1,a2);
% aDOA/norm(aDOA)