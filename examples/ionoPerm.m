%https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8895746
c=299792458;      %Speed of light
u0=4*pi*1e-7;     %Vacuum Permeability
e0=1/(c^2*u0);    %Vacuum Permittivity
me=9.1093837e-31; %Mass of electron
e=1.60217663e-19; %Electron Charge
N0=1.4e12;        %Electron Density m-3
B0=0.5e-4;        %Geomagnetic Field
v=1e3;            %Collision Frequency

wp=sqrt(N0*e^2/(me*e0));%Angular Plasma Frequency
wb=e*B0/me;        %Electron gyrofrequency
%w0=6.6e7;   
%wH=8.6e6;

f=10e6;
w=2*pi*f;

U=1-1i*v./w;
X=wp^2./w.^2;
Y=wb./w;

e1=1-X.*U./(U.^2-Y.^2);
e2=X.*Y./(U.^2-Y.^2);
e3=1-X./U;

for k=1:length(e1)
    eps(:,:,k)=[e1(k) -1i*e2(k)    0 
             1i*e2(k)     e1(k)  0
                    0         0    e3(k)];
end