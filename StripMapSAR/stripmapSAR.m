clear;clc;close all;
%%================================================================
%%Parameter--constant
C=3e8;                            %propagation speed
%%Parameter--radar characteristics
Fc=1e9;                          %carrier frequency 1GHz
lambda=C/Fc;                 %wavelength 
%%Parameter--target area
Xmin=0;                          %target area in azimuth is within[Xmin,Xmax]
Xmax=50;                  
Yc=10000;                      %center of imaged area
Y0=500;                          %target area in range is within[Yc-Y0,Yc+Y0]
                                       %imaged width 2*Y0
%%Parameter--orbital information
V=100;                            %SAR velosity 100 m/s
H=5000;                          %height 5000 m
R0=sqrt(Yc^2+H^2);
%%Parameter--antenna
D=4;                                %antenna length in azimuth direction
Lsar=lambda*R0/D;         %SAR integration length
Tsar=Lsar/V;                   %SAR integration time
%%Parameter--slow-time domain
Ka=-2*V^2/lambda/R0;    %doppler frequency modulation rate
Ba=abs(Ka*Tsar);           %doppler frequency modulation bandwidth
PRF=Ba;                         %pulse repitition frequency
PRT=1/PRF;                   %pulse repitition time
ds=PRT;                         %sample spacing in slow-time domain
Nslow=ceil((Xmax-Xmin+Lsar)/V/ds); %sample number in slow-time domain
Nslow=2^nextpow2(Nslow);              %for fft
sn=linspace((Xmin-Lsar/2)/V,(Xmax+Lsar/2)/V,Nslow);%discrete time array in slow-time domain
PRT=(Xmax-Xmin+Lsar)/V/Nslow;    %refresh
PRF=1/PRT;
ds=PRT;
%%Parameter--fast-time domain
Tr=5e-6;                         %pulse duration 10us
Br=30e6;                        %chirp frequency modulation bandwidth 30MHz
Kr=Br/Tr;                        %chirp slope
Fsr=3*Br;                        %sampling frequency in fast-time domain
dt=1/Fsr;                         %sample spacing in fast-time domain
Rmin=sqrt((Yc-Y0)^2+H^2);
Rmax=sqrt((Yc+Y0)^2+H^2+(Lsar/2)^2);
Nfast=ceil(2*(Rmax-Rmin)/C/dt+Tr/dt);%sample number in fast-time domain
Nfast=2^nextpow2(Nfast);                   %for fft
tm=linspace(2*Rmin/C,2*Rmax/C+Tr,Nfast); %discrete time array in fast-time domain
dt=(2*Rmax/C+Tr-2*Rmin/C)/Nfast;    %refresh%%%%%%%%%%%%%%%%%
Fsr=1/dt;
%%Parameter--resolution
DY=C/2/Br;                           %range resolution
DX=D/2;                                %cross-range resolution
%%Parameter--point targets
Ntarget=2;                            %number of targets
%format [x, y, reflectivity]
Ptarget=[Xmin,Yc,1               %position of targets
              Xmin,Yc+10*DY,1
              Xmin+20*DX,Yc+50*DY,1];  
disp('Parameters:')
disp('Sampling Rate in fast-time domain');disp(Fsr/Br)
disp('Sampling Number in fast-time domain');disp(Nfast)
disp('Sampling Rate in slow-time domain');disp(PRF/Ba)
disp('Sampling Number in slow-time domain');disp(Nslow)
disp('Range Resolution');disp(DY)
disp('Cross-range Resolution');disp(DX)     
disp('SAR integration length');disp(Lsar)     
disp('Position of targets');disp(Ptarget)
%%================================================================
%%Generate the raw signal data
K=Ntarget;                                %number of targets
N=Nslow;                                   %number of vector in slow-time domain
M=Nfast;                                   %number of vector in fast-time domain
T=Ptarget;                                %position of targets
Srnm=zeros(N,M);
for k=1:1:K
    sigma=T(k,3);
    Dslow=sn*V-T(k,1);
    R=sqrt(Dslow.^2+T(k,2)^2+H^2);
    tau=2*R/C;
    Dfast=ones(N,1)*tm-tau'*ones(1,M);%%%%%%%%%%%%%%%%%
    phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(R'*ones(1,M));
    Srnm=Srnm+sigma*exp(j*phase).*(0<Dfast&Dfast<Tr).*((abs(Dslow)<Lsar/2)'*ones(1,M));%%%%%%%%%%%%%%%%%range
end
%%================================================================
%%Range compression
tr=tm-2*Rmin/C;
Refr=exp(j*pi*Kr*tr.^2).*(0<tr&tr<Tr);
Sr=ifty(fty(Srnm).*(ones(N,1)*conj(fty(Refr))));
Gr=abs(Sr);
%%Azimuth compression
ta=sn-Xmin/V;
Refa=exp(j*pi*Ka*ta.^2).*(abs(ta)<Tsar/2);
Sa=iftx(ftx(Sr).*(conj(ftx(Refa)).'*ones(1,M)));
Ga=abs(Sa);
%%================================================================
%%draw the intensity image of signal
colormap(gray);
figure(1)
subplot(211);
row=tm*C/2-2008;col=sn*V-26;%%%%%%%%%%%%%%%%%
imagesc(row,col,255-Gr);           %intensity image of Sr
 axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('\rightarrow\itRange in meters'),ylabel('\itAzimuth in meters\leftarrow'),
title('Stripmap SAR after range compression'),
subplot(212);
imagesc(row,col,255-Ga);          %intensity image of Sa
 axis([Yc-Y0,Yc+Y0,Xmin-Lsar/2,Xmax+Lsar/2]);
xlabel('\rightarrow\itRange in meters'),ylabel('\itAzimuth in meters\leftarrow'),
title('Stripmap SAR after range and azimuth compression'),
%%================================================================
%%draw 3D picture
figure(2)
waterfall(real(Srnm((200:205),:)));axis tight
xlabel('Range'),ylabel('Azimuth'),
title('Real part of the raw signal'),
figure(3)
waterfall(Gr((200:205),(600:1000)));axis tight %%%%%%%%%%%%%%%%%how to choose200:205
xlabel('Range'),ylabel('Azimuth'),
title('Stripmap SAR after range compression'),
figure(4)
mesh(Ga((200:300),(750:860)));axis tight
xlabel('Range'),ylabel('Azimuth'),
title('Stripmap SAR after range and azimuth compression'),
%%================================================================
%%draw -3dB contour
figure(5)
a=max(max(Ga));
contour(row,col,Ga,[0.707*a,a],'b');grid on
 axis([9995,10050,-20,20]),
xlabel('\rightarrow\itRange in meters'),ylabel('\itAzimuth in meters\leftarrow'),
title('Resolution Demo: -3dB contour');
%%================================================================ 