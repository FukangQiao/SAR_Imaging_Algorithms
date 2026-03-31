
clear all

thetaT=0;%平台波束斜视角
thetaT=thetaT*pi/180;%rad
c=3e8;%光速
fc=1.5e9;%载频
lambda=c/fc;%波长

%%测绘带区域
X0=400;%方位向[-X0,X0]
H=4000;%高度
Y=10000;
Rb=sqrt(H^2+Y^2);
Rc=Rb/cos(thetaT);
R0=500;%测绘带宽2*R0


%%方位向(Azimuth,Cross-Range),x/u domain
v=100;%SAR 平台速度
D=4;
Lsar=lambda*Rc/D;%合成孔径长度
Na=1024;
x=linspace(-X0,X0,Na);%u域序列
u=x/v;%方位向
du=2*X0/v/Na;
fu=linspace(-1/2/du,1/2/du,Na);%fu域序列
fdc=2*v*sin(thetaT)/lambda;%Doppler调频中心频率
fdr=-2*(v*cos(thetaT))^2/lambda/Rc;%Doppler调频斜率

%%距离向(Range),r/t domain
Tr=5e-6;%LFM信号脉宽 5us 
Br=50e6; %LFM信号带宽 50MHz
Kr=Br/Tr; %调频斜率
Nr=512;
Rmin=sqrt(H^2+(Y-R0/2)^2);
Rmax=sqrt(H^2+(Y+R0/2)^2+(Lsar/2)^2);
dt=(2*(Rmax-Rmin)/c+Tr)/Nr; %sampling frequency
t=linspace(2*Rmin/c-Tr/2,2*Rmax/c+Tr/2,Nr);%块时间
r=t*c/2;
f=linspace(-1/2/dt,1/2/dt,Nr);%f域序列


%%目标位置
Ptar=[Y,0;Y-50,0;Y,50];%距离向坐标,方位向坐标
[Ntar,n]=size(Ptar);

%%产生回波
s_ut=zeros(Nr,Na);
U=ones(Nr,1)*u;%扩充为矩阵 %方位向
T=t'*ones(1,Na);%块时间
for k=1:Ntar
    rn=Ptar(k,1);
    xn=Ptar(k,2);
    R=sqrt(H^2+rn^2+(xn-v*U).^2);
    DT=T-2*R/c;
    phase=pi*Kr*DT.^2-4*pi/lambda*R;
    s_ut=s_ut+exp(j*phase).*(abs(DT)<Tr/2).*(abs(v*U-xn)<Lsar/2);
end

%加海明窗
%wr=hamming(Nr);
Nwr=floor(Tr/dt);%快时间
wr=hamming(Nwr);
N=ceil((Nr-Nwr)/2);
wr=[zeros(N,1);wr;zeros(Nr-Nwr-N,1)];
% s_ut=s_ut%.*(wr*ones(1,Na));%距离加窗
 s_ut=s_ut.*(wr*ones(1,Na));%距离加窗
wa=hamming(Na);
% s_ut=s_ut%.*(ones(Nr,1)*wa');%方位加窗
 s_ut=s_ut.*(ones(Nr,1)*wa');%方位加窗


%%距离压缩
p0_t=exp(j*pi*Kr*(t-2*Rc/c).^2).*(abs(t-2*Rc/c)<Tr/2);%距离向LFM信号
p0_f=fftshift(fft(fftshift(p0_t)));
s_uf=fftshift(fft(fftshift(s_ut)));%距离向FFT
src_uf=s_uf.*(conj(p0_f).'*ones(1,Na));%距离压缩
src_ut=fftshift(ifft(fftshift(src_uf)));%距离压缩后的信号?
src_fut=fftshift(fft(fftshift(src_ut).')).';%距离多普勒域

%%距离徙动校正原理
F=f'*ones(1,Na);%扩充为矩阵，距离向
FU=ones(Nr,1)*fu; %方位向
Hrcc=exp(-j*pi/fc/fdr*FU.^2.*F);%距离徙动校正函数
src_fuf=fftshift(fft(fftshift(src_fut)));%距离压缩后的二维频谱?
s2rc_fuf=src_fuf.*Hrcc;
s2rc_fut=fftshift(ifft(fftshift(s2rc_fuf)));%距离多普勒域

%%方位压缩
p0_2fu=exp(j*pi/fdr*(FU-fdc).^2);%方位向压缩因子%%%%%？
s2rcac_fut=s2rc_fut.*p0_2fu;%方位压缩
s2rcac_ut=fftshift(ifft(fftshift(s2rcac_fut).')).';%方位向IFFT%%%%%转置？


%
figure(1)
Sr=abs(src_ut(Nr/4:3*Nr/4,Na/2));
Srmax=max(Sr);
Sr=20*log10(Sr/Srmax);
rs=((Nr/4:3*Nr/4)*dt+2*Rmin/c-Tr/2)*c/2; %%%%%？
ys=sqrt(rs.^2-H^2);
plot(ys,Sr)
axis tight
xlabel('距离/m')
ylabel('幅度归一化/dB')
grid on

figure(2)
G=(abs(s2rcac_ut));
gm=max(max(G));
gn=min(min(G));
G=255/(gm-gn)*(G-gn); %%%%%？
yr=sqrt((r*cos(thetaT)).^2-H^2);
imagesc(x,yr,255-G),%%%%%？
colormap(gray)
axis tight,
xlabel('Azimuth/m')
ylabel('Range/m')
title('目标图象')

figure(3)
subplot(121)
G=(abs(src_fut(Nr/2-30:Nr/2+20,:)));
gm=max(max(G));
gn=min(min(G));
G=255/(gm-gn)*(G-gn);
yr=sqrt((r*cos(thetaT)).^2-H^2);
imagesc(x,yr,255-G),colormap(gray)
grid on,axis tight,
xlabel('Doppler/Hz')
ylabel('Range/m')
title('(a)距离多普勒域')

subplot(122)
G=(abs(s2rc_fut(Nr/2-30:Nr/2+20,:)));
gm=max(max(G));
gn=min(min(G));
G=255/(gm-gn)*(G-gn);
yr=sqrt((r*cos(thetaT)).^2-H^2);
imagesc(x,yr,255-G),colormap(gray)
grid on,axis tight,
xlabel('Doppler/Hz')
ylabel('Range/m')
title('(b)校正后的RD域')

figure(4)
Src_ut=abs(s2rcac_ut).';
yr=((Nr/2-50:Nr/2+50)*dt+2*Rmin/c-Tr/2)*c/2;
ys=sqrt(yr.^2-H^2);
xa=(Na/2-50:Na/2+100)*du*v-X0;
waterfall(ys,xa,Src_ut(Na/2-50:Na/2+100,Nr/2-50:Nr/2+50));
axis tight
xlabel('Range/m')
ylabel('Azimuth/m')
title('目标3D图象')
