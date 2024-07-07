% SAR_Figure_3_10
% 2016.10.12

close all;clear all;clc

T = 7.24e-6;                % 信号持续时间
B = 14.2e6;                 % 信号带宽
K = B/T;                    % 调频率
ratio = 1.25;               % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = ((0:N-1)-N/2)/N*T;      % 时间轴
f = ((0:N-1)-N/2)/N*Fs;     % 频率轴

st = exp(1i*pi*K*t.^2);     % 生成信号
Sf = fft(st);               % FFT
Hf = exp(1i*pi*f.^2/K);     % 频域匹配滤波器

Out = Sf.*Hf;               % 频域匹配滤波

window = kaiser(N,2.5)';    % 时域窗
Window = fftshift(window);  % 频域窗

st_window = st.*window;     % 时域加窗
Out_Window = Out.*Window;    % 频域匹配滤波后加窗

figure,set(gcf,'Color','w');
subplot(2,2,1),plot(t*1e6,window);axis tight
title('时域窗函数'),ylabel('幅度');
subplot(2,2,2),plot(Window);axis tight
title('频域窗函数');
subplot(2,2,3),plot(t*1e6,real(st_window));axis tight
title('信号实部'),xlabel('时间（\mus）'),ylabel('幅度');
subplot(2,2,4),plot(real(Out_Window));axis tight
title('频谱实部'),xlabel('频率（单元）');