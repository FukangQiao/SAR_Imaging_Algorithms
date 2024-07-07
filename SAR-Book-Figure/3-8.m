% SAR_Figure_3_8
% 2016.08.31
close all;clear all;clc

T = 7.24e-6;                % 信号持续时间
B = 5.8e6;                  % 信号带宽
K = B/T;                    % 调频率
ratio = 10;                 % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = ((0:N-1)-N/2)/N*T;      % 时间轴

st = exp(1i*pi*K*(t-1e-6).^2);     % 生成信号
ht = conj(fliplr(st));             % 匹配滤波器
out = fftshift(ifft(fft(st).*fft(ht)));     % 计算循环卷积

Z = abs(out);
Z = Z/max(Z);
Z = 20*log10(eps+Z);

figure,set(gcf,'Color','w');
subplot(2,2,1),plot(t*1e6,real(st));
title('(a)输入阵列信号的实部');ylabel('幅度');
subplot(2,2,2),plot(t*1e6,Z);axis([-1 1 -30 0]);
title('(c)压缩后的信号(经扩展)'),ylabel('幅度(dB)');
subplot(2,2,3),plot(t*1e6,out);
title('(b)压缩后的信号'),xlabel('相对于t_{0}时间(\mus)'),ylabel('幅度');
subplot(2,2,4),plot(t*1e6,angle(out));axis([-1 1 -5 5]);
title('(d)压缩后信号的相位(经扩展)'),xlabel('相对于t_{0}时间(\mus)'),ylabel('相位(弧度)');