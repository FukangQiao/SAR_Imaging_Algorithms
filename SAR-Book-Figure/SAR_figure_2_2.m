% SAR_Figure_2_2
% 2016.10.31
clear all;clc;close all;

%% 参数设置
M = 256;                % 矩阵高度
N = 256;                % 矩阵宽度
top = M/8+1;
bottom = M*7/8;
left = N/8+1;
right = N*7/8;
theta = pi/12;          % 扭曲或旋转角度

%% 生成信号
% 原始信号
S0 = zeros(M,N);
S0(top:bottom,left:right) = 1;
% 扭曲信号
S1 = zeros(M,N);
for ii = 1:M
    for jj = 1:N
        x = jj-N/2;
        y = (M+1-ii)-M/2;
        xx = round(x+N/2);
        yy = M+1-round(x*sin(-theta)+y*cos(-theta)+M/2);
        if(yy>=1 && yy<= M)
            S1(ii,jj) = S0(yy,xx);
        end
    end
end
% 旋转信号
S2 = zeros(M,N);
for ii = 1:M
    for jj = 1:N
        x = jj-N/2;
        y = (M+1-ii)-M/2;
        xx = round(x*cos(-theta)-y*sin(-theta)+N/2);
        yy = M+1-round(x*sin(-theta)+y*cos(-theta)+M/2);
        if(xx>=1 && xx<= N && yy>=1 && yy<=M)
            S2(ii,jj) = S0(yy,xx);
        end
    end
end
%% 二维傅里叶变换
% 原始信号的二维傅里叶变换
S0_ff = fftshift(fft2(fftshift(S0)));
S0_ff = abs(S0_ff);
S0_ff = S0_ff/max(max(S0_ff));
S0_ff = 20*log10(S0_ff+1e-4);
% 原始信号二维傅里叶变换
S1_ff = fftshift(fft2(fftshift(S1)));
S1_ff = abs(S1_ff);
S1_ff = S1_ff/max(max(S1_ff));
S1_ff = 20*log10(S1_ff+1e-4);
% 原始信号二维傅里叶变换
S2_ff = fftshift(fft2(fftshift(S2)));
S2_ff = abs(S2_ff);
S2_ff = S2_ff/max(max(S2_ff));
S2_ff = 20*log10(S2_ff+1e-4);

%% 画图
figure,set(gcf,'Color','w');colormap jet
subplot(2,3,1),imagesc(S0);axis image off
title('（a）时域，原始信号');
subplot(2,3,4),imagesc(S0_ff);axis image off
title('（b）原始信号频谱');
subplot(2,3,2),imagesc(S1);axis image off
title('（c）时域，扭曲信号');
subplot(2,3,5),imagesc(S1_ff);axis image off
title('（d）扭曲信号频谱');
subplot(2,3,3),imagesc(S2);axis image off
title('（e）时域，旋转信号');
subplot(2,3,6),imagesc(S2_ff);axis image off
title('（f）旋转信号频谱');