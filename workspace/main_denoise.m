clear all;clc;
%随机选择2s的数据点 2f=600个                 
%load('G:\叶明珠\st\original_fhr_dat.mat');
val=importdata('1.txt');
Signal=original_fhr;
f=4;  %频率
signal=Signal(2,1:1800);
[row,col]=size(signal);

%------基于小波硬阈值的循环平移心电消噪-----  m:denoise
out=zeros(1,col);original_fhr=[];
nspin=8;                                   %参数设置：平移次数
original_fhr=[];
for i=1:552
     [fhr4,n2,loss,quality]=GUIfhrPreprocess(Signal(i,:),10,4);
     original_fhr(i,:)=denoise(out,fhr4,nspin,col);
end
figure(1)
subplot(2,1,1);
plot(signal(1,:));title('Original FHR signal');ylabel('Bpm');
subplot(2,1,2);
plot(original_fhr(1,:));title('Preprocessed FHR signal');ylabel('Bpm');ylim([0,200])
