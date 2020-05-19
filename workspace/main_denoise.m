clear all;clc;
%���ѡ��2s�����ݵ� 2f=600��                 
%load('G:\Ҷ����\st\original_fhr_dat.mat');
val=importdata('1.txt');
Signal=original_fhr;
f=4;  %Ƶ��
signal=Signal(2,1:1800);
[row,col]=size(signal);

%------����С��Ӳ��ֵ��ѭ��ƽ���ĵ�����-----  m:denoise
out=zeros(1,col);original_fhr=[];
nspin=8;                                   %�������ã�ƽ�ƴ���
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
