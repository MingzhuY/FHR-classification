clear all;clc;
%load original_fhr_;
Signal=original_fhr;
f=4;  %ÆµÂÊ
signal=Signal(3,1:1800);
[row,col]=size(signal);

[fhr4,n2,loss,quality]=GUIfhrPreprocess(original_fhr,10,4);
figure(1)
subplot(2,1,1);
plot(original_fhr(1,:));title('Original FHR signal');ylabel('Bpm');xlabel('Time/s');
subplot(2,1,2);
plot(signal(1,:));title('Preprocessed FHR signal');ylabel('Bpm');xlabel('Time/s');