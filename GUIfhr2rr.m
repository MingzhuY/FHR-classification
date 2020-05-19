function [rr1,rr2]=GUIfhr2rr(fhr,fs)
% GUIfhr2rr: Convert fhr to RR interval (unit:ms)
% Input:
%       fhr:    fetal heart rate signal	(xlabel:sample,ylabel:bpm)
%       fs:     sample frequency (Hz)
% Output:
%       rr:     1*length(fhr) array     (xlabel:sample)
%       rr2:    2*ength(fhr)/fs array	(xlabel:second -Row1)
    rr1=60./fhr.*1000;   % unit: ms
    n=length(rr1);       % total samples
    rr2=zeros(2,n/fs);
    rr2(1,:)=1:1:n/fs;
    rr2(2,:)=rr1(fs:fs:end);
end