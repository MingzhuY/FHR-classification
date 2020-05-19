function [fhr4,n2,loss,quality]=GUIfhrPreprocess(fhr,time,fs)
% GUIfhrPreprocess: preprocessing of FHR signal  (更改自fhrPreprocess1.m)
% Input:
%       fhr:	fetal heart rate signal
%       time:   choose signal length    (unit:min)
%       fs:     sample frequency        (unit:Hz)
% Output:
%       fhr4:   preprocessed fetal heart rate signal
%       n:      signal length           (sample)
%       loss,quality
    time=time*60*fs;
    loss=round(length(find(fhr(1:time)<=50))/time*100*10)/10;
    %% Step 1: 选择正常点作为信号起始点
    for i=1:length(fhr)
        if fhr(i)>=50 && fhr(i)<=200
            fhr1=fhr(i:end);
            break;
        end
    end
    n1=length(fhr1);
    fhr2=fhr1;
    %% Step 2: Missing data
    point=0;
    note=0;
    for i=1:n1
        if fhr2(i)==0
            point=point+1;
            if i==n1
                note=1;
            else
                continue;
            end
        end
        if point>=15*fs
            if note==1
                fhr2(i-point+1:i)=NaN;
            else
                fhr2(i-point:i-1)=NaN;
            end
        elseif point~=0
            if note==1
                fhr2(i-point+1:i)=fhr2(i-point);
            else
                matrix=[i-point-1,i];
                fhr2(i-point:i-1)=interp1(matrix,fhr2(matrix),i-point:i-1,'linear');
            end
        end
        point=0;
    end
    fhr2(isnan(fhr2))=[];   % 去除持续15s以上信号值为0的部分
    fhr3=fhr2(1:time);
    n2=length(fhr3);
    %% Step 3: Unstable
    point1=0;
    for i=2:n2
        if abs(fhr3(i)-fhr3(i-1))>25
            for j=i+1:n2
                if j>=n2-1  % five adjacent samples which don't differ more than 10bpm
                    if abs(fhr3(j-1)-fhr3(j))<10  &&  abs(fhr3(j-2)-fhr3(j-1))<10 && ...
                            abs(fhr3(j-3)-fhr3(j-2))<10 && abs(fhr3(j-4)-fhr3(j-3))<10 && ...
                            abs(fhr3(j-5)-fhr3(j-4))<10
                        break;
                    end
                else    % 五点=当前点+左右两个点
                    if abs(fhr3(j-1)-fhr3(j-2))<10  &&  abs(fhr3(j)-fhr3(j-1))<10 && ...
                            abs(fhr3(j+1)-fhr3(j))<10 && abs(fhr3(j+1)-fhr3(j))<10 && ...
                            abs(fhr3(j+2)-fhr3(j+1))<10
                        break;
                    end
                end
            end
            matrix=[i-1,j];
            fhr3(i:j-1)=interp1(matrix,fhr3(matrix),i:j-1,'linear','extrap');
            point1=point1+(j-1-i+1);
        end
    end
    %% Step 4: Artifact
    fhr4=fhr3;
    point=0;point2=0;
    for i=1:n2
        if fhr4(i)<50 || fhr4(i)>200
            point=point+1;
            continue;
        end
        if point~=0
            matrix=[i-point-1,i];
            fhr4(i-point:i-1)=interp1(matrix,fhr4(matrix),i-point:i-1,'spline');
            point2=point2+point;
            point=0;
        end
    end
    quality=round((100-(point1+point2)/n2*100)*10)/10;
end