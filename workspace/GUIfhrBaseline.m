function baseline=GUIfhrBaseline(fhr)
% fhrBaseline: Calculate baseline of FHR signal  (更改自baseline6.m))
% Input:
%       fhr:    FHR signal after preprocessing
% Output:
%       baseline
    n=length(fhr);
    % 去除信号中的不稳定部分、加速以及减速（>15)
    k=0;
    for i=2:n
        if abs(fhr(i)-fhr(i-1))>15
            k=k+1;
            continue;
        end
        if k~=0
            matrix=[i-k-1,i];
            fhr(i-k:i-1)=interp1(matrix,fhr(matrix),i-k:i-1,'spline');
            k=0;
        end
    end
    m=5*60*4; % width of window function: 10-min
    for j=1:n
        if (j-m)<1 && (j+m)<=n
            baseline(j)=mean(fhr(1:j+m));
        elseif (j-m)>=1 && (j+m)<=n
            baseline(j)=mean(fhr(j-m:j+m));
        elseif (j-m)>=1 && (j+m)>n
            baseline(j)=mean(fhr(j-m:n));
        elseif (j-m)<1 && (j-m)>n
            baseline(j)=mean(fhr(1:n));
        end
    end
end