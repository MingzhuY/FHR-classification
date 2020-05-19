function [acc_time,acc_index]=GUIfhrAccDetection(fhr,baseline)
% GUIfhrAccDetection: Detect acceleration of FHR signl (¸ü¸Ä×ÔfhrAccDetection.m)
% Input:
%       fhr:    fetal heart rate signal
%       baseline
% Output:
%       acc: acceleration
    n=length(fhr);
    if length(baseline)==1
        baseline=ones(1,n).*baseline;
    end
    big_time=0;
    acc_index=[];acc_time=0;
    % Method 1:
    i=1;
    while i<=n
        if (fhr(i)-baseline(i))>15
            big_time=big_time+1;
            for j=i+1:n
                if(fhr(j)-baseline(j))>15
                    big_time=big_time+1;
                else
                    break;
                end
            end
        end
        if big_time>15*4
            acc_time=acc_time+1;
            acc_index=[acc_index; i,(j-1)];
            i=j+1;
        else
            i=i+1;
        end
        big_time=0;
    end
end