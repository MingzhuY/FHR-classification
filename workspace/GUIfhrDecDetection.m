function [dec_mild,dec_prolong,dec_severe,dec_index]=GUIfhrDecDetection(fhr,baseline)
% GUIfhrDecDetection: Detect deceleration of FHR signl  (¸ü¸Ä×ÔfhrDecDetection.m)
% Input:
%       fhr:    fetal heart rate signal
%       baseline
% Output:
%       dec: deceleration
    n=length(fhr);
    if length(baseline)==1
        baseline=ones(1,n).*baseline;
    end
    small_time=0;
    dec_mild=0;dec_prolong=0;dec_severe=0;
    dec_index=[];
	% Method 1:
    i=1;
    while i<=n
        if (baseline(i)-fhr(i))>15
            small_time=small_time+1;
            for j=i+1:n
                if (baseline(j)-fhr(j))>15
                    small_time=small_time+1;
                else 
                    break;
                end
            end
        end
        if small_time>=15*4
            if small_time>=300*4
                dec_severe=dec_severe+1;
            elseif 120*4<small_time &&small_time<300*4
                dec_prolong=dec_prolong+1;
            else
                dec_mild=dec_mild+1;
            end
            dec_index=[dec_index; i,(j-1)];
            i=j+1;
        else
            i=i+1;
        end
        small_time=0;
    end
end