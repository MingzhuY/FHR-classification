function output=GUIfhrTimeDomain(fhr)
% fhrTimeDomain: calculate time-domain parameters of FHR
% Input:    
%       fhr: FHR signal (xlabel: sample)
% Output
%       output: relevant indexes
    output.mean=round(mean(fhr)*10)/10;     % unit: bpm
    output.sd=round(std(fhr)*10)/10;        % unit: bpm
    output.min=round(min(fhr)*10)/10;       % unit: bpm
    output.max=round(max(fhr)*10)/10;       % unit: bpm
    fs=4;   %  samples/second
    %% Calculate sFHR: the value of FHR taken each 2.5sec
    STV1=zeros(1,length(fhr)/fs/60);
    II1=zeros(1,length(fhr)/fs/60);
    j=1;
    for i=2.5*fs:24*2.5*fs:length(fhr)
        sfhr1=fhr(i:2.5*fs:i+24*2.5*fs-1);sfhr1(end)=[];
        sfhr2=fhr(i:2.5*fs:i+24*2.5*fs-1);sfhr2(1)=[];
        STV1(j)=mean(abs(sfhr2-sfhr1));
        if STV1(j)==0
            II1(j)=0;
        else
            II1(j)=STV1(j)/std([sfhr1,fhr(i+23*2.5*fs)]); % 分母不能为0
        end
        j=j+1;
    end
    output.STV=round(mean(STV1)*10)/10;
    output.II=round(mean(II1)*10)/10;	% Interval Index
    clear STV1 II1 sfhr1 sfhr2 i j 
    sfhr=fhr(1:1*60*fs);    % one-minute window
    cha=max(sfhr)-min(sfhr);
    if cha>25
        OSC_type=1; % saltatory oscillation
    elseif cha<=25 && cha>10
        OSC_type=2; % undulatory oscillation
    elseif cha<=10 && cha>=5
        OSC_type=3; % narrow undulatory oscillation
    elseif cha<5
        OSC_type=4; % silent oscillation
    end
    output.OSC_type=OSC_type;
    clear sfhr1 sfhr2 sfhr cha
	%% Calculate distribution: m
    xp=fhr;xp(end)=[]; xm=fhr;xm(1)=[];
    m=sqrt(xp.^2+xm.^2);
%    Q1=prctile(m,[25]);Q2=prctile(m,[75]);
    % help iqr
    output.LTI=round(iqr(m)*10)/10;          % Long Term Irregularity
    clear xp xm m
	%% Calculate value of max and min within each minute
    fhr_time=fhr(fs:fs:end);    % Change xlabel to second
    m=length(fhr_time)/60;      % signal length: minute
    delta=zeros(1,m);
    j=1;
    for i=1:60:length(fhr_time)-59   % signal minutes= *60 seconds
        delta(j)=max(fhr_time(i:i+59))-min(fhr_time(i:i+59)); % each minute as a segment
        j=j+1;
    end    
    output.delta=round(mean(delta)*10)/10;
    output.delta_total=round((max(fhr)-min(fhr))*10)/10;  
end