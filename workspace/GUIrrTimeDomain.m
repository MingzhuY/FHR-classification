function output=GUIrrTimeDomain(rr,win,xx,fs)
% GUIrrTimeDomain: calculate time-domain parameters of RR interval series
% Input:    
%       rr:     1 dim RR series (ms)
%       win:    window size to use for SDNNi (s)
%       xx:     value to use for NNx and pNNx (ms)
%       fs:     sample frequency (Hz)
% Output
%       output: relevant indexes
% A part of FHRAS software.
%   timeDomainHRV.m
    if size(rr,1)>size(rr,2)
        rr=rr';             % xlabel: sample    ylabel: ms
    end
    rr2=rr(fs:fs:end);      % xlabel: second    ylabel: ms
    % Calculate and round to nearest 1 decimal point   
    % unit: ms
    output.max=round(max(rr2)*10)/10;           
    output.min=round(min(rr2)*10)/10;
    output.mean=round(mean(rr2)*10)/10;         
    output.median=round(median(rr2)*10)/10;
    output.SDNN=round(std(rr2)*10)/10;          
    output.SDANN=round(SDANN(rr2,win)*10)/10;   
    output.SDNNi=round(SDNNi(rr2,win)*10)/10;
    output.RMSSD=round(RMSSD(rr2)*10)/10;
	[p,n]=pNNx(rr2,xx);
    output.NNx=round(n*10)/10;              % unit: count
    output.pNNx=round(p*10)/10;             % unit: %
    clear p n
    
    % 必须确保xp-xm不等于0    
    for i=1:24*2.5*fs:length(rr)-24*2.5*fs
        rri=rr(i:2.5*fs:i+24*2.5*fs-1);
        xp=rri;xp(end)=[];  % value of signal taken each
        xm=rri;xm(1)=[];    % 2.5 sec space (in one minute)     
        if  (xp-xm)~=0 
            break;
        end
    end
    output.STV=round(mean(abs(xm-xp))*10)/10;           % Short Term Variability
    output.II=round(std(abs(xm-xp))/output.STV*10)/10;	% Interval Index
    output.STV_HAA=round(1000*iqr(atan(xm./xp))*10)/10;
    D=(xp-xm)./(xp+xm);
    output.STV_YEH=round(sqrt(sum((D-output.mean).^2)/(60*4-2))*10)/10;
    clear i rri xp xm D 
    for i=1:60*fs:length(rr)-60*fs      % choose 60sec=240samples
        rri=rr(i:3.75*fs:i+60*fs-1);    % value of each 3.75s space
        if min(rri)~=max(rri)
            break;
        end
    end
    output.STV_Sonicaid=round(abs(max(rri)-min(rri))*10)/10;
    clear i rri
    
    m=[];
    for i=1:60:length(rr2)
        mm=rr2(i:i+60-1);   % max and min in 1 minute
        m=[m,max(mm)-min(mm)];
    end
    output.delta=round(mean(m)*10)/10;
    output.delta_total=round((max(rr2)-min(rr2))*10)/10; 
    clear m i mm
    
    xp=rr2(1:3*60);xp(end)=[]; xm=rr2(1:3*60);xm(1)=[];         
    output.LTI_HAA=round(iqr(sqrt(xp.^2+xm.^2))*10)/10; % Long Term Irregularity
    clear xp xm
    
    % Geometric FHRV
    % calculate number of bins to use in histogram
    dt=max(rr2)-min(rr2);
    binWidth=1/128*1000; % 1/128 second
    % Reference: (1996) Heart rate variability: standards of measurement,
    %                   physiological interpretation and clinical use.
    nBins=round(dt/binWidth);
    % temp   % nBins=32;
    output.FHRVTi=round(FHRVTi(rr2,nBins)*10)/10;
    output.TINN=round(TINN(rr2,nBins)*10)/10;        % unit: ms
    clear dt binWidth
    
	% Poincare
    SD=Poincare1(rr2);
    output.SD1_method2=SD.SD1_method2;
    output.SD2_method2=SD.SD2_method2;
end

function output=SDANN(rr2,t)
% SDANN: SDANN index is the std of all the mean NN intervals from each segment of length t.
%     n=length(rr2);
%     tmp=zeros(ceil(n/t),1);
%     j=1;
%     for i=1:t:n
%         tmp(j)=mean(rr2(i:i+t-1));
%         j=j+1;
%     end
    a=0;i1=1;
    tmp=zeros(ceil(sum(rr2)/t),1);
    for i2=1:length(rr2)
        if sum(rr2(i1:i2))>=t
            a=a+1;
            tmp(a)=mean(rr2(i1:i2));
            i1=i2;
        end
    end
    output=std(tmp);
end

function output=SDNNi(rr2,t)
% SDNNi: SDNNi index is the mean of all the standard deviations of NN
%        (normal RR) intervals for all window of length t.
%     n=length(rr2);
%     tmp=zeros(ceil(n/t),1);
%     j=1;
%     for i=1:t:n
%         tmp(j)=std(rr2(i:i+t-1));
%         j=j+1;
%     end
    a=0;i1=1;
    tmp=zeros(ceil(sum(rr2)/t),1);
    for i2=1:length(rr2)
        if sum(rr2(i1:i2))>=t
            a=a+1;
            tmp(a)=std(rr2(i1:i2));
            i1=i2;
        end
    end
    output=mean(tmp);
end

function output=RMSSD(rr2)
% RMSSD: root mean square of successive RR differences
    differences=abs(diff(rr2)); % successive rr2 diffs (ms)
    output=sqrt( sum(differences.^2)/length(differences) );
end

function [p n]=pNNx(rr2,x)
% pNNx: percentage of successive/adjacement NN intereval differing by x(ms) or more. 
    differences=abs(diff(rr2));
    n=sum(differences>x);
    p=(n/length(differences))*100;
end

function output=FHRVTi(rr2,nbin)
% FHRVTi: FHRV triangular index
    % calculate samples in bin(n) and x location of bins(xout)
    [n,~]=hist(rr2,nbin);
    output=length(rr2)/max(n); % FHRV ti
end

function output=TINN(rr2,nbin)
% TINN: triangular interpolation of NN interval histogram
% Reference: Standards of Measurement, Physiological Interpretation, 
%            and Clinical Use Circulation. 1996; 93(5):1043-1065.
    % calculate histogram of rr2 using nbin bins
    [nout,xout]=hist(rr2,nbin);
%     % 补救措施
%     if nout(1)==max(nout)
%         nout=nout(2:end);
%         xout=xout(2:end);
%         nbin=nbin-1;
%     end
%     if nout(end)==max(nout)
%         nout=nout(1:end-1);
%         xout=xout(1:end-1);
%         nbin=nbin-1;
%     end
    D=nout;
    peaki=find(D==max(D));
    if length(peaki)>1
        peaki=round(mean(peaki));
    end  
    i=1;
    d=zeros((peaki-1)*(nbin-peaki),3);
    for m=(peaki-1):-1:1
        for n=(peaki+1):nbin
            % define triangle that fits the histogram
            q=zeros(1,length(D));
            q(1:m)=0;
            q(n:end)=0;
            q(m:peaki)=linspace(0,D(peaki),peaki-m+1);
            q(peaki:n)=linspace(D(peaki),0,n-peaki+1);
            % integrate squared difference
            d(i,1)=trapz((D-q).^2);
            d(i,2:3)=[m,n];
            % plot(D),hold on, plot(q,'r'),hold off
            % title([ 'd^2=' num2str(d(i,1))]);
            i=i+1;
        end
    end
    % find where minimum square diff occured
    i=find(d(:,1)==min(d(:,1)));
    i=i(1);         % make sure there is only one choise
    m=d(i,2); n=d(i,3);
    % calculate TINN in (ms)
    output=abs(xout(n)-xout(m));  
%     % plot    
%     X=xout(peaki); M=xout(m); N=xout(n); Y=nout(peaki);
%     hist(rr2,nbin)
%     xlimits=get(gca,'xlim');
%     hold on;
%     plot(xout,nout,'k')
%     line([M X N M],[0 Y 0 0],'color','r','linewidth',1.5,'LineStyle','--')
%     line([X X],[0 1000],'LineStyle','-.','color','k')
%     line([0 2000],[Y Y],'LineStyle','-.','color','k')
%     colormap white
%     xlabel('RR (ms)');
%     ylabel('Number of RR')
%     legend({'Histogram','D(t)','q(t)'})
%     set(gca,'xtick',[xout(m) xout(peaki) xout(n)],'xticklabel',{'N','X','M'}, ...
%         'ytick',Y,'yticklabel','Y')
%     set(gca,'xlim',xlimits);
end

function output=Poincare1(rr)
% poincareHRV: calculates poincare HRV
% Input:    
%       rr:     1 dim array (ylabel:ms)  
% Output:   a structure containg HRV.
    %% 2. method2
    rr=rr(1:500);
    xp=rr;xp(end)=[];
    xm=rr;xm(1)=[];
    SD1_method2=std(xp-xm)/sqrt(2);
    SD2_method2=std(xp+xm)/sqrt(2);
    SDRR_method2=sqrt(SD1_method2^2+SD2_method2^2)/sqrt(2);
    % format decimal places
    output.SD1_method2=round(SD1_method2*10)/10;    % ms
    output.SD2_method2=round(SD2_method2*10)/10;    % ms    
    output.SDRR_method2=round(SDRR_method2*10)/10;  % ms
end
        