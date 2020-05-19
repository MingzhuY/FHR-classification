function output=GUIrrFreqDomain(rr2,VLF,LF,MF,HF,AR_order,window,noverlap,nfft,fs,methods,tu)
% GUIfhrFreqDomain: Calculates frequency domain parameters of rr using FFT(Welch), 
%                 AR(Burger),and Lomb-Scargle(LS) methods For Power Spectral Density(PSD).
% Inputs:   rr2         = 2*N array of time (s) and inter-beat interval (ms)
%           VLF,LF      = Very Low Frequency,Low Frequency
%           MF,HF       = Mild Frequency (fetal movement),High Frequency
%           AR_order    = order of AR model
%           window      = number of samples in window 
%           noverlap    = number of samples to overlap (must <window)
%           nfft        = number of points in the frequency axis—Points in PSD (vector length)
%           fs          = sample rate of signal (Hz)
%           methods     = methods of calculating freqDomain. Default is all 3 methhods
%           tu          = whether to plot (1=yes,2=no)
% Outputs:  output is a structure containg all Fhrv, one field for each PSD method.
%           Output units include:
%               peakVLF,peakLF,peakMF,peakHF                        (Hz)
%               power_VLF,power_LF,power_MF,power_HF                (ms^2)
%               percent_VLF,percent_LF,percent_MF,percent_HF        (%)
%               nLF,nMf,nHF                                         ()
%               ratio                                               ()
%               psd                             (ms^2/Hz)   (nfft/2+1)*1 array onesided
%               f                               (Hz)            nfft*1 	 array twosided
% A part of FHRAS software.
%   freqDomainHRV.m
%       lomb2.m
    % Check inputs
    if nargin<10
        error('Not enough input arguments!')
    elseif nargin<11
        methods={'welch','ar','lomb'};
        tu=2;
    elseif nargin<12
        tu=2;
    elseif nargin>12
        error('Too many input argumentgs')
    end    
    % Choose methods for PSD calculation
    flagWelch=false; flagAR=false; flagLomb=false;
    for m=1:length(methods)
        if strcmpi(methods{m},'welch')
            flagWelch=true;
        elseif strcmpi(methods{m},'ar')
            flagAR=true;
        elseif strcmpi(methods{m},'lomb')
            flagLomb=true;
        end
    end   
    % Input signal/data
    t=rr2(1,:);     % time (1s interval)
    y=rr2(2,:);     % analyzed signal (ms)
    maxF=fs/2;      % used for emptyData.m (onesided 单边带频率)
    window=hamming(window); % 试验发现PSD结果与直接使用window，不加窗相同，即这行程序可以省去
    % Prepare
    y=detrend(y,'linear');
    y=y-mean(y);
    % Frequncy range (energy band) of PSD
    % Welch FFT
    if flagWelch
        [output.welch.psd,output.welch.f]=calcWelch(y,window,noverlap,nfft,fs);
        output.welch.fhrv=calcAreas(output.welch.f,output.welch.psd,VLF,LF,MF,HF,0);
    else
        output.welch=emptyData(nfft,maxF);
    end
    % AR(Burg)
    if flagAR
        [output.ar.psd,output.ar.f]=calcAR(y,AR_order,nfft,fs);
        output.ar.fhrv=calcAreas(output.ar.f,output.ar.psd,VLF,LF,MF,HF,0);
    else
        output.ar=emptyData(nfft,maxF);
    end 
    % Lomb
    if flagLomb
        [output.lomb.psd,output.lomb.f]=calcLomb(t,y,nfft,maxF);
        output.lomb.fhrv=calcAreas(output.lomb.f,output.lomb.psd,VLF,LF,MF,HF,0);
    else
        output.lomb=emptyData(nfft,maxF);
    end
    % plot all three psd
    if tu==1
        plotPSD(output.welch.f,output.welch.psd,VLF,LF,MF,HF,[0 0.6],[]);
        title('Frequency analysis of FHRV using Welch(FFT)','Fontsize',8)
    end
end

function [PSD,F]=calcWelch(y,window,noverlap,nfft,fs)
% calcWelch: Calculates the PSD using Welch method.
    % Calculate PSD—help pwelch
    [PSD,F]=pwelch(y,window,noverlap,nfft,fs,'onesided'); % uses a hamming window     
end

function output=calcAreas(F,PSD,VLF,LF,MF,HF,flagNorm)
% calcAreas: Calculates areas/energy under the PSD curve within the freq
%            bands defined by VLF, LF, MF and HF. 
% Inputs:
%       PSD:            PSD vector
%       F:              Frequency vector
%       VLF,LF,MF,HF:   array containing VLF, LF,MF and HF freq limits
%       flagNormalize:  option to normalize PSD to max(PSD)
    % Check inputs
    if nargin<6
        error('Not enough input arguments')
    elseif nargin<7
        flagNorm=false;
    elseif nargin>7
        error('Too many input arguments')
    end   
    % Normalize PSD if needed
    if flagNorm
        PSD=PSD/max(PSD);
    end
    % Find the indexes corresponding to the VLF,LF,MF and HF bands
    iVLF= (F>=VLF(1)) & (F<VLF(2));
    iLF = (F>=LF(1))  & (F<LF(2));
    iMF = (F>=MF(1))  & (F<MF(2));
    iHF = (F>=HF(1))  & (F<HF(2));   
    % Calculate percent of FHRV-PSD bandwith
    percentVLF  = length(find(iVLF==1))/length(F)*100;
    percentLF   = length(find(iLF==1)) /length(F)*100;
    percentMF   = length(find(iMF==1)) /length(F)*100;
    percentHF   = length(find(iHF==1)) /length(F)*100;
    % Find peaks
      % VLF Peak
      tmpF=F(iVLF); tmpPSD=PSD(iVLF);
      [pks,ipks]=zipeaks(tmpPSD);
      if ~isempty(pks)
        [~,i]=max(pks);    peakVLF=tmpF(ipks(i));
      else
        [~,i]=max(tmpPSD); peakVLF=tmpF(i);
      end
       % LF Peak
       tmpF=F(iLF); tmpPSD=PSD(iLF);
       [pks,ipks]=zipeaks(tmpPSD);
       if ~isempty(pks)
         [~,i]=max(pks);    peakLF=tmpF(ipks(i));
       else
         [~,i]=max(tmpPSD);	peakLF=tmpF(i);
       end
       % MF Peak
       tmpF=F(iMF);	tmpPSD=PSD(iMF);
       [pks,ipks]=zipeaks(tmpPSD);
       if ~isempty(pks)
         [~,i]=max(pks);	peakMF=tmpF(ipks(i));
       else
         [~,i]=max(tmpPSD); peakMF=tmpF(i);
       end       
       % HF Peak
       tmpF=F(iHF); tmpPSD=PSD(iHF);
       [pks,ipks]=zipeaks(tmpPSD);
       if ~isempty(pks)
         [~,i]=max(pks);    peakHF=tmpF(ipks(i));
       else
         [~,i]=max(tmpPSD); peakHF=tmpF(i);
       end   
    % Calculate raw areas (power under curve), within the freq bands (ms^2)
    aVLF   = trapz(F(iVLF),PSD(iVLF));
    aLF    = trapz(F(iLF), PSD(iLF));
    aMF    = trapz(F(iMF), PSD(iMF));
    aHF    = trapz(F(iHF), PSD(iHF));
    aTotal = aVLF+aLF+aMF+aHF;
    % Calculate areas relative to the total area (%)
    pVLF = (aVLF/aTotal)*100;
    pLF  = (aLF/aTotal) *100;
    pMF  = (aMF/aTotal) *100;
    pHF  = (aHF/aTotal) *100;
    % Calculate normalized areas 
    nLF = aLF/(aLF+aMF+aHF);
    nMF = aMF/(aLF+aMF+aHF);
    nHF = aHF/(aLF+aMF+aHF);
    % Calculate LF/HF ratio
    lfhf=aLF/aHF; 
    ratio=aLF/(aMF+aHF);
    % Create output structure
    output.percentVLF   = round(percentVLF*100)/100;    % round
    output.percentLF    = round(percentLF*100)/100;
    output.percentMF    = round(percentMF*100)/100;
    output.percentHF    = round(percentHF*100)/100;
	output.power_VLF 	= round(aVLF*100)/100;          % ms^2
	output.power_LF   	= round(aLF*100)/100;
    output.power_MF    	= round(aMF*100)/100;
	output.power_HF     = round(aHF*100)/100;
	output.power_Total 	= round(aTotal*100)/100;
    output.percent_VLF	= round(pVLF*100)/100;          % %
    output.percent_LF  	= round(pLF*100)/100;
    output.percent_MF  	= round(pMF*100)/100;
    output.percent_HF 	= round(pHF*100)/100;
    output.nLF          = round(nLF*100)/100;
    output.nMF          = round(nMF*100)/100;
    output.nHF          = round(nHF*100)/100;
    output.lfhf         = round(lfhf*100)/100;
    output.ratio        = round(ratio*100)/100;
    output.peakVLF      = round(peakVLF(1)*100)/100;
    output.peakLF       = round(peakLF(1)*100)/100;
    output.peakMF       = round(peakMF(1)*100)/100;
    output.peakHF       = round(peakHF(1)*100)/100;
end

function [pks,locs]=zipeaks(y)
% zippeaks: finds local maxima of input signal y.
% Output: 2x(number of maxima) array
%       pks     : value at maximum
%       locs	: index value for maximum
    % Check inputs
    if isempty(y)
        pks=[]; locs=[];
        return
    end
    [rows,cols]=size(y);
    if cols==1 && rows>1 % all data in 1st col
        y=y';            % input must be row vector
    elseif cols==1 && rows==1 
        pks=[]; locs=[];
        return    
    end         
    % Find locations of local maxima
    % yD=1 at maxima, yD=0 otherwise, end point maxima excluded
    N=length(y)-2;
    yD=[0 (sign(sign(y(2:N+1)-y(3:N+2))-sign(y(1:N)-y(2:N+1))-0.1)+1) 0];
    % Indices of maxima and corresponding values of y
    Y=logical(yD);
    I=1:length(Y);
    locs=I(Y);
    pks=y(Y);
end

function output=emptyData(nfft,maxF)
% emptyData: create output structure of zeros
    output.fhrv.percentVLF	=0;
    output.fhrv.percentLF  	=0;
    output.fhrv.percentMF 	=0;
    output.fhrv.percentHF 	=0;
	output.fhrv.areaVLF     =0;
	output.fhrv.areaLF      =0;
    output.fhrv.areaMF      =0;
	output.fhrv.areaHF      =0;
	output.fhrv.areaTotal	=0;
    output.fhrv.pVLF        =0;
    output.fhrv.pLF         =0;
    output.fhrv.pMF         =0;
    output.fhrv.pHF         =0;
    output.fhrv.nLF         =0;
    output.fhrv.nMF         =0;
    output.fhrv.nHF         =0;
    output.fhrv.LFHF        =0;
    output.fhrv.ratio       =0;
    output.fhrv.peakVLF     =0;
    output.fhrv.peakLF      =0;
    output.fhrv.peakMF      =0;
    output.fhrv.peakHF      =0;
    % PSD with all zeros 
    output.f                =linspace(0.0,maxF,nfft/maxF+1)';
    output.psd              =zeros(length(output.f),1);
end

function [PSD,F]=calcAR(y,AR_order,nfft,fs)
% calAR: Calculates the PSD using Auto Regression model. 
    % Calculate PSD
    % Method 1
%     [A, variance]=arburg(y,AR_order); % AR using Burg method
%     [H,F]=freqz(1,A,nfft,fs);
%     PSD=(abs(H).^2).*(variance/fs);   % malik, p.67    
    % Method 2
    [PSD,F]=pburg(y,AR_order,nfft,fs,'onesided');
    % Method 3
%      h=spectrum.burg;
%      hpsd=psd(h,y,'NFFT',nfft,'Fs',fs);
%      F=hpsd.Frequencies;
%      PSD=hpsd.Data;
end

function [PSD,F]=calcLomb(t,y,nfft,maxF)
% calLomb: Calculates the PSD using Lomb-Scargle method.        
    % Calculate PSD
    deltaF=0;       % deltaF=maxF/nfft;
    F=linspace(0.0,maxF-deltaF,nfft/maxF+1);
    PSD=lomb2(t,y,F,false); % calc lomb psd
end
function [Pn]=lomb2(t,y,f,flagNorm)
%  Uses Lomb's method to compute normalized periodogram values "Pn" as a function of
%  supplied vector of frequencies "f" for input vectors "t" (time) and "y" (observations).
%  Also returned is probability "Prob" of same length as Pn (and f) that the null hypothesis
%  is valid. If f is not supplied it assumes f =  [1/1024 : 1/1024 : 0.5/min(diff(t))];		
%  x and y must be the same length.
% See also:  
% [1] N.R. Lomb, ``Least-squares frequency analysis of 
% unequally spaced data,'' Astrophysical and Space Science, 
% (39) pp. 447--462, 1976.   ... and 
% [2] J.~D. Scargle, ``Studies in astronomical time series analysis. 
% II. Statistical aspects of spectral analysis of unevenly spaced data,''
% Astrophysical Journal, vol. 263, pp. 835--853, 1982.
% [3] T. Thong, "Lomb-Welch Periodogram for Non-uniform Sampling",
% Proceedings for the 26th anual international conference of the IEEE EMBS,
% Sept 1-5, 2004.
    if nargin<4
        flagNorm=true;
    end
    if nargin<3
        upper_freq = 0.5/min(diff(t));
        f =  [1/1024 : 1/1024 : upper_freq];
    end
    if nargin < 2
        fprintf('assuming y=diff(t)\n');
        y=diff(t); % RR tachogram?
        t=t(1:length(y)); % shorter by one time stamp now
    end
    % check inputs
    if length(t) ~= length(y); 
        error('t and y not same length');
    end
    % subtract mean, compute variance, initialize Pn
    z = y - mean(y);
    var = std(y)^2;
    N=length(f);
    Pn=zeros(size(f));
    % now do main loop for all frequencies
    for i=1:length(f)
        w=2*pi*f(i);
        if w > 0 
            twt = 2*w*t;
            tau = atan2(sum(sin(twt)),sum(cos(twt)))/2/w;
            wtmt = w*(t - tau);
            Pn(i) = (sum(z.*cos(wtmt)).^2)/sum(cos(wtmt).^2) + ...
            (sum(z.*sin(wtmt)).^2)/sum(sin(wtmt).^2);
        else
            Pn(i) = (sum(z.*t).^2)/sum(t.^2);
        end
    end
    if flagNorm     % normalize by variance
        Pn=Pn./(2*var);
    else            % return denormalized spectrum (see T. Thong)
        Pn=Pn./length(y);
    end
end

function plotPSD(F,PSD,VLF,LF,MF,HF,limX,limY)
    % Find the indexes corresponding to the VLF, LF, and HF bands
    iVLF= find( (F>=VLF(1)) & (F<VLF(2)) );
    iLF = find( (F>=LF(1))  & (F<LF(2))  );
    iMF = find( (F>=MF(1))  & (F<MF(2))  );
    iHF = find( (F>=HF(1))  & (F<HF(2))  );
    PSD=PSD/max(PSD);
    % Plot area under PSD curve
    area(F(:),PSD(:),'FaceColor',[0.80 0.80 0.80]);        
    hold on
    area(F(iVLF(1):iVLF(end)+1),PSD(iVLF(1):iVLF(end)+1),'FaceColor',[0.9882 0.902 0.902]);
    area(F(iLF(1):iLF(end)+1),PSD(iLF(1):iLF(end)+1),'FaceColor',[0.9255 0.9255 0.9882]);
    area(F(iMF(1):iMF(end)+1),PSD(iMF(1):iMF(end)+1),'FaceColor',[0.9255 0.9255 0.9882]);
    area(F(iHF(1):iHF(end)),PSD(iHF(1):iHF(end)),'FaceColor',[1 1 0.8863]);
    if isempty(limX)
        limX=[min(F) max(F)];
    else
        set(gca,'xlim',limX)
    end
    if isempty(limY)
        limY=[min(PSD) max(PSD)];
    else
        set(gca,'ylim',limY)
    end
    % Draw vertical lines around freq bands
    line1=line([VLF(2) VLF(2)],[limY(1) limY(2)]);
    set(line1,'color',[1 0 0],'parent',gca);
    line2=line([LF(2) LF(2)],[limY(1) limY(2)]);
    set(line2,'color',[1 0 0],'parent',gca);
    line3=line([MF(2) MF(2)],[limY(1) limY(2)]);
    set(line3,'color',[1 0 0],'parent',gca);
    line4=line([HF(2) HF(2)],[limY(1) limY(2)]);
    set(line4,'color',[1 0 0],'parent',gca);
    text(VLF(2)-0.04,limY(2)-0.2,{' VLF';'band'});
    text((LF(1)+LF(2))/2,limY(2)-0.2,{'  LF';'band'});
    text((MF(1)+MF(2))/2,limY(2)-0.2,{'  MF';'band'});
    text((HF(1)+HF(2))/2,limY(2)-0.2,{'  HF';'band'});
    xlabel('Frequency/Hz','Fontsize',8),ylabel('PSD/Normalized','Fontsize',8)
    set(gca,'Fontsize',8)   % 设置坐标刻度大小
    hold off
end
