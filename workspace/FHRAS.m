%% 在创建fig时,matlab自动产生的代码
function varargout=FHRAS(varargin)
% FHRAS MATLAB code for FHRAS.fig
%      FHRAS, by itself, creates a new FHRAS or raises the existing
%      singleton*.
%
%      H = FHRAS returns the handle to a new FHRAS or the handle to
%      the existing singleton*.
%
%      FHRAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FHRAS.M with the given input arguments.
%
%      FHRAS('Property','Value',...) creates a new FHRAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FHRAS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FHRAS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help fhras
% Last Modified by GUIDE v2.5 30-Dec-2017 17:08:49
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FHRAS_OpeningFcn, ...
                   'gui_OutputFcn',  @FHRAS_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% --- Executes just before FHRAS is made visible.
function FHRAS_OpeningFcn(hObject,eventdata,handles,varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FHRAS (see VARARGIN)
% Choose default command line output for FHRAS
handles.output=hObject;
mainHandle=interface();     % 启动初始界面,调用interface.m的GUI
pause(3);
close(mainHandle);
% Update handles structure
guidata(hObject,handles);
% UIWAIT makes FHRAS wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout=FHRAS_OutputFcn(hObject,eventdata,handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
javaFrame=get(gcf,'JavaFrame'); % 启动界面时，直接最大化
set(javaFrame,'Maximized',1);
varargout{1}=handles.output;

%% 界面显示设置
%% edit,text 文本框：输入,显示
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject,eventdata,handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string',30); % 设置默认显示值
function edit2_CreateFcn(hObject,eventdata,handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string',36); % 设置默认显示值
function edit3_CreateFcn(hObject,eventdata,handles)
% 输入信号长度（分）
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% % 以字符串的形式来存储数据文本框1的内容.如果字符串不是数字，则显示空白内容
% % 这段代码使得输入被严格限制，我们不能试图输入一个非数字。
% input=str2num(get(hObject,'String'));
% % 检查输入是否为空. 如果为空,则默认显示为0
% if (isempty(input))
%   set(hObject,'String','0')
% end
% guidata(hObject,handles);
set(hObject,'string',20); % 设置默认显示值
function text7_CreateFcn(hObject,eventdata,handles)
% str1=('一切准备就绪 ！');
% str2=('开始分析信号 ！');
str1=('Everything is ready ！');
str2=('Begin to analyze the signal ！');
set(hObject,'string',{str1;str2});
%% axes 图形显示
% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject,eventdata,handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% set(hObject,'visible','off');     % 清除axes的全部
set(hObject,'xTick',[]);            % 不显示axes的坐标轴
set(hObject,'ytick',[]);
function axes2_CreateFcn(hObject,eventdata,handles)
% set(hObject,'visible','off');     % 清除axes的全部
set(hObject,'xTick',[]);            % 不显示axes的坐标轴
set(hObject,'ytick',[]);
function axes3_CreateFcn(hObject,eventdata,handles)
% set(hObject,'visible','off');     % 清除axes的全部
set(hObject,'xTick',[]);            % 不显示axes的坐标轴
set(hObject,'ytick',[]);
function axes4_CreateFcn(hObject,eventdata,handles)
% set(hObject,'visible','off');     % 清除axes的全部
set(hObject,'xTick',[]);            % 不显示axes的坐标轴
set(hObject,'ytick',[]);    
%% uitable 表格显示数据
% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject,eventdata,handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% set(hObject,'ColumnName',{'特征参数','数值','单位'}) % 设置列标题
set(hObject,'ColumnName',{'Feature','Value','Unit'})    % 设置列标题
set(hObject,'ColumnWidth',{150,60,50})                  % 设置列宽
set(hObject,'ColumnEditable',logical(ones(1,3)))        % 设置列表可编辑

%% 第1步: open菜单 打开、导入信号数据
function Open_Callback(hObject,eventdata,handles)
function Open_data_Callback(hObject,eventdata,handles)
% 打开文件，导入数据 help uigetfile
[filename,pathname]=uigetfile(...  
    {'*.mat';'*.txt';'*.xls';'*.mat;*.txt;*.xls'},'Choose one FHR signal'); % 设置默认路径导入数据 
if isequal(filename,0)  
  % disp('User pressed Cancel.') 是在MATLAB命令窗中显示  
  set(handles.text7,'string',{'The user pressed Cancel.';'Data not loaded !'}); % 显示在GUI界面上
else  
  % load得到的是struct类型的数据，所以必须要转换成矩阵格式，否则不能做任何操作
  val=cell2mat(struct2cell(load([pathname,filename])));
  handles.name=filename(1:end-4);   % 取得当前信号的标号:字符串格式
  handles.data=val(1,:)/100;
  % disp(['User selected: ',fullfile(pathname,filename)])  
  str1=(['The user selected: ',fullfile(pathname,filename)]);
  % disp('Data loaded successful!')
  str2=('Data loaded successfully !');
  set(handles.text7,'string',{str1;str2});
  % 当打开、导入第二个信号数据时，清空axes中的图形，否则会与第一个信号的波形将发生重叠
  axes(handles.axes1);    cla reset
  axes(handles.axes2);    cla reset
  axes(handles.axes3);    cla reset
  axes(handles.axes4);    cla reset
end
guidata(hObject,handles);
function Open_exit_Callback(hObject,eventdata,handles)
% 退出GUI
clc;clear;close all

%% 第2步: buttton按钮 Step1-9 信号处理
% --- Executes on button press in FHR_original_plot.
function FHR_original_plot_Callback(hObject,eventdata,handles)
% hObject    handle to FHR_original_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fhr_original=handles.data;
n=length(fhr_original); % 所画信号长度
axes(handles.axes1);    cla reset
grid off
hold on
% Horizontal axis
mat1=[255 245 245]/255;
mat2=[255 230 230]/255;
for i=0:10:240 % y轴每隔10横画一条直线，平行x轴（颜色较浅）
    line([0,n],[i,i],'color',mat1);
end
for i=0:30:240 % y轴每隔30横画一条直线，平行x轴（颜色较深）
    line([0,n],[i,i],'color',mat2);
end
% Verical axis
for i=0:n/60:n   % x轴每隔time/60竖画一条直线，平行y轴（颜色较浅）
    line([i,i],[0,240],'color',mat1);
end
for i=0:n/10:n	% x轴每隔time/10竖画一条直线，平行y轴（颜色较浅）
    line([i,i],[0,240],'color',mat2);
end
set(gca,'xtick',[]);
set(gca,'ytick',[]);
for j=n/6:n/3:n
    for i=0:30:240
        text(j,i,num2str(i),'FontSize',8);
    end
    text(j+2400,0,num2str(j+2400),'Fontsize',8);
end
text(0,0,num2str(0),'Fontsize',8);
plot(1:n,fhr_original)
xlabel('Sample','Fontsize',8);ylabel('Beat per minute','Fontsize',8)
title('Fetal Heart Rate Signal','Fontsize',8)
set(handles.text7,'string','Original FHR signal plotted successfully !');
guidata(hObject,handles);
% --- Executes on button press in FHR_preprocess_plot.
function FHR_preprocess_plot_Callback(hObject,eventdata,handles)
fhr_original=handles.data;
time=str2num(get(handles.edit3,'string'));  % unit:min
[fhr4,n,loss,quality]=GUIfhrPreprocess(fhr_original,time,4);
% 画图：预处理后的FHR信号
axes(handles.axes2);    cla reset
plot((0:n-1)/4,fhr4,'b');
xlabel('Time/s','Fontsize',8);ylabel('Beat per minute','Fontsize',8)
set(gca,'Fontsize',8)   % 设置坐标刻度大小
title(['FHR preprocessed ',num2str(time),'min'],'Fontsize',8)
% text(n-1000,170,['Loss=',num2str(loss),'%']) 
% text(n-1000,180,['Quality=',num2str(quality),'%']) 
handles.fhr1=fhr_original;	% 原始信号
handles.fhr2=fhr4;          % 去噪信号
handles.loss=loss;          % 信号质量
handles.quality=quality;
set(handles.text7,'string','FHR signal preprocessed successfully !');
guidata(hObject,handles);
% --- Executes on button press in Baseline_estimation.
function Baseline_estimation_Callback(hObject,eventdata,handles)
baseline=GUIfhrBaseline(handles.fhr2);
n=length(baseline);
axes(handles.axes2);    cla reset
plot((0:n-1)/4,handles.fhr2,'b')
hold on
plot((0:n-1)/4,baseline,'r')
xlabel('Time/second','Fontsize',8);ylabel('Beat per minute','Fontsize',8)
set(gca,'Fontsize',8)   % 设置坐标刻度大小
axis([0,n/4,0,200])
legend('FHR','Baseline','Location','best')
title('Baseline calculated from preprocessed FHR','Fontsize',8)
handles.baseline=baseline;
handles.meanBL=round(mean(baseline)*10)/10;  handles.sdBL=round(std(baseline)*10)/10;
handles.minBL=round(min(baseline)*10)/10;    handles.maxBL=round(max(baseline)*10)/10;
set(handles.text7,'string','Baseline calculated successfully !');
guidata(hObject,handles);
% --- Executes on button press in ACC_DEC_Detection.
function ACC_DEC_Detection_Callback(hObject,eventdata,handles)
fhr=handles.fhr2;
baseline=handles.baseline;
n=length(fhr);
[acc_time,acc_index]=GUIfhrAccDetection(fhr,baseline);
[dec_mild,dec_prolong,dec_severe,dec_index]=GUIfhrDecDetection(fhr,baseline);
axes(handles.axes2);    cla reset
plot((0:n-1)/4,fhr,'b')
hold on         % 为了将加速和减速同时显示在图上
plot((0:n-1)/4,baseline,'r')
[x,~]=size(acc_index);
for i=1:x
    acc_mat=acc_index(i,1):acc_index(i,2);
    t1=(acc_mat(1):acc_mat(end))/4;
	plot(t1,fhr(acc_mat),'g')
end
[x,~]=size(dec_index);
for i=1:x
    dec_mat=dec_index(i,1):dec_index(i,2);
    t1=(dec_mat(1):dec_mat(end))/4;
	plot(t1,fhr(dec_mat),'k')
end
text(150,80,'Fetal heart rate','color','b','fontsize',12)
text(150,60,'Baseline','color','r','fontsize',12)
text(150,40,['Acceleration:',num2str(acc_time)],'color','g','fontsize',12)
text(150,20,['Deceleration:',num2str(dec_mild+dec_prolong+dec_severe)],'color','k','fontsize',12)
xlabel('Time/second','Fontsize',8);ylabel('Beat per minute','Fontsize',8)
set(gca,'Fontsize',8);axis([0,n/4,0,200])
title('FHR+BL+ACC+DEC 20min','Fontsize',8)   % 为了将加速和减速同时显示在图上
hold off
handles.acc=acc_time;
handles.dec_mild=dec_mild;
handles.dec_prolong=dec_prolong;
handles.dec_severe=dec_severe;
data=[{'Baseline mean',handles.meanBL,            '     bpm'};...   元胞数组：字符串+数字
      {'Baseline std',handles.sdBL,        	      '     bpm'};...
      {'Baseline min',handles.minBL,              '     bpm'};...
      {'Baseline max',handles.maxBL,              '     bpm'};...
      {'Acceleration',handles.acc,                '     count'};...
      {'Deceleration_mild',handles.dec_mild,      '     count'};...
      {'Deceleration_prolong',handles.dec_prolong,'     count'};...
      {'Deceleration_severe',handles.dec_severe,  '     count'}];
set(handles.uitable1,'data',data);                                  % 将数据显示在表格中
str1=('ACC and DEC detected successfully !');
str2=('Morphological features of FHR signal extracted successfully ！');
set(handles.text7,'string',{str1;str2});
guidata(hObject,handles);
% --- Executes on button press in FHR_time_domain.
function FHR_time_domain_Callback(hObject,eventdata,handles)
% FHR 时域参数
output=GUIfhrTimeDomain(handles.fhr2);                      % 调用函数
handles.FHRmean=output.mean;    handles.FHRstd=output.sd;   % 便于全局调用
handles.FHRmin=output.min;      handles.FHRmax=output.max;
handles.FHRSTV=output.STV;      handles.FHRII=output.II;    handles.FHRLTI=output.LTI;
handles.FHRdelta=output.delta;  handles.FHRdelta_total=output.delta_total;
handles.FHROSC_type=output.OSC_type;
axes(handles.axes3);   cla reset                                     % 画图axes4
hist(handles.fhr2,output.min:output.max);
xlabel('FHR/bpm','Fontsize',8),ylabel('Number of samples','Fontsize',8),set(gca,'Fontsize',8)  % 设置坐标刻度大小
title('FHR histogram','Fontsize',8)
data=[{'Signal loss',handles.loss,              '      %'};...
      {'Signal quality',handles.quality,        '      %'};...
      {'FHR mean',handles.FHRmean,              '     bpm'};...	元胞数组：字符串+数字
      {'FHR sd',handles.FHRstd,                 '     bpm'};...
      {'FHR min',handles.FHRmin,                '     bpm'};...
      {'FHR max',handles.FHRmax,                '     bpm'};...
      {'FHR STV',handles.FHRSTV,                '     bpm'};...
      {'FHR LTI',handles.FHRLTI,                '     bpm'};...
      {'FHR II',handles.FHRII,                  '       /'};...
      {'FHR delta',handles.FHRdelta,            '     bpm'};...
      {'FHR delta total',handles.FHRdelta_total,'     bpm'};...
      {'FHR OSC type',handles.FHROSC_type,      '       /'}];
set(handles.uitable1,'data',data);                          % 将数据显示在表格中
set(handles.text7,'string','Fetal heart rate time domain parameters extracted successfully ！');
guidata(hObject,handles);
% --- Executes on button press in FHRtoRR.
function FHRtoRR_Callback(hObject,eventdata,handles)
[rr,rr2]=GUIfhr2rr(handles.fhr2,4);
axes(handles.axes4);    cla reset
plot((1:length(rr))/4/60,rr)
xlabel('Time/min','Fontsize',8),ylabel('RR/ms','Fontsize',8),set(gca,'Fontsize',8)   % 设置坐标刻度大小
title(' Epoch-to-epoch variation')
handles.rr=rr;
handles.rr2=rr2;
str1=('RR converted successfully !');
str2=('Start analysing FHR variability signal !');
set(handles.text7,'string',{str1;str2});
guidata(hObject,handles);
% --- Executes on button press in RR_time_domain.
function RR_time_domain_Callback(hObject,eventdata,handles)
% RR 时域参数
output=GUIrrTimeDomain(handles.rr,100,50,4);                        % 调用函数
handles.RRmean=output.mean;     handles.RRmin=output.min;   handles.RRmax=output.max;
handles.RRmedian=output.median; handles.RRSDNN=output.SDNN; handles.RRSDANN=output.SDANN;
handles.RRSDNNi=output.SDNNi;   handles.RRRMSSD=output.RMSSD;       % 便于全局调用   
handles.RRNNx=output.NNx;       handles.RRpNNx=output.pNNx;     
handles.RRSTV=output.STV;       handles.RRII=output.II;     handles.RRLTI=output.LTI_HAA;
handles.RRdelta=output.delta;   handles.RRdelta_total=output.delta_total;
handles.RRFHRVTi=output.FHRVTi; handles.RRTINN=output.TINN;
% plot histogram
axes(handles.axes3);    cla reset
hist(handles.rr,20);
xlabel('RR/ms','Fontsize',8),ylabel('Number of samples','Fontsize',8),set(gca,'Fontsize',8)   % 设置坐标刻度大小
title('RR histogram','Fontsize',8)
data=[{'RR mean',handles.RRmean,              '     ms'};...   元胞数组：字符串+数字
      {'RR min',handles.RRmin,                '     ms'};...
      {'RR max',handles.RRmax,                '     ms'};...
      {'RR median',handles.RRmedian,          '     ms'};...
      {'RR SDNN',handles.RRSDNN,              '     ms'};...
      {'RR SDANN',handles.RRSDANN,            '     ms'};...
      {'RR SDNNi',handles.RRSDNNi,            '     ms'};...
      {'RR RMSSD',handles.RRRMSSD,            '     ms'};...
      {'RR NNx',handles.RRNNx,                '   sample'};...
      {'RR pNNx',handles.RRpNNx,              '      %'};...
      {'RR STV',handles.RRSTV,                '     ms'};...
      {'RR II',handles.RRII,                  '      /'};...
      {'RR LTI',handles.RRLTI,                '     ms'};...
      {'RR delta',handles.RRdelta,            '     ms'};...
      {'RR delta total',handles.RRdelta_total,'     ms'};...
      {'RR FHRVTi',handles.RRFHRVTi,          '       /'};...
      {'RR TINN',handles.RRTINN,              '       /'}];
set(handles.uitable1,'data',data);                                  % 将数据显示在表格中
set(handles.text7,'string',...
    'Fetal heart rate variability time domain parameters extracted successfully ！');
guidata(hObject,handles);
% --- Executes on button press in RR_freq_domain.
function RR_freq_domain_Callback(hObject,eventdata,handles)
% RR 频域参数
axes(handles.axes4);    cla reset
output=GUIrrFreqDomain(handles.rr2,[0,0.05],[0.05,0.15],[0.15,0.30],[0.30,0.50],16,512,256,1024,4,{'welch','',''},1);   % 调用函数
handles.RRpowerVLF=output.welch.fhrv.power_VLF;    handles.RRpowerLF=output.welch.fhrv.power_LF;
handles.RRpowerMF=output.welch.fhrv.power_MF;      handles.RRpowerHF=output.welch.fhrv.power_HF;
handles.RRpowerTotal=output.welch.fhrv.power_Total;
handles.RRpercentVLF=output.welch.fhrv.percent_VLF;    handles.RRpercentLF=output.welch.fhrv.percent_LF;
handles.RRpercentMF=output.welch.fhrv.percent_MF;      handles.RRpercentHF=output.welch.fhrv.percent_HF;
handles.RRpeakVLF=output.welch.fhrv.peakVLF;  handles.RRpeakLF=output.welch.fhrv.peakLF; 
handles.RRpeakMF=output.welch.fhrv.peakMF;    handles.RRpeakHF=output.welch.fhrv.peakHF;
handles.RRlfhf=output.welch.fhrv.lfhf;          handles.RRratio=output.welch.fhrv.ratio;
data=[{'RR power VLF',handles.RRpowerVLF,       '     ms^2'};...
      {'RR power LF',handles.RRpowerLF,         '     ms^2'};...
      {'RR power MF',handles.RRpowerMF,         '     ms^2'};...
      {'RR power HF',handles.RRpowerHF,         '     ms^2'};...
      {'RR power total',handles.RRpowerTotal,   '     ms^2'};...
      {'RR percent VLF',handles.RRpercentVLF,     '     %'};...
      {'RR percent LF',handles.RRpercentLF,       '     %'};...
      {'RR percent MF',handles.RRpercentMF,       '     %'};...
      {'RR percent HF',handles.RRpercentHF,       '     %'};...
      {'RR peak VLF',handles.RRpeakVLF,         '       Hz'};... 
      {'RR peak LF',handles.RRpeakLF,           '       Hz'};...
      {'RR peak MF',handles.RRpeakMF,           '       Hz'};...
      {'RR peak HF',handles.RRpeakHF,           '       Hz'};...
      {'RR lfhf',handles.RRlfhf,                '       /'};...
      {'RR ratio',handles.RRratio,              '        /'}];
set(handles.uitable1,'data',data);                          % 将数据显示在表格中
set(handles.text7,'string',...
    'Fetal heart rate variability frequency domain parameters extracted successfully ！');
guidata(hObject,handles);
% --- Executes on button press in RR_nonlinear.
function RR_nonlinear_Callback(hObject,eventdata,handles)
% RR 非线性参数
set(handles.text7,'string',...
    'Fetal heart rate variability nonlinear parameters calculation ！');
data=handles.rr2(2,:);
ApEn1=ApEnFast(data,2,0.1);          % 调用函数 
SampEn1=SampEn(data,2,0.1);
fdH=FD_Higuchi(data,100);
[hurst,~]=FD_Variance(data);
C=LZC(data,2,3,0);
axes(handles.axes1);    cla reset
D=DFA(data,100,400,100,1);
axes(handles.axes2);    cla reset
PR=PRSA(data,5,10,5,1);
axes(handles.axes3);    cla reset
Poin=Poincare(data,1);

handles.RRApEn=ApEn1;           handles.RRSampEn=SampEn1;   % 便于全局调用
handles.RRFD=fdH;               handles.RRhurst=hurst;
handles.RRLZC=C;                handles.RRalpha=D.alpha;
handles.RRAAC=PR.AAC;           handles.RRADC=PR.ADC;
handles.RRPoinSD1=Poin.SD1_method2; handles.RRPoinSD2=Poin.SD2_method2;
data=[{'RR Approximate Entropy',handles.RRApEn,     '      /'};...
      {'RR Sample Entropy',handles.RRSampEn,        '      /'};...      元胞数组：字符串+数字
      {'RR Fractal dimension',handles.RRFD,      	'      /'};...
      {'RR Hurst',handles.RRhurst,                  '      /'};...
      {'RR Complexity',handles.RRLZC,               '      /'};...
      {'RR DFA alpha',handles.RRalpha,              '      /'};...
      {'RR PRSA AAC',handles.RRAAC,                 '    ms'};...
      {'RR PRSA ADC',handles.RRADC,                 '	    ms'};...
      {'RR Poincare SD1',handles.RRPoinSD1,         '    ms'};...
      {'RR Poincare SD2',handles.RRPoinSD2,         '    ms'}];
set(handles.uitable1,'data',data);                          % 将数据显示在表格中
set(handles.text7,'string',...
    'Fetal heart rate variability nonlinear parameters extracted successfully ！');
guidata(hObject,handles);
% --- Executes on button press in Generate_Number.
function Generate_Number_Callback(hObject,eventdata,handles)
% 产生孕妇的标号
a=randperm(6); % 随机产生1到6的6个整数
number=a(1)*10^5+a(2)*10^4+a(3)*10^3+a(4)*10^2+a(5)*10^1+a(6)*10^0;
handles.number=number;
set(handles.text5,'String',handles.number); % 在text5(文本框)中显示结果
guidata(hObject,handles);

%% 第3步: save菜单 保存图像和数据
function Save_Callback(hObject,eventdata,handles)
function Save_axes1_Callback(hObject,eventdata,handles)
% if isempty(handles.axes1)
%    set(handles.text7,'string','The current interface(axes4) has no picture !'); return;
% end
axes(handles.axes1);            % 取得axes1的句柄
[filename,pathname]=uiputfile({ '*.png','figure type(*.png)'},'Save figure',...
  	fullfile('E:\workspace\CTG\CTG_data\physionet\1001-1552\',handles.name)); % 保存文件的默认名称
if isequal(filename,0)||isequal(pathname,0) % 如果用户选择“取消”，则退出
    set(handles.text7,'string','The user pressed cancel !');
    return;
else
    % 由于直接保存axes1上的图像有困难，所以保存在新建的figure中的谱图
    newFig=figure('Visible','off');             % 设置新建的figure为不可见
    newAxes=copyobj(handles.axes1,newFig);      % 将axes1中的图复制到新建的figure中
    set(newAxes,'Units','default','Position','default');    % 设置图显示的位置
    fpath=fullfile(pathname,filename);
    % 此行程序保存的图像没有标题和坐标
    % [f,p]=uiputfile({'*.png'},'保存图像');
    % str=strcat(p,f);
    % pix=getframe(handles.axes1);
    % saveas(gcf,str,'png');
    % imwrite(pix.cdata,str,'png');
    f=getframe(newFig);
    f=frame2im(f);
    imwrite(f,fpath);
    % close(gcf);
    str1=(['The user selected: ',fpath]);
    str2=('Picture saved successfully !');
    set(handles.text7,'string',{str1;str2});
end
function Save_axes2_Callback(hObject,eventdata,handles)
axes(handles.axes2);  
[filename,pathname]=uiputfile({ '*.png','figure type(*.png)'},'Save figure',...
  	fullfile('E:\workspace\CTG\CTG_data\physionet\1001-1552\',handles.name));
if isequal(filename,0)||isequal(pathname,0)
    set(handles.text7,'string','The user pressed cancel !');
    return;
else
    newFig=figure('Visible','off');     newAxes=copyobj(handles.axes1,newFig);
    set(newAxes,'Units','default','Position','default'); 
    f=getframe(newFig);         f=frame2im(f);
    fpath=fullfile(pathname,filename);
    imwrite(f,fpath);
    str1=(['The user selected: ',fpath]);       str2=('Picture saved successfully !');
    set(handles.text7,'string',{str1;str2});
end
function Save_axes3_Callback(hObject,eventdata,handles)
axes(handles.axes3);  
[filename,pathname]=uiputfile({ '*.png','figure type(*.png)'},'Save figure',...
  	fullfile('E:\workspace\CTG\CTG_data\physionet\1001-1552\',handles.name));
if isequal(filename,0)||isequal(pathname,0)
    set(handles.text7,'string','The user pressed cancel !');
    return;
else
    newFig=figure('Visible','off');     newAxes=copyobj(handles.axes1,newFig);
    set(newAxes,'Units','default','Position','default'); 
    f=getframe(newFig);         f=frame2im(f);
    fpath=fullfile(pathname,filename);
    imwrite(f,fpath);
    str1=(['The user selected: ',fpath]);       str2=('Picture saved successfully !');
    set(handles.text7,'string',{str1;str2});
end
function Save_axes4_Callback(hObject,eventdata,handles)
axes(handles.axes4);  
[filename,pathname]=uiputfile({ '*.png','figure type(*.png)'},'Save figure',...
  	fullfile('E:\workspace\CTG\CTG_data\physionet\1001-1552\',handles.name));
if isequal(filename,0)||isequal(pathname,0)
    set(handles.text7,'string','The user pressed cancel !');
    return;
else
    newFig=figure('Visible','off');     newAxes=copyobj(handles.axes1,newFig);
    set(newAxes,'Units','default','Position','default'); 
    f=getframe(newFig);         f=frame2im(f);
    fpath=fullfile(pathname,filename);
    imwrite(f,fpath);
    str1=(['The user selected: ',fpath]);       str2=('Picture saved successfully !');
    set(handles.text7,'string',{str1;str2});
end
function Save_data_Callback(hObject,eventdata,handles)
% 保存文件（相关参数） help uiputfile
[filename,pathname]=uiputfile({'*.txt';'*.xls'},'Save file',...
    fullfile('E:\workspace\CTG\CTG_data\physionet\1001-1552\',handles.name));
if isequal(filename,0) || isequal(pathname,0)
    str1=('User pressed Cancel.');      str2=('Info not saved !');
    set(handles.text7,'string',{str1;str2});
    return
else
    str1=(['User selected: ',fullfile(pathname,filename)]);
    %% 0、注释
    part0=      ('        Feature             Unit              Value ');
    %% 1、孕妇的相关信息
    MA=get(handles.edit1,'String');
    GA=get(handles.edit2,'String');
    Number=get(handles.text5,'String');
    part1=      ('                    孕妇的相关信息             ');
    Number=     ['Pregnant Number                    :          ',num2str(Number)];
    MA=         ['Maternal age                (year) :          ',num2str(MA)];
    GA=         ['Gestational age             (week) :          ',num2str(GA)];
    %% 2、FHR信号的形态学参数
    part2=      ('                    胎心率信号：形态学特征      ');
    BL_mean=    ['Baseline mean               (bpm)  :          ',num2str(handles.meanBL)];
    BL_sd=      ['Baseline sd                 (bpm)  :          ',num2str(handles.sdBL)];
    BL_min=     ['Baseline min                (bpm)  :          ',num2str(handles.minBL)];
    BL_max=     ['Baseline max                (bpm)  :          ',num2str(handles.maxBL)];
    ACC=        ['Acceleration                       :          ',num2str(handles.acc)];
    DEC_mild=   ['Deceleration_mild                  :          ',num2str(handles.dec_mild)];
    DEC_prolong=['Deceleration_prolong               :          ',num2str(handles.dec_prolong)];
    DEC_severe= ['Deceleration_severe                :          ',num2str(handles.dec_severe)];
    %% 3、FHR信号的时域参数
    part3=      ('                    胎心率信号：时域特征                ');
    FHR_Loss=   ['Signal loss                   (%)  :          ',num2str(handles.loss)];
    FHR_Quality=['Signal quality                (%)  :          ',num2str(handles.quality)];
    FHR_mean=   ['mean                        (bpm)  :          ',num2str(handles.FHRmean)];
    FHR_std=    ['sd                          (bpm)  :          ',num2str(handles.FHRstd)];
    FHR_min=    ['min                         (bpm)  :          ',num2str(handles.FHRmin)];
    FHR_max=    ['max                         (bpm)  :          ',num2str(handles.FHRmax)];
    FHR_STV=    ['STV                         (bpm)  :          ',num2str(handles.FHRSTV)];
    FHR_II=     ['II                                 :          ',num2str(handles.FHRII)];
    FHR_LTI=    ['LTI                         (bpm)  :          ',num2str(handles.FHRLTI)];
    FHR_delta=  ['delta                       (bpm)  :          ',num2str(handles.FHRdelta)];
FHR_delta_total=['delta_total                 (bpm)  :          ',num2str(handles.FHRdelta_total)];
FHR_OSC_type=   ['OSC_type                           :          ',num2str(handles.FHROSC_type)];
    %% 4、RR信号的时域参数
    part4=      ('                    胎心率变异信号：时域特征    ');
    RR_mean=    ['mean                         (ms)  :          ',num2str(handles.RRmean)];
	RR_min=     ['min                          (ms)  :          ',num2str(handles.RRmin)];
    RR_max=     ['max                          (ms)  :          ',num2str(handles.RRmax)];
    RR_median=  ['median                       (ms)  :          ',num2str(handles.RRmedian)];
    RR_SDNN=    ['SDNN                         (ms)  :          ',num2str(handles.RRSDNN)];
    RR_SDANN=   ['SDANN                        (ms)  :          ',num2str(handles.RRSDANN)];
    RR_SDNNi=   ['SDNNi                        (ms)  :          ',num2str(handles.RRSDNNi)];
    RR_RMSSD=   ['RMSSD                        (ms)  :          ',num2str(handles.RRRMSSD)];
    RR_NNx=     ['NNx                      (sample)  :          ',num2str(handles.RRNNx)];
    RR_pNNx=    ['pNNx                          (%)  :          ',num2str(handles.RRpNNx)];
    RR_STV=     ['STV                          (ms)  :          ',num2str(handles.RRSTV)];
    RR_II=      ['II                                 :          ',num2str(handles.RRII)];
    RR_LTI=     ['LTI                          (ms)  :          ',num2str(handles.RRLTI)];
    RR_delta=   ['delta                        (ms)  :          ',num2str(handles.RRdelta)];
RR_delta_total= ['delta_total                  (ms)  :          ',num2str(handles.RRdelta_total)];
    RR_FHRVTi=  ['FHRVTi                             :          ',num2str(handles.RRFHRVTi)];
    RR_TINN=    ['TINN                               :          ',num2str(handles.RRTINN)];
    %% 5、RR信号的频域参数
    part5=      ('                    胎心率变异信号：频域特征              ');
    RR_powerVLF=['powerVLF                   (ms^2)  :          ',num2str(handles.RRpowerVLF)];
    RR_powerLF= ['powerLF                    (ms^2)  :          ',num2str(handles.RRpowerLF)];
    RR_powerMF= ['powerMF                    (ms^2)  :          ',num2str(handles.RRpowerMF)];
    RR_powerHF= ['powerHF                    (ms^2)  :          ',num2str(handles.RRpowerHF)];
RR_powerTotal=  ['powerTotal                 (ms^2)  :          ',num2str(handles.RRpowerTotal)];
RR_percentVLF=  ['powerVLF                      (%)  :          ',num2str(handles.RRpercentVLF)];
RR_percentLF=   ['powerLF                       (%)  :          ',num2str(handles.RRpercentLF)];
RR_percentMF=   ['powerMF                       (%)  :          ',num2str(handles.RRpercentMF)];
RR_percentHF=   ['powerHF                       (%)  :          ',num2str(handles.RRpercentHF)];
   RR_peakVLF=  ['peakVLF                      (Hz)  :          ',num2str(handles.RRpeakVLF)];
   RR_peakLF=   ['peakLF                       (Hz)  :          ',num2str(handles.RRpeakLF)];
   RR_peakMF=   ['peakMF                       (Hz)  :          ',num2str(handles.RRpeakMF)];
   RR_peakHF=   ['peakHF                       (Hz)  :          ',num2str(handles.RRpeakHF)];
    RR_lfhf=    ['lfhf                          (%)  :          ',num2str(handles.RRlfhf)];
    RR_ratio=   ['ratio                         (%)  :          ',num2str(handles.RRratio)];
    %% 6、RR的非线性参数
    part6=      ('                    胎心率变异信号：非线性特征  ');
    RR_ApEn=    ['Approximate Entropy                :          ',num2str(handles.RRApEn)];
    RR_SampEn=  ['Sample Entropy                     :          ',num2str(handles.RRSampEn)];
    RR_FD=      ['Fractal dimension                  :          ',num2str(handles.RRFD)];
    RR_hurst=   ['Hurst                              :          ',num2str(handles.RRhurst)];
    RR_LZC=     ['Lempl Ziv Complexity               :          ',num2str(handles.RRLZC)];
    RR_alpha=   ['DFA alpha                          :          ',num2str(handles.RRalpha)];
    RR_AAC=     ['PRSA AAC                     (ms)  :          ',num2str(handles.RRAAC)];
    RR_ADC=     ['PRSA ADC                     (ms)  :          ',num2str(handles.RRADC)];
	RR_SD1=     ['Poincare SD1                 (ms)  :          ',num2str(handles.RRPoinSD1)];
    RR_SD2=     ['Poincare SD2                 (ms)  :          ',num2str(handles.RRPoinSD2)];
    %% 合成
    data={ part0;part1;Number;MA;GA;part2;BL_mean;BL_sd;BL_min;BL_max;ACC;DEC_mild;DEC_prolong;DEC_severe;...
        part3;FHR_Loss;FHR_Quality;FHR_mean;FHR_std;FHR_min;FHR_max;FHR_STV;FHR_II;FHR_LTI;FHR_delta;FHR_delta_total;FHR_OSC_type;...
        part4;RR_mean;RR_min;RR_max;RR_median;RR_SDNN;RR_SDANN;RR_SDNNi;RR_RMSSD;RR_NNx;RR_pNNx;RR_STV;RR_II;...
                RR_LTI;RR_delta;RR_delta_total;RR_FHRVTi;RR_TINN;...
        part5;RR_powerVLF;RR_powerLF;RR_powerMF;RR_powerHF;RR_powerTotal;RR_percentVLF;RR_percentLF;RR_percentMF;RR_percentHF;...
                RR_peakVLF;RR_peakLF;RR_peakMF;RR_peakHF;RR_lfhf;RR_ratio;...
        part6;RR_ApEn;RR_SampEn;RR_FD;RR_hurst;RR_LZC;RR_alpha;RR_AAC;RR_ADC;RR_SD1;RR_SD2;};    % 不同长度的字符串连接，用元胞矩阵
    %% 保存文件路径
    str=[pathname filename];
    fid=fopen(char(str), 'w');
    fwrite(fid,'','integer*4');
    for i=1:length(data)    % 自动换行
        fprintf(fid,'%s\n',data{i});
    end
    fclose(fid);
    % disp('Data saved successful!')
    str2=('Data saved successfully!');
    set(handles.text7,'string',{str1;str2}); 
end
guidata(hObject,handles);
