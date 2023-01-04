%% clear workspace
clc, clear all, close all

%% Use LaTeX within figures
set(0, 'DefaultTextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Some constants and user inputs
SparkFreq = 5;            % LIP formulations per second
nCh = 4;                  % NumberOfChannels
p0 = 20e-6; p0sq = p0*p0; % ref pressure
Fs = 748800;              % sampling frequency

%% DATA STRUCTURE
RefMicIdx = 4; % channel index of the reference microphone
PDid = 5;      % channel index of photodetector signal
ValidChannels = [1:4]; % these are channel indices where we have mics

%% Load data
load('SampleData_v2.mat')
MicSens = importdata('MicSensitivities.txt', ' ', 1);
MicSens = (MicSens.data(:,2)./1000);% mic sensitivities, V/Pa;

%% Find peaks in photodetector signal
nStamps = size(data,1);
MinPeakDist = Fs*(1/SparkFreq)/4; % Minimum peak-to-peak disance, # of timestamps
[PDpeaks, PDpIdx] = findpeaks(data(:,PDid),'MinPeakDistance',MinPeakDist,'MinPeakProminence',max(data(:,PDid))*0.7);

%% Plot raw data
figure(1); clf, hold on
plot(timeStamps, data(:,PDid),'linewidth', 1.2)
plot(timeStamps, data(:,1),'linewidth', 1.1)
plot(timeStamps(PDpIdx), data(PDpIdx,PDid),'ro')
xlim([0.0 0.5])
xlabel('Time, s')
ylabel('Voltage')
legend('Photodetector', 'Microphone')
grid on

%% Decide which block has plasma formation (when LIP generated sound)
RefCh = 1  % we will use this channel to decide if there was LIP formation in a given block

PkTh = 0.4; % Peak Threshold ratio wrt. max pressure value in time series
tblk = [0:MinPeakDist-1]/Fs; % timestamps
LIPform = 0.*PDpIdx; % logical array storing LIP generation blocks
RefChData = data(:,RefCh); % using a ref channel to detect LIP formation

for j = 1:length(PDpeaks)-1
    x = RefChData(PDpIdx(j):PDpIdx(j)+MinPeakDist-1);
    [mypk, mypidx] = findpeaks(x,'MinPeakDistance',MinPeakDist/2,'MinPeakProminence',max(RefChData)*PkTh);
    
    if ~isempty(mypk)
        LIPform(j) = 1;
        figure(2), clf, hold on, plot(tblk,x), plot(tblk(mypidx),x(mypidx),'ro'), xlabel('Time, s'), ylabel('Voltage, V')
    else
        %figure(2), hold on, plot(tblk,x), xlabel('Time, s'), ylabel('Pressure, Pa')
        disp('No LIP')
    end
end
clear RefChData % we dont need this data anymore
nLIP = sum(LIPform); % number of LIP formations in our dataset

%% CREATE BLOCKS
ta_Percent = 50; % percentage of pressure value wrt peak where we define the arrival time

for i = ValidChannels
    disp(['Block forulation for channel ', num2str(i)])
    k = 1;
    for j = 1:length(PDpeaks)-1 % throw out last LIP block in case it is too short for gating
        if LIPform(j)
            x = data(PDpIdx(j):PDpIdx(j)+MinPeakDist-1,i);
            [mypk, mypidx] = findpeaks(x,'MinPeakDistance',MinPeakDist/2,'MinPeakProminence',max(x)*PkTh);
            [maxi,maxi_idx] = max(x);
            %figure(31),clf,plot(tblk(1:maxi_idx),x(1:maxi_idx)),pause
            [~, ind] = unique(x(1:maxi_idx));
            %figure(32),clf,plot(tblk(1:maxi_idx),x(1:maxi_idx),'-rx'),hold on,plot(tblk(ind),x(ind),'k.'),pause(0.1)
            try
                ta = interp1(x(ind),tblk(ind)',max(x)*ta_Percent/100);
            catch
                disp('problem'), ta=nan;
            end
            
            Ch{i}.Block(k).x = x; % time series data
            Ch{i}.Block(k).pkIdx = mypidx; % peak location
            Ch{i}.Block(k).pk = mypk; % peak value
            Ch{i}.Block(k).ta = ta; % arrival time
            
            [~,min_idx] = min(x);
            Ch{i}.Block(k).tmin = min_idx/Fs; % time instant of min value
            
            % figure(33), hold on, plot(tblk,x), plot(tblk(mypidx),x(mypidx),'ro'), %xlim([0 0.01]), pause(0.01)
            figure(34), hold on, plot(tblk,x), plot(ta,max(x)*ta_Percent/100,'rx'), %xlim([0 0.01]), pause(0.01)
            k = k+1;
        end
        
    end
    clf
    
end

%% Calculate statistics of arrival times
for k = ValidChannels
    ArrivalTimes(k) = mean(vertcat(Ch{k}.Block(:).ta),'omitnan');
    ArrivalTimes_std(k)= std(vertcat(Ch{k}.Block(:).ta),'omitnan');
end
ArrivalTimes = ArrivalTimes';
ArrivalTimes_std = ArrivalTimes_std';
ArrivalTimes_std(ArrivalTimes_std==0)=nan;
ArrivalTimes(ArrivalTimes==0)=nan;

% distance correction wrt. ref mic location
Delta_dB_R = 20*log10(ArrivalTimes./ArrivalTimes(RefMicIdx));

%%  We dont need the data matrix anymore
clear data % free up workspace

%% Two-sided, Gated  Fourier transform data
% CHOOSE GATE LENGTH... this depends on the DAQ setup
GateL = round(32)*ones(nCh,1);  
GateR = round(48)*ones(nCh,1); 

LenGated = int32(1/double(SparkFreq) * Fs / 2); % Gate length to increase frequency resolution

for i = ValidChannels
    disp(['Calculating FFT for channel ', num2str(i)])
    
    for k = 1:nLIP
        n0 = Ch{i}.Block(k).pkIdx; % index of peak location
        x_ = Ch{i}.Block(k).x(n0-GateL(i):n0+GateR(i))';
        t_ = [0:length(x_)-1]/Fs;
        
        % if high-pass filtering is needed,
        % this is a good time to do:
        %x_ = highpass(x_,1000,Fs);
        
        % x_ = x_./max(x_); % Normalization, depending on analysis goals...
        
        % Two sided FFT
        [ f_Gated, myFFT, f_Gated_half, myFFT_half ] = fx_LIP_FFT(x_,Fs,LenGated);
        
        % Store data
        % Two sided FFT result
        Ch{i}.Block(k).FFT = myFFT;
        Ch{i}.f = f_Gated;
        % One sided FFT result
        Ch{i}.Block(k).FFT_half = myFFT_half;
        Ch{i}.f_half = f_Gated_half;
        % Store gated data too
        Ch{i}.Block(k).xG = x_;
        Ch{i}.tG = t_;
        
        %figure(40), subplot(2,1,1), plot(t_*1000,x_), xlabel('t, ms'), ylabel('$p$'), grid on,hold on,
        %subplot(2,1,2), semilogx(f_Gated_half,20*log10(myFFT_half)), hold on, grid on, xlim([1000 100000])
        
    end
    %figure(40), clf
end


%% Plot normalized gated pressure signals and normalized FFT results
PlotGatedSignal = 1 % set to 0 if you do not want to plot
if PlotGatedSignal
    for i = ValidChannels
        % Constructing matrixes from structs
        dtG = GateL(i)/Fs; % timespan of gated signal
        Ts = vertcat(Ch{i}.Block(:).tmin)-vertcat(Ch{i}.Block(:).ta); % Time period
        tas = vertcat(Ch{i}.Block(:).ta); % arrival times
        pkIdxs = vertcat(Ch{i}.Block(:).pkIdx); % indices of peak pressure
        xGs = vertcat(Ch{i}.Block(:).xG); % gated pressure data
        pks = vertcat(Ch{i}.Block(:).pk); % peak pressure values
        maxPk = max(pks); minPk = min(pks); % max(p_peak), min(p_peak), for colormap
        % color the lines according to the amplitude of peak pressure...
        myLvL = @(peak) interp1(linspace(minPk, maxPk, 100), parula(100), peak);
        
        figure(50+i), clf, hold on
        for k = 1:2:nLIP % plot every other only...
            x_ = Ch{i}.Block(k).xG;
            [xMax,xMaxi] = max(x_);
            % normalization: pressure
            x_ = x_/xMax;
            
            % normalization: time
            ti = Ch{i}.tG - dtG - (tas(k) - pkIdxs(k)/Fs);
            ti = ti/Ts(k);
            
            subplot(2,1,1), hold on
            plot(ti,x_, 'color', myLvL(pks(k))), hold on,
            xlabel('$ (t-t_a) / (t_{min} - t_a)$')
            ylabel('$p_G/$max$(p_G)$')
            xlim([-1 3])
            ylim([-1 1]*1.2)
            title(sprintf('Microphone \\#%d',i))
            grid on
            
            subplot(2,1,2), hold on
            toPlot = (Ch{i}.Block(k).FFT_half)/(pks(k));
            semilogx(Ch{i}.f_half*Ts(k), 20*log10((toPlot)),'-','linewidth',1,'color', myLvL(pks(k)))
            set(gca,'xscale','log')
            xlim([1e3 1e5]*mean(Ts))
            xlabel('$f (t_{min} - t_a)$')
            ylabel('$FFT(p_G)/ $max$(p_G)$, dB')
            grid on
            ylim([-50 0])
            
            
        end
    end
end


%% Correlation plot

figure(60), clf, hold on
for i = ValidChannels
    Ts = vertcat(Ch{i}.Block(:).tmin)-vertcat(Ch{i}.Block(:).ta);
    pks = vertcat(Ch{i}.Block(:).pk)./MicSens(i);
    
    plot(Ts*1000, pks, '.', 'displayname', sprintf('Mic. \\# %d', i))
    xlabel('$(t_{min} - t_a)$, ms')
    ylabel('max$(p_G)$, Pa')
end
legend('location','northwest')
grid on

%% Plot FFT results from all channels for a given LIP firing
PlotFFT = 1

if PlotFFT
    k = 1; % LIP firing index
    MyColors = parula(length(ValidChannels));
    for i = ValidChannels
        figure(70),% clf
        hold on
        if i == RefMicIdx, Corr = 0; else Corr = 0; end
        toPlot = (Ch{i}.Block(k).FFT_half)/(MicSens(i));
        semilogx(Ch{i}.f_half, 20*log10(toPlot),'-','linewidth',2,...
            'displayname', ['Ch ', num2str(i)]), hold on
        set(gca,'xscale','log')
        xlim([1e3 1e5])
        xlabel('Frequency, Hz')
        ylabel('FFT re. 20 $\mu$Pa, dB')
        
    end
    grid on
    legend('location','best')
end

%% CALCULATE TRANSFER FUNCTION(s)
RefMicIdx % = 1

for i = ValidChannels
    disp(['Transfer function calculation for channel ', num2str(i)])
    
    for k = 1:nLIP
        Ch{i}.Block(k).Txy      = (Ch{i}.Block(k).FFT./Ch{RefMicIdx}.Block(k).FFT)';
        Ch{i}.Block(k).Txy_half =  Ch{i}.Block(k).Txy(1:length(Ch{i}.f_half));
        try % calculate time lag difference for phase correction
            dt = (Ch{i}.Block(k).ta(1)-Ch{i}.Block(k).pkIdx(1)/Fs)-(Ch{RefMicIdx}.Block(k).ta(1)-Ch{RefMicIdx}.Block(k).pkIdx(1)/Fs);
        catch
            dt = 0;
        end
        Ch{i}.Block(k).dt       = dt;
    end
    
end

%% Calculate transfer function

for i = ValidChannels
    disp(['Response calculation for channel ', num2str(i)])
    
    % Average values
    M_Txy = horzcat(Ch{i}.Block(:).Txy_half); % construct matrix of Txy
    Txy_bar = mean(M_Txy,2, 'omitnan'); % ensemble average of transfer function
    
    my_dt = -vertcat(Ch{i}.Block(:).dt); % time lag between signals
    myPhi = exp(-1i.*my_dt*(Ch{i}.f_half)*2*pi)'; % calculate phase lag
    
    M_myPhase = -rad2deg(unwrap((angle( M_Txy.*myPhi)  ),1*pi)); % unwrap phase
    % calculate average value of phase at low freqs
    myPhase_b = round(mean(M_myPhase(501:100:1001,:),1)/pi,0)*pi;
    M_myPhase = M_myPhase - myPhase_b; % iron out phase at low frequencies
    myPhase_bar = mean(M_myPhase,2, 'omitnan'); % average phase
    
    % Build struct for export
    MicCal{i}.f = Ch{i}.f_half;
    MicCal{i}.Txy = Txy_bar;
    MicCal{i}.Phi = myPhase_bar;
    
end

%% Plot transfer function
MakePlots = 1
if MakePlots
    for i =  3 % Choose which channel to plot
        figure(80), clf
        for k = 1:nLIP
            my_dt = Ch{i}.Block(k).dt;
            
            subplot(2,1,1), hold on
            semilogx(Ch{i}.f_half, 20*log10(abs(Ch{i}.Block(k).Txy_half)),'-','linewidth',2), hold on
            set(gca,'xscale','log')
            xlim([1e3 1e5])
            xlabel('Frequency, Hz')
            ylabel('$T_{xy}$, dB', 'interpreter', 'latex')
            title(['Channel ',num2str(i)])
            %ylim([-1 1]*2.5)
            grid on
            
            myPhase = -rad2deg(unwrap((angle( Ch{i}.Block(k).Txy_half.*exp(-1i.*my_dt.*(Ch{i}.f_half')*2*pi))  ),1*pi));
            myPhase_b = round(mean(myPhase(501:100:1001))/pi,0)*pi;
            subplot(2,1,2)
            semilogx(Ch{i}.f_half, myPhase - myPhase_b,'-','linewidth',2), hold on, ylabel('Phase, deg')
            set(gca,'xscale','log')
            xlim([1e3 1e5])
            xlabel('Frequency, Hz')
            grid on
            ylim([-1 1]*20)
            
        end
        %% Plot ensemble averaged values
        % Average values
        M_Txy = horzcat(Ch{i}.Block(:).Txy_half); % construct matrix of Txy
        Txy_bar = mean(M_Txy,2, 'omitnan'); % ensemble average of transfer function
        
        
        figure(80)
        % Average  values
        subplot(2,1,1), hold on
        semilogx(Ch{i}.f_half, 20*log10((abs(Txy_bar))),'k-','linewidth',2), hold on
        set(gca,'xscale','log')
        xlim([1e3 1e5])
        xlabel('Frequency, Hz')
        ylabel('$T_{xy}$, dB')
        title(['Channel \#',num2str(i)])
        ylim([-2 1]*5)
        grid on,
        
        my_dt = -vertcat(Ch{i}.Block(:).dt);
        myPhi = exp(-1i.*my_dt*(Ch{i}.f_half)*2*pi)';
        M_myPhase = -rad2deg(unwrap((angle( M_Txy.*myPhi)  ),1*pi));
        myPhase_b = round(mean(M_myPhase(501:100:1001,:),1)/pi,0)*pi;
        M_myPhase = M_myPhase - myPhase_b;
        myPhase_bar = mean(M_myPhase,2, 'omitnan');
        
        subplot(2,1,2), hold on
        semilogx(Ch{i}.f_half, myPhase_bar ,'k-','linewidth',2), hold on, ylabel('Phase, deg')
        set(gca,'xscale','log')
        xlim([1e3 1e5])
        ylim([-1 1]*45)
        xlabel('Frequency, Hz')
        grid on
        
        
    end
end
%% Save data
save('Calibration_Reduced_1.mat','MicCal','Delta_dB_R','ArrivalTimes')

%%