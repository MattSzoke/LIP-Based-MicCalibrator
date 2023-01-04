%%
clc
clear all
close all

%% Use LaTeX in figures
set(0, 'DefaultTextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Constants
ValidChannels = [1:3]; % Ch4 is reference mic

%% Load data
load('Calibration_Reduced_1.mat')
% note to self: Delta_dB_R = 20*log10(ArrivalTimes./ArrivalTimes(RefMicIdx));

NominalMicSens = importdata('MicSensitivities.txt', ' ', 1);
NominalMicSens = (NominalMicSens.data(:,2)./1000)';% V/Pa;
RefSens = NominalMicSens(4); % sensitivity of reference microphone, V/Pa, found using pistonphone
RefSens_dB = 20*log10(RefSens); % sensitivity of reference microphone, dB

%% We are going to use a mask to taper off the sensitivity to 0 dB at low frequencies
% freq range where we will taper off the sensitivity linearly to zero
f0 = 1000; f1 = 10000; 
mymask = @(f) interp1([0 f0 f1 1000*f1], [0 0 1 1], f);

%Lets make a plot to visualize this
figure(1), clf, hold on
plot(MicCal{1}.f, mymask(MicCal{1}.f))
set(gca,'xscale','log')
xlabel('Frequency, Hz')
ylabel('Weight')
grid on

%% load the manufactuer provided response of our ref mic so that we can correct for it
load('SN407760.mat')
% interpolate ref mic response to the frequencies of the LIP calibration
BareRef = interp1(SN407760.f, SN407760.Resp,  MicCal{1}.f)';
BareRef(isnan(BareRef)) = 0;
BareRefC = 10.^(BareRef/20);

figure(2), clf, hold on
semilogx(SN407760.f, SN407760.Resp)
set(gca,'xscale','log')
xlabel('Frequency, Hz')
ylabel('Response, dB')
grid on
title('Ref. mic. nominal response')
xlim([1e3 1e5])
%% Calculate mic sensitivities
RefMicNominal = BareRefC;

for IDX = ValidChannels
    LIP = MicCal{IDX}; % save a few characters by making a new variable
    iii = 1000<LIP.f & LIP.f<2000; % indices of freqs within 1 kHz and 2 kHz
    
    % What is the difference in response at low frequencies
    % lets correct for propagation distance differences too
    delta_dB(IDX) = mean(20*log10(abs(LIP.Txy(iii)))) + Delta_dB_R(IDX);
    delta_dB_Resp(IDX) = mean(20*log10(abs(LIP.Txy(iii))));
    % Calculate microphone sensitivity based on response difference wrt.
    % reference mic
    Sens(IDX)     = RefSens*10^(delta_dB(IDX)/20);
    Sens_dB(IDX)  = 20*log10(Sens(IDX));
    
    % Calculate transfer function
    Txy_0(:,IDX)  = LIP.Txy/(10^(delta_dB_Resp(IDX)/20))  ;
    Txy_0(:,IDX)  = Txy_0(:,IDX).*RefMicNominal ;
    Txy_dB(:,IDX) = 20*log10(abs(Txy_0(:,IDX)));
end


%% Re-structure calibration data and apply mask
Fnew = logspace(2,5,1000);
for IDX = ValidChannels
    %IDX
    LIP = MicCal{IDX};
    Txy_ip(IDX,:) = interp1(LIP.f,Txy_0(:,IDX),Fnew).^(mymask(Fnew));
    Ms(IDX,:) = interp1(LIP.f,Txy_dB(:,IDX),Fnew).*(mymask(Fnew));
    Ps(IDX,:) = interp1(LIP.f,LIP.Phi,Fnew).*(mymask(Fnew));
end
Fs = Fnew;

%% Lets get rid of data at frequencies higher than 100 kHz, it is invalid data
[~,idx] = min(abs(Fs-100000));
Ms(idx:end,:) = nan;

%% Plot all mic responses
myXLim = [1000 100000];
figure(3), clf,
subplot(2,1,1), hold on
set(gca,'ColorOrderIndex',1)
plot(Fs,Ms(ValidChannels,:),'-', 'linewidth', 1.5)
plot(Fs,mean(Ms(ValidChannels,:),1),'k-', 'linewidth', 1.5)
set(gca,'xscale','log')
xlim(myXLim)
ylim([-15 2])
grid on
xlabel('Frequency, Hz')
ylabel('Response, dB')

subplot(2,1,2), hold on
set(gca,'ColorOrderIndex',1)
plot(Fs,Ps(ValidChannels,:),'-', 'linewidth', 1.5)
plot(Fs,mean(Ps(ValidChannels,:),1),'k-', 'linewidth', 1.5)
set(gca,'xscale','log')
xlim(myXLim)
ylim([-1 1]*30)
grid on
xlabel('Frequency, Hz')
ylabel('Phase, deg')

%% Lets compare nominal sensitivities  to calculated ones
figure(4), clf
subplot(2,1,1)
plot(ValidChannels,Sens(ValidChannels)*1000,'ro','displayname','Calculated using LIP'), hold on
plot(ValidChannels,NominalMicSens(ValidChannels)*1000,'k.','displayname','Nominal'), hold on
xlabel('Channel number')
ylabel('Mic sensitivity, mV/Pa')
grid on
ylim([1 2.2])
legend('location', 'southwest')

subplot(2,1,2)
plot(ValidChannels,20*log10(Sens(ValidChannels)),'ro','displayname','Calculated using LIP'), hold on
plot(ValidChannels,20*log10(NominalMicSens(ValidChannels)),'k.','displayname','Nominal')% Nominal
xlabel('Channel number')
ylabel('Sensitivity, dB re. 1 V')
grid on
legend('location', 'southwest')

