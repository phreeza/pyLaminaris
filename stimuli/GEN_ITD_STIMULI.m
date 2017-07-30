%% gen pure tone stim with ITD
fs = 195312.50;
dur = 200*1e-3; % 200ms sound
nSamp = round(fs*dur);
MAX = 0.8; % constrain each wave to be within [-0.8,0.8] to allow some modulation of amplitude, either increase or decrease of amplitude
% load following calibration data from 03012016
% amplitude difference between right ane left earplug, reference is the left
% left column is frequency calibrated, right column is difference in dB
ampDiff = [...  
    1000.00	0.68;
    1250.00	1.54;
    1500.00	1.78;
    1750.00	0.04;
    2000.00	-5.84;
    2250.00	-7.94;
    2500.00	-9.54;
    2750.00	-8.63;
    3000.00	-8.36;
    3250.00	-8.66;
    3500.00	-6.61;
    3750.00	-3.74;
    4000.00	-1.49;
    4250.00	-0.49;
    4500.00	-1.90;
    4750.00	-3.42;
    5000.00	-4.69;
    5250.00	-5.75;
    5500.00	-6.66;
    5750.00	-7.29;
    6000.00	-7.53;
    6250.00	-7.74;
    6500.00	-7.62;
    6750.00	-7.26;
    7000.00	-6.56;
    7250.00	-5.67;
    7500.00	-4.71;
    7750.00	-3.81;
    8000.00	-3.07]; 

% amplitude difference between right ane left earplug, reference is the left
% left column is frequency calibrated, right column is time difference in us
phaseDiff = [...
    1000.00	19.25;
    1250.00	-6.16;
    1500.00	-39.83;
    1750.00	-75.77;
    2000.00	-78.39;
    2250.00	-40.75;
    2500.00	-17.82;
    2750.00	-0.86;
    3000.00	4.21;
    3250.00	15.43;
    3500.00	26.74;
    3750.00	25.57;
    4000.00	17.98;
    4250.00	4.51;
    4500.00	-4.28;
    4750.00	-6.18;
    5000.00	-5.57;
    5250.00	-4.22;
    5500.00	-1.81;
    5750.00	1.15;
    6000.00	3.72;
    6250.00	6.06;
    6500.00	8.55;
    6750.00	10.97;
    7000.00	12.71;
    7250.00	13.47;
    7500.00	13.40;
    7750.00	12.62;
    8000.00	11.39];

path2save = 'C:\Documents and Settings\i suck\My Documents\Google Drive\Stim_cheetah';
% path2save = 'F:\Google Drive\Stim_cheetah';

% save stimulus to separate folders according to frequency
freq = 3000:250:7000;
nFreq = length(freq);

% set ITD range
ITD_array = (-300:30:300)*1e-3; % ITD in millisecond
nITD = length(ITD_array);

% generate linear ramp
t = (0:nSamp-1)/fs;
ramp = 5e-3; % rise/fall time in seconds
rmp = ones(1,nSamp);
sel = (t<=ramp);
srmp = linspace(0,1,sum(sel));
rmp(sel) = srmp;
sel = (t>=dur-ramp);
srmp = linspace(0,1,sum(sel));
rmp(sel) = srmp(end:-1:1);
%% go through each frequency and generate ITD stimuli
for iFreq = 1:nFreq
    
    folder_name = sprintf('stim_%dHz_withITD',freq(iFreq));
    if ~isdir(fullfile(path2save,folder_name))
        mkdir(fullfile(path2save,folder_name));
    end
    for iITD = 1:nITD
        
        ind = find(phaseDiff(:,1) == freq(iFreq));
        % true_ITD is the actual ITD
        true_ITD = ITD_array(iITD);
        % ITD is the ITD we want with correction for speaker delay
        ITD = - phaseDiff(ind,2)*1e-3 + true_ITD;
        % save additional gain on channel 2/ right channel
        gain2 = ampDiff(ind,2);
        % generate actual sound wave
        wave = rmp .* sin(2*pi*freq(iFreq)*t); % note the maximum of abs(wave) should be no more than 1 in order not to exceed RX6 output limit
        wave = wave / max(abs(wave)) * MAX;
        % save info to matfile
        fn = sprintf('stim_%dHz_%dus.mat',freq(iFreq),true_ITD*1e3);
        param = [];
        FREQ = freq(iFreq);
        save(fullfile(path2save,folder_name,fn),'wave','param','ITD','gain2','true_ITD','FREQ','-v7.3');
        
    end
    
end

%% generate click train 
FREQ = 5000; % 5kHz, one cycle
t1 = 0:1/fs:1/FREQ;
nSamp1Cyc = length(t1);

wave = [ sin(2*pi*t1*FREQ) zeros(1,nSamp - nSamp1Cyc)]; % want to pad zeros up to 200ms
wave = wave / max(abs(wave)) * MAX;

ind = find(phaseDiff(:,1) == FREQ);

folder_name = 'stim_click_withITD';
if ~isdir(fullfile(path2save,folder_name))
    mkdir(fullfile(path2save,folder_name))
end
for iITD = 1:nITD
        % true_ITD is the actual ITD
        true_ITD = ITD_array(iITD);
        % ITD is the ITD we want with correction for speaker delay
        ITD = - phaseDiff(ind,2)*1e-3 + true_ITD;
        % save additional gain on channel 2/ right channel
        gain2 = ampDiff(ind,2);
        % save info to matfile
        fn = sprintf('stim_click_%dus.mat',true_ITD*1e3);
        param = [];
        save(fullfile(path2save,folder_name,fn),'wave','param','ITD','gain2','true_ITD','FREQ','-v7.3');

end


%% generate noise stimuli
nNoise = 10;

if ~isdir(fullfile(path2save,'stim_noise_withITD'))
    mkdir(fullfile(path2save,'stim_noise_withITD'))
end

for iNoise = 1:nNoise
    % for now we do not correct for time delay or amplitude distortion in noise
    wave = rmp .* randn(1,nSamp);
    wave = wave / max(abs(wave)) * MAX;  % constrain abs(wave) to be within [-MAX,MAX]
    for iITD = 1:nITD
        true_ITD = ITD_array(iITD);
        ITD = ITD_array(iITD);
        fn = sprintf('stim_noise%d_%dus.mat',iNoise,ITD*1e3);
        gain2 = 0;
        param = [];
        save(fullfile(path2save,'stim_noise_withITD',fn),'wave','param','ITD','gain2','true_ITD','-v7.3');
    end
end
%% use code from Thomas to generate chirp signal, uncalibrated
nChirp = 1;

if ~isdir(fullfile(path2save,'stim_chirp_withITD'))
    mkdir(fullfile(path2save,'stim_chirp_withITD'))
end
c = 0.205; % input argument for dauChirpOAE
a = -0.706; % input argument for dauChirpOAE
for iNoise = 1:nChirp
    % for now we do not correct for time delay or amplitude distortion in
    % chirp
    wave = zeros(1,nSamp);
    chirp = dauChirpOAE(fs, c, a, [500,10000]);
    sz = length(chirp);
    if size(chirp,2) == 1
        wave(1:sz) = chirp';
    else
        wave(1:sz) = chirp;
    end
    
    wave = wave / max(abs(wave)) * MAX;  % constrain abs(wave) to be within [-MAX,MAX]
    for iITD = 1:nITD
        true_ITD = ITD_array(iITD);
        ITD = ITD_array(iITD);
        fn = sprintf('stim_chirp%d_%dus.mat',iNoise,ITD*1e3);
        gain2 = 0;
        param = [];
        save(fullfile(path2save,'stim_chirp_withITD',fn),'wave','param','ITD','gain2','true_ITD','-v7.3');
    end
end