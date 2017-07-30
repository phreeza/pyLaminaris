%% gen pure tone stim with ITD
fs = 195312.50;
dur = 200*1e-3; % 200ms sound
nSamp = round(fs*dur);
MAX = 0.8; % constrain each wave to be within [-0.8,0.8] to allow some modulation of amplitude, either increase or decrease of amplitude
% load following calibration data from 03012016
% amplitude difference between right ane left earplug, reference is the left
% left column is frequency calibrated, right column is difference in dB
ampDiff = [...
    1000.00	-7.33
    1100.00	-6.17
    1200.00	-4.80
    1300.00	-3.12
    1400.00	-1.56
    1500.00	-0.89
    1600.00	-0.42
    1700.00	-0.87
    1800.00	-2.02
    1900.00	-3.33
    2000.00	-5.15
    2100.00	-5.86
    2200.00	-5.24
    2300.00	-4.48
    2400.00	-3.97
    2500.00	-3.45
    2600.00	-2.58
    2700.00	-1.34
    2800.00	-0.39
    2900.00	0.44
    3000.00	1.15
    3100.00	1.42
    3200.00	1.20
    3300.00	0.58
    3400.00	-0.21
    3500.00	-0.82
    3600.00	-1.13
    3700.00	-1.26
    3800.00	-1.48
    3900.00	-1.60
    4000.00	-1.71
    4100.00	-1.66
    4200.00	-1.47
    4300.00	-1.36
    4400.00	-1.36
    4500.00	-1.32
    4600.00	-1.33
    4700.00	-1.39
    4800.00	-1.67
    4900.00	-2.09
    5000.00	-2.67
    5100.00	-3.26
    5200.00	-3.84
    5300.00	-4.31
    5400.00	-4.67
    5500.00	-4.93
    5600.00	-5.02
    5700.00	-4.97
    5800.00	-4.89
    5900.00	-4.81
    6000.00	-4.78
    6100.00	-4.84
    6200.00	-5.03
    6300.00	-5.26
    6400.00	-5.39
    6500.00	-5.45
    6600.00	-5.37
    6700.00	-5.32
    6800.00	-5.33
    6900.00	-5.20
    7000.00	-5.01
    7100.00	-4.93
    7200.00	-4.72
    7300.00	-4.62
    7400.00	-4.48
    7500.00	-4.24
    7600.00	-3.97
    7700.00	-3.77
    7800.00	-3.51
    7900.00	-3.24
    8000.00	-3.07];

% amplitude difference between right ane left earplug, reference is the left
% left column is frequency calibrated, right column is time difference in us
phaseDiff = [...
    1000.00	110.25
1100.00	111.84
1200.00	106.82
1300.00	97.02
1400.00	77.01
1500.00	55.18
1600.00	36.92
1700.00	17.54
1800.00	6.20
1900.00	1.35
2000.00	6.34
2100.00	21.75
2200.00	32.88
2300.00	35.46
2400.00	36.44
2500.00	36.12
2600.00	37.87
2700.00	36.50
2800.00	32.38
2900.00	27.08
3000.00	20.64
3100.00	13.23
3200.00	6.11
3300.00	1.35
3400.00	-0.64
3500.00	-0.36
3600.00	0.67
3700.00	0.89
3800.00	0.80
3900.00	0.95
4000.00	1.61
4100.00	1.88
4200.00	1.53
4300.00	0.67
4400.00	-0.05
4500.00	-0.81
4600.00	-1.98
4700.00	-3.44
4800.00	-4.98
4900.00	-6.26
5000.00	-7.04
5100.00	-7.02
5200.00	-6.43
5300.00	-5.36
5400.00	-4.08
5500.00	-2.63
5600.00	-1.28
5700.00	-0.16
5800.00	0.49
5900.00	0.90
6000.00	1.06
6100.00	0.95
6200.00	0.95
6300.00	1.38
6400.00	2.23
6500.00	3.12
6600.00	3.89
6700.00	4.41
6800.00	5.03
6900.00	5.73
7000.00	6.20
7100.00	6.47
7200.00	6.73
7300.00	7.08
7400.00	7.37
7500.00	7.52
7600.00	7.57
7700.00	7.49
7800.00	7.36
7900.00	7.11
8000.00	6.95];




path2save = 'C:\Documents and Settings\i suck\My Documents\Google Drive\Stim_cheetah';
% path2save = 'F:\Google Drive\Stim_cheetah';

% save stimulus to separate folders according to frequency
freq = 3000:100:8000;
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