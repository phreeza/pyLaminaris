%  [chirpsig] = dauChirpOAE(fs, c, a, fRange)
%                   OR
%   [chirpsig] = dauChirpOAE(fs, c, a, fRange, plotMe)
%
% This function generates the "O-chirp" from Fobel and Dau (2004).  Note
% that the equations are generalized such that users can modify the primary
% parameters, e.g. see the first footnote in Fobel and Dau (2004).
%
%
% INPUTS:
%   fs      sampling frequency (samples/sec)
%   c       constant scalar for the time vs. frequency function that
%           describes the basilar membrane group delay.
%               Fobel and Dau (2004) use c = 0.15 in the main part of the
%               paper but suggest a value of c = 0.43 in the footnote.
%   a       constant (specified as 'alpha' in the paper) which is the
%           exponent of the time vs. frequency function that describes the
%           basilar membrane group delay
%               Fobel and Dau (2004) use a = -0.5 in the main part of the
%               paper but suggest a value of a = -0.63 in the footnote.
%   fRange  can specify either [fLow] only or [fLow, fHigh].  If only
%           'fLow' is specified, fHigh = fs/2.  
%   plotMe  (OPTIONAL) if ~0, generate plots for the chirpsig.
%
%
% OUTPUTS:
%   chirpsig    output chirp signal
%
%
% USAGE:
%   [chirpsig] = dauChirpOAE(48000, 0.43, -0.63, 50);
%   [chirpsig] = dauChirpOAE(48000, 0.43, -0.63, [50 10000]);
%   [chirpsig] = dauChirpOAE(48000, 0.43, -0.63, [50 10000], 1);
%
%
% REFERENCE:
%   "Searching for the optimal stimulus eliciting auditory brainstem
%   responses in humans"
%       Fobel and Dau (2004)
%       JASA 116:2213
%   Apply general expressions for Eqs 2-7 for arbitrary values of alpha.
%
%
%   Created by Rainer Beutelmann (May 13, 2014)
%   Last updated by Darrin K. Reed (May 13, 2014)... converted from script
%   to a function and added parameter comments and checks.
%   Last updated by Darrin K. Reed (Feb 29, 2016)... added comments to aid
%   users of the function.
%--------------------------------------------------------------------------

function [chirpsig] = dauChirpOAE(fs, c, a, fRange, varargin)

%% parameter check
switch nargin
    case 4
        plotMe = 0;
    case 5
        plotMe = varargin{1};
    otherwise
        error('Invalid number of input arguments.')
end


% minimum (start) and maximum (end) frequency of the chirp
switch length(fRange)
    case 1
        f0 = fRange;
        fmax = fs/2;
    case 2
        if fRange(1) < fRange(2)
            f0 = fRange(1);
            fmax = fRange(2);
        else
            error('Must specify fRange as [fLow, fHigh].')
        end
    otherwise
        error('Invalid range of frequencies for sweep (fRange)')
end


%% preparations
% corresponding time points
t0   = c*f0^a;
tmax = t0-c*fmax^a;

% time vector
t = (0:floor(tmax*fs)).'/fs;

%% equations from the paper
% instantaneous frequency as a function of time, i.e. a generalized version
% of Eq. 3 in Fobel & Dau (2004).
%
%         / t0 - t \1/a
%  f_O =  | ------ |
%         \   c    /
%
f_O   = ((t0 - t)/c).^(1/a);
% = (c/(t0-t)).^2 for alpha = -0.5;

% instantaneous phase as a function of time, i.e. a generalized version
% of Eq. 5 in Fobel & Dau (2004).
%
%                    
%                       /       1                 1     \
%                       |       - + 1             - + 1 |
%            2 pi a c   |       a                 a     |
%  phi_O =  ----------  | / t0 \        / t0 - t \      |
%              a + 1    | | -- |      - | ------ |      |
%                       \ \ c  /        \   c    /      /
%
phi_O = 2*pi*a/((a+1)*c^(1/a))*(t0^(1/a+1)-(t0 - t).^(1/a+1));
% = 2*pi*c^2*(1./(t0-t)-1/t0) for alpha = -0.5;

% amplitude factor for a flat spectrum as a function of time, i.e. a
% generalized version of Eq. 7 in Fobel & Dau (2004).
% 
%         /             1     \1/2
%         |             - - 1 |
%         |             a     |
%         |   / t0 - t \      |
%  A_O =  |   | ------ |      |
%         |   \   c    /      |
%         | - --------------- |
%         \        a c        /
% 
A_O   = sqrt((-(t0-t).^(1/a-1))/(a*c^(1/a)));
% = sqrt(2*c^2./(t0-t).^3) for alpha = -0.5;
A_O = A_O / max(A_O);

% the chirp waveform
chirpsig = A_O.*sin(phi_O);


%% plots
if plotMe ~= 0
    
    figure(1)
    clf

    subplot(3,1,1);
    semilogy(t,f_O);
    xlim([0 max(t)]);
    ylim([f0 fmax]);
    set(gca,'ytick',[10 20 50 100 200 500 1000 2000 5000 10000 20000]);
    xlabel('time / s');
    ylabel('instantaneous frequency / Hz');

    subplot(3,1,2);
    ax = plotyy(t,phi_O/2/pi,t,A_O);
    set(ax,'xlim',[0 max(t)]);
    set(ax(1),'ylim',[0 max(phi_O/2/pi)]);
    set(ax(1),'ytick',0:3:100);
    set(ax(2),'ylim',[0 max(A_O)]);
    xlabel('time / s');
    ylabel(ax(1),'phase / cycles');
    ylabel(ax(2),'amplitude factor');

    subplot(3,1,3);
    plot(t,chirpsig);
    xlim([0 max(t)]);
    ylim([-1 1]*max(A_O));
    title('chirp waveform');
    xlabel('time / s');

    if ishandle(2)
        close(2)
    end
    figure(2)

    subplot(2,1,1);
    FFTLen = 2^16;
    win = hann(16); win = [win(1:end/2);ones(length(chirpsig)-length(win),1);win(end/2+1:end)];
    semilogx((0:FFTLen-1).'/FFTLen*fs,db(abs(fft(win.*chirpsig,FFTLen)))-10*log10(sum(chirpsig.^2)))
    xlim([f0/2 min(fmax*2,fs/2)]);
    % ylim([-1 1]*3);
    xlabel('frequency / Hz');
    ylabel('relative magnitude / dB');
    title('long-time spectrum');

    subplot(2,1,2);
    M = 128;
    [S,F,T]=spectrogram(chirpsig,gausswin(M),15/16*M,1024,fs);
    pcolor(T,F,db(abs(S)));
    set(gca,'yscale','log','ylim',[f0 fmax],'ytick',[100 200 500 1000 2000 5000 10000 20000]);
    shading flat
    caxis([-40 20])
    ylabel('frequency / Hz');
    xlabel('time / s');
    title(sprintf('short-time spectrogram (win: %1.2f ms, shift: %1.2f ms)',M/fs/1e-3,1/16*M/fs/1e-3));
end