% 521273S Biosignal Processing I 
% Lab 6. Frequency-domain Analysis of EMG Signals
% Objectives: 
%       + To study methods for frequency-domain analysis of EMG signal.
%
% Input:
%       521273S_EMGforce1.txt
%       521273S_EMGforce2.txt
%       The sampling rate is 2000 Hz per channel and the EMG sample values are in mV. 
% Output:      
%       power spectral density (PSD) of EMG signals
% 
% Useful MATLAB commands
%       diff, buffer, pwelch, meanfreq, mean, for
%
% Power spectral density function (PSD) shows the strength of the variations(energy) 
% as a function of frequency. In other words, it shows at which frequencies variations 
% are strong and at which frequencies variations are weak. The unit of PSD is 
% energy per frequency(width) and you can obtain energy within a specific frequency 
% range by integrating PSD within that frequency range. Computation of PSD is done 
% directly by the method called FFT or computing autocorrelation function and 
% then transforming it. 
% http://www.cygres.com/OcnPageE/Glosry/SpecE.html
%
% $Id: PSD_EMG,v1.0 2016/11/28 10:31:34 lhuynh Exp $
% $Id: PSD_EMG,v1.1 2016/12/05 11:47:16 lhuynh Exp $
% $Id: PSD_EMG,v1.2 2016/12/06 18:43:01 lhuynh Exp $

%% Analysize EMGforce1
%1.a.load EMGforce1 signals
fs = 2000; %Hz
threshold = 10; %  10(%MVC)
x = load('521273S_emgforce1.txt');
time = x(:, 1);
force = x(:, 2);
emgOrig = x(:, 3);
minForce = min(force);
maxForce = max(force);
fprintf('***** Analysize EMGforce1 *****\n');

%1.b,c. Normalize force signal, remove DC bias from emg signals
for i=1:length(force)
    force(i) = (force(i) - minForce) / (maxForce - minForce) * 100;
end
emg = detrend(emgOrig,0); %remove DC bias

%1.d. Identify and plot force and emg signal
figure('Name', 'Force && EMG signal 1', 'NumberTitle', 'off');
ax1 = subplot(2,1,1); %*** subplot #1
plot(ax1, time, force, 'color', [0.0 0.749 1.0]);
xlabel(ax1, 'Time(seconds)');
ylabel(ax1, 'Force(%MVC)');
title(ax1,'Force signal');
ax2 = subplot(2,1,2); %*** subplot #1
plot(ax2, time, emg, 'b');
xlabel(ax2, 'Time(seconds)');
ylabel(ax2, 'mV');
title(ax2,'EMG signal');

%Identify segment & mark
%init vector
beg10 = [];
end10 = [];
max10 = [];
beg75 = [];
end75 = [];
for i=1:length(force) - 1
    if(force(i) < threshold && force(i+1) >= threshold )
        beg10 = [beg10 i];
        text(i/fs, force(i), '*', 'parent', ax1, 'color', 'green', 'FontSize', 25);
    elseif (force(i) >= threshold && force(i+1) < threshold)
        end10 = [end10 i];
        text(i/fs, force(i), 'o', 'parent', ax1, 'color', 'red', 'HorizontalAlignment', 'center');
    end
end

for i=1:length(beg10)
    max10 = [max10 max(force(beg10(i):end10(i)))];
    flag = false;
    for j=beg10(i):end10(i)
        if(force(j) < 0.75*max10(i) && force(j+1) >= 0.75*max10(i))
            beg75 = [beg75 j];
            text(j/fs, force(j), '*', 'parent', ax1, 'color', [1.0 0.549 0.0], 'FontSize', 25);
        elseif (force(j) >= 0.75*max10(i) && force(j+1) < 0.75*max10(i))
            end75 = [end75 j];
            text(j/fs, force(j), 'o', 'parent', ax1, 'color', 'blue', 'HorizontalAlignment', 'center');
        elseif (force(j) == max10(i) && ~flag)
            text(j/fs, force(j)-0.1, '*', 'parent', ax1, 'color', [0.294 0.0 0.5098], 'FontSize', 25, 'HorizontalAlignment', 'center');
            flag = true;
        end
    end 
end

%2,3,4. Study the one-sided modified periodogram estimate of the PSD of each
%section && compute the mean frequency of the average PSD
%a. Segment the identified signal section into consecutive(partially
%overlapping, 100 sample overlap) parts of length 750 samples. 
N = length(beg75);
nfft = 256; % number of DFT points
segmentL = 750; % segment length of samples
overlapL = 100; % overlapping samples
figure('Name', 'PSD of EMGforce 1', 'NumberTitle', 'off');
for i=1:N
    y = buffer(emg(beg75(i):end75(i)), segmentL, overlapL);
    %computer the one-side PSD for each identified section 
    [pxx, f] = pwelch(y, segmentL, overlapL, nfft, fs); 
    [sz1, sz2] = size(pxx);    
    pMeanxx = zeros(1, sz1);
    for j=1:sz1
        pMeanxx(j) = mean(pxx(j, :));
    end
    
    meanFrequency = sum(f.*(pMeanxx)') / sum(pMeanxx);
    freq = meanfreq(pMeanxx, f);
    fprintf('[%d] meanFrequency = %.5f, matlab meanfreq = %.5f\n', i, meanFrequency, freq);
    
    %figure with each of the EMG signal sections and the corresponding average PSD?s.
    %one for time(s) & one for frequency(Hz)
    %*** subplot #1
    ax(2*i-1) = subplot(N,2,2*i-1); 
    plot(ax(2*i-1), time(beg75(i):end75(i)), emg(beg75(i):end75(i)), 'color', 'b');
    xlabel(ax(2*i-1), 'Time(seconds)');
    ylabel(ax(2*i-1), 'mV');
    title(ax(2*i-1),['Section [', num2str(i),'] in time domain']);
    %*** subplot #2
    ax(2*i) = subplot(N,2,2*i); 
    %plot(ax(2*i), f, pxx);
    plot(ax(2*i), f, pMeanxx, 'color', [0.0 0.749 1.0]);
    xlabel(ax(2*i), 'Hz');
    ylabel(ax(2*i), 'mV^2 / Hz'); %assume resistance = 1 ohm
    title(ax(2*i),['Section [', num2str(i),'] in frequency domain']);
    
end % end of identified sections travel
%%end EMGforce1


%% Analysize EMGforce2
%load EMGforce2 signals
x = load('521273S_emgforce2.txt');
time = x(:, 1);
force = x(:, 2);
emgOrig = x(:, 3);
minForce = min(force);
maxForce = max(force);
fprintf('***** Analysize EMGforce2 *****\n');

%Normalize force signal, remove DC bias from emg signals
for i=1:length(force)
    force(i) = (force(i) - minForce) / (maxForce - minForce) * 100;
end
emg = detrend(emgOrig,0); %remove DC bias

%Identify and plot force and emg signal
figure('Name', 'Force && EMG signal 2', 'NumberTitle', 'off');
ax1 = subplot(2,1,1); %*** subplot #1
plot(ax1, time, force, 'color', [0.0 0.749 1.0]);
xlabel(ax1, 'Time(seconds)');
ylabel(ax1, 'Force(%MVC)');
title(ax1,'Force signal');
ax2 = subplot(2,1,2); %*** subplot #1
plot(ax2, time, emg, 'b');
xlabel(ax2, 'Time(seconds)');
ylabel(ax2, 'mV');
title(ax2,'EMG signal');

%Identify segment & mark
%init vector
beg10 = [];
end10 = [];
max10 = [];
beg75 = [];
end75 = [];
for i=1:length(force) - 1
    if(force(i) < threshold && force(i+1) >= threshold )
        beg10 = [beg10 i];
        text(i/fs, force(i), '*', 'parent', ax1, 'color', 'green', 'FontSize', 25);
    elseif (force(i) >= threshold && force(i+1) < threshold)
        end10 = [end10 i];
        text(i/fs, force(i), 'o', 'parent', ax1, 'color', 'red', 'HorizontalAlignment', 'center');
    end
end

for i=1:length(beg10)
    max10 = [max10 max(force(beg10(i):end10(i)))];
    flag = false;
    for j=beg10(i):end10(i)
        if(force(j) < 0.7*max10(i) && force(j+1) >= 0.7*max10(i))
            beg75 = [beg75 j];
            text(j/fs, force(j), '*', 'parent', ax1, 'color', [1.0 0.549 0.0], 'FontSize', 25);
        elseif (force(j) >= 0.7*max10(i) && force(j+1) < 0.7*max10(i))
            end75 = [end75 j];
            text(j/fs, force(j), 'o', 'parent', ax1, 'color', 'blue', 'HorizontalAlignment', 'center');
        elseif (force(j) == max10(i) && ~flag)
            text(j/fs, force(j)-0.1, '*', 'parent', ax1, 'color', [0.294 0.0 0.5098], 'FontSize', 25, 'HorizontalAlignment', 'center');
            flag = true;
        end
    end 
end

%Study the one-sided modified periodogram estimate of the PSD of each
%section && compute the mean frequency of the average PSD
%Segment the identified signal section into consecutive(partially
%overlapping, 100 sample overlap) parts of length 750 samples. 
N = length(beg75);
nfft = 256; % number of DFT points
segmentL = 750; % segment length of samples
overlapL = 100; % overlapping samples
figure('Name', 'PSD of EMGforce 2', 'NumberTitle', 'off');
for i=1:N
    y = buffer(emg(beg75(i):end75(i)), segmentL, overlapL);
    %computer the one-side PSD for each identified section 
    % mW = mV^2 / om
    [pxx, f] = pwelch(y, segmentL, overlapL, nfft, fs); 
    [sz1, sz2] = size(pxx);    
    pMeanxx = zeros(1, sz1);
    for j=1:sz1
        pMeanxx(j) = mean(pxx(j, :));
    end
    
    meanFrequency = sum(f.*(pMeanxx)') / sum(pMeanxx);
    freq = meanfreq(pMeanxx, f);
    fprintf('[%d] meanFrequency = %.5f, matlab meanfreq = %.5f\n', i, meanFrequency, freq);
    
    %figure with each of the EMG signal sections and the corresponding average PSD?s.
    %one for time(s) & one for frequency(Hz)
    %*** subplot #1
    ax(2*i-1) = subplot(N,2,2*i-1); 
    plot(ax(2*i-1), time(beg75(i):end75(i)), emg(beg75(i):end75(i)), 'color', 'b');
    xlabel(ax(2*i-1), 'Time(seconds)');
    ylabel(ax(2*i-1), 'mV');
    title(ax(2*i-1),['Section [', num2str(i),'] in time domain']);
    %*** subplot #2
    ax(2*i) = subplot(N,2,2*i); 
    %plot(ax(2*i), f, pxx);
    plot(ax(2*i), f, pMeanxx, 'color', [0.0 0.749 1.0]);
    xlabel(ax(2*i), 'Hz');
    ylabel(ax(2*i), 'mV^2 / Hz'); %assume resistance = 1 ohm
    title(ax(2*i),['Section [', num2str(i),'] in frequency domain']);
end % end of identified sections travel

%%end EMGforce2