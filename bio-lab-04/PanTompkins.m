% 521273S Biosignal Processing I 
% Lab 4. Pan-Tompkins Algorithm for QRS Detection
% Objectives: 
%       + to detect QRS complexes in ECG signal using the Pan?Tompkins algorithm
%
% Input:
%       521273S_ecg.txt
%       The sampling rate of the signal is 200 Hz. 
% Output:      
%       + output of the Pan-Tompkins algorithm.
%       + Doing more(s)
%           1.Total number of beats detected
%           2.Average RR interval (in ms)
%           3.Standard deviation of RR intervals (in ms) 
%           4.Heart rate in beats per minute
% 
% Useful MATLAB commands
%       filter, ones
%
% $Id: PanTompkins,v1.0 2016/11/20 11:42:40 lhuynh Exp $

function PanTompkins(inputDetectQRS)
    
    %import data
    sRate = 200; %Hz 
    x_orig     = importdata(inputDetectQRS); %import the ecg input signal
    xAxis = 1/sRate:1/sRate:length(x_orig)/sRate; %xAxis for plot (time-second)
    x = zeros(size(x_orig));
    for i=1:length(x_orig)
        x(i) = x_orig(i) - x_orig(1);
    end
    
    %Plot the original ECG signal.
    figure('Name', 'QRS detection', 'NumberTitle','off'); 
    ax1  = subplot(6,1,1);
    plot(ax1,xAxis,x_orig,'b');
    xlabel(ax1,'Time(s)');
    ylabel(ax1,'AU(mV)');
    title(ax1,'Original ECG signal');
    x(1)
    
    %1.Lowpass filter, equation 4.7
    %fx = 200Hz -> lp filter delay 25(ms)
    delay = 5; % 25 * 200 / 100 (samples)
    b_lp(1)  = 1;
    b_lp(7)  = -2;
    b_lp(13) = 1;
    b_lp     = b_lp * 0.03125; % 1/32 * b_lp
    a_lp     = [1, -2, 1];
    y_lp     = filter(b_lp, a_lp, x);
    ax2      = subplot(6,1,2);
    plot(ax2,xAxis,y_lp,'b');
    xlabel(ax2,'Time(s)');
    ylabel(ax2,'AU(mV)');
    title(ax2,'Lowpass filter output');
    
    %2.Highpass filter, equation 4.11
    %fx = 200Hz -> lp filter delay 80(ms)
    %{
    %in the book
    delay = delay + 16; % 80 * 200 / 1000
    b_hp(33) = 0.03125; %1/32 = 0.03125
    b_hp(18) = -1;
    b_hp(17) = 1;
    b_hp(1)  = -0.03125;
    a_hp     = [1, -1];
    y_hp     = filter(b_hp, a_hp, y_lp);
    %}
    delay = delay + 16; % 80 * 200 / 1000
    b_hp(1) = -1;
    b_hp(17) = 32;
    b_hp(18) = -32;
    b_hp(33) = 1;
    b_hp = b_hp / 32;
    a_hp     = [1, -1];
    y_hp     = filter(b_hp, a_hp, y_lp);% - 0.03125*y_lp;
    
    ax3      = subplot(6,1,3);
    plot(ax3,xAxis,y_hp,'b');
    xlabel(ax3,'Time(s)');
    ylabel(ax3,'AU(mV)');
    title(ax3,'Highpass filter output');
    
    %3.Derivative filter, formular from dtask4 handout
    b_df = [1, 2, 0, -2, -1];
    b_df = b_df * 0.125; % 1/8 * b_df
    a_df = 1;
    y_df = filter(b_df, a_df, y_hp);
    ax4  = subplot(6,1,4);
    plot(ax4,xAxis,y_df,'b');
    xlabel(ax4,'Time(s)');
    ylabel(ax4,'AU(mV)');
    title(ax4,'Derivative filter output');
    
    %4.Squaring, equation 4.11
    %make the result positive and large differences 
    %resulting from QRS complexes
    y_sq = y_df.^2;
    ax5  = subplot(6,1,5);
    plot(ax5,xAxis,y_sq,'b');
    xlabel(ax5,'Time(s)');
    ylabel(ax5,'AU(mV)');
    title(ax5,'Squared output');
    
    %5.Intergrated filter, formular 4.15
    %Smoothing the output from squaring
    %N = 30, with fx = 200Hz
    winSize = 30;    
    b_ig    = (1/winSize)*ones(1,winSize);
    a_ig    = 1;
    y_ig    = filter(b_ig, a_ig, y_sq);
    
    
    [pks,locs] = findpeaks(y_ig,'MINPEAKDISTANCE',round(0.2*sRate));
    pks_z(length(pks)) = 0; 
    ax6        = subplot(6,1,6);
    plot(ax6,xAxis,y_ig,'b');
    xlabel(ax6,'Time(s)');
    ylabel(ax6,'AU(mV)');
    title(ax6,'Intergrated filter output');
    %text((locs/sRate), pks_z, num2str((1:numel(pks))'), 'parent', ax1, 'color', 'red');
    %text(((locs-delay)/sRate), pks_z, num2str((1:numel(pks))'), 'parent', ax6, 'color', 'red');
    
    [QRSstart, QRSEnd] = detectQRS(y_ig, 50, 500, 2500); 
    for i=1:length(QRSstart)
        text(((QRSstart(i))/sRate), y_ig(QRSstart(i)), 'x', 'parent', ax6, 'color', 'green');
        text(((QRSEnd(i))/sRate), y_ig(QRSEnd(i)), 'x', 'parent', ax6, 'color', 'red');
    end
    for i=1:length(QRSstart)
        text(((QRSstart(i)-delay)/sRate), x_orig(QRSstart(i)-delay), 'x', 'parent', ax1, 'color', 'green');
        text(((QRSEnd(i)-delay)/sRate), x_orig(QRSEnd(i)-delay), 'x', 'parent', ax1, 'color', 'red');
    end
    %text(((QRSstart)/sRate), y_ig, 'x', 'parent', ax6, 'color', 'green');
    %text(((QRSEnd)/sRate), y_ig, 'x', 'parent', ax6, 'color', 'red');
    %[QRSstart, QRSEnd] = detectQRS(y_ig, 50, 500, 1650); 
    %QRSstart
    %QRSEnd
    
end

function [QRSStart, QRSEnd] = detectQRS(data, blankingInterval, treshold1, treshold2) 
% This function determines the beginnings and endings of QRS complexes in given
% data based on the Thresholds and blankingInterval. The results are stored in
% the vectors QRSStart and QRSEnd (equal length) indicating the positions of
% QRS complex beginnings and endings, respectively.
% data: the P-T output from which you want to detect the QRS complexes 
% blankingInterval = the time after QRS start is blank, i.e. new start of
% QRS is not allowed
% treshold 1 = Q-wave begins here
% treshold 2 = S-wave, see Fig. 3 for thresholding

% The vectors are initialized. 
QRSStart = [];
QRSEnd   = [];

for i = 1:length(data)-1
    % Every position where the threshold is crossed is stored as QRS 
    % complex start or end depending on the direction.
    if data(i) <= treshold1 && data(i+1) > treshold1
        QRSStart = [QRSStart i]; 
    end
    if data(i) <= treshold2 && data(i+1) > treshold2 
        QRSEnd = [QRSEnd i];
    end
end

trueQRS = QRSStart;
for i = 2:length(QRSStart)
    % The QRS complexes occuring too early (during blanking interval) 
    % after previous QRS complex are left out.
    if QRSStart(i)-QRSStart(i-1) < blankingInterval
        trueQRS(i) = 0; 
    end
end

QRSStart = QRSStart(find(trueQRS));
QRSStart = QRSStart(1:length(QRSEnd)); %keep the vectors same length 
end

%%