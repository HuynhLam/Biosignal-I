% 521273S Biosignal Processing I 
% Lab 7. Characterization of Respiratory Patterns
% Objectives: 
%       + To compute parameters related to respiratory patterns
%
% Input:
%       + 521273S_respiratoryairflow.txt
%       + The sampling rate is 50(Hz) and the airflow signal unit is (ml/s). 
%       + Inspiration is below and expiration is above zero line.
% Output:      
%       Compute inspiratory flow and expiratory flow
% 
% Useful MATLAB commands
%       refline
%
% $Id: respiratoryPatterns,v1.0 2016/12/08 16:00:16 lhuynh Exp $

%% Analysis the respiratory airflow signal
%1. Plot the signal
fs        = 50; %(Hz)
airFlow   = load('521273S_respiratoryairflow.txt');
N         = length(airFlow);
time      = 1/fs:1/fs:N/fs;


%2,3,4,5. Compute respiratory rate [breaths per minute]
% breath = inspiration + expiration
inBegin = [];
inEnd   = [];
exBegin = [];
exEnd   = [];
mark    = 1;
Ti      = 0;
Te      = 0;
for i=1:N-1
    if (airFlow(i)<0 && airFlow(i+1)>=0) %expiration count
        Te      = Te + (i-mark+1)/fs;
        exBegin = [exBegin mark];
        exEnd   = [exEnd i];
        mark    = i + 1;
    elseif (airFlow(i)>0 && airFlow(i+1)<=0) %inspiration count
        Ti      = Ti + (i-mark+1)/fs;
        inBegin = [inBegin mark];
        inEnd   = [inEnd i];
        mark    = i + 1;
    end
end

% analyze the signal's tail
if (airFlow(N)<0) %expiration count
    Te      = Te + (N-mark+1)/fs;
    exBegin = [exBegin mark];
    exEnd   = [exEnd N];
elseif (airFlow(N)>0) %inspiration count
    Ti      = Ti + (N-mark+1)/fs;
    inBegin = [inBegin mark];
    inEnd   = [inEnd N];
end


dutyCycle = zeros(1, length(inBegin));
Tptef     = zeros(1, length(inBegin));
inPeak    = zeros(1, length(inBegin));
exPeak    = zeros(1, length(inBegin));
for i=1:length(inBegin)
    dutyCycle(i)    = (inEnd(i)-inBegin(i)+1) / ((inEnd(i)-inBegin(i))+exEnd(i)-exBegin(i)+2);
    [pks,exPeak(i)] = min(airFlow(exBegin(i):exEnd(i)));    
    Tptef(i)        = exPeak(i) / (exEnd(i)-exBegin(i));
    exPeak(i)       = exPeak(i)+exBegin(i)-1;
    [pks,inPeak(i)] = max(airFlow(inBegin(i):inEnd(i)));
    inPeak(i)       = inPeak(i)+inBegin(i)-1;
end
%testTtot
figure('Name', 'Airflow signal', 'NumberTitle', 'off');
%MarkerFaceColor','green' %fill with the color
plot(time, airFlow, 'c', time(inPeak), airFlow(inPeak), 'bv', time(exPeak), airFlow(exPeak), 'm^');
xlim([0 N/fs]);
hline       = refline([0,0]);
hline.Color = 'r';
xlabel('Time(seconds)');
ylabel('ml / s');
title('Airflow signal');


%6. Calculate 
%http://www.statisticshowto.com/average-deviation/
%http://www.statisticshowto.com/what-is-standard-deviation/#HFSSD
meanPeak = 0;
for i=1:length(inBegin)
    meanPeak = meanPeak + airFlow(inPeak(i)) + airFlow(exPeak(i));
end
meanPeak     = meanPeak / (length(inBegin) * 2);
avgDeviation = 0;
stdDeviation = 0;
for i=1:length(inBegin)
    avgDeviation = avgDeviation + abs(airFlow(inPeak(i)) - meanPeak) + abs(airFlow(exPeak(i)) - meanPeak);
    stdDeviation = stdDeviation + (airFlow(inPeak(i)) - meanPeak)^2 + (airFlow(exPeak(i)) - meanPeak)^2;
end
avgDeviation = avgDeviation / (length(inBegin) * 2);
stdDeviation = stdDeviation / ((length(inBegin) * 2) - 1);
stdDeviation = sqrt(stdDeviation);
%double check the standard deviation
%{
x = zeros(1, 24);
for i=1:12
    x(i*2-1) = airFlow(inPeak(i));
    x(i*2)   = airFlow(exPeak(i));
end
std(x)
%}

%7.Compute the average expiratory volume(in millimeter)
%use intergration by trapz. http://www.mathworks.com/help/matlab/ref/trapz.html#bua4lsr
expiratoryVolume = 0;
for i=1:length(exBegin)
    expiratoryVolume = expiratoryVolume + trapz(time(exBegin(i):exEnd(i)), airFlow(exBegin(i):exEnd(i)));
    %fprintf('begin point, [%.2f,%.2f]\n', time(exBegin(i)), airFlow(exBegin(i)));
    %fprintf('end point,   [%.2f,%.2f]\n', time(exEnd(i)), airFlow(exEnd(i)));
    %fprintf('[%d], expiratory volume = %.5f\n', i, expiratoryVolume);
end

%8. print all the results
%fprintf('Average duty cycle(Ti/Ttot) = %.5f\n', Ti/(Ti+Te));
fprintf('2.Respiration Rate                        = %.5f(bpm)\n', min(length(exBegin), length(inBegin))*60/(N/fs));
fprintf('3.Average duty cycle                      = %.5f\n', mean(dutyCycle));
fprintf('4.Average time to peak expiratory airflow = %.5f\n', mean(Tptef));
fprintf('6.Average Deviation of peak-to-peak       = %.5f\n', avgDeviation);
fprintf('6.Standard Deviation of peak-to-peak      = %.5f\n', stdDeviation);
fprintf('7.The average expiratory volume           = %.5f(ml)\n', abs(expiratoryVolume)/length(exBegin));

%%end respiratoryPatterns
