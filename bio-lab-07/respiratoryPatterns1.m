%% [Lab 7] * 521273S Biosignal Processing I 
% Characterization of Respiratory Patterns
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
% $Id: respiratoryPatterns,v1.1 2016/12/12 22:44:25 lhuynh Exp $

%% Init data
fs        = 50; %(Hz)
airFlow   = load('521273S_respiratoryairflow.txt');
N         = length(airFlow);
time      = 1/fs:1/fs:N/fs;
threshold = 20;


%% Compute respiratory rate [breaths per minute]
%breath = inspiration + expiration
inBegin = [];
inEnd   = [];
exBegin = [];
exEnd   = [];
mark    = 1;
for i=1:N-1
    if (airFlow(i)>0 && airFlow(i+1)<=0) %expiration count (ABOVE)
        if(abs(airFlow(mark))<threshold)
            exBegin = [exBegin mark];
            exEnd   = [exEnd i];
        end
        mark = i + 1;
    elseif (airFlow(i)<0 && airFlow(i+1)>=0) %inspiration count (BELOW)
        if(abs(airFlow(mark))<threshold)
            inBegin = [inBegin mark];
            inEnd   = [inEnd i];
        end
        mark = i + 1;
    end
end
%print result #2
fprintf('2.Respiration Rate                        = %.5f(bpm)\n', min(length(inBegin), length(exBegin))*60/ ((exEnd(length(exEnd))-inBegin(1))/fs));


%% Compute average duty cycle Ti/Ttot && average Tptef/Te
M = length(inEnd); %number of breath
dutyCycle = zeros(1, M);
Tptef     = zeros(1, M);
inPeak    = zeros(1, M);
exPeak    = zeros(1, M);
for i=1:M
    Ti            = inEnd(i)-inBegin(i);
    Te            = exEnd(i)-exBegin(i);
    dutyCycle(i)  = Ti / (Ti + Te);
    
    [~,exPeak(i)] = max(airFlow(exBegin(i):exEnd(i)));    
    Tptef(i)      = exPeak(i) / (Te);
    exPeak(i)     = exPeak(i) + exBegin(i);
    
    [~,inPeak(i)] = min(airFlow(inBegin(i):inEnd(i)));
    inPeak(i)     = inPeak(i) + inBegin(i);
end
%print result #3 & #4
fprintf('3.Average duty cycle                      = %.5f\n', mean(dutyCycle));
fprintf('4.Average time to peak expiratory airflow = %.5f\n', mean(Tptef));


%% Plot the figure
figure('Name', 'Airflow signal', 'NumberTitle', 'off');
plot(time, airFlow, 'c', time(inPeak), airFlow(inPeak), 'b*', time(exPeak), airFlow(exPeak), 'r*');
xlim([0 N/fs]);
hline       = refline([0,0]);
hline.Color = 'k';
xlabel('Time(seconds)');
ylabel('ml / s');
title('Airflow signal');


%% Compute average and standard deviation of peak-to-peak airflow
meanP2P = zeros(1, M);
for i=1:M
    meanP2P(i) = abs(airFlow(inPeak(i))) + airFlow(exPeak(i));
end
%print result #6
fprintf('6.Average Deviation of peak-to-peak       = %.5f\n', mean(meanP2P));
fprintf('6.Standard Deviation of peak-to-peak      = %.5f\n', std(meanP2P));


%% Compute the average expiratory volume(in millimeter)
expiratoryVolume = 0;
for i=1:M
    expiratoryVolume = expiratoryVolume + trapz(time(exBegin(i):exEnd(i)), airFlow(exBegin(i):exEnd(i)));
end
fprintf('7.The average expiratory volume           = %.5f(l)\n', abs(expiratoryVolume)/M/1000);
