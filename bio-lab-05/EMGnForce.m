% 521273S Biosignal Processing I 
% Lab 5. Analysis of the Relationship between Parameters of the EMG Signal and Muscular Force
% DR, RMS, ZCR, TCR
% Objectives: 
%       + To characterize the level of activity in EMG signals.
%       + To study the relationship between parameters of the EMG signal and muscular force.
%
% Input:
%       521273S_emgforce.txt
%       The sampling rate is 2000 Hz per channel and the EMG sample values are in mV. 
% Output:      
%       variation DR, RMS, ZCR, TCR in EMG
% 
% Useful MATLAB commands
%       polyfit, polyval, find
%
% $Id: EMGnForce,v1.0 2016/11/28 10:31:34 lhuynh Exp $

%function EMGnForce()
%1.load EMG signals
fs = 2000; %Hz
x = load('521273S_emgforce.txt');
time = x(:, 1);
force = x(:, 2);
emgOrig = x(:, 3);
minForce = min(force);
maxForce = max(force);

%2,3. Normalize force signal, remove DC bias from emg signals
for i=1:length(force)
    force(i) = (force(i) - minForce) / (maxForce - minForce) * 100;
end
emg = detrend(emgOrig,0);

%4. Plot %MVC and emg signal
Fig1 = figure('Name', 'Force and EMG signals', 'NumberTitle', 'off');
ax1 = subplot(2,1,1);
plot(Fig1, ax1, time, force, 'color', [0.0 0.749 1.0]);
xlabel(ax1, 'Time(seconds)');
ylabel(ax1, 'Force(%MVC)');
title(ax1,'Force signal');
ax2 = subplot(2,1,2);
plot(Fig1, ax2, time, emg, 'b');
xlabel(ax2, 'Time(seconds)');
ylabel(ax2, 'mV');
title(ax2,'EMG signal');

%5,6. identify segment & mark
%init vector
beg10 = [];
end10 = [];
max10 = [];
beg75 = [];
end75 = [];
for i=1:length(force) - 1
    if(force(i) < 10 && force(i+1) >= 10)
        beg10 = [beg10 i];
        text(i/fs, force(i), '*', 'parent', ax1, 'color', 'green', 'FontSize', 25);
    elseif (force(i) >= 10 && force(i+1) < 10)
        end10 = [end10 i];
        text(i/fs, force(i), 'o', 'parent', ax1, 'color', 'red', 'HorizontalAlignment', 'center');
    end
end
%{
beg10
tmp10 = [];
fg = true;
sr = binSearch(force, length(force), 10);
while (fg)
    fprintf('sr = %d\n', sr);
    if(sr == -1)
       break; 
    else
        tmp10 = [tmp10 sr];
        sr = binSearch(force(1:sr), length(force(1:sr)), 10);
    end
end
tmp10
%}

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

%7,8,9.Average force exerted
% expressed in %MVC(page 21)
N = length(beg75); % N=5 in this exercise
avgFExerted = []; 
DRemg = [];
RMSemg = [];

r_2DRSum = [0,0,0];
r_2RMSSum = [0,0,0];

for i=1:N
    avgFExerted = [avgFExerted mean( force(beg75(i):end75(i)) )];
    DRemg = [DRemg ( max(emg(beg75(i):end75(i))) - min(emg(beg75(i):end75(i))) )];
    tmpSum = 0;
    for j=beg75(i):end75(i)
        tmpSum = tmpSum + emg(j)^2;
    end
    RMSemg = [RMSemg sqrt( tmpSum / (end75(i)-beg75(i)) )];
    
    %calculate sum correlation coefficients
    r_2DRSum(1) = r_2DRSum(1) + avgFExerted(i)*DRemg(i);
    r_2DRSum(2) = r_2DRSum(2) + avgFExerted(i)^2; 
    r_2DRSum(3) = r_2DRSum(3) + DRemg(i)^2;
    
    r_2RMSSum(1) = r_2RMSSum(1) + avgFExerted(i)*RMSemg(i);
    r_2RMSSum(3) = r_2RMSSum(3) + RMSemg(i)^2;
    
end

%10,11. Plot DR & RMS
Fig2 = figure('Name', 'DR & RMS parameters', 'NumberTitle', 'off');
ax1 = subplot(2,1,1);
%https://www.mathworks.com/help/matlab/data_analysis/programmatic-fitting.html#f1-7407
%e.g modeling this data using a first-degree polynomial function(line)
p1 = polyfit(avgFExerted,DRemg,1); % return polynomial coefficients
%Evaluate the polynomial at uniformly spaced times, t1(50 points)
t1 = 21.2075:1.323164:87.3657;
f1 = polyval(p1,t1);
plot(Fig2, ax1, avgFExerted, DRemg, 'o', t1, f1, 'b');
xlabel(ax1, 'Force(%MVC)');
ylabel(ax1, 'DR value of EMG(%MVC)');
title(ax1,'Dynamic Range');

ax2 = subplot(2,1,2);
p2 = polyfit(avgFExerted,RMSemg,1);
f2 = polyval(p2,t1);
plot(Fig2, ax2, avgFExerted, RMSemg, 'o', t1, f2, 'b');
xlabel(ax2, 'Force(%MVC)');
ylabel(ax2, 'RMS value of EMG(mV)');
title(ax2,'Root Mean Square');

%12. correlation coefficient of DR and RMS
avgDash = mean(avgFExerted);
DRDash = mean(DRemg);
RMSDash = mean(RMSemg);
r_2DR = sqrt( ((r_2DRSum(1) - N*avgDash*DRDash)^2) / ( (r_2DRSum(2) - N*(avgDash^2)) * (r_2DRSum(3) - N*(DRDash^2))) );
r_2RMS = sqrt( ((r_2RMSSum(1) - N*avgDash*RMSDash)^2) / ( (r_2DRSum(2) - N*(avgDash^2)) * (r_2RMSSum(3) - N*(RMSDash^2))) );
fprintf('Correlation coefficient for Dynamic Range, r = %.5f\n', r_2DR);
fprintf('Correlation coefficient for Root Mean Square, r = %.5f\n', r_2RMS);

%Good to know that in text
% 'r^2' will print r in square
% 'r_2' will print r 2(small + under)
text(25, 2.5, ['r = ', num2str(r_2DR)], 'parent', ax1, 'color', 'black');
text(25, 0.32, ['r = ', num2str(r_2RMS)], 'parent', ax2, 'color', 'black');

%*******
%13. ZCR & TCR
%plot for check points
%{
Fig3 = figure('Name', 'EMG signals', 'NumberTitle', 'off');
ax1 = subplot(1,1,1);
plot(Fig3, ax1, time, emg, 'b');
xlabel(ax1,'Time(seconds)');
ylabel(ax1,'mV');
title(ax1,'EMG signal');
%}

ZCR = zeros([1 N]);
TCR = zeros([1 N]);

r_2ZCRSum = [0,0,0];
r_2TCRSum = [0,0,0];
for i=1:N
    cntZ = 0;
    cntT = 0;
    mark0 = [];
    for j=beg75(i):end75(i)-1
        %find cross zero point
        if( (emg(j)<0 && emg(j+1)>=0) || (emg(j)>0 && emg(j+1)<=0) ) 
            cntZ = cntZ + 1;
        end
        %find every turn points
        if( (emg(j-1)<=emg(j) && emg(j)>emg(j+1)) || (emg(j-1)>=emg(j) && emg(j)<emg(j+1)) )
            mark0 = [mark0 j]; 
            cntT = cntT + 1;
        end
    end
    
    %fprintf('Every turns point = %d\n', cntT);
    %fprintf('time[%d] = %.3f\n', i, (end75(i)-beg75(i)) / fs);
    ZCR(i) = cntZ / ((end75(i)-beg75(i)) / fs);
    TCR(i) = cntT;
    %find significant turn point 
    cntT = 0;
    for j=1:TCR(i)-1
        if( abs( emg(mark0(j)) - emg(mark0(j+1)) ) >= 0.1 )
            cntT = cntT + 1;
        end
    end
    %fprintf('Significants turns point = %d\n', cntT);
    TCR(i) = cntT / ((end75(i)-beg75(i)) / fs);                                                                                                                        
    
    %calculate sum correlation coefficients
    r_2ZCRSum(1) = r_2ZCRSum(1) + avgFExerted(i)*ZCR(i);
    r_2ZCRSum(3) = r_2ZCRSum(3) + ZCR(i)^2;
    
    r_2TCRSum(1) = r_2TCRSum(1) + avgFExerted(i)*TCR(i);                                                                                               
    r_2TCRSum(3) = r_2TCRSum(3) + TCR(i)^2;
end

%correlation coefficient of DR and RMS
ZCRDash = mean(ZCR);                                                                                                                                                            
TCRDash = mean(TCR);
r_2ZCR = sqrt( ((r_2ZCRSum(1) - N*avgDash*ZCRDash)^2) / ( (r_2DRSum(2) - N*(avgDash^2)) * (r_2ZCRSum(3) - N*(ZCRDash^2))) );
r_2TCR = sqrt( ((r_2TCRSum(1) - N*avgDash*TCRDash)^2) / ( (r_2DRSum(2) - N*(avgDash^2)) * (r_2TCRSum(3) - N*(TCRDash^2))) );
fprintf('Correlation coefficient for Zero-crossing Rate, r = %.5f\n', r_2ZCR);
fprintf('Correlation coefficient for Turns Count Rate, r = %.5f\n', r_2TCR);

% plot figure for ZCR & TCR
Fig3 = figure('Name', 'ZCR & TCR parameters', 'NumberTitle', 'off');
ax1 = subplot(2,1,1);
%https://www.mathworks.com/help/matlab/data_analysis/programmatic-fitting.html#f1-7407
%e.g modeling this data using a first-degree polynomial function(line)
p1 = polyfit(avgFExerted,ZCR,1); % return polynomial coefficients
% 0.40336 110.6836
p1
%Evaluate the polynomial at uniformly spaced times, t1(50 points)
t1 = 21.2075:1.323164:87.3657;
f1 = polyval(p1,t1);
plot(Fig3, ax1, avgFExerted, ZCR, 'o', t1, f1, 'b');
xlabel(ax1, 'Force(%MVC)');
ylabel(ax1, 'ZCR value of EMG(Hz)');
title(ax1,'Zero-crossing Rate');

ax2 = subplot(2,1,2);
p2 = polyfit(avgFExerted,TCR,1);
% 1.7632 90.6385
p2
f2 = polyval(p2,t1);
plot(Fig3, ax2, avgFExerted, TCR, 'o', t1, f2, 'b');
xlabel(ax2, 'Force(%MVC)');
ylabel(ax2, 'TCR value of EMG(Hz)');
title(ax2,'Turns Count Rate');

%{
mark0 = [];
cnt = 0;
for j=beg75(1):end75(1)-1
        %if( (emg(i)<0 && emg(i+1)>=0) || (emg(i)>0 && emg(i+1)<=0) )
        if( (emg(j-1)<=emg(j) && emg(j)>emg(j+1)) || (emg(j-1)>=emg(j) && emg(j)<emg(j+1)) )
            mark0 = [mark0 j]; 
            cnt = cnt + 1;
        end
end
cnt
mark1 = [];
cntTC = 0;
for i=1:cnt-1
    if( abs( emg(mark0(i)) - emg(mark0(i+1)) ) > 0.1 )
        mark1 = [mark1 mark0(i)];
        mark1 = [mark1 mark0(i+1)];
        cntTC = cntTC + 1;
    end
end
cntTC
for i=1:length(mark1)/5
    text(mark1(i)/fs, emg(mark1(i)), 'x', 'color', 'red');
end
%}


%%end


%{
function idx = binSearch(A, n, num)
% Syntax: [idx] = binSearch(A, n, num);            
% Complexity: O( log_2(n) )
% Inputs:       
%         + A:sorted array
%         + n:length(A)
%         + num:number I need to find               
% Outputs:      
%         + index:  position in A that A(index) == num
%                   -1 if num does not exist in A
left = 1;
right = n;
flag = 0;

while left <= right
    mid = ceil((left + right) / 2);
    if (A(mid) < num && A(mid+1) >= num)
        idx = mid;
        flag = 1;
        break;
    elseif (A(mid) > num)
        right = mid - 1;
    else    
        left = mid + 1;
    end
end

if (flag == 0)
    idx = -1;
end
end
%}
%%