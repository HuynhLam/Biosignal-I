% 521273S Biosignal Processing I 
% Lab 2. Filtering of the ECG Signal for the Removal of Noise
% Objectives: 
%       +removal of high-frequency noise
%       +removal of low-frequency noise
%       +study the use of different kind of filters
%
% Input:
%       ecg_signal.dat
%       The sampling rate of the signal is 1000 Hz. 
% Output:      
%       Design and interpret all filters with this sampling rate.
% 
% Useful MATLAB commands
%       filter, freqz, ones, conv, sgolayfilt
%
% $Id: ecgFilters,v1.0 2016/11/08 21:43:40 lhuynh Exp $

function ecgFilters(ecgSignal)
    
    %import data
    sRate         = 1000; %Hz 
    startInterval = 2; %interval start at 2nd second
    endInterval   = 3; %interval end at 3rd second
    x             = importdata(ecgSignal); %import the ecg input signal
    xAxis         = 1/sRate:1/sRate:(length(x)*1.0)/sRate; %xAxis for plot (time-second)
    oneCycle      = startInterval:1/sRate:endInterval; %xAxis for plot (one cycle)
    fprintf('length(x) = %.3f, length(x)/sRate = %.3f\n', length(x), length(x)/sRate);
    
    %1.Plot the original ECG signal.
    figure('Name', 'ECG Filtering', 'NumberTitle','off');
    ax1    = subplot(4,2,1);
    ax2    = subplot(4,2,2);
    plot(ax1,xAxis,x,'b');
    xlabel(ax1,'Time(s)');
    ylabel(ax1,'AU');
    title(ax1,'Original');
    plot(ax2,oneCycle,x(sRate*startInterval:sRate*endInterval),'b');
    xlabel(ax2,'Time(s)');
    ylabel(ax2,'AU');
    title(ax2,'Original(one cycle)');
    
    %2.Construct the moving average filter
    winSize = 10;
    a_MA    = 1;
    b_MA    = (1/winSize)*ones(1,winSize);
    y_MA    = filter(b_MA, a_MA, x);
    ax3     = subplot(4,2,3);
    ax4     = subplot(4,2,4);
    %plot MA
    plot(ax3,xAxis,y_MA,'b');
    xlabel(ax3,'Time(s)');
    ylabel(ax3,'AU');
    title(ax3,'MA filtered');
    plot(ax4,oneCycle,y_MA(sRate*startInterval:sRate*endInterval),'b');
    xlabel(ax4,'Time(s)');
    ylabel(ax4,'AU');
    title(ax4,'MA filtered(one cycle)');
    
    %3.Construct the derivative based filter
    % Derive derivative filter (DF) coefficients
    % H(z) = G * (1/T) * [1 - z^-1] / [1 - 0.995z^-1]
    % Assume that sampling interval T = 1
    % [1 - 0.995z^-1]Y(z) = GX(z)[1-z^-1]
    % Y(z) - Y(z)0.995z^-1 = GX(z) - GX(z)z^-1
    % Y(z) = GX(z) - GX(z)z^-1 + Y(z)0.995Z^-1
    % y(n) = Gx(n) - Gx(n-1) + 0.995y(n-1)
    % Apply to form a(1)y(n) + a(2)y(n-1) = b(1)x(n) + b(2)x(n-1) + b(3)x(n-2)
    a_DF(1) = 1;
    a_DF(2) = -0.995;
    b_DF(1) = 1;
    b_DF(2) = -1;

    % normalize gain at z=-1 to 1
    % H(z) = G (1 - z^-1) / (1 - 0.995z^-1)
    % H(z=-1) = 1 = G (1 - (-1) / ( 1 - 0.995(-1)
    % 1 = G * (2 / 1.9950) => G = 0.9975
    dGain = 0.9975;
    b_DF  = b_DF * dGain; % combine into b coefficients
    y_DF  = filter(b_DF, a_DF, x);
    %plot DF
    ax5   = subplot(4,2,5);
    ax6   = subplot(4,2,6);
    plot(ax5,xAxis,y_DF,'b');
    xlabel(ax5,'Time(s)');
    ylabel(ax5,'AU');
    title(ax5,'Derivative filtered');
    plot(ax6,oneCycle,y_DF(sRate*startInterval:sRate*endInterval),'b');
    xlabel(ax6,'Time(s)');
    ylabel(ax6,'AU');
    title(ax6,'Derivative filtered(one cycle)');
    
    %4.Comb fileter by convolution y_MA and y_DF
    b_comb = conv(b_MA, b_DF);
    a_comb = conv(a_MA, a_DF);
    y_comb = filter(b_comb, a_comb, x);
    %plot DF
    ax7   = subplot(4,2,7);
    ax8   = subplot(4,2,8);
    plot(ax7,xAxis,y_comb,'b');
    xlabel(ax7,'Time(s)');
    ylabel(ax7,'AU');
    title(ax7,'Comb filtered');
    plot(ax8,oneCycle,y_comb(sRate*startInterval:sRate*endInterval),'b');
    xlabel(ax8,'Time(s)');
    ylabel(ax8,'AU');
    title(ax8,'Comb filtered(one cycle)');
    
    %{
    winSize = 10; %windowSize is 10
    
    
    s  = importdata(spirometer); 
    b  = importdata(beltSignal);
    e1 = importdata(reCoefficient1);
    e2 = importdata(reCoefficient2);
    
    %resample spirometer signal
    s50 = resample(s,50,100);
    
    
    %init
    SS_err_f1 = 0; 
    SS_err_f2 = 0;  
    SS_tot    = 0; 
    RMSE_f1   = 0; 
    RMSE_f2   = 0;  
    R2_f1     = 0;
    R2_f2     = 0;
    y_dash    = mean(s50); % mean of spirometer signal
    
    %Calculate 2 predicted respiratory airflow models
    %element-wise multiplication
    f1=b(:,1).*e1(1) + b(:,2).*e1(2); %(1)
    f2=f1 + (b(:,1).^2).*e2(3) + (b(:,2).^2).*e2(4); %(2)
    
    %Evaluate the predicted respiratory airflow signals
    for i=1:length(s50)
        SS_err_f1 = SS_err_f1 + (s50(i)-f1(i))^2;
        SS_err_f2 = SS_err_f2 + (s50(i)-f2(i))^2;
        SS_tot    = SS_tot + (s50(i)-y_dash)^2;   
    end
    R2_f1   = 1 - (SS_err_f1/SS_tot);
    R2_f2   = 1 - (SS_err_f2/SS_tot);
    RMSE_f1 = sqrt(SS_err_f1/length(s50));
    RMSE_f2 = sqrt(SS_err_f2/length(s50));
    fprintf('R2 value for model 1: %.5f\n', R2_f1);
    fprintf('R2 value for model 2: %.5f\n', R2_f2);
    fprintf('RMSE value for model 1: %.5f\n', RMSE_f1);
    fprintf('RMSE value for model 2: %.5f\n', RMSE_f2);
    
    %Plot figure
    %{
    figure;
    ax1 = subplot(2,1,1);
    ax2 = subplot(2,1,2);
    x=1/50:1/50:60;
    plot(ax1,x,s50,'k',x,f1,'r',x,f2,'b');
    xlabel(ax1,'seconds');
    ylabel(ax1,'magnitude');
    title(ax1,'Spirometer, First predicted respiratory and Second predicted respiratory');
    plot(ax2,x,b(:,1),'b',x,b(:,2),'g');
    xlabel(ax2,'seconds');
    ylabel(ax2,'magnitude');
    title(ax2,'Chest and Abdomen signals');
    fprintf('abcd\n');
    %}
    
    numerator coefficients of the rational transfer function,
    a_MA = 1; ???
    DC component ???
    
    %}
    
    
    
end

%%