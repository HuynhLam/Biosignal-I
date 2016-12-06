% 521273S Biosignal Processing I 
% Lab 1. Respiration analysis
% Objectives: to become familiar with the basics of Matlab programming
%
% Input:
% spirometer.txt, beltSignals.txt (includes chest (1. column) and abdomen (2. column) 
% respiratory effort belt signals), regressionCoefficients1.txt and regressionCoefficients2.txt
% already time-synchronized
% 
% respiratory effort belt signals 50Hz
% spirometer signal 100Hz
%
% $Id: respirationAnalysis,v1.0 2016/11/01 08:43:40 lhuynh Exp $

function respirationAnalysis(spirometer,beltSignal,reCoefficient1,reCoefficient2)
    
    %import data
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
    SS_err_f1_sum = sum(( s50(i)-f1(i) )^2);
    fprintf('SS error f1 sum: %.5f\n', SS_err_f1);
    fprintf('SS error f1 sum sum: %.5f\n', SS_err_f1_sum);
    R2_f1   = 1 - (SS_err_f1/SS_tot);
    R2_f2   = 1 - (SS_err_f2/SS_tot);
    RMSE_f1 = sqrt(SS_err_f1/length(s50));
    RMSE_f2 = sqrt(SS_err_f2/length(s50));
    fprintf('R2 value for model 1: %.5f\n', R2_f1);
    fprintf('R2 value for model 2: %.5f\n', R2_f2);
    fprintf('RMSE value for model 1: %.5f\n', RMSE_f1);
    fprintf('RMSE value for model 2: %.5f\n', RMSE_f2);
    
    %Plot figure
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
   
    
end

%%