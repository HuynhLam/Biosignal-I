% 521273S Biosignal Processing I 
% Lab 3. Lab 3. Adaptive Filtering
% Objectives:
%       +Introduction to LMS adaptive filtering
%       +Example application: Separation of mother and fetus ECG signals
%
% Input:
%       Signals.mat, file contain:
%       + Fetus signal (fhb) (pure fetus signal that is used to assess analysis results)
%       +Mother?s chest signal (mhb) (pure mother?s signal that is used to assess analysis results)
%       +Abdomen signals (abd_sig1, abd_sig2 and abd_sig3) (mixed from fetus and mother
%       signals)
%       +A real abdomen signal and mother?s chest signal (abd_sig_real and mhb_real)
%       +The sampling rate of the signal is 1000 Hz. 
% Output:      
%       fetus ECG signal.
% 
% Useful MATLAB commands(8)
%       filter, for, plot, hold, adaptfilt.lms, xlim, immse, corrcoef
%
% ->>>>> 1->10. It is all about chosing the best parameters for the
% adaptive filter. ->>>>> 11->13(addition tasks). build my own adaptive filter, 
% it could produce better results. can use the suggetion from slide. 
% Required knowledge: Optimization, Signal processing(filter & adaptive filters)
% $Id: ecgFilters,v1.0 2016/11/15 16:37:40 lhuynh Exp $

function LMSFilter()
%import data
load('521273S_signals.mat');
tm = 1/Fs:1/Fs:10;


%% section 1
%1.init time vector(tm), plot mother's chest signal
Fig1 = figure('Name', 'Case I figure', 'NumberTitle','off');
ax1  = subplot(3,1,1);
plot(Fig1, ax1, tm, mhb(1:10000), 'b');
xlabel(ax1,'Time(s)');
ylabel(ax1,'AU(mV)');
title(ax1,'Mother''s chest signal');

%2.Case I: plot the first abdomen signal
ax2  = subplot(3,1,2);
plot(Fig1, ax2, tm, abd_sig1(1:10000), 'b');
xlabel(ax2,'Time(s)');
ylabel(ax2,'AU(mV)');
title(ax2,'Abdomen signals I');
%3+4.Case I (abd_sig1):
% R2016b:adaptfilt.lms will be removed in a future release. Use dsp.LMSFilter instead.
%{
% running experiments
len            = [1,5,11,15,21]; % array of adaptive lengths
optima_cs1_len = -1;
optima_cs1_mu  = -1;
mse_cs1        = 99999;
correlate_cs1  = -99999;

for i=1:5
    for j=0.1:0.1:0.9 % 0<c<1, with mu = c/energy, energy case I = 186
        ha_cs1 = adaptfilt.lms(len(i), j/186); %adaptive filter
        [y_cs1,e_cs1]  = filter(ha_cs1, mhb, abd_sig1); %output signal 
        mse    = getMSE(fhb, e_cs1); %calculate mean square error
        corre  = corrcoef(fhb, e_cs1); %calculate correlation coefficient        
        fprintf('[i,j] = [%d,%.2f], mse = %.5f, corre = %.5f\n', len(i),j, mse, corre(2));
    end
end
%}

% the result got from running experiments, 
% lenght = 11, mu = 0.9/186, mse = 0.0.00419, correlation coefficients =
% 0.86608
ha_cs1 = adaptfilt.lms(11, 0.9/186); 
[y_cs1, e_cs1]  = filter(ha_cs1, mhb, abd_sig1);
ax3 = subplot(3,1,3);
plot(Fig1, ax3, tm, fhb(1:10000), 'b', tm, e_cs1(1:10000), 'r');
xlabel(ax3,'Time(s)');
ylabel(ax3,'AU(mV)');
title(ax3,'Evaluate case I');
%fprintf('Case I, MSE = %.5f\nCase I, Correlation Coefficient = %.5f\n', getMSE(fhb, y_cs1), getCorrCoeff(fhb, y_cs1));
fprintf('Case I, MSE = %.5f\n', getMSE(fhb, e_cs1));
a_cs1 = corrcoef(fhb, e_cs1);
fprintf('Case I, Correlation Coefficient = %.5f\n', a_cs1(2));
% end of cell 1



%% section 2
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%5+6.Case II (abd_sig2):
Fig2 = figure('Name', 'Case II figure', 'NumberTitle','off');
ax1  = subplot(3,1,1);
plot(Fig2, ax1, tm, mhb(1:10000), 'b');
xlabel(ax1,'Time(s)');
ylabel(ax1,'AU(mV)');
title(ax1,'Mother''s chest signal');

ax2  = subplot(3,1,2);
plot(Fig2, ax2, tm, abd_sig2(1:10000), 'b');
xlabel(ax2,'Time(s)');
ylabel(ax2,'AU(mV)');
title(ax2,'Abdomen signals II');
%{
% running experiments
len            = [1,5,11,15,21]; % array of adaptive lengths
optima_cs2_len = -1;
optima_cs2_mu  = -1;
mse_cs2        = 99999;
correlate_cs2  = -99999;

for i=1:5
    for j=0.1:0.1:0.9 % 0<c<1, with mu = c/energy, energy case I = 186
        ha_cs2 = adaptfilt.lms(len(i), j/186); %adaptive filter
        [y_cs2, e_cs2]  = filter(ha_cs2, mhb, abd_sig2); %output signal
        mse    = getMSE(fhb, e_cs2); %calculate mean square error
        corre  = corrcoef(fhb, e_cs2); %calculate correlation coefficient        
        fprintf('[i,j] = [%d,%.2f], mse = %.5f, corre = %.5f\n', len(i),j, mse, corre(2));
    end
end
%}


% the result got from running experiments, 
% lenght = 11, mu = 0.9/186, mse = 0.00425, correlation coefficients = 0.86820
ha_cs2 = adaptfilt.lms(11, 0.9/186); 
[y_cs2, e_cs2]  = filter(ha_cs2, mhb, abd_sig2);
ax3 = subplot(3,1,3);
plot(Fig2, ax3, tm, fhb(1:10000), 'b', tm, e_cs2(1:10000), 'r');
xlabel(ax3,'Time(s)');
ylabel(ax3,'AU(mV)');
title(ax3,'Evaluate case II');
fprintf('Case II, MSE = %.5f\n', getMSE(fhb, e_cs2));
a_cs2 = corrcoef(fhb, e_cs2);
fprintf('Case II, Correlation Coefficient = %.5f\n', a_cs2(2));
% end of cell 2



%% section 3
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%7+8.Case III (abd_sig3):
Fig3 = figure('Name', 'Case III figure', 'NumberTitle','off');
ax1  = subplot(3,1,1);
plot(Fig3, ax1, tm, mhb(1:10000), 'b');
xlabel(ax1,'Time(s)');
ylabel(ax1,'AU(mV)');
title(ax1,'Mother''s chest signal');

ax2  = subplot(3,1,2);
plot(Fig3, ax2, tm, abd_sig3(1:10000), 'b');
xlabel(ax2,'Time(s)');
ylabel(ax2,'AU(mV)');
title(ax2,'Abdomen signals III');
%{
% running experiments
len            = [1,5,15,21,50]; % array of adaptive lengths
optima_cs3_len = -1;
optima_cs3_mu  = -1;
mse_cs3        = 99999;
correlate_cs3  = -99999;

for i=3:5
    for j=0.9:0.01:0.99 % 0<c<1, with mu = c/energy, energy case I = 186
        ha_cs3 = adaptfilt.lms(len(i), j/186); %adaptive filter
        [y_cs3, e_cs3]  = filter(ha_cs3, mhb, abd_sig3); %output signal
        mse    = getMSE(fhb, e_cs3); %calculate mean square error
        corre  = corrcoef(fhb, e_cs3); %calculate correlation coefficient        
        fprintf('[i,j] = [%d,%.2f], mse = %.5f, corre = %.5f\n', len(i),j, mse, corre(2));
    end
end
%}


% the result got from running experiments, 
% lenght = 21, mu = 0.99/186, mse = 0.02013, correlation coefficients = 0.61142
ha_cs3 = adaptfilt.lms(21, 0.99/186); 
[y_cs3, e_cs3]  = filter(ha_cs3, mhb, abd_sig3); %output signal
ax3 = subplot(3,1,3);
plot(Fig3, ax3, tm, fhb(1:10000), 'b', tm, e_cs3(1:10000), 'r');
xlabel(ax3,'Time(s)');
ylabel(ax3,'AU(mV)');
title(ax3,'Evaluate case III');
fprintf('Case III, MSE = %.5f\n', getMSE(fhb, e_cs3));
a_cs3 = corrcoef(fhb, e_cs3);
fprintf('Case III, Correlation Coefficient = %.5f\n', a_cs3(2));
% end of section 3



%% section 4
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%9+10.Case IV real data (abd_sig_real & mhb_real):
Fig4 = figure('Name', 'Case IV figure', 'NumberTitle','off');
ax1  = subplot(3,1,1);
plot(Fig4, ax1, tm, abd_sig_real(1:10000), 'b');
xlabel(ax1,'Time(s)');
ylabel(ax1,'AU(mV)');
title(ax1,'Real Abdomen signals');

ax2  = subplot(3,1,2);
plot(Fig4, ax2, tm, mhb_real(1:10000), 'b');
xlabel(ax2,'Time(s)');
ylabel(ax2,'AU(mV)');
title(ax2,'Real Mother''s chest signal');

% evaluate by looking at the subfigure 
ha_cs4 = adaptfilt.lms(11, 0.9/157); % len = 11, mu = 0.9/157
[y_cs4, e_cs4]  = filter(ha_cs4, mhb_real, abd_sig_real);
ax3 = subplot(3,1,3);
plot(Fig4, ax3, tm, abd_sig_real(1:10000), 'b', tm, e_cs4(1:10000), 'r');
xlabel(ax3,'Time(s)');
ylabel(ax3,'AU(mV)');
title(ax3,'Evaluate case IV');
% end of cell 4


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%11+12+13. Implement the adaptive LMS

%}


end

function MSE = getMSE(des, out)
    MSE = mean((des - out).^2);
end

%{ 
function COR = getCorrCoeff(des, out)
    y_dash = mean(des);
    COR    = 1 - ( sum( (des - out).^2 ) / sum( (out - y_dash).^2 ) );
end
%} 

%%