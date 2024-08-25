function res = processSwayAcc(acc_AP, acc_ML, sampleFreq, plotOpt)

% this function computes stabilometric parameters according to
% Martinez-Mendez R, Sekine M, Tamura T. Postural sway parameters using a triaxial accelerometer: comparing elderly and young healthy adults. Comput Methods Biomech Biomed Engin. 2012;15(9):899-910. doi:10.1080/10255842.2011.565753
% written by Pietro Picerno on August 2024 for submission to the Journal of Athletic Training

% plotOpt = 1 to plot sway path, 0 not to

testDuration = length(acc_AP)/sampleFreq; % in sec

% compute AP-ML acceleration resultant vector:
A=sqrt(power(acc_AP,2)+power(acc_ML,2)); % acc_AP e acc_ML sono array Nx1

aam = sum(sqrt(power(acc_AP,2)+power(acc_ML,2)))/length(acc_AP); % eq. A2 Martinez-Mendez et al 2012
dist = sqrt(power(acc_AP,2)+power(acc_ML,2)); % eq. 3 Prieto
mDist = sum(dist)/length(acc_AP); % eq. 4 Prieto
% calcola root mean square
rms = sqrt(sum(A.^2)/length(A));
% compute range
range = max(A)-min(A);
% compute "path" length
pathLength = sum(sqrt(power(diff(acc_AP),2)+power(diff(acc_ML),2))); % eq. 8 Prieto 1996
% note that pathLength equations are checked using raw data and results of https://nbviewer.org/github/BMClab/BMC/blob/master/notebooks/Stabilography.ipynb
% pathLength_ap = sum(abs(diff(cop_ap))) % Duarte Python code, eq. 9 Prieto 1996
% pathLength_ml = sum(abs(diff(cop_ml))) % Duarte Python code, eq. 9 Prieto 1996
% pathLength_hor = sum(sqrt(power(diff(cop_ap),2)+power(diff(cop_ml),2))) % eq. 8 Prieto 1996

for i = 1:length(acc_AP)-1
    swayArea_temp(i) = abs((acc_AP(i+1)*acc_ML(i))-(acc_AP(i)*acc_ML(i+1)));
end
swayArea = sum(swayArea_temp)/(2*testDuration);

jerk = sum(sqrt(power(diff(acc_AP),2)+power(diff(acc_ML),2)))/testDuration; % eq. 8 Prieto 1996

% The mean frequency (MFREQ) is the rotational frequency, in revolutions per second or Hz, of the COP if it had traveled
% the total excursions around a circle with a radius of the mean distance
swayFreq = pathLength/(2*pi*mDist*15); % eq. 20 Prieto 1996
% which is equal to swayFreq = jerk/dot(2*pi,aam); % eq. A12 Rigoberto


% calcola area dell'ellipse e asse principale di sway
[ellipseArea, mainSwayDir] = computePEA(acc_ML,acc_AP,0.95,1);


res = [jerk, mDist, rms, range, swayArea, ellipseArea, swayFreq]; 