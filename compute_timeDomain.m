function res = compute_timeDomain(acc_AP, acc_ML, sampleFreq)

testDuration = length(acc_AP)/sampleFreq; %  eq. A6 Martinez Mendez et al 2011

% acc intensity in the horizontal plane
%This line could be ommited as it is repeated in  the line for aam. Or use
%it in the computation of aam
A=sqrt(power(acc_AP,2)+power(acc_ML,2)); % acc_AP e acc_ML sono array Nx1

% Average acceleration magnitude
aam = sum(sqrt(power(acc_AP,2)+power(acc_ML,2)))/length(acc_AP); % eq. A2 Martinez Mendez et al 2011 (needed for computing jerk)

% output #1, mean distance from the centre of the trajectory 
dist = sqrt(power(acc_AP,2)+power(acc_ML,2)); % eq. 3 Prieto et al 1996
mDist = sum(dist)/length(dist); % eq. 4 Prieto et al 1996, in m/s^2

% output #2, compute root mean square
rms = sqrt(sum(A.^2)/length(A)); % eq. A3 Martinez Mendez et al 2011, in m/s^2

% output #3, compute range
%********
%I added the absolute computation of the substration. As in Eq. 9. 
range =abs(max(A)-min(A)); % eq. A5 Martinez Mendez et al 2011, in m/s^2

% output # ??, compute the length of the trajectory (dismissed, as it
% perfectly correlates with jerk)
pathLength = sum(sqrt(power(diff(acc_AP),2)+power(diff(acc_ML),2))); % eq. 8 Prieto et al 1996, in m/s^2

% output #5, smoothness of the trajectory (jerk)
% Rigoberto's 2012 ("The average change in acceleration", in m/s^3)
% jerk = sum(sqrt(power(diff(acc_AP),2)+power(diff(acc_ML),2)))/testDuration; % eq. A8 Martinez Mendez et al 2011, in m^2/s^5

% Mancini's 2012 ("Jerk was calculated as the time integral of the squared
% derivative of acceleration in both AP and ML directions.", in m^2/s^5)
dt = 1/sampleFreq;
daML = gradient(aML, dt);
daAP = gradient(aAP, dt);
jerk = sum(daML.^2 + daAP.^2) * dt; % no need to normalize with respect to time as test duration is fixed


% output #6, the mean frequency (MFREQ) is the rotational frequency, in revolutions per second or Hz, of the COP if it had traveled
% the total excursions around a circle with a radius of the mean distance
meanFreq = pathLength/(2*pi*mDist*testDuration); % eq. 20 Prieto 1996, in Hz
% which is equal to meanFreq = jerk/dot(2*pi,aam); % eq. A12 Martinez Mendez et al 2011

% output #7, sway area eq. 19 Prieto et al 1996
    for i = 1:length(acc_AP)-1
        swayArea_temp(i) = abs((acc_AP(i+1)*acc_ML(i))-(acc_AP(i)*acc_ML(i+1)));
    end
    swayArea = sum(swayArea_temp)/(2*testDuration); % in m^2/s^5

% output #8, 95% confidence ellipse area
eq17 = sum(acc_AP.*acc_ML)./length(acc_AP); % eq. 17 Prieto 1996
temp = power(std(acc_AP),2)*power(std(acc_ML),2);
ellipseArea = 2*pi*3*sqrt(temp-power(eq17,2)); % eq. 18 Prieto 1996, in m^2/s^4

res =  [jerk, mDist, rms, range, swayArea, ellipseArea, meanFreq];
