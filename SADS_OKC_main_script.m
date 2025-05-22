%% this is the MAIN script to run the algorithm from computing time-domain stabilometric 
%% parameters from x-y-z calibrated acceleration signals measured by a triaxial 
%% accelerometer placed on the dumbbell.
%% The script reads Movella DOT file and prepares AP and ML accelerations 
%% for subsequent computation of time- and frequency-domain stabilometric parameters as performed by
%% the functions "compute_timeDomain.m" and "compute_freqDomain.m"

%% written by Pietro Picerno on August 2024 for submission to XXXX
clear all
close all
% read data file from Movella DOT placed on the dumbbell
filePath = 'C:\Users\pietr\OneDrive\MATLAB\UNIVERSITY\SADS-OKC\GitHub\sampleData\' % change accordingly to user's path
fileName = 'subject7_day1_session1_dominantArm_35%RL.csv'
data_dumbbell = readtable([filePath, fileName]);
% extract time and acceleration vector (see Movella DOT manual for CSV file
% column arrangement)
time = table2array(data_dumbbell(1:2194,2));
% compute sample frequency
sampleFreq = round(1/mean(diff(time/1000000))); % The unit of the SampleTimeFine is microsecond
s_acc_temp = table2array(data_dumbbell(:,3:5)); % acceleration espressed in the sensor reference system
% NOTE THAT:
% s_acc_temp = [x_acc, y_acc, z_acc] is the 3-components arrangement used by Movella DOT
% Set the dumbbell so that sensor's z-axis is pseudo-vertical and points upward, the x_acc is pseudo-AP and points to the little finger, y is pseudo-ML and points 
% medially to the subject's body when the right arm is used, while y_acc point laterally when the left arm is used

% test duration
duration = 30; % seconds

% set start of the test
start = length(s_acc_temp)-(sampleFreq*duration);
s_acc = s_acc_temp(start:end,:);

% re-arrange (or simply rename) according to Millecamps et al. 2015
quasiAP_acc = s_acc(:,1); % must point anteriorly 
quasiML_acc = s_acc(:,2); % must point to the left of the subject
quasiV_acc = s_acc(:,3); % must point upward

% reset the orientation of the phone using data during the 1st second of
% recording:

% this needs to be computed just in the 1st second after start (the sensor must be as still as possible in this second)
m_quasiAP_acc = mean(quasiAP_acc(1:sampleFreq)); % average value 1st second to compute attitude
m_quasiML_acc = mean(quasiML_acc(1:sampleFreq));
m_quasiV_acc = mean(quasiV_acc(1:sampleFreq));

% use roll and pitch (atan2) instead of sin as Moe-Nilssen
roll_0 = atan2(m_quasiML_acc,m_quasiV_acc); 
pitch_0 = atan2(m_quasiAP_acc,sqrt(m_quasiV_acc.*m_quasiV_acc+m_quasiML_acc.*m_quasiML_acc)); 

% sensor's initial tilt correction (Millecamps et al, 2015)
acc_AP=((quasiAP_acc*cos(pitch_0)) - (quasiV_acc.*sin(pitch_0))); % Equation 1 Millecamps et al. 2015, points backward (=is negative when the trunk leans forward)
acc_V_temp=(quasiAP_acc*sin(pitch_0)) + (quasiV_acc.*cos(pitch_0)); % Equation 2 Millecamps et al. 2015
acc_ML=((quasiML_acc*cos(roll_0)) - (acc_V_temp.*sin(roll_0))); % Equation 3 Millecamps et al. 2015, points to the right
acc_V=((quasiML_acc.*sin(roll_0)) + (acc_V_temp.*cos(roll_0))); % Equation 4 Millecamps et al. 2015, points upward

% detrend (for time-domain parameters only)
acc_AP_det = detrend(acc_AP,1);
acc_ML_det = detrend(acc_ML,1);

% filtering
nyq = 0.5 * sampleFreq;
cutoff = 14; % 
normal_cutoff = cutoff / nyq; % normalized cutoff frequency (see "butter" doc)
order = 4; % paper use order of 1. 
[b,a] = butter(order,normal_cutoff,'low');

% this goes for time-domain parameters calculation
acc_AP_det_filt_30s =filtfilt(b,a,acc_AP_det);
acc_ML_det_filt_30s =filtfilt(b,a,acc_ML_det);

% this goes for frequency-domain parameters calculation (detrend is
% built-in power spectral density analysis of Vieira and coll.)
acc_AP_filt_30s =filtfilt(b,a,acc_AP);
acc_ML_filt_30s =filtfilt(b,a,acc_ML);

% cut til the first 20 seconds 
acc_AP_det_filt_20s =  acc_AP_det_filt_30s(1:20*sampleFreq,:);
acc_ML_det_filt_20s =  acc_ML_det_filt_30s(1:20*sampleFreq,:);
acc_AP_filt_20s =  acc_AP_filt_30s(1:20*sampleFreq,:);
acc_ML_filt_20s =  acc_ML_filt_30s(1:20*sampleFreq,:);

% compute stabilometric parameters both for 20 s and 30 s test duration:
% (time-domain)
% Jerk	The average rate of change of the 2D horizontal acceleration signal (a measure of sway jerkiness), m2/s5
% mDist	Mean distance (deviation) from the center of the sway path (trace of the spaghetti plot), m/s2
% RMS	Root mean square of the 2D horizontal acceleration time series (quantifies the magnitude of the signal), m/s2
% Range	The maximum amplitude of the 2D horizontal acceleration signal, m/s2
% swayArea	This parameter approximates the area enclosed by the envelop of the sway path, m2/s5
% ellipseArea	The area of an ellipse enclosing all points of the sway path with 95% confidence, m2/s4
% swayFreq	Mean sway frequency (the number, per second, of loops that have to be run by the dumbbell to cover a trajectory equal to the total sway path), Hz
% (frequency-domain stabilometric parameters)
% TP    total power [Hz]
% F50	Median frequency (frequency band that contains up to 50% of the total spectrum), Hz
% F95	Frequency band that contains up to 95% of the total spectrum, Hz
% CF	centroidal Frequency or Spectral centroid (frequency at which spectral mass is concentrated), Hz
% FD    frequency dispersion, unitless

% compute time and freq domain parameters for 30 s
res_freq_30s = compute_freqDomain(acc_AP_filt_30s, acc_ML_filt_30s, sampleFreq, cutoff,1);
res_time_30s  = compute_timeDomain(acc_AP_det_filt_30s, acc_ML_det_filt_30s, sampleFreq);
% [jerk, mDist, rms, range, swayArea, ellipseArea, meanFreq, TP, F50, F95, CF, FD]
res_30s = [res_time_30s, res_freq_30s];

% compute time and freq domain parameters for 20 s
res_freq_20s = compute_freqDomain(acc_AP_filt_20s, acc_ML_filt_20s, sampleFreq, cutoff,0);
res_time_20s  = compute_timeDomain(acc_AP_det_filt_20s, acc_ML_det_filt_20s, sampleFreq);
% [jerk, mDist, rms, range, swayArea, ellipseArea, meanFreq, TP, F50, F95, CF, FD]
res_20s = [res_time_20s, res_freq_20s];
