% Function to calculate spectral parameters
% written by Rigoberto Martinez Mendez

function res = compute_freqDomain(acc_AP, acc_ML, sampleFreq, cutFreq, plotOpt)

acc_hor = sqrt(power(acc_AP,2)+power(acc_ML,2)); 

% psd (Vieira et al. 2009)
nfft = 2^nextpow2(length(acc_hor));
[Pxx, f] = psdcop(acc_hor,sampleFreq,nfft,'hamming','linear'); 

% clip between 0.1 to cutoff frequency range
idx = (f >= 0.1) & (f <= cutFreq);

% Spectral moments
m0 = sum(Pxx(idx));      % Zeroth moment
m1 = sum(f(idx) .* Pxx(idx)); % First moment
m2 = sum(f(idx).^2 .* Pxx(idx)); % Second moment

% Total Power (equivalent to m0 in the paper)
TP = m0; 

% Cumulative power for median and 95% frequencies
cumpower = cumsum(Pxx(idx));

% Median Power Frequency
F50 = interp1(cumpower/cumpower(end), f(idx), 0.50, 'linear', 'extrap');

% 95% Power Frequency
F95 = interp1(cumpower/cumpower(end), f(idx), 0.95, 'linear', 'extrap');

% Centroidal Frequency (equivalent to sqrt(m2/m0) in the paper)
CF =sqrt(m2 / m0);

% Frequency Dispersion
FD = sqrt(1 - (m1^2 / (m0 * m2)));

res = [TP, F50, F95, CF, FD];

if plotOpt==1 
figure
plot(f,Pxx)
% plot(f,A)
title('Power spectrum (Welch mod.)')
xlabel('Frequency [Hz]')
ylabel('Power [(m/s^2)^2/Hz]')
xlim([0 sampleFreq/2])
hold
xline(F50,'--')
xline(F95,'-.')
xline(cutFreq) 
ylabel('Power Spectral Density [(m/s^2)^2/Hz]')
xlabel('frequenza [Hz]')
legend('frequency band','F50','F95','cut-off')
end

