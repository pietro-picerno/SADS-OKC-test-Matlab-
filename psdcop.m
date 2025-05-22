
function [Scop,freq] = psdcop(cop,varargin)

% Scop = psdcop(copdata,fs,nfft,windowtype,trendtype)
% Calculate the power density spectrum of the center of pressure (cop)
% time series according to the welch method with 50% overlapping and
% dividing the signal into 7 segments.
%
% A low variance is assured due to the spectral averaging (periodogram)
% The overlapping improves spectral resolution, being 1/15 Hz to a 60s
% record
% 
% example: [Scop,freq] = psdcop(COG_ap,33,round(length(COG_ap)/2),'hamming','linear')

% NOTE: this code has been taken "as is" from
% Vieira TMM, Oliveira LF, Nadal J. Estimation procedures affect the center of pressure frequency analysis. Brazilian Journal of Medical and Biological Research. 2009;42(7):665-673. doi:10.1590/S0100-879X2009000700012

narginchk(1,5); % Checking for the number of input arguments
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - setting parameters for Scop estimation - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%

fs = varargin{1}; % Defining the sampling frequency
% K =  10; % Defining the number of segments (Taian's)
% NB: il numero di segmenti dipende dalla durata del test e dalla frequenza
% di campionamento e definisce la risoluzione della PSD. Ad esempio, per un
% segnale di 30 secondi, k=3 vuol dire che ho 3 segmento da 10
% secondi. Quindi scelgo K in base alla risoluzione che voglio (quanto
% smussato deve essere lo spettrogramma):
K = 3; 
N = length(cop); % Defining the size of the stabilometric record
if rem(N,K+1) 
    N=N-rem(N,K+1);
end % Dimensionating N for segmentation
cop = cop(1:N);
nwindow = (2*N/(K+1)); % Defining the length of each segment for 50% overlapping
% res = fs/nwindow;
window = feval(varargin{3},nwindow); % Evaluating the windowtype
% function (‘hann’, ‘hamming’, ...)
U = sum(window.^2); % Defining the normalization factor for the window applied
if ~isempty(varargin{2}) % if nfft is specified by the user
nfft = varargin{2}; % Then, the number of points for FFT is defined by the user
else
nfft = nwindow; % Otherwise, the number of points is equal to the window length
end

if fs/nwindow >= 0.4 % Checking for estimator with low resolution
disp('warning: frequencies below 0.4 Hz will not be represented')
end
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - writing segments in matrix form - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
window = window(:,ones(1,K)); % each row corresponds to the same window type and length
xcop = reshape(cop,0.5*nwindow,K+1); % dimensioning xcop for applying the window to each segment
xcop = [xcop(:,1:K);xcop(:,2:K+1)]; % each row corresponds to a segment of the cop time series
xcop = detrend(xcop,varargin{4}).*window; % removing each segment trend ( constant  or ‘linear’) and multiplying each segment by the desired window
xcop = detrend(xcop,varargin{4}); % Detrending the COP segments after multiplication
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - computation of the Scop and frequency vector - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
if rem(nfft,2) % if nfft is odd
index = 1:(nfft+1)/2; % frequency vector length is even
else % if nfft is even
index = (1:nfft/2+1)'; % frequency vector length is odd
end
COP = sum(abs(fft(xcop,nfft)).^2,2); % Applying FFT for each segment
COP = COP(index); % Getting the one sided spectrum
Scop = COP*(1/(fs*K*U)); % Computation of the cop power
% density spectrum
if nargout == 2 % If two output arguments are defined
freq = (index-1)*(fs/nfft); % Defining the frequency vector and passing it to the second output
end