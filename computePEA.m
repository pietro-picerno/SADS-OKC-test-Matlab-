function [a, mainSwayDir] = computePEA(x,y,P,plotOpt)

% this function computes 95% confidence ellipse area 
% a = PEA(x,y,P): plots a prediction ellipse of the center of pressure (COP)
% data separated into x and y components with probability value P 
% and calculates the prediction ellipse area (PEA)



%inputs:

% x,y: data (column vectors)

% P : a value of the interval (0,1)

% plotOpt = 1 to plot ellipse, 0 not to

%output:

% a : represents the area of the ellipse

% mainSwayDir = principal sway direction is defined as the angle (ranging between 0° and
% 90°) between the AP axis and the direction of the main eigenvector produced by the PCA

% type PEA without inputs to plot 95% prediction ellipse of exemplary data

% code available from: Schubert and Kirchner, Gait & Posture Volume 39, Issue 1, January 2014, Pages 518-522

% majorAxis and devAngle added by Pietro Picerno 2024 
%%

%exemplary data illustration

if nargin==0 %case of no input arguments

x = [2,4,6,2,3,4,2,3,3,4,5,5,0,8,3,7]; %exemplary x component

y = [1,2,3,4,3,3,1,1,2,1,4,5,0,0,1,1]; %exemplary y component

a = PEA(x,y,0.95);

axis([-4 11 -3 7])

title('95% PEA of exemplary data')

text(-3,6,'blue points = data')

text(5.5,2,'major axis')

text(4.3,0,'minor axis')

text(7,-2,['PEA: ' num2str(a)])

return %end of function

end

%%

%begin of function

% centra a traiettoria togliendo offset
% x=x-mean(x);
% y=y-mean(y);

chisquare = chi2inv(P,2); %inverse of the chi-square cumulative distribution function with 2 degrees of freedom at P

x = x(isfinite(x));

y = y(isfinite(y));

mx = mean(x);

my = mean(y);

[vec,val] = eig(cov(x,y)); %calculation of eigenvalues

a = pi*chisquare*prod(sqrt(svd(val))); %area calculation




%ellipse

N = 100; %fixed value (the higher the smoother the ellipse)

t = linspace(0,2*pi,N);

elip = sqrt(chisquare)*vec*sqrt(val)*[cos(t); sin(t)] + repmat([mx; my],1,N);

elip = elip';

%major and minor axes

ax1 = sqrt(chisquare)*vec*sqrt(val)*[-1,1; 0,0] + repmat([mx; my],1,2); % le coordinate x e y dei due punti che definiscono l'asse

ax2 = sqrt(chisquare)*vec*sqrt(val)*[0,0; -1,1] + repmat([mx; my],1,2);

ax_dat = [ax1'; NaN,NaN; ax2'];

if plotOpt==1 
% figure
% hold on

%COP data

plot(x,y,'b');grid on

line(ax_dat(:,1),ax_dat(:,2),'Color',[0 0 0],'LineWidth', 1);
line(elip(:,1),elip(:,2),'Color', [0 0 0], 'LineWidth', 1);


xlabel('medio-lateral comp. [m/s^2]')
ylabel('antero-posterior comp. [m/s^2]')

axis equal
end

% compute principal sway direction (Pietro):
ax2_temp = ax2(:,2)-ax2(:,1);
ax1_temp = ax1(:,2)-ax1(:,1);
if sqrt(power(ax2_temp(1),2)+power(ax2_temp(2),2)) > sqrt(power(ax1_temp(1),2)+power(ax1_temp(2),2)) % identify major axis
    majorAxis=ax2;
else
    majorAxis=ax1;
end

tmp=majorAxis(:,2)-majorAxis(:,1); % questo è l'asse
% Rocchi et al. 2004:
mainSwayDir = (atan2(tmp(2),tmp(1))).*180/pi; % angolo in gradi tra asse maggiore dell'ellisse e asse medio-laterale
% un angolo di 90° vuol dire che la direzione di sway è completamente antero-posteriore
% un angolo di 0° (o di 180°) la direzione di sway sarebbe completamente medio-laterale
% (0° punta verso destra, 180° punta verso sinistra) 

% attenzione perchè questo angolo va da 0 a 180°

return %end of function