%% Cleanup

clear; clc;
% close all;

%% Import data

%filename = "sf_output_000210000";
%filename = "sf_output_000574000";
filename = "sf_output_000644000";

fid = fopen(filename, 'r');
data = textscan(fid, '%f %f','delimiter', ' ', 'HeaderLines', 4);
fclose(fid);

%% Order data

sf = data{2};

N = length(sf); % must be even
sf = ifftshift(sf);
sf = sf(2:floor(N/2)+1);

N = length(sf);
k = (1:N)';

%% Fit data

klo = 2;
khi = 10;
rg = klo:khi;

% least squares fit
f = @(x) x.^-4;
yp = f(k);
a = sum(sf(rg))/sum(yp(rg));
yp = a*f(k);

% % power law fit
%p = polyfit(log(k(rg)),log(sf(rg)),1)
%yp = exp(p(2))*k.^p(1);

%% Plot data

rg = 1:N;

figure
hold on
grid on
loglog(k(rg),sf(rg),'b.')
loglog(k(rg),yp(rg),'r-')
ylim([0.5 1.5*max(sf(rg))])
ylabel('covariance')
xlabel('k')
legend('data','least-squares fit')