function [a, r2] = linearRegression_khusus(x,y)
% linregr: linear regression curve fitting
% [a, r2] = linregr(x,y): Least squares fit of straight
% line to data by solving the normal equations
% input:
% x = independent variable
% y = dependent variable
% output:
% a = vector of slope, a(1), and intercept, a(2)
% r2 = coefficient of determination
n = length(x);
if length(y) ~= n, error('x and y must be same length'); end
x = x(:); y = y(:); % convert to column vectors
sx = sum(x); sy = sum(y);
sx2 = sum(x.*x); sxy = sum(x.*y); sy2 = sum(y.*y);
a(1) = (n*sxy - sx*sy)/(n*sx2 - sx^2);
a(2) = sy/n - a(1)*sx/n;
r2 = ((n*sxy - sx*sy)/sqrt(n*sx2 - sx^2)/sqrt(n*sy2 - sy^2))^2;
% create plot of data and best fit line
% xp = linspace(min(x),max(x),2);
% yp = a(1)*xp + a(2);
% plot(x,y,'o',xp,yp)
% grid on
alfa = exp(a(2));
beta = a(1);
xn = linspace(min(x),max(x),100);
yn2 = alfa.*exp(beta.*xn);
ytes = alfa.*exp(beta.*x);
sy = sqrt(sum((y - mean(ytes)).^2)/(n -1));
syx = sqrt(sum((y-ytes).^2)/(n-2));
fprintf('Eksponensial : y = a(e^bx)\n');
fprintf('Nilai Coefisien determinat (r2) = %f \n',r2);
% fprintf('Nilai Correlation coefficient (r)= %f \n',sqrt(r2));
% fprintf('Nilai Standar error (Syx) = %f \n',syx);
% fprintf('Nilai Standar deviasi(Sy) = %f \n',sy);
fprintf(' y = %fe^%fx \n',alfa,beta)