% FB 12/20/2025
tempdata=ds5strt;
    
irts=tempdata;
irts=sort(irts)
%irts(find(irts>60))=[]
% Now we will look at the steady state data and try to estimate the
% parameter of a exponential+wald distributions that best describes the
% data

xdata=sort(irts)'; % sort the data;

% Get the xdata (IRTs) and ydata (cumulative responses)
n = length( xdata );
ydata = (1:n) / n;   % This is to done to express the data as a cumulative density

% To show the fit as pdf
% We will first use the least squares method to fit the data
LSexpwaldcdf = @(p, x)   p(1)*expcdf(x,p(2)) + (1-p(1))*waldcdf(x,p(3),p(4)); % Define the mixture distribution as a cumulative density function
LSexpwaldpdf = @(p, x)   p(1)*exppdf(x,p(2)) + (1-p(1))*waldpdf(x,p(3),p(4)); % Define the mixture distribution as a probability density function

% Get bounds on parameters, and initial estimates
lower  = [0 zeros(1,3)];        % Lower bound for all of our parameters will be 0
upper  = [1 inf(1,3)];          % Upper bound for all of our parameters will be plus infinity, except for the p(1), which has an analytical upper limit of 1

x0 = [0.05  0.08   40   323];   % Stole this from the first time I called lsqcurvefit

% And first fit the cdf with least squares
lsfit   = lsqcurvefit( LSexpwaldcdf, x0,xdata,ydata,lower,upper);%,optimset('MaxIter',1e6))

% Now we will do the fitting using maximum likelihood estimates
MLexpwaldpdf = @(x, p,lambdae,mu,lambdaw)   p*exppdf(x,lambdae) + (1-p)*waldpdf(x,mu,lambdaw); 
MLexpwaldcdf = @(x, p,lambdae,mu,lambdaw)   p*expcdf(x,lambdae) + (1-p)*waldcdf(x,mu,lambdaw);

% Now refind the starting estimates (to be what the LS fitting found), and
% refit using MLE (which is the more appropriate fitting method)
x0 = lsfit;
mlfit = mle(xdata,'pdf',MLexpwaldpdf,'cdf',MLexpwaldcdf,'start',x0,'lowerbound',lower,'upperbound',upper,'optimfun','fmincon');%, 'options', statset('MaxIter',1e6, 'MaxFunEvals',1e6 ));

logL = sum(log(MLexpwaldpdf(xdata, mlfit(1), mlfit(2), mlfit(3), mlfit(4))));
BIC = -2 * logL + length(mlfit) * log(length(xdata));

% Get predicted values for the cdf
xvals = 0 : .1 : max( xdata );  
yls_cdf = LSexpwaldcdf( lsfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters of expwaldcdf (USE LSEs)
yml_cdf = LSexpwaldcdf( mlfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters of expwaldcdf (USE MLEs)



% Plot them both to see how they look...

% Plot CDFs along with data expressed as cumulative density
figure;
subplot(2,1,1)
plot( xdata, ydata, 'b-' ); % data cdf;
hold on;
plot( xvals, yls_cdf, 'g-' ); % from Least square method fits
plot( xvals, yml_cdf, 'r-' ); % From Maximum likelihood method fits

legend('Data','LS fit','ML fit');

% Plot PDF
subplot(2,1,2)
dx = .5;
xbins = 0 : dx : 60;
hx = hist( xdata, xbins );
hx = hx / sum( hx ); % normalize by the number of occurrences

% Get predicted values for the pdf
yls_pdf = dx *LSexpwaldpdf( lsfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters  (USE LSEs)
yml_pdf = dx * LSexpwaldpdf( mlfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters  (USE MLEs)

bar( xbins, hx, 'b' );
hold on;
plot( xvals, yls_pdf, 'g-', 'LineWidth', 2 );
plot( xvals, yml_pdf, 'r-', 'LineWidth', 2 );
xlim( [0 Schedule*4] ); % zoom 
legend('Data','LS fit','ML fit');

% Get the xdata (IRTs) and ydata (cumulative responses)
% xdata = x_Rat{10}; % from my data structure
xdata=sort(xdata)
n = length( xdata );
ydata = (1:n)' / n;
