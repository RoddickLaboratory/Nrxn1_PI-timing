% FB 12/20/2025
tempdata=wtstrt;
    
irts=tempdata;
irts=sort(irts)


xdata=sort(irts)'; % sort the data;

n = length( xdata );
ydata = (1:n) / n;   % This is to done to express the data as a cumulative density


LSwaldcdf = @(p, x)   waldcdf(x,p(1),p(2)); % Define the mixture distribution as a cumulative density function
LSwaldpdf = @(p, x)   waldpdf(x,p(1),p(2)); % Define the mixture distribution as a probability density function
lower  = [0 0];        % Lower bound for all of our parameters will be 0
upper  = [inf inf];          % Upper bound for all of our parameters will be plus infinity, except for the p(1), which has an analytical upper limit of 1

x0 = [40   323];   % Stole this from the first time I called lsqcurvefit

% And first fit the cdf with least squares
lsfit   = lsqcurvefit( LSwaldcdf, x0,xdata,ydata,lower,upper);%,optimset('MaxIter',1e6))

% Now we will do the fitting using maximum likelihood estimates
MLwaldpdf = @(x, mu,lambdaw)   waldpdf(x,mu,lambdaw); 
MLwaldcdf = @(x, mu,lambdaw)   waldcdf(x,mu,lambdaw);

% Now refind the starting estimates (to be what the LS fitting found), and
% refit using MLE (which is the more appropriate fitting method)
x0 = lsfit;
mlfit = mle(xdata,'pdf',MLwaldpdf,'cdf',MLwaldcdf,'start',x0,'lowerbound',lower,'upperbound',upper,'optimfun','fmincon');%, 'options', statset('MaxIter',1e6, 'MaxFunEvals',1e6 ));

logL = sum(log(MLwaldpdf(xdata, mlfit(1), mlfit(2))));
BIC = -2 * logL + length(mlfit) * log(length(xdata));

% Get predicted values for the cdf
xvals = 0 : .1 : max( xdata );  
yls_cdf = LSwaldcdf( lsfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters of expwaldcdf (USE LSEs)
yml_cdf = LSwaldcdf( mlfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters of expwaldcdf (USE MLEs)

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
yls_pdf = dx *LSwaldpdf( lsfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters  (USE LSEs)
yml_pdf = dx * LSwaldpdf( mlfit, xvals ); % Determine values estimated for different xs defined above given best fitting parameters  (USE MLEs)

bar( xbins, hx, 'b' );
hold on;
plot( xvals, yls_pdf, 'g-', 'LineWidth', 2 );
plot( xvals, yml_pdf, 'r-', 'LineWidth', 2 );

legend('Data','LS fit','ML fit');

% Get the xdata (IRTs) and ydata (cumulative responses)
% xdata = x_Rat{10}; % from my data structure
xdata=sort(xdata)
n = length( xdata );
ydata = (1:n)' / n;
