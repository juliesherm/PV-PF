%% Read in data
addpath(genpath("c:\Users\Julie\Documents\MATH\MCRN\Vortez\ESMDA Analysis\"));
%eqs = readmatrix("lameqpts_68.csv");
%eqs60 = readmatrix("lameqpts_60.csv");
%eqs40 = readmatrix("lameqpts_40.csv");

namestem = "1000mem64step\";

A = readmatrix(sprintf("%sH1lam1.gnudat",namestem), 'FileType','text');
h = A(:,1);
lam = A(:,2);
A = readmatrix(sprintf("%sH1U1.gnudat",namestem), 'FileType','text');
Upost = A(:,2);
A = readmatrix(sprintf("%sobs.gnudat",namestem), 'FileType','text');
Utrue = A(:,3);
daynum = A(:,2);

%% Fourier fit
f = fit((1:length(lam))',lam*1000,'fourier1');
%period
2*pi/f.w


%% Ruzmaikin fit
x = (1:length(lam))';
y = lam*1000;
% Define Start points, fit-function and fit curve
x0 = [0.001 0.001 0.03 264 1]; 
fitfun = fittype( @(lam0,lama,epsil,doy,hmm,x) lam0+lama*sin((x-doy)*2*pi/365.25)+epsil*lam0*sin((x-hmm)*2*pi/(365.25*11)).^2 );
[fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);
% Save the coeffiecient values for a,b,c and d in a vector
coeffvals = coeffvalues(fitted_curve);

%% Plot lambda fit
figure
scatter(x, y, 1, 'black')
hold on
plot(x,fitted_curve(x), 'LineWidth',1)
plot(x,f(x))
ylabel("$\Lambda(t)$ (m/s/km)", 'Interpreter','latex')
xticks([1, 1+365, 2+2*365, 2+3*365, 2+4*365, 2+5*365, 3+6*365, 3+7*365, 3+8*365, 3+9*365, 4+10*365,...
    4+11*365, 4+12*365, 4+13*365, 5+14*365, 5+15*365, 5+16*365, 5+17*365, 6+18*365, 6+19*365])
xticklabels({' ','Jan 1, 2000',' ','Jan 1, 2002','','Jan 1, 2004',' ','Jan 1, 2006',' ',...
'Jan 1, 2008','','Jan 1, 2010',' ','Jan 1, 2012',' ','Jan 1, 2014','','Jan 1, 2016',...
' ','Jan 1, 2018'})
xtickangle(30)
xlim([0 7305])
legend("Posterior","Ruzmaikin Fit","Fourier Fit")

hold off

%% Look at a specific example of polar vortex (Jan 2-8, 2014, see Wiki)
%Wiki says Jan 2 - April 10 though..
figure
scatter(x, y, 10, 'black','filled')
hold on
plot(x,fitted_curve(x), 'LineWidth',1)
plot(x,f(x))
title("Jan 2014")
xticks([3+15*365,4+15*365, 5+15*365, 6+15*365,7+15*365,8+15*365,9+15*365,10+15*365,11+15*365,12+15*365,13+15*365,14+15*365])
xticklabels({'', 'Dec 31, 2013','','Jan 2, 2014','','Jan 4, 2014','','Jan 6, 2014','',...
    'Jan 8, 2014','','Jan 10, 2014'})
xtickangle(30)
xline( 5+15*365 + 1, 'Color', "blue") %add a line corresponding to 2014 vortex (Jan 2, 2014)
xline( 5+15*365 + 7, 'Color', "blue") %add a line corresponding to 2014 vortex (Jan 8, 2014)
xlim([5+15*365 - 15, 5+15*365 + 50])

hold off
%% Look at SSW
%there also appears to be a wiki list of sudden stratospheric warming
%events for different reanalysis products including ERA-Interim...
%26 Feb, 1999 %57
%20 March, 2000 %445
%11 Feb, 2001 %773
%30 Dec, 2001 %1095
%18 Jan, 2003 %1479
%5 Jan, 2004 %1831
%21 Jan, 2006 %2578
%24 Feb, 2007 %2977
%22 Feb, 2008 %3340
%24 Jan, 2009 %3677
%9 Feb, 2010 %4058
%24 Mar, 2010 %4101
%7 Jan, 2013 %5121
%12 Feb, 2018 %6982
dayidxs = [57, 445, 773, 1095, 1479, 1831, 2578, 2977, 3340, 3677, 4058, 4101,5121,6982, 6982-365,6982-2*365 ]; %16 here, last two are late Jan in a non-SSW year
figure
for i=1:length(dayidxs)
    dayidx = dayidxs(i);
    subplot(4,4,i)
    scatter(x, y, 10, 'black','filled')
%    scatter(x(3:end), Utrue, 10, 'black','filled')
%    scatter(x, Upost, 10, 'black','filled')
    hold on
    plot(x,fitted_curve(x), 'LineWidth',1)
    plot(x,f(x))
    xticks([dayidx-14, dayidx-7, dayidx, dayidx+7, dayidx+14])
    xticklabels({'-2 weeks', '-1 week','SSW','+1 week','+2 weeks'})
    xtickangle(30)
    xline(dayidx)
    xlim([dayidx - 15, dayidx+15])
end


%% Look at U-lambda phase space (like ruz fig 5)
figure
scatter(eqs60(:,1),eqs60(:,2),1,'black')
hold on
plot(lam(1:600), Upost(1:600))

ylim([-40,60])

%% Look at average lambda over years
lamav = zeros(20,1);
lyc = 0;
for i=1:20
    if mod(i+1998,4)==0
        lyc = lyc+1;
    end
    lamav(i) = mean(lam((1+(i-1)*365+lyc):(i*365+lyc)));
end

figure 
scatter(1999:2018,lamav*1000, 10, "black","filled")
xlabel("Year")
ylabel("Average $\Lambda$ (m/s/km)","Interpreter","latex")
% check for trend
lm = fitlm(1:20, lamav*1000)

%% Just for fun, correlation
%-0.62 between observed wind speed and area of marginal ice zone
%-0.63 with posterior wind speed and area
%-0.47 with posterior lambda and area
%0.68 between observed wind speed and lambda
%0.63 between posterior wind speed and lambda


% hdims = readmatrix('areas_IPT2_1980_2021.csv');
% %only get hdims between 1999 and 2018
% %1 is jan 1, 1980, 
% hdims = hdims(6941:14245);
% Utrue = Utrue(~isnan(hdims));
% lam = lam(~isnan(hdims));
% Upost = Upost(~isnan(hdims));
% h = h(~isnan(hdims));
% hdims = hdims(~isnan(hdims));
% corr(hdims, Utrue)
% corr(hdims, lam)
% figure
% scatter(hdims, Utrue)