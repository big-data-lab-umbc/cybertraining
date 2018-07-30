%% MATLAB SCRIPT TO READ IN MODIS DATA and Generate 2D Surrogate cloud 
%
%   Written by Scott Hottovy May 3, 2018
%
%

% Clear and close all figures
clear;
close all;
clc;

% Constant to save .nc files. Set to 1 to save .nc files and 0 to not save.
% If Save_NC=1 and filenames are already present in the directory, then
% there will be an error. 
Save_NC = 0;

% Format the figure font size. 
figfontsize=24;


%% Read in file 

% Read in file. Change the path correct directory where MODIS data is. 
file_name = 'MODISCOT/MOD06_L2.A2000306.2325._cot_0900_0250.txt';
fileID = fopen(file_name);
C = textscan(fileID,'%f','Delimiter',',');

% Make Cloud optical thickness (COT) into a matrix
COT = reshape(C{1},50,50);

% Take transpose 
COT= COT';



% Correct for negative values to cloud optical thickness of the average of the rest. Change to
% something meaning ful in the future. 
COT(COT<0)=sum(sum(COT>0))/(size(COT,1)*size(COT,2));

% Find the mean COT.
Mean = sum(sum(COT))/(size(COT,1)*size(COT,2));

% Read in Cloud effective radius (CER)
file_name = 'MODISCOT/MOD06_L2.A2000306.2325._cer_0900_0250.txt';
fileID = fopen(file_name);
C = textscan(fileID,'%f','Delimiter',',');

% Make COT into a matrix
CER = reshape(C{1},50,50);

% Format matrix
CER= CER';


% Define the x and y axes. 
X = 1:1:size(COT,1);
Y = 1:1:size(COT,2);

% Create a plot of COT
CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);
% plot the cloud field from MODIS
%subplot(2,2,1)
h=pcolor(X,Y,COT);
colormap(gray)
axis square
caxis([0,20])
set(h, 'EdgeColor', 'none');
title('Obs. COT','Fontsize',figfontsize);
xlabel('x [km]','Fontsize',figfontsize);
ylabel('y [km]','Fontsize',figfontsize);
box on
hcolorbar=colorbar('SouthOutside');
 set(gca,'fontsize',figfontsize);

 % Save plot as .jpg
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/MODISCloud14COT.jpg')

% Create plot of CER
CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);
 %Define cloud thickness H
%subplot(2,2,2)
h=pcolor(X,Y,CER);
colormap(jet)
axis square
caxis([0,17])
set(h, 'EdgeColor', 'none');
title('Obs. CER','Fontsize',figfontsize);
xlabel('x [km]','Fontsize',figfontsize);
ylabel('y [km]','Fontsize',figfontsize);
box on
hcolorbar=colorbar('SouthOutside');
 set(gca,'fontsize',figfontsize);

% Save plot
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/MODISCloud14CER.jpg')

% Compute the cloud top height (CTH) as a function of COT. See technical
% report for the specific calculation. 
H = (COT./(0.0145*(2.0*10^(-6))^(2/3)*(0.8*200*10^6)^(1/3))).^(3/5); % m

%Create plot of Cloud top height (CTH)
CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);
%Define cloud thickness H
%subplot(2,2,2)
surf(X,Y,H);
colormap(jet)
axis square
%set(h, 'EdgeColor', 'none');
title('Cloud Thickness','Fontsize',figfontsize);
xlabel('x [km]','Fontsize',figfontsize);
ylabel('y [km]','Fontsize',figfontsize);
zlabel('z [m]','Fontsize',figfontsize);
box on
 set(gca,'fontsize',figfontsize);

% Save plot
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/MODISCloud14CTH.jpg')

% Save the COT (tau), CTH (CHT), and CER into a .nc file. 
if Save_NC == 1
    % Save to NC file
    str=strcat('MOD06_L2.A2000306.2325._cot_0900_0250COTCTHCER.nc');
        nccreate(str,'tau',...
             'Dimensions', {'x [km]', 50, 'y [km]', 50});
        ncwrite(str,'tau',COT);
        nccreate(str,'CHT',...
             'Dimensions', {'x [km]', 50, 'y [km]', 50});
        ncwrite(str,'CHT',H);
        nccreate(str,'CER',...
             'Dimensions', {'x [km]', 50, 'y [km]', 50});
        ncwrite(str,'CER',CER);
        ncwriteatt(str,'/','creation_time',datestr(now));

end




%% Compute Power spectrum 

% Define the size of the COT matrix. Assumes that it is a square matrix. 
noValues     = (length(COT));
halfNoValues = floor(noValues/2); %half the values

% Define the power by calling the periodogram2 function. 
powerObs = periodogram2(COT);

% Plot the log of the power. 
CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);
 
kx = -25:1:24;
ky = -25:1:24;
h = pcolor(kx,ky,log(powerObs));
title('PSD of Obs. Cloud Field');
colormap(jet)
xlabel('Wave Number in x [1/km]');
ylabel('Wave Number in y [1/km]');
set(h, 'EdgeColor', 'none');
box on
hcolorbar=colorbar('SouthOutside');
 caxis([2,20])
 set(gca,'fontsize',figfontsize);

% Save figure
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/MODISPower14.jpg')


%% Create surrogate Model

% malloc
spectrum  = zeros(noValues, noValues);

% Fill the Power spectrum with values. The first number, spectrum(1), is
% the mean of the time series (DC-component) and stays zero. Only the
% positive frequencies get a value here, later those values are mirrored to
% the negative frequencies.
spectrum = sqrt(powerObs); % Calculate the (magnitude of the) Fourier coefficients from the power coefficients.

phases   = rand(noValues,noValues)*2*pi-pi; % Calculate random phases for the complex Fourier coefficients with positve frequencies.

%Special Phases for halfNoValues. These should always be real. 
phases(halfNoValues+1,halfNoValues+1) = 0;
phases(1,halfNoValues+1) = 0;
phases(halfNoValues+1,1) = 0;
phases(halfNoValues+1,halfNoValues+1)=0;


spectrum = spectrum .* exp(1i * fftshift(phases)); % Multiply the coefficients with the random phases.

% Shift the spectrum to center it. 
spectrum = fftshift(spectrum);

% Set it to mean 0. 
spectrum(1,1)=0;

% symmetric conditions
spectrum(halfNoValues+2:end,1) = conj(flipud(spectrum(2:halfNoValues,1))); % 1st Col.
spectrum(1,halfNoValues+2:end) = conj(fliplr(spectrum(1,2:halfNoValues))); % 1st Row.
spectrum(halfNoValues+2:end,halfNoValues+1) = conj(flipud(spectrum(2:halfNoValues,halfNoValues+1))); % middle row
spectrum(halfNoValues+1,halfNoValues+2:end) = conj(fliplr(spectrum(halfNoValues+1,2:halfNoValues))); % middle col

spectrum(halfNoValues+2:end,halfNoValues+2:end) = conj(rot90(spectrum(2:halfNoValues,2:halfNoValues),2)); % lower right matrix
spectrum(halfNoValues+2:end,2:halfNoValues) = conj(rot90(spectrum(2:halfNoValues,halfNoValues+2:end),2)); % Lower left matrix

% Calculate the Fourier cloud. This is a mean zero cloud. 
y = ifft2(spectrum);

% include the mean
y = y+Mean;

% Throw out negative values
y(y<0) = 0;

% Compute the CTH as a function of COT. See technical report for the
% formula. 
yH = (y./(0.0145*(2.0*10^(-6))^(2/3)*(0.8*200*10^6)^(1/3))).^(3/5); % meters

% Compute the CER as a function of CTH. See technical report for the
% formula
yCER = 10^6*(0.0620)*(2.0*10^(-6))^(1/3)*(0.8*200*10^6)^(-1/3)*yH.^(1/3); %nanometers
 % m

% Save data to nc file
if Save_NC == 1
    str=strcat('MOD06_L2.A2000306.2325._cot_0900_0250SurrCOTCHTCER_4.nc');
      nccreate(str,'tau',...
          'Dimensions', {'x [km]', 50, 'y [km]', 50});
        ncwrite(str,'tau',y);
        nccreate(str,'CHT',...
             'Dimensions', {'x [km]', 50, 'y [km]', 50});
        ncwrite(str,'CHT',real(yH));
         nccreate(str,'CER',...
             'Dimensions', {'x [km]', 50, 'y [km]', 50});
        ncwrite(str,'CER',yCER);
        ncwriteatt(str,'/','creation_time',datestr(now));
end


%% plot cloud field

%plot surrogate COT
CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);
h=pcolor(X,Y,y);
colormap(gray)
axis square
%caxis([cmin,cmax])
caxis([0,20])
set(h, 'EdgeColor', 'none');
title('Surrogate COT');
xlabel('x [km]');
ylabel('y [km]');
box on
hcolorbar=colorbar('SouthOutside');
 set(gca,'fontsize',figfontsize);

% Save figure as .jpg
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/SurrCloud14COT.jpg')


%plot surrogate CER

CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);

h=pcolor(X,Y,yCER);
colormap(jet)
axis square
%caxis([cmin,cmax])
caxis([2,13])
set(h, 'EdgeColor', 'none');
title('Surrogate CER');
xlabel('x [km]');
ylabel('y [km]');
box on
hcolorbar=colorbar('SouthOutside');
 set(gca,'fontsize',figfontsize);

 % Save figure as .jpg
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/SurrCloud14CER.jpg')

%plot surrogate CTH
CF=figure;
 set(CF, 'Position', [100, 100, 1049, 895]);
%Define cloud thickness H
%subplot(2,2,2)
surf(X,Y,yH);
colormap(jet)
axis square
%set(h, 'EdgeColor', 'none');
title('Surrogate Cloud Thickness');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [m]');
box on
 set(gca,'fontsize',figfontsize);

% Save figure as .jpg
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/SurrCloud14CHT.jpg')


%% Compute and plot the Power spectrum as a check. 

 power = periodogram2(y);
 
 % Plot surrogate PSD. 
 figure
 kx = -25:1:24;
ky = -25:1:24;
 h =  pcolor(kx,ky,(log(power))); 
 title('PSD of Surrogate Cloud Field');
 colormap(jet)
 caxis([2,20])

xlabel('Wave Number in x [1/km]');
ylabel('Wave Number in y [1/km]');
set(h, 'EdgeColor', 'none');
hcolorbar=colorbar('SouthOutside');

box on

%save figure
saveas(gcf,'~/Dropbox/Math/CyberTraining_Class/Project/Code/2DClouds/SurrPower.jpg')
