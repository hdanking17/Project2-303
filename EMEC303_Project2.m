%EMEC303_Project2

clear; clc

% Inputs
Nx=25; % Grid points in x
Ny=20; % Grid points in y
dt=0.005; % Time step
Ndays=365; % Number of days to simulate
D=10; % Diffusion coefficient (miles^2/day)
plot_freq=100; % Plot frequency (steps)

% Pollutant source
yo=???
xo=???
A=???
sigma=???
source=@(x,y) ???

% Load wind data
data=csvread(`wind.csv',1.0);
wind_mean=data(:,1); % Mean wind speed (mph)
wind_dir =data(:,3); % Wind direction in degrees from north

% Create grid
lat(1)=45.5; 
lon(1)=-111.6; % Extents
lat(2)=46.0; 
lon(2)=-110.8;
x=linspace(lon(1),lon(2),Nx); % x & y values
y=linspace(lat(1),lat(2),Ny);
dx=latlon2dist(mean(lat),lon(1),mean(lat),lon(2))/Nx; % dx in miles
dy=latlon2dist(lat(1),mean(lon),lat(2),mean(lon))/Ny; % dy in miles

% Initialize variables and arrays
C=zeros(Nx,Ny); % Concentration
t=0; % Time
plot_count=0; % Counter to only plot every plot_freq steps

%% Create plot with map
figure(1); clf(1)
% Get current axes
ax1=gca;
% Plot google map
axis([x(1),x(Nx),y(1),y(Ny)]);
plot_google_map('APIKey','AIzaSyDUz4oSBuVc8LvjAqa26WARGJR9jw4-Ghk')
plot_google_map('AutoAxis',1);
% Plot source location
plot(xo,yo,'r.','Markersize',50);
% Create second axes for contour plot
ax2=axes();set(ax2,'visible','off')
% Add contour plot with colorbar
hold on
[c,hc]=contour(x,y,C');
set(hc,'XData',x);
set(hc,'YData',y);
hcb=colorbar;
ylabel(hcb,'Concentration','Fontsize',15)
% Add axes labels
xlabel(ax1,'Longitude')
ylabel(ax1,'Latitude')
set(ax1,'Fontsize',15)
% Scale axes to match
P1=get(ax1,'Position');
P2=get(ax2,'Position');
pos(1)=max(P1(1),P2(1));
pos(2)=max(P1(2),P2(2));
pos(3)=min(P1(3),P2(3));
pos(4)=min(P1(4),P2(4));
set(ax1,'Position',pos)
set(ax2,'Position',pos)
linkaxes([ax1,ax2])

%% Loop over time
Nt=(Ndays-1)/dt;
for n=1:Nt
    
    % Update time
    t=t+???;
    
    % Compute velocity interpolated to this time
    day =1+floor(t);
    frac=1+t-day;
    % Wind speed and direction
    mywind_spd=(1-frac)*wind_mean(day)+frac*wind_mean(day+1);
    mywind_dir=(1-frac)*wind_dir (day)+frac*wind_dir (day+1);
    % x and y components of wind velocity (miles/day)
    u=-24*mywind_spd*sind(mywind_dir);
    v=-24*mywind_spd*cosd(mywind_dir);
    
    % Update concentration by solving Eq. 2
    ???
    ???
    ???
    
    % Apply Neumann boundary conditions (zero slope)
    C( 1,:)=C( 2,:);
    C(Nx,:)=C(Nx-1,:);
    C(:, 1)=C(:, 2);
    C(:,Ny)=C(:,Ny-1);
    
    % Update concentration plot
    plot_count=plot_count+1;
    if plot_count==plot_freq
        set(hc,'ZData',C');
        titlestr=sprintf('Time=%5.f days; Wind=%5.2f mph, %5.2f deg.',...
            t,mywind_spd,mywind_dir);
        title(ax1,titlestr)
        drawnow
        plot_count=0;
    end
end