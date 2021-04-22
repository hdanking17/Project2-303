%EMEC 303
%Project 2
%Diego Armstrong, Hannah King, Carter Storrusten

clear; clc

% Inputs
Nx=25; % Grid points in x
Ny=20; % Grid points in y
dt=0.005; % Time step
Ndays=365; % Number of days to simulate
D=10; % Diffusion coefficient (miles^2/day)
plot_freq=100; % Plot frequency (steps)

% Pollutant source
yo=45.817315;
xo=-111.073837;
A=1e4;
sigma=.01;
source=@(x,y) A*exp(-((x-xo)^2+(y-yo)^2)/(2*sigma^2));


% Load wind data
data=csvread('wind.csv',1.0);
wind_mean=data(:,1); % Mean wind speed (mph)
wind_dir =data(:,3); % Wind direction in degrees from nort

% Create grid
lat(1)=45.5; lon(1)=-111.6; % Extents
lat(2)=46.0; lon(2)=-110.8;
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
plot_google_map('APIKey','AIzaSyDUz4oSBuVc8LvjAqa26WARGJR9jw4-Ghk');
plot_google_map('AutoAxis',1);

% Plot source location
plot(xo,yo,'r.','Markersize',50);

% Create second axes for contour plot
ax2=axes();
set(ax2,'visible','off')

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
<<<<<<< HEAD
%update T for interior points
% Cnew=C;
% for i=2:Nx-1
%     Cnew(i)=C(i)+dt*(-Nx*C(i+1)-C(i-1))/(2*dx);
% end
%C=rand(size(C));
=======
>>>>>>> 4125fd4da8bb00ea6eaaca023a6528f56eca7ad1
Nt=(Ndays-1)/dt;
for n=1:Nt
    % Update time
    t=t+dt;
    
    % Compute velocity interpolated to this time
    day =1+floor(t);
    frac=1+t-day;
    
    % Wind speed and direction
    mywind_spd=(1-frac)*wind_mean(day)+frac*wind_mean(day+1);
    mywind_dir=(1-frac)*wind_dir (day)+frac*wind_dir (day+1);
    
    % x and y components of wind velocity (miles/day)
    u=-24*mywind_spd*sind(mywind_dir);
    v=-24*mywind_spd*cosd(mywind_dir);
    %%
    % Update concentration by solving Eq. 2
    %C=rand(size(C));
    C_star=C;
    for j=2:Ny-1
        %ADVECTION
        for i=2:Nx-1
            %dCdx & u
            if     u>=0
                dCdx=((C(i,j)-C(i-1,j))/dx);
            elseif u<0
                dCdx=((C(i+1,j)-C(i,j))/dx);
            end
            
            %dCdy & v
            if     v>=0
                dCdy=((C(i,j)-C(i,j-1))/dy);
            elseif v<0
                dCdy=((C(i,j+1)-C(i,j))/dy);
            end
        end
        %DIFFUSION
                dCdx2  = D*(C(i+1,j)-2*C(i,j)+C(i-1,j))/(dx^2);
                dCdy2  = D*(C(i,j+1)-2*C(i,j)+C(i,j-1))/(dy^2);
        %Right Hand Side of equation
                RHS   = -v*dCdy-u*dCdx+dCdx2+dCdy2+source(x(i),y(j));
        %TIME DERIVATIVE
                C_star(i,j)= C(i,j)+dt*RHS;

    end
        %update C
                C=C_star;
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