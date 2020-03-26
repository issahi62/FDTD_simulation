%******************* 
%% Matlab codes for 1 slit eperiement
% *******************

%% INITIALISE 
clear  
close all; 
clc


%% DASHBOARD 
fig1 = figure('color' ,'w'); 
set(fig1,'NumberTitle', 'off',  'Name', 'FDTD analysis');
% Constants 
e0 = 8.852e-12; 
slitn = 1; %% Make changes of slits here 1-diffraction 2- young's double slit, 3-grating like
%% Source Parameters
c = 3e8;  % Speed of EM wave
freq = 3e9; % Frequency of EM wave = 3GHz
lambda = c/freq; % Wavelength of EM wave 

%% Grid 
a =2; % xlength of the box = 2m 
b =4; % ylength of the box = 4m

dx = lambda/5; %Mesh size along X-direction
dy = lambda/5; %Mesh size along Y-direction

x = 0:dx:a; 
Nx = length(x); 

y = 0:dy:b; 
Ny = length(y);

R = 0.5; %CFL stability condition 0<x<.5-.7

time_duration = 200; 

dt= R*dx/c; 
tsteps = time_duration; 

%% Compute PML Parameters
% -----------------------------|---|
% PML along the x-axis   30    | P |
% -----------------------------|---|
%                              | M | 
% 0                            | L |
% source                       |30 |
%                              | Y |  
%                              |   |
% -----------------------------|---|
% PML along the x-axis  30     |   |
% -----------------------------|---|

Nx2 = 2*Nx; 
Ny2 = 2*Ny; 
NPML = [0, 20, 20, 20]; 

%%%******* Forming Matrix ********
sigx = zeros(Nx2, Ny2); 
for i = 1:2*NPML(1)
    i1 = 2*NPML(1) -i+1;
    sigx(i1, :) = (.5*e0/dt)*(i/2/NPML(1))^3; 
end

for i = 1:2*NPML(2)
    i1 = Nx2 - 2*NPML(2) +i;
    sigx(i1, :) = (.5*e0/dt)*(i/2/NPML(2))^3; 
end

sigy = zeros(Nx2, Ny2); 
for j = 1:2*NPML(3)
    j1 = 2*NPML(3) -j+1;
    sigy(:, j1) = (.5*e0/dt)*(j/2/NPML(3))^3; 
end

for j = 1:2*NPML(4)
    j1 = Ny2 - 2*NPML(4) +j;
    sigy(:, j1) = (.5*e0/dt)*(j/2/NPML(4))^3; 
end
%%% ********************** end **********%%%
%pcolor((sigx+sigy)')
URxx =1; 
URyy =1; 
ERzz =1; 
%%%************************************
%% COMPUTE UPDATE COEFFICIENTS
%%%************************************
sigHx = sigx(1:2:Nx2, 2:2:Ny2); 
sigHy = sigy(1:2:Nx2, 2:2:Ny2); 
mHx0  = (1/dt) + sigHy/(2*e0); 
mHx1  = ((1/dt) - sigHy/(2*e0))./mHx0; 
mHx2  = -c./URxx./mHx0; 
mHx3  = -(c*dt/e0) * sigHx./URxx./mHx0; 
sigHx = sigx(2:2:Nx2, 1:2:Ny2); 
sigHy = sigy(2:2:Nx2, 1:2:Ny2);
mHy0  = (1/dt) + sigHx/(2*e0); 
mHy1  = ((1/dt) - sigHx/(2*e0))./mHy0; 
mHy2  = -c./URyy./mHy0; 
mHy3  = -(c*dt/e0) * sigHy./URyy./mHy0; 

sigDx = sigx(1:2:Nx2, 1:2:Ny2); 
sigDy = sigy(1:2:Nx2, 1:2:Ny2);
mDz0  = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2); 
mDz1  = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2); 
mDz1  = mDz1./mDz0; 
mDz2  = c./mDz0; 
mDz4  = -(dt/e0^2) * sigDx.*sigDy./mDz0; 
mEz1  = 1/ERzz; 


%% Initialize field matrices
Dz = zeros(Nx, Ny); 
Ez = zeros(Nx, Ny); 
Hx = zeros(Nx, Ny); 
Hy = zeros(Nx, Ny); 

%% Curl matrices 
 CEx = zeros(Nx, Ny); 
 CEy = zeros(Nx, Ny); 
 CHz = zeros(Nx, Ny); 
 
 %% Integration matricss 
 ICEx = zeros(Nx, Ny); 
 ICEy = zeros(Nx, Ny); 
 IDz = zeros(Nx, Ny); 

 for t = 1:tsteps
     % Defining slit
     if slitn == 1
     Ez(2, 1:(floor(Ny/2)-2)) = 0; 
     Ez(2,(floor(Ny/2)+2):Ny) = 0; 
     end
     if slitn == 2
         Ez(2, 1:(floor(Ny/2)-8)) = 0; 
         Ez(2, (floor(Ny/2)-4):(floor(Ny/2)+4)) =0;
         Ez(2,(floor(Ny/2)+8):Ny) = 0;
     end 
     
     if slitn == 3
         Ez(2, 1:(floor(Ny/2)-9)) = 0; 
         Ez(2, (floor(Ny/2)-5):(floor(Ny/2)-2)) =0;
         Ez(2, (floor(Ny/2)+2):(floor(Ny/2)+5)) = 0;
         Ez(2,(floor(Ny/2)+9):Ny) = 0;
     end 
     % calculating CEx
     for i =1:Nx 
         for j =1:Ny-1
             CEx(i,j) = (Ez(i,j+1)-Ez(i,j))/dy; 
         end 
         CEx(i, Ny) = (0-Ez(i,j))/dy; % For tackling Y-high side
     end 
     
      % calculating CEy
     for j =1:Ny 
         for i =1:Nx-1
             CEy(i,j) = -(Ez(i+1,j)-Ez(i,j))/dx; 
         end 
         CEx(Nx, j) = -(0-Ez(Nx,j))/dx; % For tackling X-high side
     end 
     
     %Update H integration 
     ICEx = ICEx+ CEx; 
     ICEy = ICEy + CEy; 
     
     %Update H fields 
     Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx; 
     Hy = mHx1.*Hy + mHy2.*CEy + mHy3.*ICEy; 
     
     % Compute CHz 
     % Curl equations automatically include PEC BC
     CHz(1,1) = ((Hy(1,1)-0)/dx)-((Hx(1,1)-0)/dy); 
     
     for i = 2:Nx
         CHz(i, 1) = ((Hy(i,1)-Hy(i-1,1))/dx)-((Hx(i,1)-0)); %Note
     end
     
     for j = 2:Ny
         CHz(1, j) = ((Hy(1,j)-0)/dx)-((Hx(1,j)-Hx(1,j-1)));
         for i=2:Nx
           CHz(i, j) = ((Hy(i,j)-Hy(i-1,j))/dx)-((Hx(i,j)-Hx(i,j-1))/dy); 
         end 
     end 
     
     % Update D integration 
     IDz = IDz+Dz; 
     Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz; 
     
     % for sine wave mode 
     source = sin(((2*pi*(freq)*t*dt)));
     
     % Assigning source on y-low edge
     for j= 1:Ny 
         Dz(1, j) = sin((pi*(j-1)*dy)/b)*source; % source from -y to y
     end
     
     % Update Ez field 
     Ez = mEz1.*Dz; 
     Ez(:, 1) =0; 
     
     % Plotting Ez-wave
     [yy, xx] = meshgrid(y,x); 
     pcolor(x, y, real(Ez)'); 
     shading interp; 
     xlabel('X \rightarrow'); 
     ylabel('\leftarrow Y'); 
     titlestring =['fontsize{20} Slit at timestep = ', num2str(t)]; 
     title(titlestring, 'color', 'k'); 
     axis([0 a 0 b -1 1]); 
     view(0,90); 
     caxis([-1, 1]); 
     colormap(hsv); 
     colorbar; 
     getframe; 
 end
 
