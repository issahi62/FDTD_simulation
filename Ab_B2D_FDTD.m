%% Absorption boundary condition formulation 


%%%%%%%%%%% Ibrahim Issah %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  UEF school   %%%%%%%%%%%%%%%%%%


%% INITIALIZE 
clear; 
close all; 
clc; 
tsteps = 100;
%% DASHBOARD 

c    = 3e8; % speed of light
freq = 3e7; % frequency in Hz
epsr = 1;   % relative permittivity
lambda = c/freq/sqrt(epsr); % wavelength

xdim = 100; 
Rx = .5;
dx = lambda/10;  % x-position steps 
x = 0:dx:xdim; 
xsteps = length(x); 

ydim = xdim;
Ry =.5; 
dy = dx;  % y-position steps 
y = 0:dy:ydim; 
ysteps = length(y); 


%%%%%%%%%%% Total simulation time %%%%%%%%%%%%%%%%%%
R = 0.5; % CFL condition R<1
dt = R*dx/c; 
 

%% Source position 
xsource = floor(xsteps/2); 
ysource = floor(ysteps/2); 

%% BUILDING FDTD simulation 

%%% Field vectors 
Ez = zeros(ysteps, xsteps);
Hx = zeros(ysteps, xsteps);
Hy = zeros(ysteps, xsteps);

Ex2 = zeros(tsteps, xsteps); 
Exlast_1 = zeros(tsteps, xsteps);
Ey2 = zeros(tsteps, ysteps);
Eylast_1 = zeros(tsteps, ysteps);



for n = 1+ceil(1/min(Rx,Ry)):tsteps
    
    for l = 1:xsteps
        for m = 1:ysteps-1
            Hx(m, l) = Hx(m,l)- Ry*(Ez(m+1, l)-Ez(m,l));
        end
    end
    
     for m1 = 1:ysteps
        for l1 = 1:xsteps-1
            Hy(m1, l1) = Hy(m1,l1)+ Rx*(Ez(m1,l1+1)-Ez(m1,l1));
        end
     end
    
     for m2 = 2:ysteps
        for l2 = 2:xsteps
            Ez(m2, l2) = Ez(m2,l2) + (Rx*(Hy(m2,l2)-Hy(m2,l2-1))-Ry*(Hx(m2,l2)-Hx(m2-1,l2)));
        end
     end
    
     %% Absorbing boundary condition 
     
   % in x-direction
    Ex2(n,:) = Ez(:,2);  
    Ez(:,1)= Ex2(n-1/Rx,:);
    Exlast_1(n,:) = Ez(:,xsteps-1);
    Ez(:,xsteps) = Exlast_1(n-1/Rx,:);
    
% in y-direction
    Ey2(n,:) = Ez(2,:);  
    Ez(1,:)= Ey2(n-1/Ry,:);
    Eylast_1(n,:) = Ez(ysteps-1,:);
    Ez(ysteps,:) = Eylast_1(n-1/Ry,:);
     
     %% source input 
      source=sin(((2*pi*(freq)*n*dt)));
     
     %%% Guassian source 
      pulse = (10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*pi/20))/32;
      
      %%% Assigning source to the Ez component
      Ez(ysource, xsource) = pulse +source; %hard source like metal
   % Ez(ysource, xsource) = Ez(ysource, xsource) +source; % soft source 
   
   %% Plotting Ez-wave 
   
   surf(x, y, Ez); 
   shading interp 
   xlabel('X\rightarrow'); 
   ylabel('\leftarrow Y'); 
   caxis([-1, 1]); 
   %zlabel('E_z \rightarrow'); 
   titlestring =['Ez wave at timestep = ', num2str(n)]; 
   title(titlestring, 'color', 'k');  
   colorbar; 
   axis([0 xdim 0 ydim -1 1]);
   getframe;
end    




