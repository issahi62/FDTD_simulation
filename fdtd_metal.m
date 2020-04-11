% Matlab script of 2D-FDTD (TE-mode) surrounded
% with PML boundary conditions and plane wave source incident
% on metallic grating structure
%************************************************************************
%
%************************************************************************
% This MATLAB M-file models a thin metal film with single aperture
% surrounded on either sides by grooves in Input Configuration (IC)
% if "ho" is set to zero, input and output corrugations (IOC) if "hi"
% and "ho" are non-zero and no corrugations (NC) if hi=ho=0.
%
% This program is inspired by a FDTD code written by Susan C. Hagness
% and found in the supplementary CD-ROM of the book:
% Taflove, A. and Hagness, S. C., ?Computational Electrodynamics: the
% finite-difference time-domain method, 2nd Ed.?, Boston, Artech
% House (2000), 852 pages.
%
%************************************************************************
clear
%************************************************************************
% Physical constants
%************************************************************************
cc = 2.99792458e8; % speed of light in free space
muz = 4.0*pi*1.0e-7; % permeability of free space
epsz = 1.0/(cc*cc*muz); % permittivity of free space
h = 6.6262e-34; % Planck's constant
hbar = h/2/pi; % Dirac's constant
e = 1.60219e-19; % electron's charge
% Type of source excitation (CW = 1; Gaussian modulated pulse = 2) :
SourceType = 1;
% Wavelenght of source excitation (m) :
lambda = 500e-9;

% Frequency of source excitation (Hz)
freq = cc/lambda;
% Angular frequency of source excitation (rad/s)
omegalight = 2.0*pi*freq;
%************************************************************************
% Material parameters
%************************************************************************
% The grid can be modeled with 3 types of material with their optical
% properties defined via vectors:
% vector = [(vacuum) (dielectric) (metal)];
eps = [1.0 2.25 1.0]; % relative permittivity (at infinite frequency)
sig = [0.0 0.0 0.0]; % conductivity [Siemens/m]
mur = [1.0 1.0 1.0]; % coefficient of permeability
sim = [0.0 0.0 0.0]; % coefficient of magnetic loss
% Select the material composition of the structure to be modeled later on
% where media = 1 is free space, 2 is dielectric and 3 is metal.
media = 3;
%************************************************************************
% Lorentz-Drude model parameters for the metal (Rakic et al. 1999.)
%************************************************************************
% LEGEND:
% material = 1=Ag ; 2=Al ; 3=Au ; 4=Cu ; 5=Cr ; 6=Ni ; 7=W
% omegap = Plasma frequency
% f = Oscillators' strenght
% Gamma = Damping frequency (1/s)
% omega = Oscillators' resonant frequency (rad/s)
% order = Number of resonant frequencies
material = 1; % Select the metal to be used in the structure
if material == 1
% LD parameters for Ag:
omegap = 9.01*e/hbar;
f = [0.845 0.065 0.124 0.011 0.840 5.646];
Gamma = [0.048 3.886 0.452 0.065 0.916 2.419]*e/hbar;
omega = [0.000 0.816 4.481 8.185 9.083 20.29]*e/hbar;
order = length(omega);
elseif material == 2
% LD parameters for Al:
omegap = 14.98*e/hbar;
f = [0.523 0.227 0.050 0.166 0.030];
Gamma = [0.047 0.333 0.312 1.351 3.382]*e/hbar;
omega = [0.000 0.162 1.544 1.808 3.473]*e/hbar;
order = length(omega);

elseif material == 3
% LD parameters for Au:
omegap = 9.03*e/hbar;
f = [0.760 0.024 0.010 0.071 0.601 4.384];
Gamma = [0.053 0.241 0.345 0.870 2.494 2.214]*e/hbar;
omega = [0.000 0.415 0.830 2.969 4.304 13.32]*e/hbar;
order = length(omega);
elseif material == 4
% LD parameters for Cu:
omegap = 10.83*e/hbar;
f = [0.575 0.061 0.104 0.723 0.638];
Gamma = [0.030 0.378 1.056 3.213 4.305]*e/hbar;
omega = [0.000 0.291 2.957 5.300 11.18]*e/hbar;
order = length(omega);
elseif material == 5
% LD parameters for Cr:
omegap = 10.75*e/hbar;
f = [0.168 0.151 0.150 1.149 0.825];
Gamma = [0.047 3.175 1.305 2.676 1.335]*e/hbar;
omega = [0.000 0.121 0.543 1.970 8.775]*e/hbar;
order = length(omega);
elseif material == 6
% LD parameters for Ni:
omegap = 15.92*e/hbar;
f = [0.096 0.100 0.135 0.106 0.729];
Gamma = [0.048 4.511 1.334 2.178 6.292]*e/hbar;
omega = [0.000 0.174 0.582 1.597 6.089]*e/hbar;
order = length(omega);
elseif material == 7
% LD parameters for W:
omegap = 13.22*e/hbar;
f = [0.206 0.054 0.166 0.706 2.590];
Gamma = [0.064 0.530 1.281 3.332 5.836]*e/hbar;
omega = [0.000 1.004 1.917 3.580 7.498]*e/hbar;
order = length(omega);
else
disp('ERROR!! Not a valid choice of material. CORRECT: 1<= material <=7');
end
%**********************************************************************
% FDTD simulation parameters
%************************************************************************
% Size of the square grid cell''s for spatial discretization (m)
dx = 10e-9;
% Courant number
Courant = 0.5;
% Timestep increment (seconds)
dt = (Courant*dx)/cc;
% Number of grid cells per lambda
N = round(lambda/dx);
% Time of one full oscillation cycle in # of timesteps
period = round((1/freq)/dt);
% relative permittivity of 1D-linear source grid
% (the value must match the background permittivity of the main grid)
epslin = 1.0;
% Size of the Scattered Field (SF) region in # cells
SF_width = 20;
%************************************************************************
% Structure's geometrical parameters
%************************************************************************
% N.B.
% In this section we define the parameters to model a metallic periodic
% grating on either sides of the film (input + output) and symmetrical
% about a central slit milled through the film.
% Number of input grating periods on either side (top + bottom) of the slit
Ni = 2;
% Slit width (# of cells)
a = 25;
% Grating parameters on the input side
di = 25; % Grating''s ridge size in (# of cells)
ei = 25; % Grating''s corrugation size in (# of cells)
LambdaGi = di+ei; % Grating''s period (# of cells)
hi = 0; % Depth of the corrugations (# of cells)
% Grating parameters on the output side
do = di; % Grating''s ridge size in (# of cells)
eo = ei; % Grating''s corrugation size in (# of cells)
LambdaGo = do+eo; % Grating''s period (# of cells)
ho = 0; % Depth of the corrugations (# of cells)
% Number of output grating periods on either side of the slit (may not be
% an integer number)
No = (Ni*LambdaGi+di)/LambdaGo;
FullNo = floor(No);
% Core thickness of the metal film in (# of cells)
t = 30;
% Total thickness of the metal film in (# of cells)
W = hi+t+ho;
% Add a dielectric substrate to the structure? (1= YES; 0= NO) :
Substrate = 0;
if Substrate == 1
    
% Thickness in # cells of the substrate behind the metal film:
    SubstrateThickness = input('Enter substrate'' thickness in # cells:');
% Choose the medium of the substrate where 1=(vacuum); 2=(dielectric)
% 3=(metal)
    SubstrateMedium = input('Choose the medium of the substrate :');
else
     SubstrateThickness = 0;
     SubstrateMedium = 1;
end
% Remove the bottom grating part? (1= YES; 0= NO) :
BottomGrooves = 0;
if BottomGrooves == 1
    BottomGratingLength = 0;
else
    BottomGratingLength = Ni*LambdaGi;
end
TopGratingLength = Ni*LambdaGi;
% Enter the size of the BOTTOM space in # cells between the bottom SF
% interface and the first bottom groove:
BottomSpace = round(1.5*47); % (47 is approximately the SPP wavelength)
% Enter the size of the TOP space in # cells between the top SF interface
% and the first top groove:
TopSpace = BottomSpace;
% The variable "gap1" defines the width of the empty space between the
% left-PML and the structure (in # cells)
if SourceType == 2
gap1 = round(7*rtau/(4*dt));
else
gap1 = N+10;
end
% The variable "gap2" defines the width of the empty space between the
% structure and the right-PML (in # cells)
gap2 = 30;
% Block all scattered waves from the plate top and bottom edges from
% disturbing the output EM fields by prolonging the metal plate extremeties
% extremities inside the SF (and possibly the PML) zone? (1=YES; 0= NO)
EdgeWavesBlockers = 0;
% Size of the computation grid along the x-axis
ie = 2*SF_width + gap1 + W + SubstrateThickness + gap2;
% Size of the computation grid along the y-axis
je = 2*SF_width + BottomSpace + BottomGratingLength + a...
+ TopGratingLength + TopSpace;
% Define various grid points position
ib = ie+1;
jb = je+1;
i0 = SF_width; % define the scattered-field zones around Region1
i1 = ie-SF_width; % ''
j0 = SF_width; % ''
j1 = je-SF_width; % ''
m0 = i0; % Origin of the linear source grid
% Width of the PML Region in # cells
PML_width = 12;
% Define some PML grid points
iebc = PML_width; % thickness of left and right PML region
jebc = PML_width; % thickness of front and back PML region
rmax = 0.00001;
orderbc = 2; % Polynomial's order of the absorbing PML
ibbc = iebc+1;
jbbc = jebc+1;
iefbc = ie+2*iebc;
jefbc = je+2*jebc;
ibfbc = iefbc+1;
jbfbc = jefbc+1;
%************************************************************************
% Incident beam(s) parameters
%************************************************************************
% N.B. This section of the program allows to model up to 2 incident plane
% waves independently of each other
% Width of the first incident beam in # cells
FirstBeamWidth = je-(2*SF_width);
if FirstBeamWidth ~= 0
FirstBeamPosition = j0 + BottomSpace...
+ BottomGratingLength+round((a-FirstBeamWidth)/2);
end
% Width of the second incident beam in # cells
SecondBeamWidth = 0;
% SecondBeamWidth = ei+4;
if SecondBeamWidth ~= 0
SecondBeamPosition = j0 + BottomSpace + BottomGratingLength...
+ a + di + round((ei-SecondBeamWidth)/2);
end
% The right TF/SF interface produces noisy back reflections if not disabled
% when there very low transmitted fields propagating toward this interface
% Disable the Right TF/SF interface ? (1= YES; 0= NO, default= 1) :
RightSF_interface = 1;
% Set the total number of time steps (duration of simulaton)
nmax = 24*period;
% Position in time (# timesteps) for starting the recording of intensity
% values at the OUTPUT detector (transmission):
IntensityTr_start = 0;
% Position in time (# timesteps) for starting the recording of intensity
% values at the INPUT detector (incident energy):
IntensityIn_start = round((1/Courant)*(sqrt(eps(1))*(gap1)));
%************************************************************************
% Visualization
%************************************************************************
% Variables visualization
% Choose which physical quantities(s) you want to visualize the evolution
% in real-time. Enter chosen components in vector form:
% Example: ['ex ey hz S ex2 ey2 hz2 Jpx Jpy']
VizComponents = ['ex ey hz S '];
if mod(length(VizComponents),3) ~= 0
disp('ERROR!! Each field component must be specified using 3 characters (and blanks to fill empty spaces if necessary)')
%break
end
% Real-time visualization start
% Position in time (# timesteps) for starting the real-time visualization
% screen. Otherwise, the program's computation runs without visualization.
VizStart = input('VizStart: Enter wave vizualisation start offset in # timesteps (ex: nmax-2*period): ');
% Set the # frames skipped to determine vizualisation's refresh rate
VizRate = 4;
%************************************************************************
% Initialization of the update coefficients (in vacuum)
%************************************************************************
% Calculate the values of the coefficients
cb = zeros(length(media), 1); 
db = zeros(length(media), 1); 
for i=1:media
ca = 1;
cb(i)=dt/epsz/eps(i)/dx;
da = 1;
db(i)= dt/muz/mur(i)/dx;
end

% Coefficients for the 1D linear source grid
caeyinc(1:ib)=ca(1);
cbeyinc(1:ib)=cb(1);
dahzinc=da(1);
dbhzinc=db(1);
% Coefficients for the 2D main (computation) grid
cbex(1:ie,1:jb)=cb(1);
cbey(1:ib,1:je)=cb(1);
dbhz(1:ie,1:je)=db(1);
%************************************************************************
% Initialization of matrix buffers for the polarization vectors and
% current density vectors
%************************************************************************
% Current density vector components
Jpx = zeros(ie,jb,order);
Jpy = zeros(ib,je,order);
% Polarization vector components
Px = zeros(ie,jb,order);
Py = zeros(ib,je,order);
% Calculate the values of the updating coefficients
for i = 1:order
c1(i) = (2 - dt*Gamma(i))/(2 + dt*Gamma(i));
c2(i) = 2*dt/(2 + dt*Gamma(i));
end
% Coefficients for the current density vector components
c1Jpx = ones(ie,jb,order);
c1Jpy = ones(ib,je,order);
% Coefficients for the polarization vector components
c2Px = zeros(ie,jb,order);
c2Py = zeros(ib,je,order);
%************************************************************************
% Initialization of Incident and Transmitted intensity buffers
%************************************************************************
Io = zeros(floor((nmax-IntensityTr_start)/(period)),1); % output buffer
%************************************************************************
% Poynting vector initialization
%************************************************************************
% Poynting vector matrix buffer
S = zeros(ie,je);

%************************************************************************
% Add the scattering structure inside the main grid
%************************************************************************
% We assume that the basic structure to be modeled is a thin metal film
% with a single flanked by periodic corrugations on either side.
% (1) INPUT SURFACE
% N.B. The structure is built starting from the center and up/down
% (1.1) Bottom part (Building from the bottom
% edge of center slit and down)
if BottomGrooves == 0
for x=1:Ni
for j=1:di
cbex(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGi) +1 -j) = cb(3);
cbey(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGi) +1 -j) = cb(3);
for k=1:order
c1Jpx(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGi) +1 -j, k) = c1(k);
c2Px(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGi) +1 -j, k) = c2(k);
c1Jpy(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGi) +1 -j, k) = c1(k);
c2Py(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGi) +1 -j, k) = c2(k);
end
end
end
end
% (1.1) Bottom part (last bottom groove)
cbex(i0+gap1 +1:i0+gap1 +hi, j0+1:j0+BottomSpace) = cb(3);
cbey(i0+gap1 +1:i0+gap1 +hi, j0+1:j0+BottomSpace) = cb(3);
for k=1:order
c1Jpx(i0+gap1 +1:i0+gap1 +hi, j0+1:j0+BottomSpace, k) = c1(k);
c2Px(i0+gap1 +1:i0+gap1 +hi, j0+1:j0+BottomSpace, k) = c2(k);
c1Jpy(i0+gap1 +1:i0+gap1 +hi, j0+1:j0+BottomSpace, k) = c1(k);
c2Py(i0+gap1 +1:i0+gap1 +hi, j0+1:j0+BottomSpace, k) = c2(k);
end
% (1.2) Top part (building from top edge of center slit and up)
for x=1:Ni
for j=1:di
cbex(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGi) +j) = cb(3);
cbey(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGi) +j) = cb(3);
for k=1:order
c1Jpx(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGi) +j, k) = c1(k);
c2Px(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGi) +j, k) = c2(k);
c1Jpy(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
    +a +(x-1)*(LambdaGi) +j, k) = c1(k);
c2Py(i0+gap1 +1:i0+gap1 +hi, j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGi) +j, k) = c2(k);
end
end
end
% (1.2) Top part (last top groove)
cbex(i0+gap1 +1:i0+gap1 +hi, j1-TopSpace:j1) = cb(3);
cbey(i0+gap1 +1:i0+gap1 +hi, j1-TopSpace:j1) = cb(3);
for k=1:order
c1Jpx(i0+gap1 +1:i0+gap1 +hi, j1-TopSpace:j1, k) = c1(k);
c2Px(i0+gap1 +1:i0+gap1 +hi, j1-TopSpace:j1, k) = c2(k);
c1Jpy(i0+gap1 +1:i0+gap1 +hi, j1-TopSpace:j1, k) = c1(k);
c2Py(i0+gap1 +1:i0+gap1 +hi, j1-TopSpace:j1, k) = c2(k);
end
% (2) OUTPUT SURFACE BU030205
% N.B. Building from the center and up/down
% (2.1) Bottom part (Building from the bottom edge
% of center slit and down)
if BottomGrooves == 0
for x=1:FullNo
for j=1:do
cbex(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGo) +1 -j) = cb(3);
cbey(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGo) +1 -j) = cb(3);
for k=1:order
c1Jpx(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGo) +1 -j, k) = c1(k);
c2Px(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGo) +1 -j, k) = c2(k);
c1Jpy(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGo) +1 -j, k) = c1(k);
c2Py(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
-(x-1)*(LambdaGo) +1 -j, k) = c2(k);
end
end
end
else % (2.1) Bottom part (last bottom groove)
cbex(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+BottomSpace) = cb(3);
cbey(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+BottomSpace) = cb(3);
for k=1:order
c1Jpx(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+BottomSpace, k) = c1(k);
c2Px(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+BottomSpace, k) = c2(k);
c1Jpy(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+BottomSpace, k) = c1(k);
c2Py(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+BottomSpace, k) = c2(k);
end
end
% % (2.1) Bottom part (filling last bottom groove)
cbex(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+do) = cb(3);
cbey(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+do) = cb(3);
for k=1:order
c1Jpx(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+do, k) = c1(k);
c2Px(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+do, k) = c2(k);
c1Jpy(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+do, k) = c1(k);
c2Py(i0+gap1+hi+t +1:i0+gap1 +W, j0+1:j0+do, k) = c2(k);
end
% (2.2) Top part (building from top edge of center slit and up)
for x=1:FullNo
for j=1:do
cbex(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGo) +j) = cb(3);
cbey(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGo) +j) = cb(3);
for k=1:order
c1Jpx(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGo) +j, k) = c1(k);
c2Px(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGo) +j, k) = c2(k);
c1Jpy(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGo) +j, k) = c1(k);
c2Py(i0+gap1+hi+t +1:i0+gap1 +W,j0+BottomSpace+BottomGratingLength...
+a +(x-1)*(LambdaGo) +j, k) = c2(k);
end
end
end
% (2.1) Top part (filling last top groove)
cbex(i0+gap1+hi+t +1:i0+gap1 +W, j1-do:j1) = cb(3);
cbey(i0+gap1+hi+t +1:i0+gap1 +W, j1-do:j1) = cb(3);
for k=1:order
c1Jpx(i0+gap1+hi+t +1:i0+gap1 +W, j1-do:j1, k) = c1(k);
c2Px(i0+gap1+hi+t +1:i0+gap1 +W, j1-do:j1, k) = c2(k);
c1Jpy(i0+gap1+hi+t +1:i0+gap1 +W, j1-do:j1, k) = c1(k);
c2Py(i0+gap1+hi+t +1:i0+gap1 +W, j1-do:j1, k) = c2(k);
end
% (3) Adding metallic core base plate of the grooves
% (3.1) Bottom part of base (relative to hole)
cbex(i0+gap1+hi +1:i0+gap1+hi +t, j0 +1:j0+BottomSpace...
+ BottomGratingLength) = cb(3);
cbey(i0+gap1+hi +1:i0+gap1+hi +t, j0 +1:j0+BottomSpace...
+ BottomGratingLength) = cb(3);
for k=1:order
    c1Jpx(i0+gap1+hi +1:i0+gap1+hi +t, j0 +1:j0+BottomSpace...
+ BottomGratingLength, k) = c1(k);
c2Px(i0+gap1+hi +1:i0+gap1+hi +t, j0 +1:j0+BottomSpace...
+ BottomGratingLength, k) = c2(k);
c1Jpy(i0+gap1+hi +1:i0+gap1+hi +t, j0 +1:j0+BottomSpace...
+ BottomGratingLength, k) = c1(k);
c2Py(i0+gap1+hi +1:i0+gap1+hi +t, j0 +1:j0+BottomSpace...
+ BottomGratingLength, k) = c2(k);
end
% (3.2) Upper part of base plate (relative to hole)
cbex(i0+gap1+hi+1:i0+gap1+hi+t,j0+BottomSpace+BottomGratingLength+a+1:...
j0+a+BottomGratingLength+TopGratingLength+BottomSpace+TopSpace)=cb(3);
cbey(i0+gap1+hi+1:i0+gap1+hi+t,j0+BottomSpace+BottomGratingLength+a+1:...
j0+a+BottomGratingLength+TopGratingLength+BottomSpace+TopSpace)=cb(3);
for k=1:order
c1Jpx(i0+gap1+hi+1:i0+gap1+hi+t,j0+BottomSpace+BottomGratingLength+a+1:...
j0+a+BottomGratingLength+TopGratingLength+BottomSpace+TopSpace,k)=c1(k);
c2Px(i0+gap1+hi+1:i0+gap1+hi+t,j0+BottomSpace+BottomGratingLength+a+1:...
j0+a+BottomGratingLength+TopGratingLength+BottomSpace+TopSpace,k)=c2(k);
c1Jpy(i0+gap1+hi+1:i0+gap1+hi+t,j0+BottomSpace+BottomGratingLength+a+1:...
j0+a+BottomGratingLength+TopGratingLength+BottomSpace+TopSpace,k)=c1(k);
c2Py(i0+gap1+hi+1:i0+gap1+hi+t,j0+BottomSpace+BottomGratingLength+a+1:...
j0+a+BottomGratingLength+TopGratingLength+BottomSpace+TopSpace,k)=c2(k);
end
% (4) Adding dielectric substrate in main grid BU030804
if SubstrateThickness ~= 0
cbex(i0+gap1+W +1 : i0+gap1+W +SubstrateThickness, j0 +1: j1) =cb(2);
cbey(i0+gap1+W +1 : i0+gap1+W +SubstrateThickness, j0 +1: j1) =cb(2);
end
% (5) Block all scattered waves from the plate top and bottom edges from
% disturbing the output EM fields by prolonging the metal plate extremeties
% inside the SF (and possibly the PML) zone
if EdgeWavesBlockers == 1
cbex(i0+gap1+1:i0+gap1+W, 1:j0) = cb(3);
cbey(i0+gap1+1:i0+gap1+W, 1:j0) = cb(3);
cbex(i0+gap1+1:i0+gap1+W, j1+1:jb) = cb(3);
cbey(i0+gap1+1:i0+gap1+W, j1+1:je) = cb(3);
for k=1:order
c1Jpx(i0+gap1+1:i0+gap1+W, 1:j0, k) = c1(k);
c2Px(i0+gap1+1:i0+gap1+W, 1:j0, k) = c2(k);
c1Jpy(i0+gap1+1:i0+gap1+W, 1:j0, k) = c1(k);
c2Py(i0+gap1+1:i0+gap1+W, 1:j0, k) = c2(k);
c1Jpx(i0+gap1+1:i0+gap1+W, j1+1:jb, k) = c1(k);
c2Px(i0+gap1+1:i0+gap1+W, j1+1:jb, k) = c2(k);
c1Jpy(i0+gap1+1:i0+gap1+W, j1+1:je, k) = c1(k);
c2Py(i0+gap1+1:i0+gap1+W, j1+1:je, k) = c2(k);
end
end
%************************************************************************
% Creating intensity calculation index matrix
% Coordinates of the OUTPUT intensity detector
%************************************************************************
% (1) Setting the radius R of the semi-circular detector (collector)
R = a; % (# cells)
% (2) Setting the coordinates of the origin of the semi-circular detector
% Manually input coordinate origins of output detector or keep default
% position? (Default=1; User defined=2)
Question = 1;
if Question == 1
    OriginX = i0+gap1+W+1;
    OriginY = SF_width + BottomSpace + BottomGratingLength + round(a/2)+1;
else
    OriginX = input('OriginX: Enter x-coordinate origin of output intensity detector :');
    OriginY = input('OriginY: Enter y-coordinate origin of output intensity detector :');
end
% (3) Values at theta = pi/4 and -pi/4
IndexX = round(R*cos(pi/4));
IndexY = round(R*sin(pi/4));
% (3.1) Introduce origin coordinates offset
IndexX = [OriginX+IndexX;OriginX+IndexX];
IndexY = [OriginY+IndexY;OriginY-IndexY];
% (4.1) Values for theta = 0:pi/4 and 0:-pi/4
for j = 0:round(R*sin(pi/4))-1
tmp = OriginX + round(sqrt(R^2 - j^2));
% range: theta = 0:pi/4
IndexX = [IndexX;tmp];
IndexY = [IndexY;OriginY+j];
% range: theta = 0:-pi/4
IndexX = [IndexX;tmp];
IndexY = [IndexY;OriginY-j];
end

% (4.2) Values for theta = pi/4:pi/2 and -pi/4:-pi/2
for i = 0:round(R*cos(pi/4))-1
tmp = round(sqrt(R^2 - i^2));
% range: theta = pi/4:pi/2
IndexX = [IndexX;OriginX+i];
IndexY = [IndexY;OriginY+tmp];
%range: theta = -pi/4:-pi/2
IndexX = [IndexX;OriginX+i];
IndexY = [IndexY;OriginY-tmp];
end
%***********************************************************************
% Grating preview plot
%***********************************************************************
% (1) Initializing matrix buffer for grating and detector visualization
SemiCircDetector = c1Jpx(:,:,2); % same size as "c1Jpx" matrix
% (2) Filling Grating position coordinates with zeros
for i = 1:size(SemiCircDetector,1)
for j = 1:size(SemiCircDetector,2)
if SemiCircDetector(i,j) ~= 1
SemiCircDetector(i,j) = 0;
end
end
end
% (2) Filling Input and Output Intensity Detector position coordinates with
% zeros to identify their position (N.B. smaller value = darker line)
% (2.1) Output detector
for k=1:length(IndexX)
SemiCircDetector(IndexX(k),IndexY(k)) = 0.3;
end
% (3) Creating the plot preview
subplot(1,1,1),pcolor(SemiCircDetector');
shading flat;
axis([1 ie 1 je]);
colormap(gray); % Grayscale figure
title('Grating and Intensity detector preview');
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
Question = input('Continue simulation ? (1= YES; 0= NO) :');
if Question == 1
close % closes the figure
else
    %break
end
% (4) Erasing the zeros which were added to draw the semi-circular detector

% in order to only draw the metallic structure later on
for k=1:length(IndexX)
SemiCircDetector(IndexX(k),IndexY(k)) = 1;
end
%************************************************************************
% Source excitation
%************************************************************************
% Amplitude of source
Eo = 1;
% PLANE WAVE SOURCE
if SourceType == 1
for n = 1:nmax
source(n) = Eo*sin(omegalight*n*dt);
end
end
% GAUSSIAN PULSE SOURCE
if SourceType == 2
lambda0 = lambda; % central wavelength of pulse
freq0 = cc/lambda0; % central frequency of pulse
omega0 = 2.0*pi*freq0;
rtau = input('Enter the half-width value at 1/e of the Gaussian pulse [sec] (default: 1.5e-15):');
tau = rtau/dt; % half-width value at 1/e in # timesteps
delay = 3*tau; % time delay in # timesteps
source = zeros(1,nmax);
for n = 1:nmax
source(n) = Eo*sin(omega0*(n-delay)*dt)*exp(-((n-delay)^2/tau^2));
end
end
% PLOT THE SOURCE TIME EVOLUTION
figure
TimeAxis = dt*(1:length(source));
plot(TimeAxis,source, '-.')
xlabel('Time [sec]')
ylabel('Amplitude at the hard source point')
title(['Source Amplitude as a function of time']);
% clear TimeAxis;
% PLOT THE SOURCE POWER SPECTRAL DENSITY
NPoints = 10*length(source); % Sets spectral resolution of the FFT
source_fft = fft(source, NPoints);
PowerSpectrum = source_fft.*conj(source_fft)/NPoints;
FreqAxis = (dt^(-1))/NPoints*(0:(NPoints/2-1));
% position of frequencies to zoom-in
F1Position = round(0.0020*NPoints);
F2Position = round(0.0190*NPoints);
% Frequencies at 1/e :'half-width frequencies'
HWFreqs = zeros(1,1);

PowerSpectrumMax = max(PowerSpectrum(1:NPoints/2));
for k = 2:F2Position
MiddleDifference = abs(PowerSpectrum(k) - PowerSpectrumMax/exp(1));
LeftDifference = abs(PowerSpectrum(k-1) - PowerSpectrumMax/exp(1));
RightDifference = abs(PowerSpectrum(k+1) - PowerSpectrumMax/exp(1));
if (MiddleDifference <= LeftDifference &&...
        MiddleDifference <= RightDifference)
if HWFreqs == 0
   HWFreqs = k;
else
    HWFreqs = [HWFreqs k];
end
end
end
% Corresponding (free-space) RMS wavelengths
HWlambda = cc./FreqAxis(HWFreqs);
% Find central frequency
CentFreq = find(PowerSpectrum(1:F2Position)...
==max(PowerSpectrum(1:F2Position)));
Centlambda = cc./FreqAxis(CentFreq);
% Plot full spectrum
figure
plot(FreqAxis,PowerSpectrum(1:(NPoints/2)), 'LineWidth',1.5)
title('Power spectral density (full spectrum)')
xlabel('Frequency (Hz)')
% Zoom-in on the interesting left-part of the spectrum
figure
plot(FreqAxis(F1Position:F2Position),...
PowerSpectrum(F1Position:F2Position),'-o')
title('Power spectral density')
xlabel('Frequency [Hz]')
hold on
% Left HWFreq
line([FreqAxis(HWFreqs(1)) FreqAxis(HWFreqs(1))],...
[0 PowerSpectrum(HWFreqs(1))], 'Color', 'k','LineStyle',':')
% Right HWFreq
line([FreqAxis(HWFreqs(2)) FreqAxis(HWFreqs(2))],...
[0 PowerSpectrum(HWFreqs(2))], 'Color', 'k','LineStyle',':')
% Center HWFreq
line([FreqAxis(CentFreq) FreqAxis(CentFreq)],...
[0 PowerSpectrum(CentFreq)], 'Color', 'k','LineStyle',':')
% Left lambda
text(FreqAxis(HWFreqs(1)+4), PowerSpectrum(HWFreqs(1)),...
['\lambda_{HW 1} = ', num2str(HWlambda(1)), ' m'], 'FontSize', 8)
% Right lambda
text(FreqAxis(HWFreqs(2)+4), PowerSpectrum(HWFreqs(2)),...
['\lambda_{HW 2} = ', num2str(HWlambda(2)), ' m'], 'FontSize', 8)
% Center lambda
text(FreqAxis(CentFreq+4), PowerSpectrum(CentFreq),...
['\lambda_{HW 2} = ', num2str(Centlambda), ' m'], 'FontSize', 8)
hold off

%************************************************************************
% Field components initialization
%************************************************************************
% Fields in linear source grid
eyinc(1:ib)=0.0;
hzinc(1:ie)=0.0;
% ABC buffers needed for Ez_inc
eyinc_lbc(1:nmax)=0.0;
eyinc_rbc(1:nmax)=0.0;
% Fields in main grid
ex = zeros(ie,jb);
ey = zeros(ib,je);
hz = zeros(ie,je);
% Steady-state fields
ex_2_dist = zeros(ie,jb);
ey_2_dist = zeros(ib,je);
hz_2_dist = zeros(ie,je);
S_dist = zeros(ie,je);
% Current density vectors
Jpx_2_dist = zeros(ie,jb);
Jpy_2_dist = zeros(ib,je);
% Fields in FRONT PML region
exbcf = zeros(iefbc,jebc);
eybcf = zeros(ibfbc,jebc);
hzxbcf = zeros(iefbc,jebc);
hzybcf = zeros(iefbc,jebc);
% Fields in BACK PML region
exbcb=zeros(iefbc,jbbc);
eybcb=zeros(ibfbc,jebc);
hzxbcb=zeros(iefbc,jebc);
hzybcb=zeros(iefbc,jebc);
% Fields in LEFT PML region
exbcl=zeros(iebc,jb);
eybcl=zeros(iebc,je);
hzxbcl=zeros(iebc,je);
hzybcl=zeros(iebc,je);
% Fields in RIGHT PML region
exbcr=zeros(iebc,jb);
eybcr=zeros(ibbc,je);
hzxbcr=zeros(iebc,je);
hzybcr=zeros(iebc,je);
%************************************************************************
% Fill the PML regions
%************************************************************************

delbc = iebc*dx;
sigmam = -log(rmax/100.0)*epsz*cc*(orderbc+1)/(2*delbc);
bcfactor = eps(1)*sigmam/(dx*(delbc^orderbc)*(orderbc+1));
% FRONT region
caexbcf(1:iefbc,1)=1.0;
cbexbcf(1:iefbc,1)=0.0;
for j=2:jebc
y1=(jebc-j+1.5)*dx;
y2=(jebc-j+0.5)*dx;
sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
ca1=exp(-sigmay*dt/(epsz*eps(1)));
cb1=(1.0-ca1)/(sigmay*dx);
caexbcf(1:iefbc,j)=ca1;
cbexbcf(1:iefbc,j)=cb1;
end
sigmay = bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmay*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmay*dx);
caex(1:ie,1)=ca1;
cbex(1:ie,1)=cb1;
caexbcl(1:iebc,1)=ca1;
cbexbcl(1:iebc,1)=cb1;
caexbcr(1:iebc,1)=ca1;
cbexbcr(1:iebc,1)=cb1;
for j=1:jebc
y1=(jebc-j+1)*dx;
y2=(jebc-j)*dx;
sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
sigmays=sigmay*(muz/(epsz*eps(1)));
da1=exp(-sigmays*dt/muz);
db1=(1-da1)/(sigmays*dx);
dahzybcf(1:iefbc,j)=da1;
dbhzybcf(1:iefbc,j)=db1;
caeybcf(1:ibfbc,j)=ca(1);
cbeybcf(1:ibfbc,j)=cb(1);
dahzxbcf(1:iefbc,j)=da(1);
dbhzxbcf(1:iefbc,j)=db(1);
end
% BACK region
caexbcb(1:iefbc,jbbc)=1.0;
cbexbcb(1:iefbc,jbbc)=0.0;
for j=2:jebc
y1=(j-0.5)*dx;
y2=(j-1.5)*dx;
sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
ca1=exp(-sigmay*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmay*dx);
caexbcb(1:iefbc,j)=ca1;
cbexbcb(1:iefbc,j)=cb1;
end
sigmay = bcfactor*(0.5*dx)^(orderbc+1);

ca1=exp(-sigmay*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmay*dx);
caex(1:ie,jb)=ca1;
cbex(1:ie,jb)=cb1;
caexbcl(1:iebc,jb)=ca1;
cbexbcl(1:iebc,jb)=cb1;
caexbcr(1:iebc,jb)=ca1;
cbexbcr(1:iebc,jb)=cb1;
for j=1:jebc
y1=j*dx;
y2=(j-1)*dx;
sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
sigmays=sigmay*(muz/(epsz*eps(1)));
da1=exp(-sigmays*dt/muz);
db1=(1-da1)/(sigmays*dx);
dahzybcb(1:iefbc,j)=da1;
dbhzybcb(1:iefbc,j)=db1;
caeybcb(1:ibfbc,j)=ca(1);
cbeybcb(1:ibfbc,j)=cb(1);
dahzxbcb(1:iefbc,j)=da(1);
dbhzxbcb(1:iefbc,j)=db(1);
end
% LEFT region
caeybcl(1,1:je)=1.0;
cbeybcl(1,1:je)=0.0;
for i=2:iebc
x1=(iebc-i+1.5)*dx;
x2=(iebc-i+0.5)*dx;
sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
ca1=exp(-sigmax*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmax*dx);
caeybcl(i,1:je)=ca1;
cbeybcl(i,1:je)=cb1;
caeybcf(i,1:jebc)=ca1;
cbeybcf(i,1:jebc)=cb1;
caeybcb(i,1:jebc)=ca1;
cbeybcb(i,1:jebc)=cb1;
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmax*dx);
caey(1,1:je)=ca1;
cbey(1,1:je)=cb1;
caeybcf(iebc+1,1:jebc)=ca1;
cbeybcf(iebc+1,1:jebc)=cb1;
caeybcb(iebc+1,1:jebc)=ca1;
cbeybcb(iebc+1,1:jebc)=cb1;
for i=1:iebc
x1=(iebc-i+1)*dx;
x2=(iebc-i)*dx;
sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
sigmaxs=sigmax*(muz/(epsz*eps(1)));

da1=exp(-sigmaxs*dt/muz);
db1=(1-da1)/(sigmaxs*dx);
dahzxbcl(i,1:je)=da1;
dbhzxbcl(i,1:je)=db1;
dahzxbcf(i,1:jebc)=da1;
dbhzxbcf(i,1:jebc)=db1;
dahzxbcb(i,1:jebc)=da1;
dbhzxbcb(i,1:jebc)=db1;
caexbcl(i,2:je)=ca(1);
cbexbcl(i,2:je)=cb(1);
dahzybcl(i,1:je)=da(1);
dbhzybcl(i,1:je)=db(1);
end
% RIGHT region
caeybcr(ibbc,1:je)=1.0;
cbeybcr(ibbc,1:je)=0.0;
for i=2:iebc
x1=(i-0.5)*dx;
x2=(i-1.5)*dx;
sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
ca1=exp(-sigmax*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmax*dx);
caeybcr(i,1:je)=ca1;
cbeybcr(i,1:je)=cb1;
caeybcf(i+iebc+ie,1:jebc)=ca1;
cbeybcf(i+iebc+ie,1:jebc)=cb1;
caeybcb(i+iebc+ie,1:jebc)=ca1;
cbeybcb(i+iebc+ie,1:jebc)=cb1;
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz*eps(1)));
cb1=(1-ca1)/(sigmax*dx);
caey(ib,1:je)=ca1;
cbey(ib,1:je)=cb1;
caeybcf(iebc+ib,1:jebc)=ca1;
cbeybcf(iebc+ib,1:jebc)=cb1;
caeybcb(iebc+ib,1:jebc)=ca1;
cbeybcb(iebc+ib,1:jebc)=cb1;
for i=1:iebc
x1=i*dx;
x2=(i-1)*dx;
sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
sigmaxs=sigmax*(muz/(epsz*eps(1)));
da1=exp(-sigmaxs*dt/muz);
db1=(1-da1)/(sigmaxs*dx);
dahzxbcr(i,1:je) = da1;
dbhzxbcr(i,1:je) = db1;
dahzxbcf(i+ie+iebc,1:jebc)=da1;
dbhzxbcf(i+ie+iebc,1:jebc)=db1;
dahzxbcb(i+ie+iebc,1:jebc)=da1;
dbhzxbcb(i+ie+iebc,1:jebc)=db1;
caexbcr(i,2:je)=ca(1);
cbexbcr(i,2:je)=cb(1);
dahzybcr(i,1:je)=da(1);
dbhzybcr(i,1:je)=db(1);
end
%************************************************************************
% Movie initialization
%************************************************************************
% MOVIE buffers
ex_MOVIE = zeros(ie,jb,round((nmax-VizStart)/VizRate));
ey_MOVIE = zeros(ib,je,round((nmax-VizStart)/VizRate));
hz_MOVIE = zeros(ie,je,round((nmax-VizStart)/VizRate));
% Movie frame initialization
rect = get(gcf,'Position');
rect(1:2) = [0 0];
%************************************************************************
% BEGIN TIME-STEPPING LOOP
%************************************************************************
% Start clock
tic
% Start loop
for n = 1:nmax
%************************************************************************
% Update electric and magnetic fields (E inc and H inc)
% inside the auxilliary linear source grid
%************************************************************************
% source excitation at hard source point: i=iL-2
eyinc(i0-1) = source(n);
% Ey_inc 1D FDTD equation
eyinc(2:ie) = caeyinc(2:ie).*eyinc(2:ie) +...
cbeyinc(2:ie).*(hzinc(1:ie-1)-hzinc(2:ie));
% ABC buffers
eyinc_lbc(n)=eyinc(1+1);
eyinc_rbc(n)=eyinc(ib-1);
% ABC (continuation):
if (n >= 1+2 && n <= nmax)
eyinc(1)=eyinc_lbc(n-2);
eyinc(ib)=eyinc_rbc(n-2);
else
    eyinc(1)=0;
    eyinc(ib)=0;
end
% Hz_inc 1D FDTD equation
hzinc(1:ie)=dahzinc*hzinc(1:ie) + dbhzinc*(eyinc(1:ie)-eyinc(2:ib));

%************************************************************************
% Update electric fields (EX and EY) in main grid
%************************************************************************
% Ex-FDTD equation
ex(:,2:je) = ca*ex(:,2:je)+...
cbex(:,2:je).*(hz(:,2:je) - hz(:,1:je-1) -...
dx*(sum(Jpx(:,2:je,:), 3)));
% Central beam TF/SF
if FirstBeamWidth ~= 0
for k = i0+1:i0+gap1
% Back TF/SF interface correction term:
ex(k,FirstBeamPosition+FirstBeamWidth+1) = ex(k,FirstBeamPosition+...
FirstBeamWidth+1) +...
cbex(k,FirstBeamPosition+FirstBeamWidth+1).*hzinc(k);
% Front TF/SF interface correction term:
ex(k,FirstBeamPosition+1) = ex(k,FirstBeamPosition+1) -...
cbex(k,FirstBeamPosition+1).*hzinc(k);
end
end
if SecondBeamWidth ~= 0
% Upper beam TF/SF
for k = i0+1:i0+gap1
% Back TF/SF interface correction term:
ex(k,SecondBeamPosition+SecondBeamWidth+1)=ex(k,SecondBeamPosition+...
SecondBeamWidth+1) +...
cbex(k,SecondBeamPosition+SecondBeamWidth+1).*hzinc(k);
% Front TF/SF interface correction term:
ex(k,SecondBeamPosition+1) = ex(k,SecondBeamPosition+1) -...
cbex(k,SecondBeamPosition+1).*hzinc(k);
end
end
% Ey-FDTD equation
ey(2:ie,:) = ca*ey(2:ie,:)+...
cbey(2:ie,:).*(hz(1:ie-1,:)-hz(2:ie,:) -...
dx*(sum(Jpy(2:ie,:,:), 3)));
% Central beam TF/SF
% Outside left TF/SF interface correction term:
if FirstBeamWidth ~= 0
ey(i0+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth)=ey(i0+...
1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth)...
+ cbey(i0+1,FirstBeamPosition+1:FirstBeamPosition+...
FirstBeamWidth).*hzinc(i0+1);
%Outside right TF/SF interface correction term:
if RightSF_interface == 0
ey(i1+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth) =...
ey(i1+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth)...
- cbey(i1+1,FirstBeamPosition+1:FirstBeamPosition+...
FirstBeamWidth).*hzinc(i1+1);
end
end
if SecondBeamWidth ~= 0
% Upper beam TF/SF BU040305
% Outside left TF/SF interface correction term:
ey(i0+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth) =...
ey(i0+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth)...
+ cbey(i0+1,SecondBeamPosition+1:SecondBeamPosition+...
SecondBeamWidth).*hzinc(i0+1);
%Outside right TF/SF interface correction term:
if RightSF_interface == 0
ey(i1+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth)=...
ey(i1+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth)...
- cbey(i1+1,SecondBeamPosition+1:SecondBeamPosition+...
SecondBeamWidth).*hzinc(i1+1);
end
end
%************************************************************************
% Update Px and Py values for each resonant freq in main grid
%************************************************************************
Px = Px + dt*Jpx;
Py = Py + dt*Jpy;
%************************************************************************
% Update EX in PML regions
%************************************************************************
% FRONT
exbcf(:,2:jebc)=caexbcf(:,2:jebc).*exbcf(:,2:jebc)-...
cbexbcf(:,2:jebc).*(hzxbcf(:,1:jebc-1)+hzybcf(:,1:jebc-1)-...
hzxbcf(:,2:jebc)-hzybcf(:,2:jebc));
ex(1:ie,1)=caex(1:ie,1).*ex(1:ie,1)-...
cbex(1:ie,1).*(hzxbcf(ibbc:iebc+ie,jebc)+...
hzybcf(ibbc:iebc+ie,jebc)-hz(1:ie,1));
% BACK
exbcb(:,2:jebc-1)=caexbcb(:,2:jebc-1).*exbcb(:,2:jebc-1)-...
cbexbcb(:,2:jebc-1).*(hzxbcb(:,1:jebc-2)+hzybcb(:,1:jebc-2)-...
hzxbcb(:,2:jebc-1)-hzybcb(:,2:jebc-1));
ex(1:ie,jb)=caex(1:ie,jb).*ex(1:ie,jb)-...
cbex(1:ie,jb).*(hz(1:ie,jb-1)-hzxbcb(ibbc:iebc+ie,1)-...
hzybcb(ibbc:iebc+ie,1));
% LEFT
exbcl(:,2:je)=caexbcl(:,2:je).*exbcl(:,2:je)-...
cbexbcl(:,2:je).*(hzxbcl(:,1:je-1)+hzybcl(:,1:je-1)-...
hzxbcl(:,2:je)-hzybcl(:,2:je));
exbcl(:,1)=caexbcl(:,1).*exbcl(:,1)-...
cbexbcl(:,1).*(hzxbcf(1:iebc,jebc)+hzybcf(1:iebc,jebc)-...
hzxbcl(:,1)-hzybcl(:,1));
exbcl(:,jb)=caexbcl(:,jb).*exbcl(:,jb)-...
cbexbcl(:,jb).*(hzxbcl(:,je)+hzybcl(:,je)-...
hzxbcb(1:iebc,1)-hzybcb(1:iebc,1));
% RIGHT
exbcr(:,2:je)=caexbcr(:,2:je).*exbcr(:,2:je)-...
cbexbcr(:,2:je).*(hzxbcr(:,1:je-1)+hzybcr(:,1:je-1)-...
hzxbcr(:,2:je)-hzybcr(:,2:je));
exbcr(:,1)=caexbcr(:,1).*exbcr(:,1)-...
cbexbcr(:,1).*(hzxbcf(1+iebc+ie:iefbc,jebc)+...
hzybcf(1+iebc+ie:iefbc,jebc)-...
hzxbcr(:,1)-hzybcr(:,1));
exbcr(:,jb)=caexbcr(:,jb).*exbcr(:,jb)-...
cbexbcr(:,jb).*(hzxbcr(:,je)+hzybcr(:,je)-...
hzxbcb(1+iebc+ie:iefbc,1)-...
hzybcb(1+iebc+ie:iefbc,1));
%************************************************************************
% Update EY in PML regions
%************************************************************************
% FRONT
eybcf(2:iefbc,:)=caeybcf(2:iefbc,:).*eybcf(2:iefbc,:)-...
cbeybcf(2:iefbc,:).*(hzxbcf(2:iefbc,:)+hzybcf(2:iefbc,:)-...
hzxbcf(1:iefbc-1,:)-hzybcf(1:iefbc-1,:));
% BACK
eybcb(2:iefbc,:)=caeybcb(2:iefbc,:).*eybcb(2:iefbc,:)-...
cbeybcb(2:iefbc,:).*(hzxbcb(2:iefbc,:)+hzybcb(2:iefbc,:)-...
hzxbcb(1:iefbc-1,:)-hzybcb(1:iefbc-1,:));
% LEFT
eybcl(2:iebc,:)=caeybcl(2:iebc,:).*eybcl(2:iebc,:)-...
cbeybcl(2:iebc,:).*(hzxbcl(2:iebc,:)+hzybcl(2:iebc,:)-...
hzxbcl(1:iebc-1,:)-hzybcl(1:iebc-1,:));
ey(1,:)=caey(1,:).*ey(1,:)-...
cbey(1,:).*(hz(1,:)-hzxbcl(iebc,:)-hzybcl(iebc,:));
% RIGHT
eybcr(2:iebc,:)=caeybcr(2:iebc,:).*eybcr(2:iebc,:)-...
cbeybcr(2:iebc,:).*(hzxbcr(2:iebc,:)+hzybcr(2:iebc,:)-...
hzxbcr(1:iebc-1,:)-hzybcr(1:iebc-1,:));
ey(ib,:)=caey(ib,:).*ey(ib,:)-...
cbey(ib,:).*(hzxbcr(1,:)+hzybcr(1,:)- hz(ie,:));
%************************************************************************
% Update magnetic fields (HZ) in main grid
%************************************************************************

hz(1:ie,1:je) = da*hz(1:ie,1:je)+...
dbhz(1:ie,1:je).*(ex(1:ie,2:jb)-ex(1:ie,1:je)+...
ey(1:ie,1:je)-ey(2:ib,1:je));
% Central beam TF/SF
% Left TF/SF interface correction term:
if FirstBeamWidth ~= 0
hz(i0+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth) =...
hz(i0+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth)...
+ dbhz(i0+1,FirstBeamPosition+1:FirstBeamPosition+...
FirstBeamWidth).*eyinc(i0+1);
% Right TF/SF interface correction term:
if RightSF_interface == 0
hz(i1+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth)=...
    hz(i1+1,FirstBeamPosition+1:FirstBeamPosition+FirstBeamWidth)...
    - dbhz(i1+1,FirstBeamPosition+1:FirstBeamPosition+...
    FirstBeamWidth).*eyinc(i1+2);
end
end
% Upper beam TF/SF BU040305
% Left TF/SF interface correction term:
if SecondBeamWidth ~= 0
hz(i0+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth) =...
hz(i0+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth)...
+ dbhz(i0+1,SecondBeamPosition+1:SecondBeamPosition+...
SecondBeamWidth).*eyinc(i0+1);
% Right TF/SF interface correction term:
if RightSF_interface == 0
hz(i1+1,SecondBeamPosition+1:SecondBeamPosition+SecondBeamWidth)=...
hz(i1+1,SecondBeamPosition+1:SecondBeamPosition+...
SecondBeamWidth)-dbhz(i1+1,SecondBeamPosition+1:...
SecondBeamPosition+SecondBeamWidth).*eyinc(i1+2);
end
end
%************************************************************************
% Update Jpx and Jpy current density values
%************************************************************************
for k = 1:order
Jpx(:,:,k) = c1Jpx(:,:,k).*Jpx(:,:,k) + c2Px(:,:,k).*...
((omegap^2)*epsz*f(k)*ex - ((omega(k))^2)*Px(:,:,k));
Jpy(:,:,k) = c1Jpy(:,:,k).*Jpy(:,:,k) + c2Py(:,:,k).*...
((omegap^2)*epsz*f(k)*ey - ((omega(k))^2)*Py(:,:,k));
end
%************************************************************************
% Update HZX in PML regions
%************************************************************************
% FRONT
hzxbcf(1:iefbc,:)=dahzxbcf(1:iefbc,:).*hzxbcf(1:iefbc,:)-...
dbhzxbcf(1:iefbc,:).*(eybcf(2:ibfbc,:)-eybcf(1:iefbc,:));
% BACK
hzxbcb(1:iefbc,:)=dahzxbcb(1:iefbc,:).*hzxbcb(1:iefbc,:)-...
dbhzxbcb(1:iefbc,:).*(eybcb(2:ibfbc,:)-eybcb(1:iefbc,:));
% LEFT
hzxbcl(1:iebc-1,:)=dahzxbcl(1:iebc-1,:).*hzxbcl(1:iebc-1,:)-...
dbhzxbcl(1:iebc-1,:).*(eybcl(2:iebc,:)-eybcl(1:iebc-1,:));
hzxbcl(iebc,:)=dahzxbcl(iebc,:).*hzxbcl(iebc,:)-...
dbhzxbcl(iebc,:).*(ey(1,:)-eybcl(iebc,:));
% RIGHT
hzxbcr(2:iebc,:)=dahzxbcr(2:iebc,:).*hzxbcr(2:iebc,:)-...
dbhzxbcr(2:iebc,:).*(eybcr(3:ibbc,:)-eybcr(2:iebc,:));
hzxbcr(1,:)=dahzxbcr(1,:).*hzxbcr(1,:)-...
dbhzxbcr(1,:).*(eybcr(2,:)-ey(ib,:));
%************************************************************************
% Update HZY in PML regions
%************************************************************************
% FRONT
hzybcf(:,1:jebc-1)=dahzybcf(:,1:jebc-1).*hzybcf(:,1:jebc-1)-...
dbhzybcf(:,1:jebc-1).*(exbcf(:,1:jebc-1)-exbcf(:,2:jebc));
hzybcf(1:iebc,jebc)=dahzybcf(1:iebc,jebc).*hzybcf(1:iebc,jebc)-...
dbhzybcf(1:iebc,jebc).*(exbcf(1:iebc,jebc)-exbcl(1:iebc,1));
hzybcf(iebc+1:iebc+ie,jebc)=...
dahzybcf(iebc+1:iebc+ie,jebc).*hzybcf(iebc+1:iebc+ie,jebc)-...
dbhzybcf(iebc+1:iebc+ie,jebc).*(exbcf(iebc+1:iebc+ie,jebc)-...
ex(1:ie,1));
hzybcf(iebc+ie+1:iefbc,jebc)=...
dahzybcf(iebc+ie+1:iefbc,jebc).*hzybcf(iebc+ie+1:iefbc,jebc)-...
dbhzybcf(iebc+ie+1:iefbc,jebc).*(exbcf(iebc+ie+1:iefbc,jebc)-...
exbcr(1:iebc,1));
% BACK
hzybcb(1:iefbc,2:jebc)=dahzybcb(1:iefbc,2:jebc).*hzybcb(1:iefbc,2:jebc)-...
dbhzybcb(1:iefbc,2:jebc).*(exbcb(1:iefbc,2:jebc)-exbcb(1:iefbc,3:jbbc));
hzybcb(1:iebc,1)=dahzybcb(1:iebc,1).*hzybcb(1:iebc,1)-...
dbhzybcb(1:iebc,1).*(exbcl(1:iebc,jb)-exbcb(1:iebc,2));
hzybcb(iebc+1:iebc+ie,1)=...
dahzybcb(iebc+1:iebc+ie,1).*hzybcb(iebc+1:iebc+ie,1)-...
dbhzybcb(iebc+1:iebc+ie,1).*(ex(1:ie,jb)-exbcb(iebc+1:iebc+ie,2));
hzybcb(iebc+ie+1:iefbc,1)=...
dahzybcb(iebc+ie+1:iefbc,1).*hzybcb(iebc+ie+1:iefbc,1)-...
dbhzybcb(iebc+ie+1:iefbc,1).*(exbcr(1:iebc,jb)-...
exbcb(iebc+ie+1:iefbc,2));

% LEFT
hzybcl(:,1:je)=dahzybcl(:,1:je).*hzybcl(:,1:je)-...
dbhzybcl(:,1:je).*(exbcl(:,1:je)-exbcl(:,2:jb));
% RIGHT
hzybcr(:,1:je)=dahzybcr(:,1:je).*hzybcr(:,1:je)-...
dbhzybcr(:,1:je).*(exbcr(:,1:je)-exbcr(:,2:jb));
%************************************************************************
% Total absolute Incident and transmitted intensity calculation
%************************************************************************
% Output intensity detector calculation:
% Total energy detected at the output detector during 1 period (cycle)
if (n >= IntensityTr_start && n <= nmax)
for k = 1:floor((nmax - IntensityTr_start)/(period))
if (n >= IntensityTr_start+(k-1)*round(period) &&...
n <= IntensityTr_start+k*round(period))
for j = 1:length(IndexX)
Io(k) = Io(k) + S(IndexX(j), IndexY(j));
end
end
end
end
if SourceType == 2
% Incident and Transmitted intensity buffers for FFT analysis
if (n >= 1 && n < round(7*tau))
    Ii_tot(n) = sum(S(i0 + 1, j0+1:j1));
end
% Total energy outputted from the hole at exit
for k = 1:length(IndexX)
    Io_tot(n) = S(IndexX(k), IndexY(k));
end
end
% Ex, Ey, Hz, Jpx and Jpy steady-state distributions
if (n >= nmax-(2*period) && n <= nmax)
ex_2_dist = ex_2_dist + (abs(ex).^2);
ey_2_dist = ey_2_dist + (abs(ey).^2);
hz_2_dist = hz_2_dist + (abs(hz).^2);
S_dist = S_dist + S;
Jpx_2_dist = Jpx_2_dist + (abs(sum(Jpx(:,:,:), 3)).^2);
Jpy_2_dist = Jpy_2_dist + (abs(sum(Jpy(:,:,:), 3)).^2);
end
%************************************************************************
% Update instantaneous Intensity module
%************************************************************************
S = (1/muz)*sqrt((ey(2:ib,:).*hz(:,:)).^2 + (ex(:,2:jb).*hz(:,:)).^2);

%************************************************************************
% Real-time fields visualization
%************************************************************************
if (n >= VizStart && n <= nmax)
if mod(n,VizRate)==0
timestep = int2str(n); % string tag
% Plot Ex, Ey, Hz and S simultaneously
if length(VizComponents)/3 >= 3
for u = 1:length(VizComponents)/3
subplot(2,2,u)
if VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ex '
colormap jet
pcolor(ex');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ex at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ey '
colormap jet
pcolor(ey');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ey at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'hz '
colormap jet
pcolor(hz');
shading flat;
ymin = -1e-3;
ymax = 1e-3;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Hz at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'S '
colormap jet
pcolor(S');
shading flat;
ymin = -1.0e+2;
ymax = 1.0e+2;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['|S| at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'Jpx'
colormap jet
pcolor(sum(Jpx(:,:,:), 3)');
shading flat;
ymin = -1.0e+5;
ymax = 1.0e+5;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Jpx at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'Jpy'
colormap jet
pcolor(sum(Jpy(:,:,:), 3)');
shading flat;
ymin = -1.0e+5;
ymax = 1.0e+5;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Jpy at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ex2'
colormap jet
pcolor((abs(ex).^2)');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ex at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ey2'
colormap jet
pcolor((abs(ey).^2)');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ey at time step = ',timestep]);

elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'hz2'
colormap jet
pcolor((abs(hz).^2)');
shading flat;
ymin = -1e-3;
ymax = 1e-3;
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Hz at time step = ',timestep]);
end
end
end
if length(VizComponents)/3 == 2
for u = 1:length(VizComponents)/3
subplot(2,1,u)
if VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ex '
colormap jet
pcolor(ex');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ex at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ey '
colormap jet
pcolor(ey');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ey at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'hz '
colormap jet
pcolor(hz');
shading flat;
ymin = -1e-3;
ymax = 1e-3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])

ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Hz at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'S '
colormap jet
pcolor(S');
shading flat;
ymin = -1.0e+2;
ymax = 1.0e+2;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['|S| at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'Jpx'
colormap jet
pcolor(sum(Jpx(:,:,:), 3)');
shading flat;
ymin = -1.0e+5;
ymax = 1.0e+5;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Jpx at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'Jpy'
colormap jet
pcolor(sum(Jpy(:,:,:), 3)');
shading flat;
ymin = -1.0e+5;
ymax = 1.0e+5;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Jpy at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ex2'
colormap jet
pcolor((abs(ex).^2)');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ex at time step = ',timestep]);

elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'ey2'
colormap jet
pcolor((abs(ey).^2)');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ey at time step = ',timestep]);
elseif VizComponents(u+(u-1)*2:3+(u-1)*3) == 'hz2'
colormap jet
pcolor((abs(hz).^2)');
shading flat;
ymin = -1e-3;
ymax = 1e-3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Hz at time step = ',timestep]);
end
end
end
if length(VizComponents)/3 == 1
subplot(1,1,1)
if VizComponents(1:3) == 'ex '
colormap jet
pcolor(ex');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Ex at time step = ',timestep]);
elseif VizComponents(1:3) == 'ey '
colormap jet
pcolor(ey');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
99
title(['Ey at time step = ',timestep]);
elseif VizComponents(1:3) == 'hz '
colormap jet
pcolor(hz');
shading flat;
ymin = -1e-3;
ymax = 1e-3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Hz at time step = ',timestep]);
elseif VizComponents(1:3) == 'S '
colormap jet
pcolor(S');
shading flat;
ymin = -1.0e+2;
ymax = 1.0e+2;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['|S| at time step = ',timestep]);
elseif VizComponents(1:3) == 'Jpx'
colormap jet
pcolor(sum(Jpx(:,:,:), 3)');
shading flat;
ymin = -5.0e+4;
ymax = 5.0e+4;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Jp_x at time step = ',timestep]);
elseif VizComponents(1:3) == 'Jpy'
colormap jet
pcolor(sum(Jpy(:,:,:), 3)');
shading flat;
ymin = -5.0e+4;
ymax = 5.0e+4;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['Jp_y at time step = ',timestep]);
elseif VizComponents(1:3) == 'ex2 '

colormap jet
pcolor((abs(ex).^2)');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['|Ex|^2 at time step = ',timestep]);
elseif VizComponents(1:3) == 'ey2'
colormap jet
pcolor((abs(ey).^2)');
shading flat;
ymin = -0.3;
ymax = 0.3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['|Ey|^2 at time step = ',timestep]);
elseif VizComponents(1:3) == 'hz2 '
colormap jet
pcolor((abs(hz).^2)');
shading flat;
ymin = -1e-3;
ymax = 1e-3;
caxis([ymin ymax]);
axis equal tight
colorbar;
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title(['|Hz|^2 at time step = ',timestep]);
end
end
% SAVE MOVIE FRAMES
nn = n/VizRate - VizStart/VizRate;
if nn == 0
nn = nn+2;
else
nn = round(nn);
end
ex_MOVIE(:,:,nn) = ex;
ey_MOVIE(:,:,nn) = ey;
hz_MOVIE(:,:,nn) = hz;
end
end
%************************************************************************
% Display final Transmission and time evolution graphics
%************************************************************************
if n == nmax
format long
% Transmittted optical intensity evolution with time
figure
plot(Io, '-o')
xlabel(['Number of periods (time) (where one period = ',...
num2str((1/freq)/1e-15),' fs)'])
ylabel('Transmitted optical intensity (Io)')
title(['Transmitted optical intensity (Io) as a function of'...
' time for \lambda = ', num2str(lambda/1e-6),' \mum']);
text(1.5,3*max(Io)/4,['Io_{final} = ', num2str(Io(end))]);
% Ex field profile along the metal surface
figure
plot(ex(i0+gap1-round(N/10),j0+1:j0+BottomSpace+BottomGratingLength+...
a+di), '-o')
line([j0+BottomSpace+BottomGratingLength+a+di-round(N/2) j0+...
BottomSpace+BottomGratingLength+a+di-round(N/2)],...
[min(ex(i0+gap1-round(N/10),j0+1:j0+BottomSpace+...
BottomGratingLength+a+di)) max(ex(i0+gap1-round(N/10),...
j0+1:j0+BottomSpace+BottomGratingLength+a+...
di))], 'Color', 'k','LineStyle',':')
legend('ex','\lambda_0 /2 distance from the groove ridge');
xlabel('cells')
ylabel('ex')
title('Ex profile along the metal surface from the bottom to the top groove ridge at the last timestep')
%************************************************************************
% Steady-State Field distributions
%************************************************************************
% Adjust the field distribution contrast by multiplying by a factor
DistributionContrast = 1;
figure
subplot(2,2,1)
contour (SemiCircDetector', ':w');
hold on
pcolor((ex_2_dist/(2*period))');
shading flat;
ymin = 0;
ymax = 0.3;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])

ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|Ex|^2 distribution');
hold off
subplot(2,2,2)
contour (SemiCircDetector', ':w');
hold on
pcolor((ey_2_dist/(2*period))');
shading flat;
ymin = 0;
ymax = 0.3;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|Ey|^2 distribution');
hold off
subplot(2,2,3)
contour (SemiCircDetector', ':w');
hold on
ContrastFactor = 1.5e5;
pcolor((ContrastFactor*hz_2_dist/(2*period))');
shading flat;
ymin = 0;
ymax = 1;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|Hz|^2 distribution');
hold off
subplot(2,2,4)
contour (SemiCircDetector', ':w');
hold on
ContrastFactor = 1.0e1;
pcolor((ContrastFactor*S/(2*period))');
shading flat;
ymin = 0;
ymax = 9.0e-2;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|S| distribution');
hold off
% OUTPUT squared field distributions (Zoom-in behind the output part
% and increasing the visualization constrast)

% Adjust the output field distribution contrast by multiplying by a
% factor
OutputDistributionContrast = 1;
figure
subplot(2,2,1),pcolor((ex_2_dist(i0+gap1+W-2:ie,:)/(2*period))');
shading flat;
ymin = 0;
ymax = 0.01;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|Ex|^2 Output distribution');
subplot(2,2,2)
pcolor((ey_2_dist(i0+gap1+W-2:ib,:)/(2*period))');
shading flat;
ymin = 0;
ymax = 0.01;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|Ey|^2 Output distribution');
ContrastFactor = 5.0e6;
subplot(2,2,3)
pcolor((ContrastFactor*hz_2_dist(i0+gap1+W-2:ie,:)/(2*period))');
shading flat;
ymin = 0;
ymax = 1;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|Hz|^2 Output distribution');
ContrastFactor = 1.0e1;
subplot(2,2,4)
pcolor((ContrastFactor*S_dist(i0+gap1+W-2:ie,:)/(2*period))');
shading flat;
ymin = 0;
ymax = 5.5e-5;
colormap hot
caxis([ymin ymax]);
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'])
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'])
title('|S| Output distribution');

% Poynting vectors plot %BU150705
figure
subplot(1,2,1)
[DX,DY] = gradient((ContrastFactor*S/(2*period))',1,1);
contour((ContrastFactor*S/(2*period))');
hold on
ArrowLengthScale = 4;
quiver(DX,DY,ArrowLengthScale,'b')
colormap jet
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
title('Poynting vectors','FontSize',12,'FontWeight','bold');
hold off
tmp = subplot(1,2,2);
contour (SemiCircDetector(1:ie,1:je)', ':w');
hold on
ContrastFactor = 1.0e1;
pcolor((ContrastFactor*S_dist(i0+gap1+W-2:ie,:)/(2*period))');
shading flat;
ymin = 0;
ymax = DistributionContrast*9.0e-2;
colormap hot
caxis([ymin ymax]);
set(tmp,'FontSize',11,'FontWeight','normal')
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
title('|S| distribution','FontSize',12,'FontWeight','bold');
hold off
% Current density distributions
if gap1 <= gap2
tmp2 = gap1;
else
    tmp2 = gap2;
end
figure
tmp = subplot(1,2,1);
pcolor((Jpx_2_dist(i0+gap1-tmp2:i0+gap1+W+tmp2,:)/(2*period))');
shading flat;
ymin = 0;
ymax = 1.0e+8;
colormap hot
caxis([ymin ymax]);
set(tmp,'FontSize',11,'FontWeight','normal')
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
title('|Jpx|^2 distribution','FontSize',12,'FontWeight','bold');
tmp = subplot(1,2,2);
pcolor((Jpy_2_dist(i0+gap1-tmp2:i0+gap1+W+tmp2,:)/(2*period))');
shading flat;
ymin = 0;
ymax = 1.0e+8;
colormap hot
caxis([ymin ymax]);
set(tmp,'FontSize',11,'FontWeight','normal')
axis equal tight
xlabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
ylabel(['x ', num2str(dx/1e-6), ' [\mum]'],...
'FontSize',11,'FontWeight','normal')
title('|Jpy|^2 distribution','FontSize',12,'FontWeight','bold');
end
%***********************************************************************
% END TIME-STEPPING LOOP
%***********************************************************************
end
format short
CPU_Time = toc; % Stop clock
disp(['Total computing time : ',num2str(CPU_Time),...
' sec or ',num2str(CPU_Time/60),' min'])