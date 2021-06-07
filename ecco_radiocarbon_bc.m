% Script to calculate steady-state D14C and age using TMM for ECCO
% By Tatsuro Tanioka, Dec 2020
% For Mathematical information refer to 
%   1. Khatiwala (2005), Ocean Modelling, 9(1), GB3001
%   2. Khatiwala (2007), GBC, 21, GB3001
% Uses: calc_steadystate_tracer.m in "~TMM2/tmm_matlab_code/TMM"
% Uses Surface D14C boundary condition from Key et al. (2004) Data
% Need to use ecco_radiocarbon.sh to process on MSI

%base_path='~/TMM2/MITgcm_2.8deg';
base_path='~/TMM2/MITgcm_ECCO';
addpath(genpath('~/TMM2/tmm_matlab_code'));
load(fullfile(base_path,'config_data'))
bgcDataPath=fullfile(base_path,'BiogeochemData');

rearrangeProfiles=0;
bigMat=0;
%dt = 3600.0*24.0/16.0;  % Tracer Timestep = 5400 s = 1/16 d   (bgc time step)
%dt = 3600.0*24.0/2.0;    % Tracer Timestep = 43,200 s = 1/2 d  (physics time step)
dt = 1

matrixPath=fullfile(base_path,matrixPath);
explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dznom','dz','da','x','y','z','deltaT','gridType')
load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

Ib=find(izBox==1);
I_interior = find(izBox~=1);

windFile=fullfile(bgcDataPath,'wind_speed');
load(windFile,'windspeed');
iceFile = fullfile(bgcDataPath,'ice_fraction');
load(iceFile,'Fice');
Ficeb=gridToMatrix(Fice,Ib,boxFile,gridFile,1);
windb=gridToMatrix(windspeed,Ib,boxFile,gridFile,1);

% Explicit transport matrix
I=speye(nb,nb);
disp('loading annual mean explicit TM')	
load(explicitAnnualMeanMatrixFile,'Aexpms')	
Aexpms2=Aexpms(I_interior,I_interior); % Aexpms with Interior Points only
Bexpms=Aexpms(I_interior,Ib);          % Explicit Boundary Matrix

% Implicit transport matrix
disp('loading annual mean implicit TM')		
load(implicitAnnualMeanMatrixFile,'Aimpms')
Aimpms2=Aimpms(I_interior,I_interior); % Aimpms with Interior Points only
Bimpms=Aimpms(I_interior,Ib);          % Implicit Boundary Matrix

seconds_yr = 3600.0*24.0*365.0;
seconds_day = 3600.0*24.0;

% -----   1   ----- %   Computing D14C
%           - Surface D14C is forced to the observation
%           - D14C [permil] = (C14_units - 100)*10 from Equation (9) of Toggweiler et al. (1989)
%           - Background D14 is from Key et al. (2004) GLODAP data (Key et al. (2004), GBC, 18(4), GB4031)

%load bkgc14_cdiac_regrid.mat               % Load BKGC14 data
%bkgc14b=gridToMatrix(bkgc14_cdiac(:,:,1),Ib,boxFile,gridFile,1);

load bkgc14_cdiac_filled_regrid.mat
ocean=(bkgc14_filled>-1000)+convert2nan(double(bkgc14_filled<-1000),1);  % ocean mask (=1), land (NaN in ferret) is assigned -1e33 in matlab

bkgc14b=gridToMatrix(bkgc14_filled(:,:,1),Ib,boxFile,gridFile,1);

lambdaDecay = 1.2097d-4./seconds_yr;       % C14 Decay timescale in sec^-1

Cbc = bkgc14b./10.0+100.0;                      % Boundary Condition = Forced with GLODAP data in model C14_units
Cbc(isnan(Cbc))=0;                              % Replace NaN with 0
BC = Bexpms*Cbc+inv(Aimpms2)*Bimpms*Cbc/dt;        % Boundary Condition Forcing
C14_units = calc_steadystate_tracer('dst3',Aexpms2,Aimpms2,dt,[],[],lambdaDecay,BC);          % Calculate 14C in model C14_units
D14C = (C14_units-100.0).*10.0d0;              % Converting model unit to D14C following Eqn(10) of Toggweiler et al. (1989)
D14C_grid = matrixToGrid(D14C,I_interior,boxFile,gridFile);
D14C_grid = cat(3,bkgc14_filled(:,:,1).*ocean(:,:,1),D14C_grid);        % Combining Surface and Interior D14C to a single matrix
%varName='D14C_bc';
%write2netcdf('d14c_inv_filled_ecco.nc',D14C_grid,x,y,z,[],'D14C_bc');            % Save as .nc file
write2netcdf('d14c_inv_filled_ecco_dt_dt1.nc',D14C_grid,x,y,z,[],'D14C_bc');            % Save as .nc file
%imagesc(D14C_grid(:,:,7));colorbar;caxis([-220 -40]);                              % Display Output

%% -----   2   -----%   Computing 14C age from D14C
%D14C_grid(D14C_grid>0.0) = NaN;
%D14C_grid(D14C_grid<-1000.0) = NaN;
%C14age_grid = -8033.0.*log(1.0+D14C_grid./1000.0);   % From Stuiver and Polach [1977]
%varName='C14Age_ECCO_Years_bc';
%write2netcdf([varName '.nc'],C14age_grid,x,y,z,[],upper(varName));                  % Save as .nc file
%%imagesc(C14age_grid(:,:,7));colorbar;caxis([0 1600]);                              % Display Output

