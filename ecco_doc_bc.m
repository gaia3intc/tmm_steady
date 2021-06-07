% Script to calculate steady-state DOCt using TMM for ECCO
% KM following TT (Dec 2020)
% For Mathematical information refer to 
%   1. Khatiwala (2005), Ocean Modelling, 9(1), GB3001
%   2. Khatiwala (2007), GBC, 21, GB3001

% Uses: calc_steadystate_tracer.m in "~TMM2/tmm_matlab_code/TMM"
% Uses Surface DOCt boundary condition from the Roshan & Devries interpolated data
% Use ecco_doc_bc.sh to process on MSI

clear all

base_path='~/TMM2/MITgcm_ECCO';
addpath(genpath('~/TMM2/tmm_matlab_code'));
load(fullfile(base_path,'config_data'))
bgcDataPath=fullfile(base_path,'BiogeochemData');

matrixPath=fullfile(base_path,matrixPath);
explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dznom','dz','da','x','y','z','deltaT','gridType')
load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

% time step & decay timescales 
dt = 3600.0*24.0/16.0                          % Tracer timestep = 5400 s = 1/16 d (std=1/16 day for BGC)
%dt = 1                                        % Effectively remove time step
%dt = 3600.0*24.0/2.0                          % Tracer timestep = 43,200 sec = 1/2 day (ocean physics time step)

seconds_yr = 3600.0*24.0*365.0;                % conversion factor
seconds_day = 3600.0*24.0;
tau_docr  = 16;                                % kilo years (default=16)
tau_docsl = 5;                                 % years (default=1.5)
lambda_docr  = (1/(tau_docr*1e3))/seconds_yr;  % inverse of tau_docr  in sec^-1
lambda_docsl = (1/tau_docsl)/seconds_yr;       % inverse of tau_docsl in sec^-1

% various switches and parameter values
k=1;                                           % k=6 --> 100 m
meank = 0;                                     % use mean (1:k) for BC for k (=1) else 0
fDOCr=0.9;                                     % spatially uniform fraction of DOCr: between 0 and 1.0

RD = 0;                                        % spatially variable fDOCr using Roshan & DeVries DOCsl map (=1) else 0
const_sfDOCr = 0;                           % when RD=0: if const_sfDOR=0, use fDOCr; else, set sfDOCr to const value (default=39.4)														       
fNDP2ex = 0.0;                                 % modifying additive factor of NDP2ex; 0.0-->no modification

mask_natl = 0;                                 % use mask (=1) or not (=0)
mask_npac = 0;                                 % use mask (=1) or not (=0)
mask_so = 0;                                   % use mask (=1) or not (=0)
docr_natl = 40;                                % regional sf value (50) 40 seems good
docr_npac =45;                                 % regional sf value (45=Hansell&Carlson's endmember but not shown in Fig 2)
docr_so = 40;                                  % regional sf value (40)

% find indices for surface boundary/ML data vs. interior data
Ib1=find(izBox==1);              % just the surface
I_interior1 = find(izBox~=1);    % interior, excluding just the surface

if k==1
  Ib=Ib1;
  I_interior=I_interior1;
else
  Ibx=Ib1;
  for n=2:k, Ibx=cat(1,Ibx,find(izBox==n)); end    % Ibx has all boxes for n=1:k
  izBox1 = ones(size(izBox));
  izBox1(Ibx) = 0;
  I_interior = find(izBox1==1);  % for any value of k

  Ib=find(izBox==k);             % just the 2D slice
end

% Explicit transport matrix
I=speye(nb,nb);
disp('loading annual mean explicit TM')	
load(explicitAnnualMeanMatrixFile,'Aexpms')	
Aexpms2=Aexpms(I_interior,I_interior);         % Aexpms with Interior Points only
Bexpms=Aexpms(I_interior,Ib);                  % Explicit Boundary Matrix

% Implicit transport matrix
disp('loading annual mean implicit TM')		
load(implicitAnnualMeanMatrixFile,'Aimpms')
Aimpms2=Aimpms(I_interior,I_interior);         % Aimpms with Interior Points only
Bimpms=Aimpms(I_interior,Ib);                  % Implicit Boundary Matrix

% Split Roshan & DeVries' interpolated DOCt into DOCsl and DOCr with their decay rates
load roshan_doct_eccogrid.mat                  % var="doct"; file is a symbolic link...may or may not work
ocean=(doct>0)+convert2nan(double(doct<0),1);  % ocean (=1), land (=NaN; note that NaN in ferret is read by matlab as -1e34)
mask = doct(:,:,1);

if mask_natl==1
  natl = doct(:,:,1);
  natl(290:end,135:end)=999;
  natl(1:15,125:end)=999;
  natl = natl.*ocean(:,:,1);
  natl_index = find(natl==999);
  mask(natl_index)=100;
end
if mask_npac==1
  npac = doct(:,:,1);
  npac(100:250,130:145)=999;
  npac = npac.*ocean(:,:,1);
  npac_index = find(npac==999);
  mask(npac_index)=110;
end
if mask_so==1
  so = doct(:,:,1);
  so(:,1:20)=999;
  so = so.*ocean(:,:,1);
  so_index = find(so==999);
  mask(so_index)=120;
end

%figure(1), clf
%subplot(2,2,1), imagesc(doct(:,:,1)), title('DOCt'),colormap(jet), caxis([0 130]), colorbar
%subplot(2,2,2), imagesc(ocean(:,:,1)), title('Ocean'),colormap(jet), colorbar
%%subplot(2,2,3), imagesc(natl), title('NAtl'),colormap(jet), caxis([0 130]), colorbar
%subplot(2,2,4), imagesc(mask), title('Mask'),colormap(jet), caxis([0 130]), colorbar

%return

%figure(1), clf
%subplot(2,2,1), plot(izBox), hold on, plot(izBox(Ib),'r'), title('izBox and izBox(k)')
%subplot(2,2,2), plot(izBox(Ib)), hold on, plot(izBox(Ib1),'r'), title('izBox(Ib and k=1)')
%subplot(2,2,3), plot(izBox(I_interior1),'r'), hold on, plot(izBox(I_interior),'b'), title('izBox interior and k=1')
%subplot(2,2,4), pcolor(ocean(:,:,1)), shading flat
%%whos; return

% spatially variable NDP:Cexport ratio according to Roshan & DeVries
% use Hirata et al. 2011 to get the fpico

if RD==1
  load woa18no3_eccogrid.mat
  load modis_chla_eccogrid.mat
  
  a0 =  0.1529;
  a1 =  1.0306;
  a2 = -1.5576;
  a3 = -1.8597;
  a4 =  2.9954;
  xc = real( log10(chla(:,:,1).*ocean(:,:,1)) );
  fpico = -1./(a0+exp(a1*xc + a2)) + a3*xc + a4;
  
  no3_k5 = mean(no3(:,:,1:5),3);
  NDP2ex0 = 0.710*fpico - 0.101*real(log10(no3_k5));
  
  NDP2ex = NDP2ex0 + fNDP2ex;
  NDP2ex(NDP2ex<0) = 0;
  NDP2ex(NDP2ex>1) = 1;

%  figure(1), clf
%  subplot(2,2,1), imagesc(fpico),  colormap(jet), caxis([0 1]), colorbar
%  subplot(2,2,2), imagesc(no3_k5), colormap(jet), caxis([0 30]), colorbar
%  subplot(2,2,3), imagesc(NDP2ex), colormap(jet), caxis([0 1]), colorbar
end

%whos; return

% doesn't work...gridToMatrix expecting full 3D grid (23 levels)...not a partial grid
%  sfdocr_matrix =gridToMatrix(   fDOCr *doct(:,:,1:k),Ib,boxFile,gridFile,0);
%  sfdocsl_matrix=gridToMatrix((1-fDOCr)*doct(:,:,1:k),Ib,boxFile,gridFile,0);

% create surface BC and filename
if meank==0

  if RD==0
    if const_sfDOCr > 0

      for i=1:size(doct,1)
	for j=1:size(doct,2)
	  if doct(i,j,k) < const_sfDOCr
	    sfdocr(i,j,k) = doct(i,j,k);
	  else
	    sfdocr(i,j,k) = const_sfDOCr;
	  end
	end
      end

      sfdocsl = doct(:,:,k)-sfdocr;
      fileName=['doc_sfdocr' num2str(const_sfDOCr) '_k' num2str(k) '_rtau' num2str(tau_docr) 'ky_sltau' num2str(tau_docsl) 'y_ecco.nc']
    else
      sfdocr  =    fDOCr *doct(:,:,k);
      sfdocsl = (1-fDOCr)*doct(:,:,k);	 
      fileName=['doc_fdocr' num2str(fDOCr) '_k' num2str(k) '_rtau' num2str(tau_docr) 'ky_sltau' num2str(tau_docsl) 'y_ecco.nc']
    end
  elseif RD==1
    sfdocr  = (1-NDP2ex).*doct(:,:,k);
    sfdocsl =    NDP2ex .*doct(:,:,k);

    if fNDP2ex < 0.01
      fileName=['doc_NDP2exSL_k' num2str(k) '_rtau' num2str(tau_docr) 'ky_sltau' num2str(tau_docsl) 'y_ecco.nc']
    else
      fileName=['doc_NDP2exSL' num2str(fNDP2ex) '_k' num2str(k) '_rtau' num2str(tau_docr) 'ky_sltau' num2str(tau_docsl) 'y_ecco.nc']
    end
    
%    sfdocr  =    NDP2ex .*doct(:,:,k);
%    sfdocsl = (1-NDP2ex).*doct(:,:,k);	 
%    fileName=['doc_NDP2exR_k' num2str(k) '_rtau' num2str(tau_docr) 'ky_sltau' num2str(tau_docsl) 'y_ecco.nc']
  end

  if mask_npac == 1
    sfdocr(npac_index) = docr_npac;
    sfdocsl = doct(:,:,k) - sfdocr;
    fileName=[fileName(1:4) 'npac' num2str(docr_npac) '_' fileName(5:end)]
  end
  if mask_so == 1
    sfdocr(so_index) = docr_so;
    sfdocsl = doct(:,:,k) - sfdocr;
    fileName=[fileName(1:4) 'so' num2str(docr_so) '_' fileName(5:end)]
  end
  if mask_natl == 1
    sfdocr(natl_index) = docr_natl;
    sfdocsl = doct(:,:,k) - sfdocr;
    fileName=[fileName(1:4) 'natl' num2str(docr_natl) '_' fileName(5:end)]
  end

elseif meank==1
  doct_mean=nanmean(doct(:,:,1:k),3);           % ensure mean has same land mask as level k 
  sfdocr  =    fDOCr *doct_mean;
  sfdocsl = (1-fDOCr)*doct_mean;
  fileName=['doc_fdocr' num2str(fDOCr) '_kmean' num2str(k) '_rtau' num2str(tau_docr) 'ky_ecco.nc']
end

%whos; return
%figure(1),clf,subplot(1,2,1),image(sfdocr),colorbar, subplot(1,2,2),image(sfdocsl),colorbar

sfdocr_matrix =gridToMatrix(sfdocr, Ib,boxFile,gridFile,1);
sfdocsl_matrix=gridToMatrix(sfdocsl,Ib,boxFile,gridFile,1);

%figure(2), clf
%subplot(2,2,1), imagesc(doct(:,:,1)), title('DOCt'),colormap(jet), caxis([0 130]), colorbar
%%subplot(2,2,2), imagesc(NDP2ex),  title('NDP2ex'),  colormap(jet), caxis([0 1]), colorbar
%subplot(2,2,3), imagesc(sfdocr), title('sfDOCr'),   colormap(jet), caxis([0 130]), colorbar
%subplot(2,2,4), imagesc(sfdocsl), title('sfDOCsl'), colormap(jet), caxis([0 130]), colorbar

%figure(3), clf
%subplot(2,2,2), imagesc(NDP2ex0), title('NDP2ex0'), colormap(jet),   caxis([0 1]), colorbar
%subplot(2,2,3), imagesc((1-NDP2ex0).*doct(:,:,k)), title('sfDOCr0'), colormap(jet), caxis([0 130]), colorbar
%subplot(2,2,4), imagesc(   NDP2ex0 .*doct(:,:,k)), title('sfDOCsl0'),colormap(jet), caxis([0 130]), colorbar

%return

% refractory
Cbc = sfdocr_matrix;                            % Boundary Condition: surface [DOCr] in umol/kg
Cbc(isnan(Cbc))=0;                              % Replace NaN with 0 - this is for C14
Cbc(Cbc<0)=0;                                   % Land, which is NaN in ferret, is assigned -1e33...set to zero; why does land matter?
BC = Bexpms*Cbc+inv(Aimpms2)*Bimpms*Cbc/dt;     % Boundary Condition Forcing (645804x1)

BCflux = Aimpms2 * BC;
BCr_interior = matrixToGrid(BCflux,I_interior,boxFile,gridFile);     % 360x160x22...short one vertical layer
BCr_3d = cat(3,BCr_interior(:,:,1),BCr_interior);

docr_matrix  = calc_steadystate_tracer('dst3',Aexpms2,Aimpms2,dt,[],[],lambda_docr,BC);      % Calculate interior docr (matrix form)
docr_in_grid = matrixToGrid(docr_matrix,I_interior,boxFile,gridFile);                        % Convert interior docr to grid form   

% combine surface (land=nan), interior DOCr
if RD==0
  if const_sfDOCr > 0
    docr_grid = cat(3,sfdocr.*ocean(:,:,1),docr_in_grid);         
  else
    docr_grid = cat(3,fDOCr*doct(:,:,1:k).*ocean(:,:,1:k),docr_in_grid);
  end
elseif RD==1
  docr_grid = cat(3,sfdocr.*ocean(:,:,1:k),docr_in_grid);
end

% semilabile
Cbc = sfdocsl_matrix;                           % Boundary Condition: surface [DOCr] in umol/kg
Cbc(isnan(Cbc))=0;
Cbc(Cbc<0)=0;                                   % Replace NaN with 0
BC = Bexpms*Cbc+inv(Aimpms2)*Bimpms*Cbc/dt;     % Boundary Condition Forcing

BCflux = Aimpms2 * BC;
BCsl_interior = matrixToGrid(BCflux,I_interior,boxFile,gridFile);
BCsl_3d = cat(3,BCsl_interior(:,:,1),BCsl_interior);
  
docsl_matrix  = calc_steadystate_tracer('dst3',Aexpms2,Aimpms2,dt,[],[],lambda_docsl,BC);    % Calculate interior docr (matrix form)
docsl_in_grid = matrixToGrid(docsl_matrix,I_interior,boxFile,gridFile);                      % Convert interior docr to grid form   

% combine surface (land=nan), interior DOCsl
if RD==0
  if const_sfDOCr > 0
    docsl_grid = cat(3,sfdocsl.*ocean(:,:,1),docsl_in_grid);
  else
    docsl_grid = cat(3,(1-fDOCr)*doct(:,:,1:k).*ocean(:,:,1:k),docsl_in_grid);
  end
elseif RD==1
  docsl_grid = cat(3,sfdocsl.*ocean(:,:,1:k),docsl_in_grid);
end

% write netcdf file
ncid=netcdf.create(fileName,'CLOBBER');

%netcdf.reDef(ncid)
dimid_lon=netcdf.defDim(ncid,'Longitude',nx);
dimid_lat=netcdf.defDim(ncid,'Latitude',ny);
dimid_dep=netcdf.defDim(ncid,'Depth',nz);

lon_varid   = netcdf.defVar(ncid,'Longitude','double',dimid_lon);
lat_varid   = netcdf.defVar(ncid,'Latitude','double',dimid_lat);
dep_varid   = netcdf.defVar(ncid,'Depth','double',dimid_dep);

docr_varid  = netcdf.defVar(ncid,'docr','double',[dimid_lon dimid_lat dimid_dep]);
docsl_varid = netcdf.defVar(ncid,'docsl','double',[dimid_lon dimid_lat dimid_dep]);

bcr_varid  = netcdf.defVar(ncid,'bcr','double',[dimid_lon dimid_lat dimid_dep]);
bcsl_varid = netcdf.defVar(ncid,'bcsl','double',[dimid_lon dimid_lat dimid_dep]);

xunit='degrees_east';
yunit='degrees_north';
zunit='meter';

netcdf.putAtt(ncid,lon_varid,'units',xunit);
netcdf.putAtt(ncid,lat_varid,'units',yunit);
netcdf.putAtt(ncid,dep_varid,'units',zunit);

netcdf.endDef(ncid)

netcdf.putVar(ncid,lon_varid,x);
netcdf.putVar(ncid,lat_varid,y);
netcdf.putVar(ncid,dep_varid,z);

netcdf.putVar(ncid,docr_varid,docr_grid);
netcdf.putVar(ncid,docsl_varid,docsl_grid);

netcdf.putVar(ncid,bcr_varid,BCr_3d);
netcdf.putVar(ncid,bcsl_varid,BCsl_3d);

netcdf.close(ncid)


return
