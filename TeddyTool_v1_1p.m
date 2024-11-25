% Teddy-Tool v1.1p: Temporal Disaggregation of Daily Climate Model Data (parallelized version)
% The Teddy-Tool disaggregates daily climate model data to hourly values for climate impact and attribution studies
% Changes in v1.1p: Parallel processing of samples
% By Florian Zabel, University of Basel, Switzerland (2024)
% Contact: florian.zabel@unibas.ch

%compile TeddyTool_v1_1p
%mcc -m TeddyTool_v1_1p.m -a WFDE5Factors.m WFDE5Factors_LST.m write_netcdf.m write_netcdf_dailyval.m parsave.m

%function TeddyTool_v1_1p

clc;
clear;
warning off;
%dbstop if error
%dbclear all

disp('Starting Teddy-Tool v1.1p for temporal empirical disaggregation of daily climate model data');
disp('Reading ini file');
%read variables from teddy_cmd.ini
filename=['TeddyTool.ini'];
fid = fopen(filename);
line=fgetl(fid);
lat=str2num(strtok(line,'!'));        %latitude
line=fgetl(fid);
lon=str2num(strtok(line,'!'));        %longitude
line=fgetl(fid);
htimestep=str2num(strtok(line,'!'));  %select hourly time step for temporal disaggregation (1=1-hourly,2=2-hourly,etc.)
line=fgetl(fid);
outdir=strtrim(strtok(line,'!'));     %output directory
line=fgetl(fid);
factordir=strtrim(strtok(line,'!'));  %directory containing preprocessed subdaily factors derived from wfde5 data
line=fgetl(fid);
wfde5dir=strtrim(strtok(line,'!'));   %directory containing hourly WFDE5 reference data
line=fgetl(fid);
climdir=strtrim(strtok(line,'!'));    %directory containing ISIMIP climate model data
line=fgetl(fid);
models=strtrim(strtok(line,'!'));     %select ISIMIP climate model(s): GFDL-ESM4,IPSL-CM6A-LR,MPI-ESM1-2-HR,MRI-ESM2-0,UKESM1-0-LL
models=deblank(models);
models=strsplit(models,',');
line=fgetl(fid);
scenarios=strtrim(strtok(line,'!'));  %select ISIMIP scenario(s): ssp126,ssp585,ssp370,historical,picontrol
scenarios=deblank(scenarios);
scenarios=strsplit(scenarios,',');
line=fgetl(fid);
period=strtrim(strtok(line,'!'));
period=strsplit(period,'-');
startyear=str2num(cell2mat(period(1)));       %select startyear for ISIMIP climate data
endyear=str2num(cell2mat(period(2)));         %select endyear for ISIMIP climate data
line=fgetl(fid);
parameters_isimip=strtrim(strtok(line,'!'));  %select variable(s): tas,hurs,rsds,rlds,ps,sfcwind,pr
parameters_isimip=deblank(parameters_isimip);
parameters_isimip=strsplit(parameters_isimip,',');
line=fgetl(fid);
period=strtrim(strtok(line,'!'));
period=strsplit(period,'-');
startyear_wfde5=str2num(cell2mat(period(1))); %select startyear for hourly wfde5 data that is used to determine the diurnal profile
endyear_wfde5=str2num(cell2mat(period(2)));   %select endyear for hourly wfde5 data that is used to determine the diurnal profile
years_wfde5=endyear_wfde5-startyear_wfde5+1;
line=fgetl(fid);
timewindow=str2num(strtok(line,'!'));         %select time window of +- n days around DOY to search for similar conditions in the past (default=11)
line=fgetl(fid);
lstflag=str2num(strtok(line,'!'));            %use LST (local solar time) (1=default) or UTC (0) for data processing?
line=fgetl(fid);
rainfallflag=str2num(strtok(line,'!'));       %consider precipitation on consecutive days: 1=yes, 0=no; doy window flag (ln 13) must be >= 1
line=fgetl(fid);
precipnanflag=str2num(strtok(line,'!'));      %write NaN values for precipitation (0=no=default) in case that no precipitation event can be found in the historical hourly reference dataset (1=yes: in this case mass and energy are not preserved!)
line=fgetl(fid);
parallel_cpus=str2num(strtok(line,'!'));      %number of parallel workers
fclose(fid);

if(floor(24/htimestep)~=24/htimestep)
  error('timestep cannot be divided by 24');
end

%load all global coordinates
if(lat(1)==-9999)
  load('coordinates_global_land.mat'); %added Florian Zabel, 27.03.2024
end

%open parallel pool
delete(gcp('nocreate'));
parpool('Threads',parallel_cpus);

%read land-sea-mask (WFDE5 Land-Sea-Mask)
disp('Reading land-sea mask');
mask=load('wfde5_mask.mat').mask;

%convert lat/lon to row,col
res=0.5;
yrow=floor((90-lat)/res)+1;
xcol=floor((180+lon)/res)+1;

%calculate local solar time offset to UTC
lst_offsets=local_solar_time(res);

%calculate local midnight in UTC
midnights=lst_offsets*(-1)+1;
midnights(midnights<=0)=midnights(midnights<=0)+24;
numdays=daysact(['1-jan-',num2str(startyear)],['31-dec-',num2str(endyear)])+1;

potRadAll=zeros(366,24,length(xcol),'single');
%read potential radiation
if(find(contains(parameters_isimip,'rsds','IgnoreCase',true)>=1))
  disp('Reading potential radiation');
  filename=['potential_radiation.nc4'];
  data=ncread([filename],'prad');
  for xloop=1:length(xcol)
    x=xcol(xloop);
    y=yrow(xloop);
    if(mask(y,x)==0)
      continue
    end
    potRadAll(:,:,xloop)=data(x,y,:,:);
  end
  clear data
end

for scenarioloop=1:length(scenarios)
  scenario=scenarios{scenarioloop};
  disp(['Processing scenario ',scenario]);
  
  for modelloop=1:length(models)
    model=models{modelloop};
    disp(['Processing model ',model]);
    
    parameters_isimip_read={'tas','hurs','rsds','rlds','ps','sfcwind','pr','tasmax','tasmin'}; %for method climate analogies, all variables are required
    parameters_wfde5_factors={'tas','hurs','rsds','rlds','ps','sfcwind','pr'};
    parameters_wfde5={'Tair','Qair','SWdown','LWdown','PSurf','Wind','Rainf'};
    
    % check if precalculated factors exist
    % for varloop=parameters_isimip
    %   varloop=char(varloop);
    %   read precalculated multi-year mean subdaily factors
    %   if(strcmpi(varloop,'pr')==1)
    %     filename_wfde5_fact=[factordir,filesep,varloop,'_WFDE5_sum_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
    %   else
    %     filename_wfde5_fact=[factordir,filesep,varloop,'_WFDE5_mean_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
    %   end
    %   if(~exist(filename_wfde5_fact,'file'))
    %     prompt = ['Precompiled daily mean/sum does not exist. Compile from WFDE5 data? (y|n)?'];
    %     CalcWFDE5FactorsFlag = input(prompt,'s');
    %     if isempty(CalcWFDE5FactorsFlag)
    %       CalcWFDE5FactorsFlag = 'y';
    %     end
    %     if(strcmpi(CalcWFDE5FactorsFlag,'y')==1)
    %       if(lstflag==1)
    %         WFDE5Factors_LST(parameters_isimip,wfde5dir,factordir,startyear_wfde5,endyear_wfde5);
    %       else
    %         WFDE5Factors(parameters_isimip,wfde5dir,factordir,startyear_wfde5,endyear_wfde5);
    %       end
    %     else
    %       disp('Stopped processing due to missing subdaily factors');
    %       return
    %     end
    %   end
    % end
    % 
    % read daily ISIMIP climate model data
    % vloop=1;
    % data_cm=zeros(numdays,length(xcol),length(parameters_isimip_read),'single');
    % time_cm=zeros(numdays,1,'double');
    % for varloop=parameters_isimip_read
    %   varloop=char(varloop);
    % 
    %   path_cm=[climdir,filesep,scenario,filesep,model];
    %   filelist = dir([path_cm,filesep,model,'*',varloop,'_global_daily*.nc']);
    %   nofileflag=0;
    %   if(isempty(filelist)==1)
    %     disp(['ERROR: No files in folder ',path_cm]);
    %     nofileflag=1;
    %     break
    %   end
    %   n=1;
    %   for fileloop=1:length(filelist)
    %     open ISIMIP climate data
    %     filename_cm=[path_cm,filesep,filelist(fileloop).name];
    %     ncid = netcdf.open([filename_cm],'NOWRITE');
    %     finfo = ncinfo([filename_cm]);
    %     xcm=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lon')==1)).Length;
    %     ycm=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lat')==1)).Length;
    %     zcm=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
    %     yearoffset=filename_cm(length(filename_cm)-11:length(filename_cm)-8);
    %     if(str2double(yearoffset)>endyear)
    %       continue
    %     end
    %     timeoffset=datenum(['01.01.',num2str(yearoffset),' 00:00:00'],'dd.mm.yyyy HH:MM:SS');
    %     timeread=ncread([filename_cm],'time');
    %     timeread=double(timeread);
    %     if(timeread(1)~=0)
    %       modelspec_offset=timeread(1);
    %     else
    %       modelspec_offset=0;
    %     end
    %     timeoffset=timeoffset+(timeread(1)-modelspec_offset);
    %     time_cm_read_all=timeread+timeoffset;
    %     time_cm_read=timeread+timeoffset;
    %     time_cm_read(time_cm_read>datenum(endyear,12,31))=[];
    %     time_cm_read(time_cm_read<datenum(startyear,1,1))=[];
    %     index=find(ismember(time_cm_read_all,time_cm_read)==1);
    %     if(isempty(index))
    %       disp(['skip ',filename_cm])
    %       continue
    %     end
    %     acm=index(1);
    %     bcm=index(end);
    %     zcm=bcm-acm+1;
    %     time_cm_read_all=timeread+timeoffset;
    %     time_cm(n:n+zcm-1,1)=time_cm_read_all(acm:bcm);
    %     disp(['reading ',filename_cm]);
    %     read ISIMIP climate data for lat/lon
    %     data=ncread([filename_cm],varloop,[1 1 acm],[720 360 zcm]);
    %     netcdf.close(ncid)
    %     for xloop=1:length(xcol)
    %       x=xcol(xloop);
    %       y=yrow(xloop);
    %       if(mask(y,x)==0)
    %         continue
    %       end
    %       data_cm(n:n+zcm-1,xloop,vloop)=data(x,y,:,:);
    %     end
    %     clear data
    %     n=n+zcm;
    %   end %fileloop
    %   convert units
    %   if(strcmpi(varloop,'pr')==1)
    %     data_cm(:,:,vloop)=data_cm(:,:,vloop)*86400; %convert mass flux density to mm/day
    %   elseif(strcmpi(varloop,'ps')==1)
    %     data_cm(:,:,vloop)=data_cm(:,:,vloop)/100;   %Pa to hPa
    %   elseif(strcmpi(varloop,'clt')==1)
    %     data_cm(:,:,vloop)=data_cm(:,:,vloop)/100;   %percent to 0-1
    %   end
    %   vloop=vloop+1;
    % end %varloop
    % 
    % if(nofileflag==1)
    %   continue
    % end
    % 
    % read preprocessed wfde5 daily mean values for all variables over the climate analogy time period (standard=1980-2018)
    % wfde5_mean=zeros(years_wfde5,366,length(xcol),length(parameters_isimip_read),'single');
    % vloop=1;
    % for varloop=parameters_isimip_read
    %   varloop=char(varloop);
    %   if(strcmpi(varloop,'pr')==1)
    %     filename_wfde5_mean=[factordir,filesep,varloop,'_WFDE5_sum_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
    %     varname='sum';
    %   elseif(strcmpi(varloop,'tasmax')==1)
    %     filename_wfde5_mean=[factordir,filesep,'tas','_WFDE5_max_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
    %     varname='max';
    %   elseif(strcmpi(varloop,'tasmin')==1)
    %     filename_wfde5_mean=[factordir,filesep,'tas','_WFDE5_min_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
    %     varname='min';
    %   else
    %     filename_wfde5_mean=[factordir,filesep,varloop,'_WFDE5_mean_doy_',num2str(startyear_wfde5),'-',num2str(endyear_wfde5),'.nc4'];
    %     varname='mean';
    %   end
    %   disp(['reading ',filename_wfde5_mean])
    %   ncid = netcdf.open([filename_wfde5_mean],'NOWRITE');
    %   data=ncread([filename_wfde5_mean],varname,[1 1 1 1],[720 300 years_wfde5 366]);
    %   netcdf.close(ncid);
    %   for xloop=1:length(xcol)
    %     x=xcol(xloop);
    %     y=yrow(xloop);
    %     if(mask(y,x)==0)
    %       continue
    %     end
    %     wfde5_mean(:,:,xloop,vloop)=data(x,y,:,:);
    %   end
    %   clear data
    %   vloop=vloop+1;
    % end

    %save data to file, added 27.03.2024
    %save(['workspace_0','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat']);
    %load(['workspace_0','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat']);
    
    % bestrank=zeros(numdays,2,length(xcol),'int32'); %int16 not sufficient for 65797 global samples, changed to int32, Florian Zabel, 27.03.2024
    % ttimesteps=size(data_cm,1);  %length(data_cm) changed, Florian Zabel, 27.03.2024
    % %parfor (xloop=1:length(xcol),parallel_cpus)
    % %for xloop=1:length(xcol)
    % for xloop=1:100
    %   x=xcol(xloop);
    %   y=yrow(xloop);
    %   if(mask(y,x)==0)
    %     continue
    %   end
    %   disp(['Processing sample ',num2str(xloop),' of ',num2str(length(xcol)),' to find most similar day']);
    % 
    %   parfor t=1:ttimesteps
    %     time=time_cm(t);
    %     doy=day(datetime(time,'ConvertFrom','datenum'),'dayofyear');
    %     year=str2double(datestr(time,'yyyy'));
    %     if(doy>=60 && leap_year(year)==0) %doy 60 is skipped if no leap year
    %       doy=doy+1;
    %     end
    % 
    %     %find most similar day in past (best fit/fingerprint)
    %     abserr=NaN(years_wfde5,timewindow*2+1,length(parameters_isimip_read),'single');
    %     abserr_reshaped=NaN(years_wfde5*(timewindow*2+1),length(parameters_isimip_read),'single');
    %     timewindowsteps_all = zeros(years_wfde5, timewindow*2+1,'int16');
    %     valuesstd = NaN(years_wfde5, timewindow*2+1, 'single');
    %     stdthreshold = NaN(length(parameters_isimip_read), 'single');
    % 
    %     for vloop=1:length(parameters_isimip_read)
    %       wfde5_years=startyear_wfde5:endyear_wfde5;
    %       for ywfde5=1:years_wfde5
    %         timewindowsteps=(doy-timewindow:doy+timewindow);
    %         if(leap_year(wfde5_years(ywfde5))==0) %wfde5 leap day is at position 60 (29.02)
    %           timewindowsteps(timewindowsteps>=60)=timewindowsteps(timewindowsteps>=60)+1;
    %         end
    %         timewindowsteps(timewindowsteps<1)=timewindowsteps(timewindowsteps<1)+366;
    %         timewindowsteps(timewindowsteps>366)=timewindowsteps(timewindowsteps>366)-366;
    %         timewindowsteps_all(ywfde5,:)=timewindowsteps;
    %         abserr(ywfde5,:,vloop)=abs(squeeze(wfde5_mean(ywfde5,timewindowsteps,xloop,vloop))-data_cm(t,xloop,vloop)');
    %         valuesstd(ywfde5,:)=squeeze(wfde5_mean(ywfde5,timewindowsteps,xloop,vloop));
    %       end
    %       stdthreshold(vloop)=std(valuesstd,0,'all'); %get standard deviation
    % 
    %       %consider rainfall classes, added 22.10.2022
    %       pr_clim=NaN(1,3);
    %       pr_ref=NaN(1,3);
    %       consecutive_pr_states=NaN(years_wfde5,3,timewindow*2+1);
    %       if(strcmpi(parameters_isimip_read(vloop),'pr')==1)
    %         abserror_precip=abserr(:,:,vloop);
    %         if(t==1)
    %           pr_clim(1,1)=data_cm(t,xloop,vloop);       %climate model precipitation
    %           pr_clim(1,2:3)=data_cm(1:t+1,xloop,vloop); %climate model precipitation
    %         elseif(t==size(data_cm,1))
    %           pr_clim(1,1:2)=data_cm(t-1:t,xloop,vloop); %climate model precipitation
    %           pr_clim(1,3)=data_cm(t,xloop,vloop);       %climate model precipitation
    %         else
    %           pr_clim(1,:)=data_cm(t-1:t+1,xloop,vloop); %climate model precipitation
    %         end
    %         pr_clim(pr_clim<1)=0; %reduce drizzle effect for precipitation in clim
    %         pr_clim(pr_clim>0)=1; %convert clim to boolean
    %         timewindowsteps=(doy-timewindow:doy+timewindow);
    %         for mw=1:length(timewindowsteps)    %moving window +- 1 day over timewindow
    %           for ywfde5=1:years_wfde5
    %             mwsteps=[timewindowsteps(mw)-1,timewindowsteps(mw),timewindowsteps(mw)+1];
    %             if(leap_year(wfde5_years(ywfde5))==0)
    %               mwsteps(mwsteps>=60)=mwsteps(mwsteps>=60)+1;
    %             end
    %             mwsteps(mwsteps<1)=mwsteps(mwsteps<1)+366;
    %             mwsteps(mwsteps>366)=mwsteps(mwsteps>366)-366;
    %             pr_ref(ywfde5,:)=wfde5_mean(ywfde5,mwsteps,xloop,vloop); %reference wfde5 precipitation
    %           end
    %           pr_ref(pr_ref<1)=0;   %reduce drizzle effect for precipitation in ref
    %           pr_ref(pr_ref>0)=1;   %convert ref to boolean
    %           consecutive_pr_states(:,:,mw)=pr_ref(:,:)==pr_clim(1,:); %check if consecutive days with or without precipitation is the same in both datasets
    %         end
    %         %same day must have the same precipitation status
    %         proof_dayspecific_pr_state=squeeze(consecutive_pr_states(:,2,:));
    %         if(sum(proof_dayspecific_pr_state,'all')>=1)
    %           abserror_precip(proof_dayspecific_pr_state==0)=NaN; %mark all days that don't have the same precipitation status on same day
    %         end
    %         if(rainfallflag==1)
    %           %consecutive days (yesterday, today, tomorrow) must have same precipitation status
    %           proof_consecutive_pr_states=squeeze(sum(consecutive_pr_states,2));
    %           if(sum(proof_consecutive_pr_states==3,'all')>=1)
    %             abserror_precip(proof_consecutive_pr_states<3)=NaN;  %mark uneven sequence of precipitation days
    %           elseif(sum(proof_consecutive_pr_states==2,'all')>=1)
    %             abserror_precip(proof_consecutive_pr_states<2)=NaN;  %if three consecutive days do not exist, take at least two secutive days
    %           end
    %         end
    %         abserr(:,:,vloop)=abserror_precip;
    %       end
    %     end
    %     abserr_reshaped(:,:)=reshape(abserr,(timewindow*2+1)*years_wfde5,length(parameters_isimip_read));
    %     yearn=repmat([1:years_wfde5]',timewindow*2+1,1);
    %     doyn=reshape(timewindowsteps_all,[(timewindow*2+1)*years_wfde5,1]);
    % 
    %     %remove nan values for sorting
    %     masknan=sum(isnan(abserr_reshaped),2);
    %     test=sum(masknan==0);
    %     if(test>=1) %only if at least one sample remains after exclusion
    %       %remove samples from population that do not fit the criterias (e.g. precipitation or standard deviation)
    %       yearn(masknan>0)=[];
    %       doyn(masknan>0)=[];
    %       abserr_reshaped(masknan>0,:)=[];
    %     end
    %     [srt,index]=sort(abserr_reshaped(:,:),1);
    %     rank=zeros(length(yearn),length(parameters_isimip_read),'single');
    %     for vloop=1:length(parameters_isimip_read)
    %       pos=index(:,vloop);
    %       idxrepeat=[false diff(srt(:,vloop))'==0];
    %       rnkNoSkip=cumsum(~idxrepeat);
    %       rnk=1:numel(pos);
    %       rnk(idxrepeat)=rnkNoSkip(idxrepeat);
    %       rank(pos,vloop)=rnk;
    %     end
    %     rank(isnan(abserr_reshaped)==1)=NaN; %in case no matching precipitation day is found, ignore precipitation for finding the most similar day
    %     ranksum=sum(rank,2,'omitnan');
    %     [mn, idx]=min(ranksum);
    %     bestrank(t,:,xloop)=[yearn(idx) doyn(idx)]; %year of best fit, doy of best fit
    %   end
    % end
    % 
    % %save data to file, added 27.03.2024
    % save(['workspace_1','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat']);
    % save(['bestrank','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat'],'bestrank');
    % load(['workspace_1','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat']);
    
    % %read hourly wfde5 reanalysis data from netcdf file
    % wfde5_all=zeros(366*24,years_wfde5,length(xcol),length(parameters_wfde5),'single');
    % wfde5_snow=zeros(366*24,years_wfde5,length(xcol),'single');
    % vloop=1;
    % for varloop=parameters_wfde5
    %   varloop=char(varloop);
    %   path=[wfde5dir,filesep,varloop];
    %   for yy=1:years_wfde5
    %     for mm=1:12
    %       year=yy+startyear_wfde5-1;
    %       filelist=dir([path,filesep,varloop,'*_',num2str(year),num2str(mm,'%.2i'),'*.nc']);
    %       filename=[filelist.name];
    %       cdir=[filelist.folder];
    %       finfo = ncinfo([cdir,filesep,filename]);
    %       varname=finfo.Variables(1,4).Name;
    %       x_nc=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lon')==1)).Length;
    %       y_nc=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'lat')==1)).Length;
    %       z_nc=finfo.Dimensions(1,find(strcmpi({finfo.Dimensions.Name},'time')==1)).Length;
    %       if(leap_year(year)==1)
    %         dayspermonth = [31 29 31 30 31 30 31 31 30 31 30 31];
    %       else
    %         dayspermonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    %       end
    %       b=sum(dayspermonth(1:mm))*24;
    %       a=b-dayspermonth(mm)*24+1;
    %       disp(['read hourly WFDE5 file ',filename]);
    %       data=ncread([cdir,filesep,filename],varname,[1 1 1],[720 360 z_nc]);
    %       for xloop=1:length(xcol)
    %         x=xcol(xloop);
    %         y=yrow(xloop);
    %         if(mask(y,x)==0)
    %           continue
    %         end
    %         %x and y for rotated wfde5 data
    %         xrot=360-y+1;
    %         yrot=x;
    %         wfde5_all(a:b,yy,xloop,vloop)=data(yrot,xrot,:);
    %       end
    %       clear data
    %       %add snowfall to rainfall
    %       if strcmpi(varloop,'Rainf')==1
    %         filename=strrep(filename,'Rainf','Snowf');
    %         cdir=strrep(cdir,'Rainf','Snowf');
    %         finfo = ncinfo([cdir,filesep,filename]);
    %         varname=finfo.Variables(1,4).Name;
    %         disp(['read hourly WFDE5 file ',filename]);
    %         data=ncread([cdir,filesep,filename],varname,[1 1 1],[720 360 z_nc]);
    %         for xloop=1:length(xcol)
    %           x=xcol(xloop);
    %           y=yrow(xloop);
    %           if(mask(y,x)==0)
    %             continue
    %           end
    %           %x and y for rotated wfde5 data
    %           xrot=360-y+1;
    %           yrot=x;
    %           wfde5_snow(a:b,yy,xloop)=data(yrot,xrot,:);
    %         end
    %         clear data
    %         wfde5_all(a:b,yy,:,vloop)=wfde5_all(a:b,yy,:,vloop)+wfde5_snow(a:b,yy,:);
    %       end
    %     end%months
    %   end%years
    %   vloop=vloop+1;
    % end%varloop

    % save(['workspace_2','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat']);
    load(['workspace_2','_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat']);    
    
    %convert wfde5 units
    for xloop=1:length(xcol)
      wfde5_all(:,:,xloop,5)=wfde5_all(:,:,xloop,5)/100; %convert surface pressure [Pa] to [hPa]
      wfde5_all(:,:,xloop,7)=wfde5_all(:,:,xloop,7)*3600; %convert mass flux density to [mm h-1]
      %convert specific humidity to relative humidity [%]
      wfde5_t=wfde5_all(:,:,xloop,1)-273.16;%[K] to [°C]
      %es according to Buck equation; t in °C; es in hPa
      a=wfde5_t;
      b=wfde5_t;
      c=wfde5_t;
      d=wfde5_t;
      %over liquid water, T>0°C
      a(wfde5_t>=0)=6.1121;
      b(wfde5_t>=0)=18.678;
      c(wfde5_t>=0)=234.5;
      d(wfde5_t>=0)=257.14;
      %over ice water, T<0°C
      a(wfde5_t<0)=6.1115;
      b(wfde5_t<0)=23.036;
      c(wfde5_t<0)=333.7;
      d(wfde5_t<0)=279.82;
      es_buck=a.*exp((b-(wfde5_t./c)).*(wfde5_t./(d+wfde5_t))); %saturation vapour pressure [hPa]
      e = wfde5_all(:,:,xloop,2) .* wfde5_all(:,:,xloop,5) ./ (0.378 .* wfde5_all(:,:,xloop,2) + 0.622);  %vapour pressure [hPa]
      rh = e ./ es_buck *100; %relative humidity [%]
      rh(rh>100)=100;
      rh(rh<0)=0;
      wfde5_all(:,:,xloop,2)=rh;
    end
    
    
    for xloop=1:length(xcol)
      x=xcol(xloop);
      y=yrow(xloop);
      if(mask(y,x)==0)
        continue
      end
      
      disp(['Assign hourly climate for sample ',num2str(xloop),' of ',num2str(length(xcol))]);
      
      %aoutdir=[outdir,filesep,'lat_',num2str(lat(xloop)),'_lon_',num2str(lon(xloop)),filesep];
      %aoutdir=strrep(aoutdir,'.','p');
      aoutdir=[outdir,filesep,'sample_',num2str(xloop),filesep];
      if(~exist(aoutdir,'dir'))
        mkdir(aoutdir)
      elseif(exist(aoutdir,'dir'))
        disp(['output directory ',aoutdir,' already exists']);
        %continue
      end
      
      lst_offset=lst_offsets(x);
      if(lstflag==1)
        midnight=24; %set to 11 pm
      else
        midnight=midnights(x);
      end
      
      %parfor (vloop=1:length(parameters_wfde5_factors),min(7,parallel_cpus))
      for (vloop=1:length(parameters_wfde5_factors))
        varloop=parameters_wfde5_factors{vloop};
        wfde5=wfde5_all(:,:,xloop,vloop);
        aa=1;
        bb=24;
        valh=zeros(length(time_cm)*24,1,'single');
        timestepsh=zeros(length(time_cm)*24,1,'double');
        for tloop=1:length(time_cm)
          disp(['processing subdaily climate ',num2str(tloop),' of ',num2str(length(time_cm)),' for ',varloop])
          time=time_cm(tloop);
          doy=day(datetime(time,'ConvertFrom','datenum'),'dayofyear');
          year=str2double(datestr(time,'yyyy'));
          yy=bestrank(tloop,1,xloop);
          year_wfde5=yy+startyear_wfde5-1;
          doy_wfde5=bestrank(tloop,2,xloop);
          if(leap_year(year_wfde5)==0 && doy_wfde5>=60)  %leap day is at position 366 (31.12)
            doy_wfde5=doy_wfde5-1;
          end
          
          if(lstflag==1) %in case of using local solar time
            a=(doy_wfde5-1)*24+1-lst_offset;
            b=a+24-1;
            if(leap_year(year_wfde5)==1)
              leap_offset=0;
            else
              leap_offset=24;
            end
            if(leap_year(year_wfde5-1)==1)
              leap_offset_m=0;
            else
              leap_offset_m=24;
            end
            if(a<1)
              wfde5_day=wfde5(1:b,yy);
              if(yy-1>=1)
                wfde5_day=[wfde5(size(wfde5,1)-leap_offset_m+a:size(wfde5,1)-leap_offset_m,yy-1);wfde5_day;];
              else %if next year does not exist, take values from next year
                wfde5_day=[wfde5(size(wfde5,1)-leap_offset+a:size(wfde5,1)-leap_offset,yy);wfde5_day;];
              end
            elseif(b>size(wfde5,1)-leap_offset)
              wfde5_day=wfde5(a:size(wfde5,1)-leap_offset,yy);
              if(yy+1<=size(wfde5,2))
                wfde5_day=[wfde5_day;wfde5(1:24-length(wfde5_day),yy+1)];
              else %if next year does not exist, take values from previous year
                wfde5_day=[wfde5_day;wfde5(1:24-length(wfde5_day),yy)];
              end
            else
              wfde5_day=wfde5(a:b,yy); %hourly data of best fit
            end
            %adjust potential radiation to local solar time
            a=1-lst_offset;
            b=24-lst_offset;
            if(a<1)
              potRad=squeeze(potRadAll(max(doy-1,1),24+a:24,xloop))';
              potRad=[potRad;potRadAll(doy,1:b,xloop)'];
            elseif(b>24)
              potRad=squeeze(potRadAll(doy,a:24,xloop))';
              potRad=[potRad;potRadAll(min(doy+1,366),1:b-24,xloop)'];
            else
              potRad=squeeze(potRadAll(doy,1:24,xloop))';
            end
          else %otherwise process data in utc
            a=(doy_wfde5-1)*24+1;
            b=a+24-1;
            wfde5_day=wfde5(a:b,yy); %hourly data of best fit
            potRad=squeeze(potRadAll(doy,1:24,xloop))'; %set potential radiation for doy
          end
          
          if(strcmpi(varloop,'pr'))
            factor=wfde5_day./sum(wfde5_day);
            if(sum(wfde5_day)==0 && data_cm(tloop,xloop,vloop)>1)
              %if no precipitation in wfde5 data, distribute precipitation > 1mm/day
              %this can be the case under very rare environments with no historical precipitation in the doy-window, e.g. in deserts
              %expand doy window to +-50 days around doy
              doysequence_precip=(doy-50)*24:(doy+50)*24;
              doysequence_precip(doysequence_precip>8784)=doysequence_precip(doysequence_precip>8784)-8784;
              doysequence_precip(doysequence_precip<1)=doysequence_precip(doysequence_precip<1)+8784;
              if(sum(wfde5(doysequence_precip,:),'all')>0)
                wfde5_pr=wfde5(doysequence_precip,:);
                wfde5_pr=wfde5_pr(:);
              else
                wfde5_pr=wfde5(:);
              end
              if(sum(wfde5_pr)>0) %if historical precipitation exists
                %linear regression between precipitation amount and duration
                wfde5_pr_ds=zeros(length(1:24:length(wfde5_pr)-24),1);
                wfde5_pr_dh=zeros(length(1:24:length(wfde5_pr)-24),1);
                n=1;
                for m=1:24:length(wfde5_pr)-24
                  wfde5_pr_ds(n,1)=sum(wfde5_pr(m:m+24-1));     %precipitation sum per day
                  wfde5_pr_dh(n,1)=sum((wfde5_pr(m:m+24-1)>0)); %hours of precipitation per day
                  n=n+1;
                end
                f = fittype('poly1');
                [myfit gof] = fit(wfde5_pr_ds,wfde5_pr_dh,f);
                coefficients=coeffvalues(myfit);
                %plot(myfit,wfde5_pr_ds,wfde5_pr_dh);
                precip_duration=min(24,round(coefficients(1)*data_cm(tloop,xloop,vloop)+coefficients(2)));
                nighttime=find(potRad==0);
                if(precip_duration<length(nighttime))
                  random_nighttime=sort(randsample(nighttime,precip_duration));  %randomly select hours
                else
                  random_nighttime=sort(randsample(1:24,precip_duration));
                end
                intensities=rand(length(random_nighttime),1);
                intensities=intensities/sum(intensities);
                wfde5_day(random_nighttime)=intensities;
                factor=wfde5_day./sum(wfde5_day);
              else %if no historical precipitation exists
                if (precipnanflag==1)
                  factor=wfde5_day./sum(wfde5_day);  %if precipitation nan flag == 1, write NaN: note that mass- and energy in not preserved in this case!
                else
                  nighttime=find(potRad==0);
                  precip_duration=randsample(1:length(nighttime),1); %randomly select duration
                  random_nighttime=sort(randsample(nighttime,precip_duration));  %randomly select hours
                  intensities=rand(length(random_nighttime),1);
                  intensities=intensities/sum(intensities);
                  wfde5_day(random_nighttime)=intensities;
                  factor=wfde5_day./sum(wfde5_day);
                end
              end
            elseif(sum(wfde5_day)==0)
              %if precipitation nan flag == 1, write NaN: note that mass- and energy in not preserved in this case!
              %if no precipitation in wfde5 data, distribute precipitation <1mm/day randomly to hour without radiation.
              %if no hour without radiation occurs, distribute precipitation to midnight
              %this should only be the case for drizzling in climate model data vs. no precipitation in wfde5
              nighttime=find(potRad==0); %all hours without sunlight
              if(isempty(nighttime))
                wfde5_day(midnight)=1; %set factor to 1 at local midnight
              else
                random_nighttime=randsample(nighttime,1);  %randomly select one hour
                wfde5_day(random_nighttime)=1;
              end
              factor=wfde5_day./sum(wfde5_day);
            end
            wfde5_pr=[];
            nighttime=[];
            random_nighttime=[];
            intensities=[];
            %clear wfde5_pr nighttime random_nighttime intensities
          else
            factor=wfde5_day./mean(wfde5_day);
          end
          
          %check for NaN-values
          if(strcmpi(varloop,'pr')==1 && precipnanflag==1)
            disp('crosscheck skipped for precipitation due to NaN flag');
          else
            if (sum(isnan(factor))>0)
              disp(['crosscheck failed, NaN values occured for ',varloop,' on ', datestr(time,'dd.mm.yyyy')]);
              %disp('execution paused, click return to continue');
              %pause
            end
          end
          
          %apply hourly fraction
          valh(aa:bb,1)=data_cm(tloop,xloop,vloop).*factor;
          timestepsh(aa:bb,1)=time:1/24:time+1-(1/24);
          
          %for temperature: scale between min and max
          if(strcmpi(varloop,'tas')==1)
            data_max=data_cm(tloop,xloop,length(parameters_isimip_read)-1);
            data_min=data_cm(tloop,xloop,length(parameters_isimip_read));
            data=valh(aa:bb,1);
            valh(aa:bb,1)=((data-min(data))/(max(data)-min(data)))*(data_max-data_min)+data_min;
          end
          
          %for temperature: correct diurnal data to conserve climate model mean while maintaining min and max
          if(strcmpi(varloop,'tas')==1)
            deltamean=data_cm(tloop,xloop,vloop)-mean(valh(aa:bb));
            dist=(min(valh(aa:bb)-data_min,abs(valh(aa:bb)-data_max))./((data_max-data_min)/2));
            valhcor=min(valh(aa:bb)+(dist./sum(dist).*deltamean*24),data_max);
            while (round(data_cm(tloop,xloop,vloop)-mean(valhcor),4)>0) %round to 4 decimals
              deltamean=data_cm(tloop,xloop,vloop)-mean(valhcor);
              dist=(min(valhcor-data_min,abs(valhcor-data_max))./((data_max-data_min)/2));
              valhcor=min(valhcor+(dist./sum(dist).*deltamean*24),data_max);
            end
            valh(aa:bb)=valhcor;
          end
          
          %for relative humidity and cloud cover: limit to max value of 100 and 1 respectively and correct mean value correspondigly
          if(strcmpi(varloop,'clt')==1 || strcmpi(varloop,'hurs')==1)
            if(strcmpi(varloop,'clt')==1)
              maxthreshold=1;
            elseif(strcmpi(varloop,'hurs')==1)
              maxthreshold=100;
            end
            if(max(valh(aa:bb))>maxthreshold)
              valhcor=min(valh(aa:bb),maxthreshold);
              while(round(data_cm(tloop,xloop,vloop)-mean(valhcor),4)>0) %round to 4 decimals
                deltamean=data_cm(tloop,xloop,vloop)-mean(valhcor);
                valhcor=min(valhcor+(deltamean*24./sum(valhcor<maxthreshold)),maxthreshold);
              end
              valh(aa:bb)=valhcor;
            end
          end
          
          %for solar radiation, check if values exceed potential maximum radiation for each hour
          if(strcmpi(varloop,'rsds')==1)
            valhor=valh(aa:bb);
            valhcor=valh(aa:bb);
            valhcor=min(valhcor,potRad);
            valhcor(valhcor<120)=valhor(valhcor<120); %WMO defined sunshine duration as direct solar irradiance > 120 W/m2. Equivalent to level of solar irradiance shortly after sunrise or before sunset in cloud-free conditions
            valhcor(valhcor<0)=0;
            while(round(data_cm(tloop,xloop,vloop)-mean(valhcor),4)>0.1)
              deltamean=data_cm(tloop,xloop,vloop)-mean(valhcor);
              n=find(valhcor>0);
              valhcor(n)=max(valhcor(n)+deltamean,0);
              valhcor=min(valhcor,potRad);
            end
            valh(aa:bb)=valhcor;
          end
          aa=aa+24;
          bb=bb+24;
        end
        
        %set negative values to 0
        valh(valh<0)=0;
        
        %aggregate hourly to n-hourly
        val=zeros(length(valh)/htimestep,1,'single');
        timesteps=zeros(length(valh)/htimestep,1,'single');
        if(htimestep>1)
          counter=1;
          for nstep=1:htimestep:length(valh)
            if(strcmpi(varloop,'pr')==1)
              val(counter)=sum(valh(nstep:nstep+htimestep-1)); %sum for precipitation
              timesteps(counter)=timestepsh(nstep);
              counter=counter+1;
            else
              val(counter)=mean(valh(nstep:nstep+htimestep-1)); %mean
              timesteps(counter,1)=timestepsh(nstep);
              counter=counter+1;
            end
          end
        else
          val=valh;
          timesteps=timestepsh;
        end
        
        %check for NaN values in the dataset
        if(strcmpi(varloop,'pr')==1 && precipnanflag==1)
          %skip
        else
          if(sum(isnan(val)==1)>0)
            test = find(isnan(val));
            teststr=datestr(double(timesteps(test(:))));
            disp(['NaN orruced in dataset ',varloop,' at timestep ',teststr]);
            %disp('execution paused, click return to continue');
            %pause
          end
        end
        
        %write output as csv file
        %formatOut='dd.mm.yyyy HH:MM:SS';
        %datestr(double(timestepsh(:,1)),formatOut)
        output=zeros(length(timesteps),5);
        datev=zeros(length(timesteps),1);
        datev=datevec(timesteps(:,1));
        output(:,1:4)=datev(:,1:4);
        output(:,5)=val;
        
        writematrix(output,[aoutdir,varloop,'_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.csv'],'Delimiter','tab');
        data_cm_var=data_cm(:,xloop,vloop);
        parsave([aoutdir,varloop,'_',model,'_',scenario,'_',num2str(startyear),'-',num2str(endyear),'.mat'],val,output,timesteps,data_cm_var);
        
      end %variable
    end %xloop
  end %model
end %scenario

%end
