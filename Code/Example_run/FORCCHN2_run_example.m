%%%%%%%%%a example of FORCCHN2model run
%%%%%%%%%this DLL package is a 64-bit program, need the support of C++ compiler
%%%%%%%%%this example is in the Harvard Forest
%%%%%%%%%you can change the different input data to calculate the forests which you need
clear;
clc;
name1=('XXX\');%%load path of the FORCCHN2 DLL package
name2=[name1,'FORCCHN2_64.dll'];%%we use the 64-bit
name3=[name1,'FORCCHN2.h'];%%load header file
loadlibrary(name2,name3);%%load the DLL package

%%%%%%%input year information
ny0=1991;%%initial year
ny1=2012;%%end year
ndays=8036;%%total days, including the days of leap years
ny=ny1-ny0+1;%%total years

%%%%%%%input geo information
lat=42.538;%%latitude,[degree]
lon=-72.171;%%longtitude,[degree]
ele=340.0;%%elevation,[meter]

%%%%%%all initialization data are the data in first year
%%%%%%input vegetation information, initialization data
ntrees=0;%%initial trees' number, set as 0 to input
lai0=7.0;%%maximum leaf area index,[m2/m2]
evergr0=0.139;%%evergreen trees percent
deci0=1.0-evergr0;%%deciduous trees percent
class1=104;%%forest type,101:DNF,102:ENF,103:MF,104:DBF,105:EBF
fj=zeros(ny,3);%%leaf phenology dates,1.year,2.SOS,3.EOS
fj(1,1)=ny0;%%the first year
fj(1,2)=127;%%leaf growth start time (SOS) in the first year,[DOY]
fj(1,3)=257;%%leaf growth end time (EOS) in the first year,[DOY]

%%%%%%%input soil data, initialization data
sfc=28.434;%%soil field capacity,[cm]
pwp=7.244;%%permanent wilting point,[cm]
vw=1274.32;%%soil volume weight,[kg/m3] 
sc0=18.067;%%soil total organic carbon,[kg C/m2]
sn0=0.775;%%soil total nitrogen,[kg C/m2]
silt=15.0;%%soil silt percent
sand=100-silt;%%soil san percent 
yxc=[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.88]';%%initial parameters of the soil carbon pool

%%%%%%%input climate data, driven data
name11=('XXX\US-Ha1_climate.txt');%%load the path of climate data
cli=importdata(name11);%%climate data(ndays,11),1.year,2.DOY,3.mean temperature [oC],4.maximum temperature [oC],5.minimum temperature [oC],6.air pressure [hPa],7.wind [m/s],8.relative humidity [%],9.precipitation [mm],10.shortwave radiation [W/m2],11.carbon dioxide concentration [ppm]
cli1=zeros(ny,366,9);%%calculate climate format to ny*366
day=0;
for i=1:ny
    iy=ny0+i-1;
    if mod(iy,4)==0
        day1=366;
    else
        day1=365;
    end
    day=day+day1;
    cli1(i,1:day1,:)=cli(day-day1+1:day,3:11);    
end
tmean=cli1(:,:,1);
tmax=cli1(:,:,2);
tmin=cli1(:,:,3);
wind=cli1(:,:,5);
rh=cli1(:,:,6);
prec=cli1(:,:,7);
ra=cli1(:,:,8);
co2=cli1(:,:,9);
pho=zeros(ny,366);%%use the latitude to calculate photoperiod
for i=1:366
    pho1=photoperiod(lat,i);
    pho(:,i)=pho1;
end

%%%%%%%output 
%%dayout:daily dynamics,1.aboveground biomass [kg C/m2],2.belowground biomass [kg C/m2],3.GPP [kg C/m2],4.aboveground respiration [kg C/m2],5.belowground respiration [kg C/m2],6.soil heterotrophic respiration [kg C/m2],7.litter-fall [kg C/m2],8.soil carbon pool [kg C/m2]
%%yearout:yearly dynamics,variables are consistent with daily dynamics
dayout=zeros(ndays,8);
yearout=zeros(ny,8);

%%%%%%%run model with DLL file
%%%%%%%can choose four results,fj(ny,3):phenology dates,
%%%%%%%yxc(11,1):parameters of the soil carbon pools, you can input this variable to the next input to get the balance statement of the soil
%%%%%%%dayout(ndays,8):daily dynamics, yearout(ny,8):yearly dynamics,
[fj,yxc,dayout,yearout]=calllib('FORCCHN2_64','forcchn2',fj,yxc,dayout,yearout,ntrees,ny0,ny,ndays,lat,lon,ele,tmax,tmin,tmean,pho,prec,ra,rh,wind,sfc,pwp,vw,sc0,sn0,silt,sand,class1,evergr0,deci0,lai0,co2);
unloadlibrary FORCCHN2_64;