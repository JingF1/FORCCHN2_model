function [h]=photoperiod(lat,doy)
pi=3.1415;
lat1=lat*pi/180.0;
darta=0.409*sin(2.0*pi*doy/365.0-1.37);
tmp1=-tan(darta)*tan(lat1);
if tmp1<= -0.99999
    tmp2=pi;
elseif tmp1>= 0.99999
    tmp2=0.0;
else
    tmp2=acos(tmp1);
end
omiga=tmp2;
h=7.639437*omiga;
end
