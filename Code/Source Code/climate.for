c===================================================================
c-------------------------------------------------------------------
	subroutine env(ju,lat,ele,lai,tmax,tmin,tmean,pre,
     ^	ra,rh,wind,w_soil,par_soil,
     ^    parmax,dl,tday,t_photo,tave,tp_resp,ts_resp,tl_resp,
     ^	dry_photo,dry_lit,dry_hum,par_tree,runoff,aet_pet,rn)!!!!!!!!!
	real lat,ele,tmax,tmin,pre,rh,wind,w_soil(4),par_soil(8),lai,ra
	integer ju,i1
	real parmax,dl,tday,tnight,t_photo(9),tp_resp,ts_resp, tl_resp,
     ^     par_tree(28,9),tmean,
     ^	 dry_photo(9),dry_lit,dry_hum,aet_pet,aetd,offd
	real tmp1,tmp2,tmp3,tave,tde,lat1,pi,daita,lhr,dr,ws,dl1,so,so1
	real alb,ea,es,rnl,tt,pres1,lamda1,psych1,darta1,u2,petd,tw,raa
	pi=3.1415926
	yxd=2.0*pi/365.0*ju
	raa=ra
	ra=raa*0.0864
	lat1=lat*pi/180.0
	daita=0.4093*sin(yxd-1.405)
	ws=acos(-tan(lat1)*tan(daita))
	lhr=24.0/pi*ws
	ra=ra*lhr/24.0
	dr=1.0+0.033*cos(yxd)
	so=38.48*dr*(ws*sin(lat1)*sin(daita)+cos(lat1)*cos(daita)*sin(ws))
	so1=so*11.574
c========06_04_23 把下列过程修改为后附
c	dl1=lhr*max((ra/so-0.25),0.3)/0.77
c	dl=min(dl1,lhr)
c	parmax=raa*lhr*pi/dl*0.55
c========06_04_23 把下列过程修改为后附
c========06_04_23
	dl1=lhr
	dl=lhr
	parmax=0.275*raa*pi
      !parmax=raa!!!!!!!!!!!!!!!
c========06_04_23
	alb=0.3-lai**0.5*0.066
	es=0.5*(satvp(tmax)+satvp(tmin))
	ea=rh*es/100.0
	rnl=(4.903e-9)*((tmax+273.16)**4+(tmin+273.16)**4)*
     ^    (0.17-0.07*ea**0.5)*(0.1+0.9*dl1/lhr)
	rn=ra*(1.0-alb)-rnl
c===============================
	tt=tmean!0.5*tmin+0.5*tmax
	if (tt.ge. -1.0) then
c	
		pres1=airpres(ele)
c
		lamda1=flamda(tt)
c
		psych1=psych(pres1,lamda1)
c
		es=0.5*(satvp(tmax)+satvp(tmin))
c
		ea=rh*es/100.0
c	
		darta1=slsat(tt)
c
		u2=u2m(wind,10.0)
		petd=(0.408*darta1*rn+psych1*900.0*u2*(es-ea)/
     ^		(tt+273.0))/(darta1+psych1*(1.0+0.34*u2))
	else
		petd=0.0
	end if
c================================
	tw=pre+w_soil(4)+w_soil(1)+w_soil(2)
	w_soil(4)=0.0
	if (tw.gt.petd) then
		aetd=petd
		w_soil(1)=min(par_soil(1),tw-aetd)
          
c
		runoff=max(0.0,tw-par_soil(1))
c		
		w_soil(2)=min(runoff,par_soil(2))
	else
c		
		aetd=max((tw/petd)*tw,0.0)
		w_soil(1)=max(tw-aetd,0.0)
		w_soil(2)=0.0
	end if
c================================
	tave=0.5*(tmax+tmin)
	tde=0.5*(tmax-tmin)
	tmp3=3.1415926/24
c----
	if (dl.eq.24.0 .or. dl.eq. 0.0) then
		tday=tave
		tnight=tave
	else
		tday=tmax-tde*(24-dl)/24
		tnight=tmin+tde*dl/24
	end if
c	
	do i1=1,9
		t_photo(i1)=min(tef_photo(tday,par_tree(22,i1),par_tree(23,i1),
     ^						par_tree(24,i1)),1.0)	
	end do
c	
	tp_resp=tef_resp(tday,tnight,dl)
	ts_resp=tef_soil(tday,tnight,dl)
c	
	tl_resp=tef_resp(tday,tnight,dl)  !*(1.-dl/24.0)	
c------------------------------------------
c
	if ((petd.le.1.0e-10) .or. (aetd.gt.(petd*0.01))) then
		aet_pet=1.0
	else
c	
		aet_pet=min((w_soil(1)/par_soil(1)+max(rh*.01-.5, .01)),1.)
	end if
c			
	do i1=1,9
		if (aet_pet .ge. 0.90) then
			dry_photo(i1)=1.0
		else
c			dry_photo(i1)=min((aet_pet/par_tree(25,i1))**0.5,1.0)
		dry_photo(i1)=min((aet_pet/par_tree(25,i1))**0.67,1.0)
		end if
	end do 
c	
	dry_hum=min(aet_pet,1.-(w_soil(1)/par_soil(1)/0.6-1.0)**2)
	dry_lit=aet_pet
c    在06－04－25前
c    dry_hum=max(aet_pet,1.-(w_soil(1)/par_soil(1)/0.6-1.0)**2)
c    在06－04－25前
	return
	end
c====================================================================
C     GAUSS RETURNS A NORMAL RANDOM NUMBER, X,S=(0,1).
C     THIS IS ROUTINE GASDEV(IDUM) FROM PRESS ET AL., ABOVE.
C
      SUBROUTINE GAU(GAUSS)
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
   21   CALL RANDOM(RAN1)
	CALL RANDOM(RAN2)
	V1=2.0*RAN1-1
	V2=2.0*RAN2-1
	R=V1*V1+V2*V2
	IF (R.GE.1.) GO TO 21
	FAC=SQRT(-2.*ALOG(R)/R)
	GSET=V1*FAC
	SS=V2*FAC
	ISET=1
      ELSE
	SS=GSET
	ISET=0
      ENDIF
	GAUSS=MIN(MAX(SS,-1.0),1.0)
      RETURN
      END
C
c     ============================================================
c	air pressure 
c     -------------
c	by Yan Xiaodong 2002
c
	real function airpres(ele)
c	output air pressure [kP]
c	ele: altitude [m]
c
	real ele
	airpres=101.3*(1.0-0.0065*ele/293.0)**5.26
	return
	end
c
c	=============================================================
c	psychmeter constant
c	-------------------
c	by Yan Xiaodong 2002
c	
	real function psych(airpres, lamda)
c	output psychrometric constant [kP/oC]
c	lamda: latent heat of vaporization [MJ/kg],2.45
c	airpres: air pressure [kP]
c
	real lamda,airpres
	psych=0.001013/0.622*airpres/lamda
c	psych=0.00665*airpres
	return
	end
c
c	==============================================================
c	from daily temperature to lamda(mj/kg)
c	-------------------------------------- 
c	Yan xiaodong  1997      	
c
      real function flamda(temd)
	real temd
	flamda=(2501-2.3601*temd)/1000.0
	return
	end
c
c	==============================================================
c	saturation vapor pressure
c	--------------------------
c	By Yan Xiaodong 2002
c	
	 real function satvp(t)
c	output saturation vapor press [kP]
c	t: temperature [oC]
c
	real t
	satvp=0.6108*exp(17.27*t/(237.3+t))!0.6108*exp(7.45*t/(241.9+t))
	return
	end
c
c	=============================================================
c	actual vapor pressure
c	----------------------
c	By Yan Xiaodong 2002
c
	real function actvp(satvp,rh)
c	output actual vapor pressure [kP]
c	satvp: saturation vapor pressure [kP]
c	rh: relative humidity [%]
	real rh, satvp
	actvp=rh*satvp/100.0
	return
	end
c
c	=============================================================
c	slope of saturation-temperature relation
c	----------------------------------------
c	By Yan Xiaodong
c
	real function slsat(t)
c	output slope of saturation vapor pressure curve [kP/oC]
c	t: temperature [oC]
	real t
	slsat=4098.0*satvp(t)/(t+237.3)**2
	return
	end
c
c	============================================================
c	day length, noon-radiation, net radiation
c	------------------------------------------
c	By Yan Xiaodong 2002
c
	subroutine rad(lat,ele,lai,julia,ra,tmin,tmax,ea,rh,dl,rn)
c	output radiation relative variables
c	lat: latitude [ o]
c 	julia: number of day
c	ele: elevation [m]
c	lai:leaf area index [-]
c	ra: extraterrestrial radiation [MJ/m2/day]
c	dl: day light length [hours]
c	ramax: maxmum radiation at noon (mj/m2/day)
c	rn: net radiation (mj/m2/day)
c	1367[W/m2]=4.9212[MJ/m2/h]
c	1[W/m2]=0.0684[MJ/m2/day]
	real lat, ele,ea,lai,tmax,tmin, darta,sigma,pai,omiga,rh,rn
	real ra,dr,tmp1,tmp2,rs0,dl,ra1,lat1,raday
	integer julia
c	from [w/m2] to Mj/m2/day
	raday=ra*0.0684
c=====
	es=0.5*(satvp(tmax)+satvp(tmin))
c	actual vapor pressure
	ea=rh*es/100.0

	ea=0.5
c=====
	pai=3.1415926
	lat1=lat*pai/180.0
c	sun-earth relative distance
	dr=1.0+0.033*cos(2*pai/365.0*float(julia))
c	solar declination
	darta=0.409*sin(2*pai*float(julia)/365.0-1.37)
c	day length 
	tmp1=-tan(darta)*tan(lat1)
	if (tmp1.le. -0.99999) then
		tmp2=pai
	elseif (tmp1.ge. 0.99999) then
		tmp2=0.0
	else	
		tmp2=acos(tmp1)
	end if	
	omiga=tmp2
c	extraterrestrial radiation (mj/m2/day)
	ra1=37.58603*cos(lat1)*cos(darta)*(sin(omiga)-omiga*cos(omiga))
c	daylength (hours)
	dl=7.639437*omiga
c	radiation at noon (mj/m2/min)
	ramax=0.0820*dr*cos(lat1-darta)
c	surface radiation related to cloudness
c	rs=ra1*(0.85-0.47*clou)
	rs=raday*dl
c	surface radiation also be calculated from actual day 
c	light hours (adl) according Angstrom formula
c	rs=ra*(0.25+0.5*adl/dl)
c	ramax=ramax*ra/ra1
c	clearday surface radiation: related to elevation
	rs0=(0.75+0.00002*ele)*ra1
c	sigma: reflect ratio related to leaf area index
	sigma=0.3-lai**0.5*0.066
c	net short wave radiation [mj/m2/day]
	rns=(1-sigma)*raday
c	net long wave radiation 
	rnl=4.903e-9*((tmax+273.16)**4+(tmin+273.16)**4)*
     ^    (0.17-0.07*ea**0.5)*(1.35*raday/rs0-0.35)
c	net radiation
	rn=rns-rnl
c
	return
	end
c
c	===============================================================
c	wind speed at 2 meters
c	------------------------
c	By Yan Xiaodong 2002
c
	real function u2m(uzm,z)
c	output wind speed [m/s]
c	uzm: at z height, wind speed [m/s]
c	z: height [m]
	real uzm,z
	u2m=uzm*4.87/log(67.8*z-5.42)
	return
	end
c
c	=============================================================
c	potential evaporation
c	-----------------------
c	By Yan Xiaodong 2002
c 
	real function petd(lat,ele,tmax,tmin,humi,uzm,rn)
c	output daily potential evaporation [mm/day]
	real lat,ele,tmax,tmin,humi,rn
	real ea, es, lamda1, psych1, pres1, darta1, tt, u2
c
c	daily temperature
	tt=0.4*tmin+0.6*tmax
	if (tt.ge. -1.0) then
c	air pressure
	pres1=airpres(ele)
c	lamda
	lamda1=flamda(tt)
c	psychmeter 
	psych1=psych(pres1,lamda1)
c	saturation vapor press
	es=0.5*(satvp(tmax)+satvp(tmin))
c	actual vapor pressure
	ea=humi*es/100.0
c	darta
	darta1=slsat(tt)
c	wind speed at 2m
	u2=u2m(uzm,10.0)
	petd=(0.408*darta1*rn+psych1*900.0*u2*(es-ea)/
     ^	(tt+273.0))/(darta1+psych1*(1.0+0.34*u2))
	else
	petd=0.0
	end if
	return
	end
c===========================================================
c
	real function tef_photo(t,tmin,topt,tmax)
	real t,tmin,topt,tmax
	if(t.ge. tmax .or. t.le.tmin) then
		tef_photo=0.0
	else
		tef_photo=((tmax-t)/(tmax-topt))**((tmax-topt)/(tmax-tmin))*
     ^			  ((t-tmin)/(topt-tmin))**((topt-tmin)/(tmax-tmin))
	end if
	if (t.le.5.0) then
		tef_photo=0.0
	end if
	return
	end
c-----
	real function tef_resp(tday,tnight,dl)
	real tday,tnight,dl
	real tg1,tg2
	if (tday.gt. 0.0) then
		tg1=2.0*exp(-0.009*(tday-15.0))
	    tef_resp1=dl/24.0*exp(log(tg1)/10.0*(tday-15.0))
	else 
		tef_resp1=0.0
	end if
	if	(tnight.gt. -2.5) then
		tg2=2.0*exp(-0.009*(tnight-15.0))
	    tef_resp2=(24.0-dl)/24.0*exp(log(tg2)/10.0*(tnight-15.0))
	else
		tef_resp2=0.0
	end if
	tef_resp=tef_resp1+tef_resp2
	return
	end
c-----
	real function tef_soil(tday,tnight,dl)
	real tday,tnight,dl,tef1,tef2
	if (tday.le.1.5) then
		tef1=0.0
	else
		tef1=dl/24.0*exp(3.36*(tday-40.0)/(tday+31.79))
	end if
	if (tnight.le. 1.5) then
		tef2=0.0
	else
		tef2=(24.0-dl)/24.0*exp(3.36*(tnight-40.0)/(tnight+31.79))
	end if
		tef_soil=tef1+tef2
	return
	end
c-----
	subroutine pheno(t,phenology)
	real t(365),phenology(366)
	integer i,n,j
	j=0
	n=0
	do i=1,187
		if (t(i).gt.5.0) then
			j=1+j
		else
			j=0
		end if
		if (j.eq. 10) then
			n=i-9
	    end if
	end do
	do i=1,187
		if (i.le. n) then
			phenology(i)=0.0
		else
			phenology(i)=1.0
		end if
	end do
c	
	j=0
	n=366
	do i=188,365
		if (t(i).lt.5.0) then
			j=1+j
		else
			j=0
		end if
		if (j.eq. 10) then
			n=i-9
	    end if
	end do
	do i=188,365
		if (i.lt. n) then
			phenology(i)=1.0
		else
			phenology(i)=0.0
		end if
	end do
	phenology(366)=phenology(355)
c
	return
	end
c------	
