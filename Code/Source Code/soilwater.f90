!!!!!!!!!!!soil initialization
	subroutine newenv(ju,lat,ele,lai,tmax,tmin,tmean,pre,&
     	ra,rh,wind,w_soil,par_soil,&
         parmax,dl,tday,t_photo,tave,tp_resp,ts_resp,tl_resp,&
     	dry_photo,dry_lit,dry_hum,par_tree,runoff,aet_pet,rn,htmax,pwp,petd)!!!!!!!!!
	real lat,ele,tmax,tmin,pre,rh,wind,w_soil(4),par_soil(8),lai,ra
	integer ju,i1
	real parmax,dl,tday,tnight,t_photo(9),tp_resp,ts_resp, tl_resp,&
          par_tree(28,9),tmean,&
     	 dry_photo(9),dry_lit,dry_hum,aet_pet,aetd,offd
	real tmp1,tmp2,tmp3,tave,tde,lat1,pi,daita,lhr,dr,ws,dl1,so,so1
	real alb,ea,es,rnl,tt,pres1,lamda1,psych1,darta1,u2,petd,tw,raa
    real f,incep,pnet,inf,npet,htmax,kcb,ke,kcmax,kr,tran
	pi=3.1415926
	yxd=2.0*pi/365.0*ju
	raa=ra
	ra=raa*0.0864
	lat1=lat*pi/180.0
	daita=0.4093*sin(yxd-1.405)
    ws=min(1.0,-tan(lat1)*tan(daita))
	ws=acos(ws)
	lhr=24.0/pi*ws
	!ra=ra*lhr/24.0
	dr=1.0+0.033*cos(yxd)
	so=38.48*dr*(ws*sin(lat1)*sin(daita)+cos(lat1)*cos(daita)*sin(ws))
	so1=so*11.574
	dl1=lhr
	dl=lhr
	parmax=0.275*raa*pi
	alb=max(0.0,0.3+lai**0.5*0.066)
	es=0.5*(satvp(tmax)+satvp(tmin))
	ea=rh*es/100.0
	rnl=(4.903e-9)*((tmax+273.16)**4+(tmin+273.16)**4)*&
        (0.17-0.07*ea**0.5)*(0.1+0.9*dl1/lhr)
	rn=max(0.0,ra*(1.0-alb)-rnl)

	tt=tmean!0.5*tmin+0.5*tmax
    if (rn .ge. 0.0)then
	    !if (tt.ge. -1.0) then
	
		    pres1=airpres(ele)

		    lamda1=flamda(tt)

		    psych1=psych(pres1,lamda1)

		    es=0.5*(satvp(tmax)+satvp(tmin))

		    ea=rh*es/100.0
	
		    darta1=slsat(tt)

		    u2=u2m(wind,10.0)
		    petd=(0.408*darta1*(rn-0.0)+psych1*900.0*u2*(es-ea)/&
     		    (tt+273.0))/(darta1+psych1*(1.0+0.34*u2))
	    !else
		   ! petd=0.0
     !   end if
    else
        petd=0.0
    endif
    
    
    f=1-exp(-0.6*lai)
    incep=(0.25*lai*(1.0-1.0/(1.0+f*pre*0.001/(0.25*lai))))*1000.0
    pnet=max(0.0,pre-incep)
    pnet=pnet/10.0
    inf=min(pnet,(par_soil(1)*0.01-w_soil(1)*0.01)*0.2*100)
    
    
    u2=u2m(wind,10.0)
    kcb=0.5
    kcmax=max(kcb+0.05,1.2+(0.04*(u2-2.0)-0.004*(rh-45.0))*((htmax/3)**0.3))
    ke=min(1.0,kcmax-kcb)
    kr=(w_soil(1)-(1.0/3.0)*pwp)/(25-(1.0/3.0)*pwp)
    aetd=petd*0.1*ke*kr
    aetd=min(w_soil(1)+w_soil(2),aetd)
    
    npet=inf-aetd
    if (npet>0.0)then
        if (pnet>inf)then
            runoff=pnet-inf
        else
            runoff=0.0
        endif
        if (w_soil(1)<20.0)then
            w_soil(1)=w_soil(1)+npet
            if (w_soil(1)>20.0)then
                w_soil(2)=w_soil(1)-20.0+w_soil(2)
                w_soil(1)=20.0
            else
                w_soil(1)=w_soil(1)
                w_soil(2)=w_soil(2)
            endif
        else
            w_soil(1)=20.0
            w_soil(2)=min(par_soil(1)-20.0,npet+w_soil(2))!
        endif
    else
        if (pnet>inf)then
            runoff=pnet-inf
        else
            runoff=0.0
        endif
        w_soil(1)=w_soil(1)+npet
        if (w_soil(1)<0.0)then
            w_soil(2)=max(0.0,w_soil(2)+w_soil(1))
            w_soil(1)=0.0
        else
            w_soil(1)=w_soil(1)
            w_soil(2)=w_soil(2)
        endif
    endif
    w_soil(3)=w_soil(1)+w_soil(2)
    
    
	!tw=pnet+w_soil(1)+w_soil(2)
	!w_soil(4)=0.0
	!if (tw.gt.petd) then
	!	aetd=petd
	!	w_soil(1)=min(par_soil(1),tw-aetd)
 !         
 !
	!	runoff=max(0.0,tw-par_soil(1))
	!
	!	w_soil(2)=min(runoff,par_soil(2))
	!else
	!	
	!	aetd=max((tw/petd)*tw,0.0)
	!	w_soil(1)=max(tw-aetd,0.0)
	!	w_soil(2)=0.0
	!end if

	tave=0.5*(tmax+tmin)
	tde=0.5*(tmax-tmin)
	tmp3=3.1415926/24

	if (dl.eq.24.0 .or. dl.eq. 0.0) then
		tday=tave
		tnight=tave
	else
		tday=tmax-tde*(24-dl)/24
		tnight=tmin+tde*dl/24
	end if
	
	do i1=1,9
		t_photo(i1)=min(tef_photo(tday,par_tree(22,i1),par_tree(23,i1),&
     						par_tree(24,i1)),1.0)	
	end do
	
	tp_resp=tef_resp(tday,tnight,dl)
	ts_resp=tef_soil(tday,tnight,dl)
	
	tl_resp=tef_resp(tday,tnight,dl)  !*(1.-dl/24.0)	

    tw=min(w_soil(1)+w_soil(4),100.0)
	if ((petd.le.1.0e-10) .or. (aetd.gt.(petd*0.01))) then
		aet_pet=1.0
	else
	
		aet_pet=min((tw/par_soil(1)+max(rh*.01-.5, .01)),1.)
	end if
			
	do i1=1,9
		if (aet_pet .ge. 0.90) then
			dry_photo(i1)=1.0
		else
			!dry_photo(i1)=min((aet_pet/par_tree(25,i1))**0.5,1.0)
		    dry_photo(i1)=min((aet_pet/par_tree(25,i1))**0.67,1.0)
		end if
    end do 
	
    tw=min(tw,par_soil(1))
	dry_hum=min(aet_pet,1.-(tw/par_soil(1)/0.6-1.0)**2)
	dry_lit=aet_pet
	return
    end    
    
    
    !!!!!!transpiration
   	subroutine tran1(ju,lat,ele,ntree,ht,parcano,lai,tmax,tmin,tmean,pre,rh,wind,tran,htmax,pwp,w_soil)!!!!!!!!!
	real lat,ele,tmax,tmin,pre,rh,wind,w_soil(4),par_soil(8),lai,ra,pwp
	integer ju,i1,ntree,ht2
	real parmax,dl,tday,tnight,t_photo(9),tp_resp,ts_resp, tl_resp,&
          par_tree(28,9),tmean,&
     	 dry_photo(9),dry_lit,dry_hum,aet_pet,aetd,offd,fj
	real tmp1,tmp2,tmp3,tave,tde,lat1,pi,daita,lhr,dr,ws,dl1,so,so1
	real alb,ea,es,rnl,tt,pres1,lamda1,psych1,darta1,u2,petd,tw,raa
    real f,incep,pnet,inf,npet,htmax,kcb,ke,kcmax,kr,tran(ntree),ht(ntree,1),parcano(120)
    
    tran=0.0
    do i1=1,ntree
        fj=max(1.0,ht(i1,1))
        ht2=int(fj)
        ra=parcano(ht2)
	    pi=3.1415926
	    yxd=2.0*pi/365.0*ju
	    raa=ra
	    ra=raa*0.0864
	    lat1=lat*pi/180.0
	    daita=0.4093*sin(yxd-1.405)
	    ws=acos(-tan(lat1)*tan(daita))
	    lhr=24.0/pi*ws
	    !ra=ra*lhr/24.0
	    dr=1.0+0.033*cos(yxd)
	    so=38.48*dr*(ws*sin(lat1)*sin(daita)+cos(lat1)*cos(daita)*sin(ws))
	    so1=so*11.574
	    dl1=lhr
	    dl=lhr
	    parmax=0.275*raa*pi
	    alb=max(0.0,0.3+lai**0.5*0.066)
	    es=0.5*(satvp(tmax)+satvp(tmin))
	    ea=rh*es/100.0
	    rnl=(4.903e-9)*((tmax+273.16)**4+(tmin+273.16)**4)*&
            (0.17-0.07*ea**0.5)*(0.1+0.9*dl1/lhr)
	    rn=max(0.0,ra*(1.0-alb)-rnl)

	    tt=tmean!0.5*tmin+0.5*tmax
        if (rn .gt. 0.0)then
	        !if (tt.ge. -1.0) then
	
		        pres1=airpres(ele)

		        lamda1=flamda(tt)

		        psych1=psych(pres1,lamda1)

		        es=0.5*(satvp(tmax)+satvp(tmin))

		        ea=rh*es/100.0
	
		        darta1=slsat(tt)

		        u2=u2m(wind,10.0)
		        petd=(0.408*darta1*(rn-0.0)+psych1*900.0*u2*(es-ea)/&
     		        (tt+273.0))/(darta1+psych1*(1.0+0.34*u2))
	        !else
		       ! petd=0.0
         !   end if
        else
            petd=0.0
        endif
        f=1-exp(-0.6*lai)
        incep=(0.25*lai*(1.0-1.0/(1.0+f*pre*0.001/(0.25*lai))))*1000.0
        pnet=max(0.0,pre-incep)
        pnet=pnet/10.0
        inf=min(pnet,(par_soil(1)*0.01-w_soil(1)*0.01)*0.2*100)
    
    
        u2=u2m(wind,10.0)
        kcb=0.5
        !tran=kcb*petd*0.1
        kcmax=max(kcb+0.05,1.2+(0.04*(u2-2.0)-0.004*(rh-45.0))*((htmax/3)**0.3))
        ke=min(1.0,kcmax-kcb)
        kr=(w_soil(1)-(1.0/3.0)*pwp)/(25-(1.0/3.0)*pwp)
        aetd=petd*0.1*ke*kr
        aetd=min(w_soil(1)+w_soil(2),aetd)
        tran(i1)=0.5*aetd
    enddo
    end
    