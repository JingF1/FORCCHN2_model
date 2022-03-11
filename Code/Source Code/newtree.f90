    !!!!!!tree height
    real function ht(dbh)!d是胸径,wang and fang的树高经验公式
    real dbh
    ht=1./(0.82*((dbh*100.)**1.25))+1./37.26!树径单位为米时，化为厘米
    ht=1./ht
    return
    end
    
    !!!!!!sapwood
    real function sapwood(dbh,h)
    real dbh,h
    sapwood=0.000102*(100.*dbh**0.750383)*(h**1.673065)
    return
    end
    
    !!!!sapwood carbon  
    real function  sapwoodc(sw,ck)
    real sw,ck
    sapwoodc=sw*ck
    return
    end
    
    
    !!!!!twig length
    real function stemlen(bd)
    real bd
    stemlen=1.054+0.854*log(100.*bd)
    if (bd.eq.0.)then
        stemlen=0.
    endif
    if (stemlen<0.) then
        stemlen=0.
    endif
    return
    end
    
    
    !!!!!twig diameter
    real function stembd(ht,dbh,dinc,bh)
    real ht,dbh,hd,dinc,rdinc,bh
    dinc=min(dinc,ht-bh)
    rdinc=dinc/(ht-bh)
    hd=ht/(dbh*100.0)
    stembd=max(0.01,0.791+0.009*ht*dinc+0.364*log(rdinc)+0.021*dbh*(100.0))!max(0.0,0.0306*ht-0.9495*hd+1.2732-(0.7663*log(dbh*100.0))*(exp(-0.3733*dinc)))!
    return
    end
    
    
    !!!!!
    real function sleaf(bd,dinc)
    real bd,dinc
    if (dinc<0.7)then
        sleaf=1.904+22.428*log(bd*100.)
    else
        sleaf=1.877+23.631*log(bd*100.)
    endif
    if (sleaf<0.)then
        sleaf=0.
    endif
    return
    end
    
    !!!!!
    real function ccf(forcn,ntrees)
    real forcn(600,19),lcw,smca
    integer i,ntrees
    ccf=0.
    do i=1,ntrees
       lcw=1.306+0.148*forcn(i,2)*100
       smca=((3.1415/4.)*(lcw**2.))
       ccf=ccf+smca
    enddo
    return
    end
    
    !!!!
    real function ccfl(forcn,ntrees,dbh)
    real forcn(600,19),lcw,smca
    integer i,ntrees
    ccfl=0.
    do i=1,ntrees
       if(forcn(i,2)>dbh)then
           lcw=1.306+0.148*forcn(i,2)*100
           smca=((3.1415/4.)*(lcw**2.))
       else
           smca=0.0
       endif
       ccfl=ccfl+smca
    enddo
    return
    end
    
    !!!!basal area in 1 hectare
    real function ba(forcn,ntrees,gaparea)
    real forcn(600,19),gaparea,bal
    integer i,ntrees
    ba=0
    do i=1,ntrees
       bal=3.1415*((forcn(i,2)/2.)**2.)*10000.0/gaparea
       ba=ba+bal
    enddo
    return
    end
    
    !!!!!
    real function bal(forcn,dbh,ntrees,gaparea)
    real forcn(600,19),dbh,gaparea,bat
    integer i,ntrees
    do i=1,ntrees
        if(forcn(i,2)>dbh)then
            bat=3.1415*((forcn(i,2)/2.)**2.)*10000.0/gaparea
        else
            bat=0.0
        endif
        bal=bal+bat
    enddo
    return
    end
    
    !!!!!twig height
    real function hcb(ba,ht,ccf,dbh,bal)
    real ba,ht,ccf,dbh,bal
    hcb=-7.8655+0.5217*ht+0.0896*bal+0.3999*ba-0.0094*ccf-0.0735*(dbh*100)
    if (hcb>ht)then
        hcb=0.5*ht
    endif
    if (hcb<=0.)then
        hcb=0.5*ht
    endif
    return
    end
    
    
    !!!!!mamxium leaf biomass
    real function wl(dbh,hcb)
    real dbh,hcb
    wl=0.0199*((dbh*100)**2.1321)*(hcb**(-0.3890))!*0.5
    return
    end
    
    !!!!!twig carbon
    real function wb1(bd)
    real bd
    wb1=0.039+0.012*log(bd*100.)!*0.5
    return
    end
    
    !!!!!
    real function nb(dinc,ht,hcb)
    real rdinc,a,b
    a=24.235
    b=1.574
    rdinc=dinc/(ht-hcb)
    nb=a*exp(-b*rdinc)
    return
    end
    
    !!!!!
    real function wb(dbh,ht)
    real dbh,ht
    wb=0.0189*((dbh*100.)**3.2361)*(ht**(-0.9692))!*0.5
    return
    end
    
    !!!!!!
    real function can(dbh,ccfl,ba,n)
    real dbh,ccfl,ba,n
    can=max(0.0,4.758+0.1567*(dbh*100.0)-0.0029*ccfl-0.1635*ba+0.0003*n)
    end
    
    !!!!
    real function blm(dinc,ba,cw,h,d)
    real dinc,cw,cr,hd,h,d
    hd=h/(d*100.0)
    cr=max(0.0,0.9522+0.4206*log(ba)-0.4509*log(cw)-0.0348*hd+0.0003*((d*100.0)**2.0)&
       -0.0753*log(h))
    blm=max(0.0,-0.112+0.962*log(dinc)+0.086*cw-0.809*cr-0.06*dinc-0.518*hd)
    return
    end
    
    !!!!!!!!!Light competition (radiation distribution)
    subroutine newphenology(forcn,ntrees,tday,par_tree,parmax,&
     				gaparea,lai,parcano,key,sn,ss,sd,sr,snc,ssc,sdc,src,phenfac)
	real forcn(600,19),dry_photo,tday,par_tree(28,9),parmax,gaparea
	real lai,laicano(120),parcano(120),key(900,3),yxd(900),&
         phenfac,yy
	integer ntrees
	integer it,i1,i2,i3,i5
	real tmp1,tmp2,tmp5(900),con(120),bro(120),con1(120),bro1(120)
    real sn,ss,sd,sr,snc,ssc,sdc,src

	do i1=1,120
		laicano(i1)=0.0
		parcano(i1)=parmax
		con(i1)=0.0
		bro(i1)=0.0
	end do

	do it=1,ntrees
		forcn(it,16)=phenfac
	end do


	tmp1=0.0
	do it=1,ntrees
		if (forcn(it,2).gt. 0.01) then
			i1=int(forcn(it,1))
			!if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then
			
				!if ((sn.ge.0.0 .and. sn.le.snc ).and. &
    !                (forcn(it,6).le. forcn(it,12)*0.01)) then
				!	yy=max(forcn(it,9)*0.9,forcn(it,12)*0.025)
				!	forcn(it,6)=yy/1.5+forcn(it,6)
				!	forcn(it,8)=yy/3.0+forcn(it,8)
				!	forcn(it,9)=forcn(it,9)-yy
				!	forcn(it,10)=forcn(it,10)-&
    ! 							yy/par_tree(9,i1)
				!	forcn(it,16)=0.0
				!end if
			!end if
	
			forcn(it,7)=max(flai(forcn(it,6),par_tree(17,i1)),1.0e-15)
			yxd(it)=(forcn(it,7)+flai(forcn(it,12),&
     			par_tree(17,i1)))*0.5
			tmp5(it)=max((yxd(it))**0.5*forcn(it,2),0.000001)
			tmp1=tmp1+tmp5(it)
			i2=min(int(forcn(it,3)),120)
			i3=max(min(int(forcn(it,5)),int(forcn(it,3))),1)
			tmp2=forcn(it,7)/gaparea/float(i2-i3+1)
			do i5=i3,i2
				if (i1.eq.7 .or. i1.eq.8 .or. i1.eq.9) then
					con(i5)=con(i5)+tmp2
				else
					bro(i5)=bro(i5)+tmp2
				end if
			end do
		else
			tmp5(it)=0.0
		end if
	end do

	if (tmp1.gt. 1.0e-15) then
		do it=1,ntrees
			forcn(it,13)=tmp5(it)/tmp1*gaparea
		end do
	else 
		do it=1,ntrees
		    forcn(it,13)=1.0/float(ntrees)*gaparea
		end do
	end if

	lai=0.0

	do i5=1,119
		lai=lai+con(i5)+bro(i5)
		con(120-i5)=con(120-i5+1)+con(120-i5)
		bro(120-i5)=bro(120-i5+1)+bro(120-i5)
	end do

	do i5=1,119
	    parcano(i5)=max(parmax*exp(-0.35*con(i5+1)-0.40*bro(i5+1)),10.)
    end do
    
	do it=1,ntrees
		i2=max(min(int(forcn(it,3)),119),1)
		i3=max(int(forcn(it,5)),1)

		key(it,1)=key(it,1)+parcano(i3)*phenfac
				key(it,2)=key(it,2)+phenfac*&
          (parcano(i2+1)-parcano(i2))/(parcano(i2+1))
	end do
	return
    end            