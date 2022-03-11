!=================================================================
!-----------------------------------------------------------------
	subroutine initialer2(trh,sc,sn,sfc,pwp,vw,silt,sand,class,evergr,deci,
     ^  lai,lat,ele,gaparea,par_soil,c_soil,n_soil,l_soil,w_soil,
     ^  hold0,wata,forcn,ntrees,par_tree,yxc,tever)
c
	real par_soil(8),c_soil(15),n_soil(15),w_soil(4),hold0,wata,
     ^	 forcn(600,19),par_tree(28,9),l_soil(15)
	real sc,sn,sfc,pwp,silt,sand,evergr,deci,lai,lat,ele,gaparea,
     ^above1,above2,below
	integer ntrees,class
	integer i1,i2,i3,i4,itrees
	real tmp1,tmp2,tmp3,pro(11),pdeci,pever,tever,gao1,yxd,lai1
	real yxc(11),yxn(11),sn0
	!data yxc/0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.88/
	data yxn/50.0,100.0,50,100.0,200.0,200.0,200.0,9.0,9.0,9.0,9.0/
!    forcn (ntrees,19)
!	   1: species
!	   2: base diameter
!	   3: tree height
!    4: twig diameter
!    5: twig height
!    6: leaf biomass
!    7: leaf area
!    8: fine root biomass
!    9: NSC pool
!   10: nitrogen
!   11: wood increment
!   12: maximum possible leaf biomass
!   13: tree area
!   14: carbon balance
!   15: nitrogen balance
!   16: none
!   17: stem biomass
!   18: twig biomass
!   19: coarse root biomass				
!
	hr=2.0
c
	gaparea=400.00
c
	ntrees=int(gaparea)
c
	par_soil(1)=(sfc-pwp)*15.0    
	par_soil(2)=50.0
c
	par_soil(4)=max(silt/100.0,0.00001)
	par_soil(5)=max(sand/100.0,0.00001)
	par_soil(6)=max(1.0-par_soil(4)-par_soil(5),0.00001)
c
	par_soil(7)=1.0
c
	sc1=max(sc,0.1)
	sn1=sc1*max(0.05,sn/sc)
	sc=sc1
	sn=sn1
c
	par_soil(8)=max(sc/sn-2.0,5.0)
c-----
	do i1=1,11
		c_soil(i1)=yxc(i1)*sc
		l_soil(i1)=0.5*c_soil(i1)
	end do
	sn0=0.0
	do i1=1,7
	n_soil(i1)=c_soil(i1)/yxn(i1)
	sn0=sn0+n_soil(i1)
	end do			
	n_soil(8)=c_soil(8)/par_soil(8)
	n_soil(9)=c_soil(9)/par_soil(8)
	sn0=sn0+n_soil(8)+n_soil(9)
	sn0=sn-sn0
	n_soil(10)=sn0*c_soil(10)/(c_soil(10)+c_soil(11))
	n_soil(11)=sn0-n_soil(10)
c
		w_soil(1)=par_soil(i1)*0.5
		w_soil(2)=0.0
		w_soil(3)=0.0
	hold0=5.0
	wata=5.0
c	
	w_soil(4)=max(lat+ele/100.0-35.0,0.0)*3.0
c
	do i1=1,11
		pro(i1)=0.0
	end do
cc
c	if (class.eq.101) then
c		evergr=0.0
c		deci=1.0
c		lai=max(2.0,lai)
c	end if
c	if (class.eq.102) then
c		evergr=1.0
c		deci=0.0
c	end if
c	if (class.eq.103) then
c		evergr=0.5
c		deci=0.5
c	end if
c	if (class.eq.104) then
c		evergr=0.0
c		deci=1.0
c		lai=max(lai,2.0)
c	end if
c	if (class.eq.105) then
c		evergr=1.0
c		deci=0.0
c	end if
c	if (class.eq.106) then
c		evergr=1.0
c		deci=0.0
c	end if
c	if (class.eq.107) then
c		evergr=0.7
c		deci=0.3
c		lai=max(2.0,lai)
c	end if
c	if (class.eq.108) then
c		evergr=1.0
c		deci=0.0
c	end if
c	if (class.eq.109) then
c		evergr=1.0
c		deci=0.0
c	end if
c
	if (class.eq.102.or.class.eq.105.or.class.eq.108.or.class.eq.109)then
          evergr=1.0
          deci=0.0
      elseif (class.eq.101.or.class.eq.104)then
          evergr=0.0
          deci=1.0
      endif
      
      evergr=lai*evergr
      deci=lai*deci
      tmp1=evergr+deci
c
	lai1=deci+evergr
c
	pever=evergr/lai1
	pdeci=1.0-pever
c
	tmp1=lat+ele/200.0
	if (class.eq. 0) then
		if (tmp1.le. 15.0) then
			class=109
		end if
		if (tmp1.le.30.0 .and. tmp1.gt.20.0) then
			class=105
		end if
          if (tmp1.le.30.0 .and. tmp1.gt.40.0) then
			class=104
		end if
		if (tmp1.le.45.0 .and. tmp1.gt. 40.0) then
			class=103
		end if
		if (tmp1.le.55.0 .and. tmp1.gt.45.0) then
			class=102
		end if
		if (tmp1.gt. 55.0) then
			class=101
		end if
	end if  
c
	if (class.eq.109) then 
c
		if (lai1 .le. 4.0) then
			pro(1)=pever*0.35
			pro(2)=pever*0.35
			pro(3)=pever*0.15
			pro(4)=pever*0.15
			pro(5)=pdeci*0.5
			pro(6)=pdeci*0.5
		else
			pro(1)=pever*0.65
			pro(2)=pever*0.05
			pro(3)=pever*0.25
			pro(4)=pever*0.05
			pro(5)=pdeci*0.90
			pro(6)=pdeci*0.10
		end if
	end if
c
	if (class.eq.108 .and. lai1.le.4.0) then
		if (lai1.le. 4.0) then
			pro(2)=pever*0.4
			pro(1)=pever*0.4
			pro(3)=pever*0.1
			pro(4)=pever*0.1
			pro(5)=pdeci*0.5
			pro(6)=pdeci*0.5
		else
			pro(2)=pever*0.10
			pro(1)=pever*0.70
			pro(3)=pever*0.15
			pro(4)=pever*0.05
			pro(5)=pdeci*0.8
			pro(6)=pdeci*0.2
		end if
	end if
c
	if (class.eq.107) then
		if (lai1.le. 4.0) then
			pro(3)=pever*0.5
			pro(4)=pever*0.5
			pro(5)=pdeci*0.5
			pro(6)=pdeci*0.5
		else
			pro(3)=pever*0.8
			pro(4)=pever*0.2
			pro(5)=pdeci*0.8
			pro(6)=pdeci*0.2
		end if
	end if
c
	if (class.eq.106) then
		if (lai1.le. 4.0) then
			pro(4)=pever*0.7
			pro(3)=pever*0.3
			pro(6)=pdeci*0.5
			pro(5)=pdeci*0.5
		else
			pro(4)=pever*0.3
			pro(3)=pever*0.7
			pro(6)=pdeci*0.3
			pro(5)=pdeci*0.7
		end if
	end if
c
	if (class.eq.105) then
		if (lai1.le. 4.0) then
			pro(3)=pever*0.5
			pro(4)=pever*0.5
			pro(6)=pdeci*0.5
			pro(5)=pdeci*0.5
		else
			pro(3)=pever*0.5
			pro(4)=pever*0.5
			pro(6)=pdeci*0.2
			pro(5)=pdeci*0.8
		end if
	end if
c
	if (class.eq.104) then
c
		pro(6)=min(0.2+lat/100.0,1.0)
		pro(5)=1.0-pro(6)
	end if
c
	if (class.eq.103) then
		if (lat.le. 40.0) then
			if (lai1.le. 4.0) then
				pro(7)=pever*0.0
				pro(8)=pever*1.0
				pro(5)=pdeci*0.9
				pro(6)=pdeci*0.1
			else
				pro(7)=pever*0.0
				pro(8)=pever*1.0
				pro(5)=pdeci*0.9
				pro(6)=pdeci*0.1
			end if
		else
			if (lai1.le. 4.0) then
				pro(7)=pever*0.5
				pro(8)=pever*0.5
				pro(5)=pdeci*0.5
				pro(6)=pdeci*0.5
			else
				pro(7)=pever*0.8
				pro(8)=pever*0.2
				pro(5)=pdeci*0.8
				pro(6)=pdeci*0.2
			end if
		end if
	end if
c
	if (class.eq.102) then
c
		if (lat .gt. 40.0) then
			pro(7)=min(0.2+lat/100.0,1.0)
			pro(8)=1.-pro(7)
		else
			pro(8)=1.0
			pro(7)=0.0
		end if
	end if
c
	if (class.eq.101) then
		pro(9)=1.0
c
	end if
c
	do i1=2,11
		pro(i1)=pro(i1-1)+pro(i1)
	end do
c------------------------------------------------------------------
c
	itrees=0
	tever=0.0
	yydd=0.0
c
	do while ((tever.le. lai1).and. (itrees.lt.555))
		itrees=itrees+1
c
		call random(tmp2)
		tmp2=min(tmp2,0.9999)
		tmp2=max(tmp2,0.0001)
c
		i2=1
		do while ( tmp2 .ge. pro(i2) )
			i2=i2+1
		end do
		forcn(itrees,1)=float(i2)
c
		tmp1=lai1/45.0+0.02
		forcn(itrees,2)=tmp1
c
		tmp2=forcn(itrees,2)*80.0+1.5*lai1
          !tmp2=trh!!!!!!!!!!!!!!!!!!!!!!!!!
		forcn(itrees,3)=tmp2+1.0
c
		forcn(itrees,5)=forcn(itrees,3)*(-lai1/50.0+0.53)
c
		forcn(itrees,4)=fh_d3(tmp1,tmp2+1.0,forcn(itrees,5))
c
		forcn(itrees,12)=fleafc(forcn(itrees,4),par_tree(15,i2))
c
		yxd=flai(forcn(itrees,12),par_tree(17,i2))
c
		tever=tever+yxd/gaparea

		forcn(itrees,11)=forcn(itrees,12)
c
		forcn(itrees,14)=0.0
	    forcn(itrees,15)=0.0
          above1=fstem3(forcn(itrees,2),forcn(itrees,3),
     ^par_tree(18,i2))
          above2=ftwig3(forcn(itrees,2),
     ^forcn(itrees,3),
     ^forcn(itrees,5),par_tree(18,i2))	
		below=froot3(forcn(itrees,2),forcn(itrees,3),2.0,par_tree(18,i2))
          
          forcn(itrees,17)=above1!!!!stem biomass
          forcn(itrees,18)=above2!!!!twig biomass
          forcn(itrees,19)=below!!!!coarse root biomass
c
		if (i2.ne.5 .and. i2.ne.6 .and. i2.ne.9) then
c	
		forcn(itrees,6)=forcn(itrees,12)
		forcn(itrees,8)=forcn(itrees,6)/2.0
		forcn(itrees,9)=forcn(itrees,17)*0.015
     ^+forcn(itrees,18)*0.04
		forcn(itrees,10)=forcn(itrees,9)/par_tree(9,i2)
c
		forcn(itrees,7)=yxd
c
		else
c
		forcn(itrees,6)=0.0
	    forcn(itrees,7)=0.0
		forcn(itrees,8)=0.0
		forcn(itrees,9)=forcn(itrees,17)*0.015
     ^+forcn(itrees,18)*0.025
		forcn(itrees,10)=forcn(itrees,9)/par_tree(9,i2)
c
		end if
		forcn(itrees,13)=1.0
	    forcn(itrees,16)=0.0
		yydd=yydd+fwood3(tmp1,tmp2,forcn(itrees,5),hr,par_tree(18,i2))
c
	end do
	ntrees=itrees
c-------------------------------------------------------------------	
      !open(111,file='C:\Users\lenovo\Desktop\lai.txt',status='new')
      !write(111,'(16f10.5)') ((forcn(i3,i4),i4=1,16),i3=1,600)
      !close(111)
	return
	end