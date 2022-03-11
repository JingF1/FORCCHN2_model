c=====================================================================
c---------------------------------------------------------------------
	real function phototree(par,dl,tlai,c,tem,dry,ns,am,sl,kl)
	real par,dl,tlai,tem,dry,ns,am,sl,kl,photo,photo1
	photo1=am*dl/kl*LOG((1.0+(1.0+kl*sl*par/am)**0.5)/
     ^			(1.0+(1.0+kl*sl*par*exp(-kl*tlai)/am)**0.5))
      photo=2.0*photo1*c*tem*dry
      !phototree=photo*(ns/(ns+2.4e-6))!!!!!!!!!!!!!!!!!
	!phototree=min(photo,ns*150.0)!!!!!!!
      phototree=photo!!!!!!!!!!!!!!!!!!!
c
	return
	end
c
c---------------------------------------------------------------------
	subroutine phenology(forcn,ntrees,tday,par_tree,parmax,
     ^					gaparea,lai,parcano,key,phenfac)
	real forcn(600,16),dry_photo,tday,par_tree(28,9),parmax,gaparea
	real lai,laicano(120),parcano(120),key(900,3),yxd(900),
     ^     phenfac,yy
	integer ntrees
	integer it,i1,i2,i3,i5
	real tmp1,tmp2,tmp5(900),con(120),bro(120),con1(120),bro1(120)
c
	do i1=1,120
		laicano(i1)=0.0
		parcano(i1)=parmax
		con(i1)=0.0
		bro(i1)=0.0
	end do
c----
	do it=1,ntrees
		forcn(it,16)=phenfac
	end do

c-----
	tmp1=0.0
	do it=1,ntrees
		if (forcn(it,2).gt. 0.01) then
			i1=int(forcn(it,1))
			if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then
c			
				if ((phenfac.eq.1.0 ).and. 
     ^                (forcn(it,6).le. forcn(it,12)*0.01)) then
					yy=max(forcn(it,9)*0.9,forcn(it,12)*0.025)
					forcn(it,6)=yy/1.5+forcn(it,6)
					forcn(it,8)=yy/3.0+forcn(it,8)
					forcn(it,9)=forcn(it,9)-yy
					forcn(it,10)=forcn(it,10)-
     ^							yy/par_tree(9,i1)
					forcn(it,16)=0.0
				end if
			end if
c--------	
			forcn(it,7)=max(flai(forcn(it,6),par_tree(17,i1)),1.0e-15)
			yxd(it)=(forcn(it,7)+flai(forcn(it,12),
     ^			par_tree(17,i1)))*0.5
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
c------
	if (tmp1.gt. 1.0e-15) then
		do it=1,ntrees
			forcn(it,13)=tmp5(it)/tmp1*gaparea
		end do
	else 
		do it=1,ntrees
		forcn(it,13)=1.0/float(ntrees)*gaparea
		end do
	end if
c-------
	lai=0.0
c
	do i5=1,119
		lai=lai+con(i5)+bro(i5)
		con(120-i5)=con(120-i5+1)+con(120-i5)
		bro(120-i5)=bro(120-i5+1)+bro(120-i5)
c
	end do
c------
	do i5=1,119
	parcano(i5)=max(parmax*exp(-0.35*con(i5+1)-0.40*bro(i5+1)),10.)
c
	end do
	do it=1,ntrees
		i2=max(min(int(forcn(it,3)),119),1)
		i3=max(int(forcn(it,5)),1)
c
		key(it,1)=key(it,1)+parcano(i3)*phenfac
				key(it,2)=key(it,2)+phenfac*
     ^      (parcano(i2+1)-parcano(i2))/(parcano(i2+1))
	end do
	return
	end
c========================================================
c	tree's geometry
c--------------------------------------------------------	
	real function stemc(diameter,height,kstem)
	real diameter,height,kstem
	real tmp1,tmp2
c
	tmp1=0.235619
	stemc=kstem*diameter**2*height*tmp1
	return
	end
c
	real function twigc(baseh,twigd,twigh,kstem)
	real baseh,twigd,twigh,kstem
c
	twigc=stemc(twigd,baseh-twigh,kstem)*2.0/3.0
	return
	end
c	
	real function h_d(d0,h0,h)
	real h0,d0,h
c
	h_d=d0*(1.0-h/h0)**0.333333333
	return
	end

	real function rootc(d0,h0,h,kstem)
	real h,d0,h0,kstem
	real tmp1,d
c	
	d=h_d(d0,h0,-h)
	tmp1=0.392699
	rootc=stemc(d,h0+h,kstem)-stemc(h0,h0,kstem)
	return
	end
c
	real function fleafc(d,kleaf)
	real kleaf,d
c
	fleafc=kleaf*d**2
	return
	end
c	
	real function flai(leafc,c_leaf)
	real leafc,c_leaf
c
	flai=c_leaf*leafc
	return
	end
c=========================================
c
	real function fh_d1(d,h,h1)
	real d,h,h1
	fh_d1=d*(1.-(h1/h))
	return
	end
c
	real function fstem1(d,h,astem)
	real d,h,astem
	fstem1=1.0471975*astem*d*d*h
	return
	end
c
	real function ftwig1(d,h,hb,astem)
	real d,h,hb,astem
	ftwig1=2.094395067*astem*d*d*h*(1.0-hb/h)**3
	return
	end
c
	real function froot1(d,h,hr,astem)
	real d,h,hr,astem
	froot1=1.0471975*astem*d*d*h*((1.0+hr/h)**3-1.0)
	return
	end
c
	real function fwood1(d,h,hb,hr,astem)
	real d,h,hb,hr,astem
	fwood1=1.0471975*astem*d*d*h*((1.+hr/h)**3+2.*(1.-hb/h)**3)
	return
	end
c-----------------------------------------------------------
c
	real function fh_d3(d,h,h1)
	real d,h,h1
      !d=0.5*d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	fh_d3=max(d*(1.-(h1/h)**3),0.0)
      !fh_d3=2*fh_d3!!!!!!!!!!!!!!!!!
	return
	end
c
	real function fstem3(d,h,astem)
	real d,h,astem
      !d=0.5*d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	fstem3=max(2.019595243*astem*d*d*h,0.0)
	return
	end
c
	real function ftwig3(d,h,hb,astem)
	real d,h,hb,astem,b
      !d=0.5*d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	b=hb/h
	ftwig3=max(1.121997357*astem*d*d*h*((1.0-b**3)**2)*(1.0-b),
     ^0.0)
	return
	end
c
	real function froot3(d,h,hr,astem)
	real d,h,hr,astem,b
      !d=0.5*d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	b=hr/h
	froot3=
     ^ 1.5*max(2.019595243*astem*d*d*h*((1.0+b**3)**2*((1.0+b))-1.0),
     ^0.0)
	return
	end
c
	real function fwood3(d,h,hb,hr,astem)
	real d,h,hb,hr,astem,b,c
c
      !d=0.5*d!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	fwood3=fstem3(d,h,astem)+ftwig3(d,h,hb,astem)+froot3(d,h,hr,astem)
c
	return
	end
c===============================================================
	real function woodc(based,baseh,twigd,twigh,rooth,kstem)
	real based,baseh,twigh,rooth,twigd,rootd,kstem
c
	rootd=based*(1+rooth/baseh)**(0.33333)
	woodc=max(stemc(rootd,baseh+rooth,kstem)+
     ^	  twigc(baseh,twigd,twigh,kstem),0.0)
	return
	end
c
	subroutine resptree(i1,leaf,above1,below1,tp_resp,tl_resp,
     ^           par_tree,resp_above,resp_below,ns)
	integer i1
	real leaf,above1,below1,tp_resp,tl_resp,par_tree(28,9),ns
	real resp_above,resp_below
	    resp_above=tp_resp*leaf*par_tree(5,i1)+!!!!!!!!!!!*1
     ^	tp_resp*above1*par_tree(6,i1)!!!!!!!!!!*2.5
	    resp_below=tp_resp*(leaf*par_tree(5,i1)*0.5
     ^    +below1*par_tree(7,i1))!!!!!!!!!!!!*2.5
 !     resp_above=(tp_resp*leaf*par_tree(5,i1)+
 !    ^	tp_resp*above1*par_tree(6,i1))*(ns/(0.01+ns))
	!resp_below=(tp_resp*(leaf*par_tree(5,i1)*0.5+below1*par_tree(7,i1)))
 !    ^*(ns/(0.01+ns))
	return
	end
c
	subroutine fall(i1,phen,leaf,firoot,par_tree,ll,frl,lln,frln)
	integer i1
	real leaf,par_tree(28,9),phen,firoot
	real ll,frl,lln,frln,kk
      kk=0.2!!0.2
	if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then
		if (phen.lt.0.8) then
			ll=leaf*kk
			frl=firoot*kk
		else
			ll=leaf*par_tree(26,i1)
			frl=firoot*par_tree(26,i1)
		end if
	end if
	if (i1.eq.1 .or. i1.eq.2 .or. i1.eq.4 .or. i1.eq. 3) then
c
			ll=leaf*par_tree(26,i1)
			frl=firoot*par_tree(26,i1)
c
	end if
	if (i1.eq.7 .or. i1.eq.8 ) then
c
			ll=leaf*par_tree(26,i1)
			frl=firoot*par_tree(26,i1)
c
	end if
	ll=max(ll,0.0)
	frl=max(frl,0.0)
c
      frln=frl/par_tree(9,i1)
      lln=ll/par_tree(9,i1)
	return
	end
c
	subroutine fallw(i1,above1,below1,par_tree,abc,bec,abn,ben)
	integer i1
	real above1,below1,par_tree(28,9)
	real abc,abn,bec,ben,xx
c
	abc=above1*par_tree(28,i1)!!!!!!!!!!!!!!!!!!!!!
	bec=below1*par_tree(28,i1)!!!!!!!!!!!!!!!!!!!
	abn=abc/par_tree(10,i1)
	ben=bec/par_tree(11,i1)
	return
	end
c=====================================================================
c-------------------------------------------------------
	subroutine param(par_tree)
	integer i,j
	real par_tree(28,9),tmp(9,28)	
	data tmp/
     ^ 5.5,   11.0,   5.5,   11.0, 5.5,  11.0, 5.5,  11.0,   11.0,		 
     ^ 5.5E-4,5.5E-4,5.5E-4,5.5E-4, 5.0E-4,5.0E-4,5.0E-4,5.0E-4,5.0E-4,
     ^ 1.3E-5,1.3E-5,1.3E-5,1.3E-5, 1.3E-5,1.3E-5,1.3E-5,1.3E-5,1.3E-5,
     ^ 4.5E-1,4.5E-1,4.5E-1,4.5E-1, 4.0E-1,4.0E-1,4.0E-1,4.0E-1,3.5E-1,
     ^ 2.0E-3,2.0E-3,2.0E-3,2.0E-3, 6.0E-3,3.0E-3,3.5E-3,3.5E-3,1.2E-2,
     ^ 1.0E-3,1.0E-3,1.0E-3,1.0E-3, 2.0E-3,2.0E-3,2.0E-3,2.0E-3,2.0E-3,
     ^ 1.5E-3,1.5E-3,1.5E-3,1.5E-3, 2.5E-3,2.5E-3,2.5E-3,2.5E-3,2.5E-3,
     ^ 0.50,	0.50,	0.40,  0.40,  0.40,	 0.40,	0.50,  0.50,  0.50,
     ^ 40.,	40.,	45.,   45.,	  40.,	 40.,   60.,   60.,   50.,
     ^ 200.,	200.,	200.,  200.,  200.,	 200.,	200.,  200.,  200.,
     ^ 40.,   40.,	45.,   45.,	  40.,	 40.,	60.,   60.,   50.,
     ^ 40.,   60.,   	50.,   40.,   40.,   40.,   60.,   60.,   50.,
     ^ 2.,	3.,    	2.,	   1.5,	  2.,    1.5,   2.,	   2.,    2.,
     ^ 200.,  100.,  	400.,  200.,  400.,  200.,	1000., 300.,  500.,
     ^ 600.,  600.,	600.,  600.,  200,	 700.,	700.,  700.,  300.,
     ^ 20.,	20.,	20.,   20.,	  30.,   30.,	15.,   15.,   28.,	
     ^ 15.,   15.,   	15.,   15.,   45.,   20.,   18.,   18.,   40.,!!!!!!µÚÎåÁÐ45. 
     ^ 350.,	350.,  	350.,  350.,  350.,	 350.,  350.,  350.,  350.,!200.,	200.,  	200.,  200.,  200.,	 200.,  200.,  200.,  200.,!!!
     ^ 1.,    3.,	    1.,	   3.,	  1.,    9.,	3.,	   9.,    9.,
     ^ 9., 	3.,	    9.,	   3.,	  9.,	 3.,	9.,	   3.,    9.,
     ^ 3.,	3.,	    3.,	   3.,	  3.,	 3.,	0.,	   0.,    0.,  
     ^ 5.0,   5.0,    3.0,   1.0,   -1.0,  -5.5,  -5.5,  -2.5, -5.5,
     ^ 27.0,  29.0,  	27.0,  25.0,  23.0,  20.0,  18.0,  23.0,  16.0,
     ^ 50.0,  50.0, 	50.0,  50.0,  45.0,  45.0,  40.0,  40.0,  35.0,
     ^ 1.0,   0.8,   	0.9,   0.8,   0.8,   0.6,   0.9,   0.7,   0.5,
     ^ 2.0E-3,2.0E-3,2.0E-3,2.0E-3, 1.1E-4,1.1E-4,2.0E-3,2.0E-3,1.1E-4,
     ^ 40.0, 	40.0,  40.0,  40.0,   30.0,  50.0, 	80.0,  80.0,  50.0,
     ^ 5.0e-5,5.0e-5,5.0e-5,5.0e-5, 4.0e-5,4.0e-5,8.0e-5,8.0e-5,8.0e-5/
c     3^ 1.3E-5,1.3E-5,1.3E-5,1.3E-5, 1.3E-5,1.3E-5,1.3E-5,1.3E-5,1.3E-5,
c     4^ 4.5E-1,4.5E-1,4.5E-1,4.5E-1, 4.0E-1,4.0E-1,4.0E-1,4.0E-1,3.5E-1,
	do i=1,28
		do j=1,9
			par_tree(i,j)=tmp(j,i)
		end do
	end do
	return
	end
