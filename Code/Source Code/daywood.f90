subroutine daywood(ii,par_tree,forcn,l_soil,c_soil,n_soil,key,gaparea,nsc2,nsc1,tmp,woodinc,aa)
	real par_tree(28,9),forcn(600,19),woodinc,dinc(900,1)
	real l_soil(15),c_soil(15),n_soil(15),nsc2,nsc1
	real cout(8),key(900,3),ttmmpp,yxd1,gaparea
	integer ntrees,i0,ii
	integer it,ih,i1,i2,i3
	real tmp,tmp1,tmp2,tmp3,fallca,fallna,fallcb,fallnb
	real d,h,b,cp,a1(1000),b1(1000),ff,aa,ck,yxd,hr,bh,fx,aacp,coutn
    
    hr=2.0
    
	fallca=0.0
	fallna=0.0
	fallcb=0.0
	fallnb=0.0
    
	do it=1,8
		cout(it)=0.0
    end do

    
    it=ii
    yxdi0=-30000.0
	i1=int(forcn(it,1))
	if (forcn(it,9).gt.yxdi0) then
		i0=i1
	end if
    if (i1.eq.3 .or. i1.eq.5 .or. i1.eq.7.or. i1.eq. 9) then
		!tmp2=min(forcn(it,12)*0.10,tmp1)
        tmp=(nsc2-nsc1)*0.75
	else
		!tmp2=min(forcn(it,12)*0.05,tmp1)
        tmp=(nsc2-nsc1)*1.5
    end	if
    if (tmp.le.0.0)then
        aa=0.0
        woodinc=0.0
    endif
    
    if (tmp.gt.0.0)then
        cp=key(it,2)*100.0
		if (cp.gt.100.0) then
			cp=100.0
		end if
		if (cp.lt. 10.0) then
			cp=10.0
		end if
        d=forcn(it,2)
		h=forcn(it,3)
		b=forcn(it,5)

		cp=key(it,2)*100.0
		if (cp.gt.100.0) then
			cp=100.0
		end if
		if (cp.lt. 10.0) then
			cp=10.0
		end if
		ck=par_tree(18,i1)
		yxd=fwood3(d,h,b,hr,ck)
		i3=1
		ff=1.0
		a1(1)=-d+0.005
		b1(1)=0.2
        fx=0.0
		do while ((abs(ff).gt. 1.0E-7) .and. (i3 .le. 555))
			aa=(a1(i3)+b1(i3))*0.5
			aacp=aa*cp
			ff=fwood3(d+aa,h+aacp,b+fx*aacp,hr,ck)-yxd-tmp

			if (ff.gt. 0.0) then
				a1(i3+1)=a1(i3)
				b1(i3+1)=aa
				i3=i3+1
			else
				a1(i3+1)=aa
				b1(i3+1)=b1(i3)
				i3=i3+1
			end if
        end do
            woodinc=fwood3(d+aa,h+aacp,b+fx*aacp,hr,ck)-yxd
			forcn(it,2)=d+aa

			forcn(it,3)=h+aacp

			forcn(it,5)=b+fx*aacp
			forcn(it,4)=fh_d3(forcn(it,2),forcn(it,3),b+fx*aacp)

			!forcn(it,12)=wl(forcn(it,2),forcn(it,5))
			!if (tmp .gt. 0.00001) then
   !
			!	forcn(it,11)=fwood3(forcn(it,2),forcn(it,3),forcn(it,5),hr,ck)-yxd
   !
			!	tmp5=forcn(it,11)/par_tree(10,i1)
   !
			!	forcn(it,9)=tmp1-ttmmpp-forcn(it,11)		
   !
			!	tmp4=forcn(it,10)-forcn(it,9)/par_tree(9,i1)-tmp5
   !         else
   !                 tmp4=0.0
			!end if
	
	 		!forcn(it,10)=forcn(it,9)/par_tree(9,i1)

			d=d+aa
			h=h+cp*aa 

			!yxd1=ttmmpp/gaparea
			!cout(7)=cout(7)+yxd1
			!c_soil(6)=c_soil(6)+yxd1
			!n_soil(6)=n_soil(6)+tmp4/gaparea
			!l_soil(6)=l_soil(6)+yxd1*par_tree(27,i1)*0.01

			!if (key(it,1).lt.par_tree(1,i1) .and. tmp.gt.0.0) then
	  !
			!	forcn(it,5)=min(b+0.5,h)
			!	bh=forcn(it,5)
	  !
			!	yxd1=ftwig3(d,h,b,ck)-ftwig3(d,h,bh,ck)
			!	yxd1=yxd1/gaparea
   !				cout(7)=cout(7)+yxd1
			!	c_soil(6)=c_soil(6)+yxd1
			!	n_soil(6)=n_soil(6)+yxd1/par_tree(10,i1)
			!	l_soil(6)=l_soil(6)+yxd1*par_tree(27,i1)*0.01
			!	
			!	forcn(it,4)=fh_d3(d,h,b+1.0)
			!
			!	forcn(it,12)=wl(forcn(it,2),forcn(it,5))
   !         end if
    endif
  
 !   cout(1)=0.0
	!cout(2)=0.0
 !
 !
	!if (forcn(it,2) .gt. 0.01) then
	!    i1=int(forcn(it,1))
	!    d=forcn(it,2)
	!    h=forcn(it,3)
	!    b=forcn(it,5)
	!    ck=par_tree(18,	i1)
	!    cout(1)=cout(1)+fstem3(d,h,ck)+ftwig3(d,h,b,ck)+forcn(it,6)!+forcn(it,9)           
	!    cout(2)=cout(2)+froot3(d,h,hr,ck)+forcn(it,8)
 !   end if
        
	cout(8)=0.0
    coutn=0.0!!!!!!!!!!!!!!!!!!!!!!!!!!
	do ip=1,11
		cout(8)=cout(8)+c_soil(ip)
        coutn=coutn+n_soil(ip)!!!!!!!!!!!!!!!!!!!!!!!!
    end do
    
    return
    end