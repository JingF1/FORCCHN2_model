	subroutine grow2(ntrees,par_tree,forcn,l_soil,c_soil,n_soil,cout,key,gaparea,coutn)
	real par_tree(28,9),forcn(600,19)
    real nsc2,nsc1,tmp11,woodinc,dinc  
	real l_soil(15),c_soil(15),n_soil(15)
	real cout(8),key(900,3),ttmmpp,yxd1,gaparea
	integer ntrees,i0
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
 
	yxdi0=-30000.0
	do it=1,ntrees
		i1=int(forcn(it,1))
	    if (forcn(it,9).gt.yxdi0) then
			i0=i1
		end if
	end do

	do it=1,ntrees
		if (forcn(it,2).gt. 0.01) then
			i1=int(forcn(it,1))

			tmp1=forcn(it,9)
			
			if (tmp1.gt.0.00001) then

				if (i1.eq.5 .or. i1.eq.6 .or. i1.eq. 9) then
					tmp2=min(forcn(it,12)*0.10,tmp1)
				else
					tmp2=min(forcn(it,12)*0.05,tmp1)
				end	if

				forcn(it,9)=tmp2
				tmp33=tmp1-tmp2
				tmp01=forcn(it,12)*par_tree(8,i1)
				tmp02=tmp33*par_tree(8,i1)
				tmp0=max(tmp01,tmp02)
				tmpo=min(tmp01,tmp02)
				if (tmp33.gt.tmp0) then					

					tmp3=tmp0
				else
					tmp3=tmp33*0.9 
				end if

				ttmmpp=tmp3

				tmp=tmp33-tmp3	
			else
				if (i1.eq.5 .or. i1.eq.6 .or. i1.eq. 9) then
					if (i0.eq.5 .or. i0.eq.6 .or. i0.eq. 9) then
						forcn(it,1)=float(i0)
						tmp=forcn(it,9)-forcn(it,12)*0.10
						forcn(it,9)=forcn(it,12)*0.10
					else
						forcn(it,1)=float(i0)
						tmp=forcn(it,9)-forcn(it,12)*1.25+forcn(it,6)+forcn(it,8)
						forcn(it,6)=forcn(it,12)*0.8
						forcn(it,8)=forcn(it,12)*0.4
						forcn(it,9)=forcn(it,12)*0.05
					end if
				else
					if (i0.eq.5 .or. i0.eq.6 .or. i0.eq. 9) then
						forcn(it,1)=float(i0)
						tmp=forcn(it,9)+forcn(it,6)+forcn(it,8)-forcn(it,12)*0.13
						forcn(it,9)=forcn(it,12)*0.10
						forcn(it,6)=forcn(it,12)*0.02
						forcn(it,6)=forcn(it,12)*0.01
					else
						forcn(it,1)=float(i0)
						tmp=forcn(it,9)-forcn(it,12)*0.05
						forcn(it,9)=forcn(it,12)*0.05
					end if
				end if
			    forcn(it,11)=tmp
                  ttmmpp=0.0!
			end if

			if (tmp .gt. 0.0) then 
				fx=0.0
			else
				fx=-1.0
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
       
            nsc2=forcn(it,9)!+prod-resp_above-resp_below-grwresp-nscstart(it)!forcn(it,9)+prod-resp_above-resp_below-grwresp-dleaf-cs-cr!!!!!!光和作用有问题，太小了
            nsc1=0.0
            if (nsc2<=0.0)then
                tmp11=0.0
            else
                call daywood(it,par_tree,forcn,l_soil,c_soil,n_soil,key,gaparea,nsc2,nsc1,tmp11,woodinc,dinc)
                !dinc1(id,it)=forcn(it,2)!dinc,日输出
                !dinc3(iy,it)=dinc
                if (i1.eq.3 .or. i1.eq.5 .or. i1.eq.7.or. i1.eq. 9) then
                    tmp11=tmp11
                else
                    !tmp2=min(forcn(it,12)*0.05,tmp1)
                    tmp11=tmp11/1.5
                end	if
            end if
            forcn(it,9)=forcn(it,9)-tmp11-tmp11*0.3
            forcn(it,17)=fstem3(forcn(it,2),forcn(it,3),par_tree(18,i1))
            forcn(it,18)=ftwig3(forcn(it,2),forcn(it,3),forcn(it,5),par_tree(18,i1))	
		    forcn(it,19)=froot3(forcn(it,2),forcn(it,3),2.0,par_tree(18,i1))

			forcn(it,12)=fleafc(forcn(it,4),par_tree(15,i1))
			if (tmp .gt. 0.00001) then

				forcn(it,11)=fwood3(forcn(it,2),forcn(it,3),forcn(it,5),hr,ck)-yxd

				tmp5=forcn(it,11)/par_tree(10,i1)

				forcn(it,9)=tmp1-ttmmpp-forcn(it,11)		

				tmp4=forcn(it,10)-forcn(it,9)/par_tree(9,i1)-tmp5
            else
                  tmp4=0.0
			end if
	
	 		forcn(it,10)=forcn(it,9)/par_tree(9,i1)

			!d=d+aa
			!h=h+cp*aa 

			yxd1=ttmmpp/gaparea
			cout(7)=cout(7)+yxd1
			c_soil(6)=c_soil(6)+yxd1
			n_soil(6)=n_soil(6)+tmp4/gaparea
			l_soil(6)=l_soil(6)+yxd1*par_tree(27,i1)*0.01

			if (key(it,1).lt.par_tree(1,i1) .and. tmp.gt.0.0) then
	
				forcn(it,5)=min(b+0.5,h)
				bh=forcn(it,5)
	
				yxd1=ftwig3(d,h,b,ck)-ftwig3(d,h,bh,ck)
				yxd1=yxd1/gaparea
   				cout(7)=cout(7)+yxd1
				c_soil(6)=c_soil(6)+yxd1
				n_soil(6)=n_soil(6)+yxd1/par_tree(10,i1)
				l_soil(6)=l_soil(6)+yxd1*par_tree(27,i1)*0.01
				
				forcn(it,4)=fh_d3(d,h,b+1.0)
			
				forcn(it,12)=fleafc(forcn(it,4),par_tree(15,i1))
			end if

		else
			i1=int(forcn(it,1))
			d=forcn(it,2)
			h=forcn(it,3)
			hb=forcn(it,5)
			ck=par_tree(18,i1)
			yxdw=fwood3(d,h,hb,hr,ck)+forcn(it,9)
			yxdl=forcn(it,6)+forcn(it,8)
			yxdwl=(yxdw+yxdl)/gaparea
			cout(7)=cout(7)+yxdwl
			c_soil(6)=c_soil(6)+yxdwl
			n_soil(6)=n_soil(6)+(yxdl/par_tree(9,i1)+yxdw/par_tree(10,i1)+forcn(it,10))/gaparea
			l_soil(6)=l_soil(6)+yxdwl*par_tree(27,i1)*0.01
			forcn(it,2)=0.0
			forcn(it,3)=0.0
			forcn(it,4)=0.0
			forcn(it,5)=0.0
			forcn(it,6)=0.0
			forcn(it,8)=0.0
			forcn(it,9)=0.0
			forcn(it,10)=0.0
		end if
	forcn(it,12)=fleafc(forcn(it,4),par_tree(15,i1))
	end do
	cout(1)=0.0
	cout(2)=0.0

	do it=1,ntrees
		if (forcn(it,2) .gt. 0.01) then
		i1=int(forcn(it,1))
		d=forcn(it,2)
		h=forcn(it,3)
		b=forcn(it,5)
		ck=par_tree(18,	i1)
		cout(1)=cout(1)+fstem3(d,h,ck)+ftwig3(d,h,b,ck)+forcn(it,6)!+forcn(it,9)
		cout(2)=cout(2)+froot3(d,h,hr,ck)+forcn(it,8)
		end if
	end do
	cout(8)=0.0
      coutn=0.0!!!!!!!!!!!!!!!!!!!!!!!!!!
	do ip=1,11
		cout(8)=cout(8)+c_soil(ip)
          coutn=coutn+n_soil(ip)!!!!!!!!!!!!!!!!!!!!!!!!
	end do

	return
	end	