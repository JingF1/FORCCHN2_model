c	===================================================================
c	soil_c_n: daily soil N,C cycle	
c	
	subroutine soil_cn1(tmean,sw,iy,id,sand,runoff,c_soil,n_soil,ph,
     ^l_soil,par_soil,ts_resp,dry_lit,dry_hum,resp,c_out,n_out,nim,no3)!!!!!!
	real c_soil(15),n_soil(15),l_soil(15),par_res(15),
     ^       ts_resp,dry_lit,dry_hum,par_soil(8),tmean,sw
	real c_n(15),tmp1,tmp2,tmp3,tmp4,tmp5,lig_c,lc
	integer ipool,i
      real soil1,soil2,resp, c_out,  n_out,decomp,decomp2,nl,nim,no3!!!!!!!!!!!
c
	data par_res/0.021,  0.10,   0.027,
     ^		     0.13,  0.01,   0.0020,
     ^			 0.0020,  0.042,   0.042,
     ^	         0.001,  3.5E-05,   1.0,
     ^             1.0,   1.0,     1.0/
c
	resp1=0.0
      tmp1=0.0!!!!!!!!
      tmp2=0.0!!!!!!!
      tmp3=0.0!!!!!!
      tmp4=0.0!!!!!
      tmp5=0.0!!!!!
      nl=0.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do ipool=1,11
	if (n_soil(ipool).gt. 1.0e-15) then
		c_n(ipool)=c_soil(ipool)/n_soil(ipool)
	else
		c_n(ipool)=20.0
	end if
	if  (n_soil(ipool).gt. c_soil(ipool)*0.5 ) then
		n_soil(ipool)=c_soil(ipool)*0.5
	end if 
      end do
c
	soil1=0.0
	do ipool=1,11
		if (l_soil(ipool).lt. 0.000001) then
			l_soil(ipool)=c_soil(ipool)*0.5
		end if
		soil1=soil1+c_soil(ipool)
	end do
c
c 1
	if (c_soil(1).gt.1.0e-10) then
		if ((c_soil(1)-l_soil(1)) .lt. 0.1  ) then
			l_soil(1)=c_soil(1)*0.5
              lig_c=0.5
          else
              lig_c=l_soil(1)/c_soil(1)
		end if
c

		lc=exp(-5.0*lig_c)
c
		decomp=lc*par_res(1)*(c_soil(1)-l_soil(1))*ts_resp*dry_lit
		resp1=resp1+decomp*0.6
		tmp3=decomp*0.4
		c_soil(8)=c_soil(8)+tmp3
		n_soil(8)=n_soil(8)+decomp/c_n(1)!c_n(1)
         !     n_soil(8)=n_soil(8)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^decomp*0.6/c_n(1)-(tmp3/c_n(8)-tmp3/c_n(1))!!!!!!!!!!!!!!
c
		decomp2=decomp/(c_soil(1)-l_soil(1))*l_soil(1)
		resp1=resp1+decomp2*0.3
		tmp4=decomp2*0.7
		c_soil(10)=c_soil(10)+tmp4
		n_soil(10)=n_soil(10)+decomp2/c_n(1)
c
		c_soil(1)=c_soil(1)-decomp-decomp2
		l_soil(1)=l_soil(1)-decomp2
		n_soil(1)=n_soil(1)-(decomp+decomp2)/c_n(1)
	else
		c_soil(1)=1.0e-10
	end if
c
c 3
	if (c_soil(3).gt.1.0e-10) then
		if ((c_soil(3)-l_soil(3)) .lt. 0.1 ) then
			l_soil(3)=c_soil(3)*0.5
		end if
c
		lig_c=l_soil(3)/c_soil(3)
		lc=exp(-5.0*lig_c)
		tmp1=lc*par_res(3)*(c_soil(3)-l_soil(3))*ts_resp*dry_hum
		resp1=resp1+tmp1*0.55
		tmp3=tmp1*0.45
		c_soil(9)=c_soil(9)+tmp3
		n_soil(9)=n_soil(9)+tmp1/c_n(3)!!!!!!!!!!
         !     n_soil(9)=n_soil(9)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^tmp1*0.55/c_n(3)-max(0.0,(tmp3/c_n(9)-tmp3/c_n(3)))!!!!!!!!!!!!!!
c
		tmp2=tmp1/(c_soil(3)-l_soil(3))*l_soil(3)
		resp1=resp1+tmp2*0.3
		tmp4=tmp2*0.7
		c_soil(10)=c_soil(10)+tmp4
		n_soil(10)=n_soil(10)+tmp2/c_n(3)
c
		c_soil(3)=c_soil(3)-tmp1-tmp2
		l_soil(3)=l_soil(3)-tmp2
		n_soil(3)=n_soil(3)-(tmp1+tmp2)/c_n(3)
	else
		c_soil(3)=1.0e-10
	end if
c
c 2
	if (c_soil(2).gt.1.0e-10) then
		tmp1=par_res(2)*c_soil(2)*ts_resp*dry_lit
		resp1=resp1+tmp1*0.6
		tmp3=tmp1*0.4
c
		c_soil(8)=c_soil(8)+tmp3
		n_soil(8)=n_soil(8)+tmp1/c_n(2)!!!!!!!!!!tmp1/c_n(2)
         !     n_soil(8)=n_soil(8)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^tmp1*0.6/c_n(2)-max(0.0,(tmp3/c_n(8)-tmp3/c_n(2)))!!!!!!!!!!!!!!
		c_soil(2)=c_soil(2)-tmp1
		n_soil(2)=n_soil(2)-tmp1/c_n(2)!!!!!!
	else
		c_soil(2)=1.0e-10
	end if
c
c 4
	if (c_soil(4).gt.1.0e-10) then
		tmp1=par_res(4)*c_soil(4)*ts_resp*dry_hum
		resp1=resp1+tmp1*0.55
		tmp3=tmp1*0.45
c
		c_soil(9)=c_soil(9)+tmp3
		n_soil(9)=n_soil(9)+tmp1/c_n(4)!!!!!!!!tmp1/c_n(4)
         !     n_soil(9)=n_soil(9)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^tmp1*0.55/c_n(4)-max(0.0,(tmp3/c_n(8)-tmp3/c_n(4)))!!!!!!!!!!!!!!
		c_soil(4)=c_soil(4)-tmp1
		n_soil(4)=n_soil(4)-tmp1/c_n(4)
	else
		c_soil(4)=1.0e-10
	end if
c
c	5
	if (c_soil(5).gt.1.0e-10) then
		if ((c_soil(5)-l_soil(5)) .lt. 0.1 ) then
			l_soil(5)=c_soil(5)*0.5
		end if
c
		lig_c=l_soil(5)/c_soil(5)
		lc=exp(-5*lig_c)
		tmp1=lc*par_res(5)*(c_soil(5)-l_soil(5))*ts_resp*dry_hum
		resp1=resp1+tmp1*0.6
		tmp3=tmp1*0.4
		c_soil(8)=c_soil(8)+tmp3
		n_soil(8)=n_soil(8)+tmp1/c_n(5)
         !     n_soil(8)=n_soil(8)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^tmp1*0.6/c_n(5)-max(0.0,(tmp3/c_n(8)-tmp3/c_n(5)))!!!!!!!!!!!!!!
c
		tmp2=tmp1/(c_soil(5)-l_soil(5))*l_soil(5)
		resp1=resp1+tmp2*0.3
		tmp4=tmp2*0.7
		c_soil(10)=c_soil(10)+tmp4
		n_soil(10)=n_soil(10)+tmp3/c_n(5)
c
		c_soil(5)=c_soil(5)-tmp1-tmp2
		l_soil(5)=l_soil(5)-tmp2
		n_soil(5)=n_soil(5)-(tmp1+tmp2)/c_n(5)
	else
		c_soil(5)=1.0e-10
	end if
c
c 6
	if (c_soil(6).gt.1.0e-10) then
		if ((c_soil(6)-l_soil(6)) .lt. 0.1  ) then
			l_soil(6)=c_soil(6)*0.5
		end if
c
		lig_c=l_soil(6)/c_soil(6)
		lc=exp(-5*lig_c)
		tmp1=lc*par_res(6)*(c_soil(6)-l_soil(6))*ts_resp*dry_hum
		resp1=resp1+tmp1*0.6
		tmp3=tmp1*0.4
		c_soil(8)=c_soil(8)+tmp3
		n_soil(8)=n_soil(8)+tmp1/c_n(6)!!!!!!!!tmp1/c_n(6)
         !     n_soil(8)=n_soil(8)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^tmp1*0.6/c_n(6)-max(0.0,(tmp3/c_n(8)-tmp3/c_n(6)))!!!!!!!!!!!!!!
c
		tmp2=tmp1/(c_soil(6)-l_soil(6))*l_soil(6)
		resp1=resp1+tmp2*0.3
		tmp4=tmp2*0.7
		c_soil(10)=c_soil(10)+tmp4
		n_soil(10)=n_soil(10)+tmp2/c_n(6)
c
		c_soil(6)=c_soil(6)-tmp1-tmp2
		l_soil(6)=l_soil(6)-tmp2
		n_soil(6)=n_soil(6)-(tmp1+tmp2)/c_n(6)
	else
		c_soil(6)=1.0e-10
	end if
c
c 7.
	if (c_soil(7).gt.1.0e-10) then
		if ((c_soil(7)-l_soil(7)) .lt. 0.1  ) then
			l_soil(7)=c_soil(7)*0.5
		end if
c
		lig_c=l_soil(7)/c_soil(7)
		lc=exp(-5*lig_c)
		tmp1=lc*par_res(7)*(c_soil(7)-l_soil(7))*ts_resp*dry_hum
		resp1=resp1+tmp1*0.55
		tmp3=tmp1*0.45
		c_soil(9)=c_soil(9)+tmp3
		n_soil(9)=n_soil(9)+tmp1/c_n(7)!c_n(7)
         !     n_soil(9)=n_soil(9)+!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !^tmp1*0.55/c_n(7)-max(0.0,(tmp3/c_n(9)-tmp3/c_n(7)))!!!!!!!!!!!!!!
c
		tmp2=tmp1/(c_soil(7)-l_soil(7))*l_soil(7)
		resp1=resp1+tmp2*0.3
		tmp4=tmp2*0.7
		c_soil(10)=c_soil(10)+tmp4
		n_soil(10)=n_soil(10)+tmp2/c_n(7)
c
		c_soil(7)=c_soil(7)-tmp1-tmp2
		l_soil(7)=l_soil(7)-tmp2
		n_soil(7)=n_soil(7)-(tmp1+tmp2)/c_n(7)
	else
		c_soil(7)=1.0e-10
	end if
c
c 8
	if (c_soil(8).gt.1.0e-10) then
		tmp1=par_res(8)*c_soil(8)*ts_resp*max(dry_hum,0.10)
		resp1=resp1+tmp1*0.6
		tmp3=tmp1*0.4
c
		c_soil(10)=c_soil(10)+tmp3
		n_soil(10)=n_soil(10)+tmp1/c_n(8)
		c_soil(8)=c_soil(8)-tmp1
		n_soil(8)=n_soil(8)-tmp1/c_n(8)
	else
	c_soil(8)=1.0e-10
	end if
c
c 9
	n_soil(9)=n_soil(9) !+ 3.0e-6
c--------
	if (c_soil(9).gt.1.0e-10) then
		tmp1=(1.0-0.75*(1.0-par_soil(4)))*
     ^            par_res(9)*c_soil(9)*ts_resp*max(dry_hum,0.10)
		tmp2=tmp1*(0.85-0.68*(1.0-par_soil(4)))
		resp1=resp1+tmp2
		tmp3=tmp1*(0.003+0.032*par_soil(6))
c
		tmp4=0.0
c
		tmp5=tmp1-tmp2-tmp3-tmp4
c
		c_soil(10)=c_soil(10)+tmp5
		n_soil(10)=n_soil(10)+tmp5/c_n(10)
		c_soil(11)=c_soil(11)+tmp3
		n_soil(11)=n_soil(11)+tmp3/c_n(11)
		c_out=tmp4
		n_out=tmp4/c_n(9)
		c_soil(9)=c_soil(9)-tmp1
		n_soil(9)=n_soil(9)-tmp5/c_n(10)-tmp3/c_n(11)-tmp4/c_n(9)
	else
		c_soil(9)=1.0e-10
      end if

c
c 10
	if (c_soil(10).gt.1.0e-10) then
		tmp1=par_res(10)*c_soil(10)*ts_resp*max(dry_hum,0.10)
		tmp2=tmp1*0.55
		resp1=resp1+tmp2
		tmp3=tmp1*(0.003+0.009*par_soil(6))
		tmp5=tmp1-tmp2-tmp3
c
		c_soil(11)=c_soil(11)+tmp3
		n_soil(11)=n_soil(11)+tmp3/c_n(10)
		c_soil(9)=c_soil(9)+tmp5
		n_soil(9)=n_soil(9)+tmp5/c_n(10)
		c_soil(10)=c_soil(10)-tmp1
		n_soil(10)=n_soil(10)-tmp3/c_n(10)-tmp5/c_n(10)
	else
		c_soil(10)=1.0e-10
	end if
c
c 11
	if (c_soil(11).gt.1.0e-10) then
		tmp1=par_res(11)*c_soil(11)*ts_resp*max(dry_hum,0.10)
		tmp2=tmp1*0.55
		resp1=resp1+tmp2
		tmp3=tmp1*0.45
c
		c_soil(9)=c_soil(9)+tmp3
		n_soil(9)=n_soil(9)+tmp1/c_n(11)
		c_soil(11)=c_soil(11)-tmp1
		n_soil(11)=n_soil(11)-tmp1/c_n(11)
	else
		c_soil(11)=1.0e-10
      end if
c
      !!!!!!!!!!!!!calculate mineralized nitrogen!!!!!!!!!!!!!!!
      nim=n_soil(13)+n_soil(14)
      !tmp5=n_soil(9)-c_soil(9)/par_soil(8)
      !if (tmp5.ge.0.0)then
      !    n_soil(9)=n_soil(9)-tmp5
      !    n_soil(13)=n_soil(13)+tmp5
      !else
      !    tmp6=n_soil(9)-(-tmp5)
      !    if(tmp6.ge.0.0)then
      !        n_soil(9)=tmp6
      !        n_soil(9)=n_soil(9)+(-tmp5)
      !    endif
      !end if
      !tmp5=n_soil(8)-c_soil(8)/par_soil(8)
      !if (tmp5.ge.0.0)then
      !    n_soil(8)=n_soil(8)-tmp5
      !    n_soil(13)=n_soil(13)+tmp5
      !else
      !    tmp6=n_soil(8)-(-tmp5)
      !    if(tmp6.ge.0.0)then
      !       n_soil(8)=tmp6
      !       n_soil(9)=n_soil(9)+(-tmp5)
      !    endif
      !end if 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmp5=c_soil(9)/n_soil(9)
	  if (tmp5.lt.par_soil(8)) then
		tmp6=n_soil(9)-c_soil(9)/par_soil(8)!)*1.0e-3!*1.0e-3
		n_soil(13)=n_soil(13)+tmp6
		n_soil(9)=n_soil(9)-tmp6
	  else
		tmp6=c_soil(9)/par_soil(8)-n_soil(9)!)*1.0e-3!*1.0e-3
		tmp7=n_soil(13)-tmp6
		if (tmp7.gt.0.0) then 
		   n_soil(13)=tmp7
			n_soil(9)=n_soil(9)+tmp6!n_soil(9)
		else
			n_soil(9)=n_soil(9)+n_soil(9)!n_soil(9)
			n_soil(13)=0.0
		end if
      end if
c	
	tmp5=c_soil(8)/n_soil(8)
	if (tmp5.lt.par_soil(8)) then
		tmp6=n_soil(8)-c_soil(8)/par_soil(8)!)*1.0e-3!*1.0e-3
		n_soil(13)=n_soil(13)+tmp6
		n_soil(8)=n_soil(8)-tmp6!n_soil(8)
	else
		tmp6=c_soil(8)/par_soil(8)-n_soil(8)!)*1.0e-3!*1.0e-3
		tmp7=n_soil(13)-tmp6
		if (tmp7.gt.0.0) then
			n_soil(13)=tmp7
			n_soil(8)=n_soil(8)+tmp6
		else
			n_soil(8)=n_soil(8)+n_soil(13)!n_soil(8)
			n_soil(13)=0.0
		end if
      end if
        !ph=7
        tmp2=5.8*(10**(ph-10))*n_soil(13)!ph-10
        n_soil(13)=n_soil(13)-tmp2
      no3=n_soil(14)
      tmp1=2**((tmean-20.)/10.)
      if (sw.lt.par_soil(1))then
          tmp2=1.17*sw/par_soil(1)+0.165
      else
          tmp2=1.0-(0.1*sw/par_soil(1))
      endif
      tmp3=(min(1.0,0.003*tmp1*tmp2))*n_soil(13)!0.003
      if(n_soil(13).ge.tmp3)then
          n_soil(13)=n_soil(13)-tmp3
          tmp3=(1-0.001)*tmp3
          n_soil(14)=n_soil(14)+tmp3
            tmp1=3**((tmean-25.)/10.)
          if (sw.le.par_soil(1))then
              tmp2=0.0
          else
              tmp2=sw/par_soil(1)
          endif
          tmp4=min(1.0,0.0015*tmp1*tmp2)*n_soil(14)
          n_soil(14)=max(0.0,n_soil(14)-tmp4)
      endif
      tmp4=runoff*0.1*n_soil(14)
     ^*(0.2+0.7*sand)*0.2*0.001!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tmp4=max(tmp4,0.0)
      n_soil(14)=max(0.,(n_soil(14)-tmp4))!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      no3=n_soil(14)-no3!!!net nitrification nitrogen
      nim=n_soil(13)+n_soil(14)-nim!!!net mineralized nitrogen
      tmp1=0.9*((sw/par_soil(1))**3.0)+0.1
      tmp5=0.5*tmp1*(n_soil(13)+n_soil(14))/(2.7+tmp1*(n_soil(13)
     ^+n_soil(14)))*exp(0.0693*tmean)
      if(n_soil(13)+n_soil(14).le.tmp5)then
          n_soil(13)=0.0
          n_soil(14)=0.0
          n_soil(15)=n_soil(13)+n_soil(14)
      elseif (n_soil(13)+n_soil(14).ne.0.0)then
          n_soil(13)=max(0.0,n_soil(13)-tmp5*n_soil(13)
     ^/(n_soil(13)+n_soil(14)))
          n_soil(14)=max(0.0,n_soil(14)-tmp5*n_soil(14)
     ^/(n_soil(13)+n_soil(14)))
          n_soil(15)=tmp5
      else
          n_soil(13)=0.0
          n_soil(14)=0.0         
      endif
          !n_soil(15)=max(0.0,n_soil(15)-tmp5)
c
	resp=max(0.0001,resp1)
c
      return
      end