
    subroutine ly(ny0,ny,iyear,gpp,ndays1,li0,liy)!
    real gpp(ndays1,1),li0,liy
    real ff1,ff,ff2(ny,1),ff3
    integer i4,iy,ny,ny0,day,ndays(ny,1),year,idays1,idays2,day1,day2,day11(ny,1),day22(ny,1)

    do i4=1,ny
       iy=ny0+i4-1
       if (mod(iy,4).eq.0) then
           day=366
           day1=183
           day2=245
       else
           day=365
           day1=182
           day2=244
       endif
       ndays(i4,1)=day
       day11(i4,1)=day1
       day22(i4,1)=day2
    enddo
    year=ny0+iyear-1
    if (iyear.eq.1)then
        idays1=day11(iyear,1)
        idays2=day22(iyear,1)
    else
        idays1=sum(ndays(1:(iyear-1),1))+day11(iyear,1)-1
        idays2=sum(ndays(1:(iyear-1),1))+day22(iyear,1)
    endif
    ff1=sum(gpp(idays1:idays2,1))
    if (ny.eq.1)then
        ff3=ff1
    else
        do i4=1,ny
           idays1=sum(ndays(1:(i4-1),1))+day11(i4,1)
           idays2=sum(ndays(1:(i4-1),1))+day22(i4,1)
           ff2(i4,1)=sum(gpp(idays1:idays2,1))
        enddo
        ff3=sum(ff2)
    endif
    liy=li0*ff1/((1.0/real(ny))*ff3)
    return
    end
    
    
    subroutine ldyk(ny0,ny,gpp,ndays1,ld0,ldy)
    real gpp(ndays1,1),ld0,ldy(ndays1,1)
    real ff1,ff2,ldy1(366,1),ldy2(365,1)
    integer ny0,ny,iyear,ndays1,id,ndays(ny,1),day,i,j,k,ndays21,ndays22
    do i4=1,ny
       iy=ny0+i4-1
       if (mod(iy,4).eq.0) then
           day=366
       else
           day=365
       endif
       ndays(i4,1)=day
    enddo
    ndays21=ndays(1,1)
    ndays22=ndays(2,1)
    if (ny.eq.1)then
        do i=1,ndays1
            ff1=sum(gpp(1:i,1))
            ff2=ff1
            ff2=max(0.00001,ff2)
            ldy(i,1)=ld0*ff1/ff2
        enddo
    else
        k=0
        do i4=1,ny
            iy=ny0+i4-1
            if (mod(iy,4).eq.0) then
                day=366
                do j=1,day
                   ff1=sum(gpp(k+1:k+j,1))
                   ff2=sum(gpp((ndays21+1):(ndays21+j),1))
                   ff2=max(0.00001,ff2)
                   ldy1(j,1)=ld0*ff1/ff2
                enddo
                ldy(k+1:k+day,1)=ldy1(:,1)
            else
                day=365
                do j=1,day
                   ff1=sum(gpp(k+1:k+j,1))
                   ff2=sum(gpp((ndays21+1):(ndays21+j),1))
                   ff2=max(0.00001,ff2)
                   ldy2(j,1)=ld0*ff1/ff2
                enddo
                ldy(k+1:k+day,1)=ldy2(:,1)
            endif
            k=k+day
        enddo
    endif
    return
    end
    
    
    real function kt(wt,wcrit)
    real a,wt,wcrit
    a=0.1
    kt=min(1.0,(1-exp(-a*wt))/(1-exp(-a*wcrit)))
    return
    end
    
    
    subroutine phenoleaf(tt,sn,snc,fn,gn)
    real sn,snc,fn,gn,a,b,tt
    a=0.185
    b=18.4
    if (tt<0.0)then
        gn=0.0
    else
        gn=1.0/(1.0+exp(-a*(tt-b)))
    endif
    
    if (sn<0.0)then
        fn=0.0
    elseif(sn>=0.0.and.sn<=snc)then
        fn=4.0*(((snc*sn)**0.5)-sn)/snc
    else
        fn=0.0
    endif 
    end
    
    !dimensinal growth: daily growth
    subroutine gf(ny,iy,id,t,ts,mt,par_soil,sfc,gn,gs,gd,gr,sn,ss,sd,sr,fn,fs,fd,fr,snc,ssc,sdc,src,fj)
    real a,b,ff,pi,mt,mt1,par_soil(8)
    real gn,gs,gd,gr
    real snd,ssd,sdd,srd,sn,ss,sd,sr,snc,ssc,sdc,src
    real fn,fs,fd,fr
    integer id,iy,fj(ny,3)

    pi=3.1415
    a=0.185
    b=18.4
    mt1=mt/par_soil(1)!par_soil(1)
    if (t<0.0)then
       gn=0.0
    else
        gn=(1.0+exp(-a*(t-b)))**(-1.0)
    endif
    gs=gn
    gd=gn
    
    ff=10.0
    if (ts<0.0)then!ts为土壤平均温
        gr=0.0
        gr=(1-exp(-ff*mt1))*((1+exp(-a*(t-b)))**(-1.0))
    endif
    
    if (id>79) then
        sdd=gd
    else
        sdd=0.0
    endif
    sd=sd+sdd!sd0=-3.724
    
    snd=gn
    sn=sn+snd!sn0=-8.376
    ssd=gs
    ss=ss+ssd!ss0=-1.359
    srd=gr
    sr=sr+srd!sr0=-3
    
    !sdc=25.76
    !snc=27.708
    !ssc=14.596
    !src=25.0
    
    if (ss<=0.0)then
        fs=0.0
    elseif(ss<ssc.and.ss>0.0)then
        fs=0.5*(sin((2*pi/ssc)*(ss-(ssc/4)))+1)
    else
        fs=0.0
    endif
    
    if (sr<=0.0)then
        fr=0.0
    elseif(sr<src.and.sr>0.0)then
        fr=0.5*(sin((2*pi/src)*(sr-(src/4)))+1)
    else
        fr=0.0
    endif
    
    if (sd<0.0)then
        fd=0.0
    elseif(sd>=0.0.and.sd<=sdc)then
        fd=4*(((sdc*sd)**0.5)-sd)/sdc
    else
        fd=0.0
    endif
    
    if (id<fj(iy,2))then
        fn=0.0
    elseif(id>=fj(iy,2).and.id<=fj(iy,3))then
        fn=4.0*(((snc*sn)**0.5)-sn)/snc
    else
        fn=0.0
    endif 
    end
   

    real function te(sd)
    real sd,q,sdc,teearly,telate
    q=9.9
    sdc=25.76
    teearly=10.69
    telate=8.79
    if (sd<(sdc/q))then
        te=teearly
    else
        te=telate
    endif
    return
    end
    

    real function twa(sd)
    real sd,q,sdc,twaearly,twalate
    q=9.9
    sdc=25.76
    twaearly=10.69
    twalate=8.79
    if (sd<(sdc/q))then
        twa=twaearly
    else
        twa=twalate
    endif
    return
    end   
    

    subroutine nrows(dbh,ht,sd,nrows0,lcell,dcell)
    real dbh,ht,q,qq,pi,lcell,dcell,sd,sdc,nrows0
    qq=0.6
    q=9.9
    pi=3.1415
    sdc=25.76
    nrows0=qq*(ht/lcell)*pi*(dbh/dcell)
    return
    end
    
    real function cet(te,t,sd,lcell,dcell)!早晚材物候(细胞增大),空气温度,物候
    real te,ii,v,msuc,r,t,mc,sd,sdc,q,lcell,dcell
    t=t+273.15
    q=9.9
    sdc=25.76
    ii=2.0
    msuc=342.3
    r=8.314
    mc=12.01
    if (sd<(sdc/q))then
        lcell=2.59
        !dcell=35.7*(10.0**(-3.0))
    else
        lcell=2.73
        !dcell=24.2*(10.0**(-3.0))
    endif
    v=lcell*(10.0**-1.0)*((dcell*(10.0**-1.0))**2.0)
    cet=ii*v*(msuc/(1000.0*r*t))*12.0/(msuc*te)
    return
    end
    
    real function an(fleac,spe)!可能的最大叶量，树种
    real fleac
    integer i1,spe
    i1=int(spe)
    if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then!落叶
        zn=1.0
        ln=35
    else
        zn=1.5
        ln=34.2
    endif
    an=fleac/(zn*ln)
    return
    end
    
    real function pmax(i2,age)
    real am,age,lf
    integer k,i2
    k=int(age/30.0)
    if (i2.eq.1.or.i2.eq.2.or.i2.eq.3.or.i2.eq.4)then
        if (k>12)then
           k=12
        elseif(k.eq.0)then
           k=1
        endif
    else
        if (k>12)then
          k=12
        elseif(k.eq.0)then
          k=1
        endif
    endif
    
    lf=real(k)
    if (i2.ne.5 .and. i2.ne.6 .and. i2.ne.9) then
        pmax=exp(4.3-0.65*log(lf))!lf
    else
        pmax=exp(4.0-0.50*log(lf))!lf
    endif
    
    !!!!
    if (i2.eq.1.or.i2.eq.2.or.i2.eq.3.or.i2.eq.4)then
        pmax=-0.96*(lf-12.0)**2.0+121.0!lf
        !yy2=lf**2.0+20.0;
    endif
    !!!!
    pmax=pmax*(10.0**-9.0)*(10.0**-3.0)*12.0*3600.0*500.0
    if (i2.eq.1 .or. i2.eq.2 .or. i2.eq.3 .or. i2.eq.4) then
        pmax=pmax
    else
        pmax=pmax/1.1
    endif
    
    return
    end
    
    
    real function treeage(dbh,h)
    real y,dbh,a,b,c,d
    y=dbh!cm
    a=210.115
    b=5.3523
    c=0.2655
    d=0.0064
    treeage=-a+exp(b+c*dbh+d*h);!(log(1.0-(y/60.0)**(1.0/1.165)))/(-0.07)
    return
    end
    
    

    subroutine treeld(par,parmax,treeld0)
    real par,parmax,a0,a3
    a0=0.43
    a3=0.2348
    treeld0=a0*exp(-a3*(par/parmax))
    !treed=exp(-a3*(par/parmax))
    end
    
    
   
    real function tens(sfc,sw,silt,sand)
    real sfc,sw,silt,sand
    real q0,k,qs,p
    p=sw/sfc
    if (silt<50.0)then
        q0=-3.68;
        k=19.92;
        qs=(q0+k*log(1.0+p))/(1.0+p);
    else
        q0=-4.49;
        k=26.74;
        qs=(q0+k*log(1.0+p))/(1.0+p);
    endif
    tens=qs
    if (tens>0.0)then
        tens=0.0
    endif
    end


    subroutine daysin2(tmax,tmin,tt)
    real tmp,prod,tran,rleaf,rabove,stem,xx,rstem,twig,rtw,leaf,rl
    real k,pi,fj,fj1,fj2,x1,x2,c1,c2,c3,c4,c33,c44,dl1,dl11,dl2,dl22,kk
    real a,b,fj11,fj22,a11,a22,b11,b22,c11,c22,d11,d22
    real fjj1,fjj2,fjj3,fjj4
    real tt(86400)
    integer i
    
    pi=3.1415

    d11=tmin+(tmax-tmin)/2.0
    c11=(tmax-tmin)/2.0
    a11=pi/(12.0*3600.0)
    b11=5.0*pi/6.0!!!8.5
    
    do i=1,86400
       !kk=real(i)/3600.0
       kk=real(i)
       tt(i)=c11*sin(a11*kk-b11)+d11
    enddo
    end
    
    