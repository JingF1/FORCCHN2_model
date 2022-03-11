    !!!!srping phenology model
    subroutine pheno11(tt1,ta,pheno)
    real tt,a,b,pheno1,tt1(366,1),f,pheno(366,1)
    a=0.185
    b=18.4
    pheno1=0.0
    do i=1,366
        tt=tt1(i,1)
        if (tt<ta)then
            f=0.0
        else
            f=1.0/(1.0+exp(-a*(tt-b)))
        endif
        pheno1=pheno1+f
        pheno(i,1)=pheno1
    enddo
    return 
    end
    
    !!!!autumn phenology model
    subroutine phe(pheno,tt1,dayl1,tb,ds) 
    real tt1(366),dayl1(366),pheno(366),tb,cdd,tt,ds,dayl
    integer i
    pheno(1)=0
    do i=2,366
        tt=tt1(i)
        dayl=dayl1(i)
        if (dayl>=ds)then
            pheno(i)=0.0
        else
            if (tt<tb)then
                cdd=(tb-tt)*dayl/ds
            else
                cdd=0
            endif
            pheno(i)=pheno(i-1)+cdd
        endif
    enddo
    end
    
    !!!!calculate photoperiod
    subroutine period2(h,lat,doy)
    real lat1,lat,pi,darta,tmp1,tmp2,omiga,h
    integer doy
    pi=3.1415
    lat1=lat*pi/180.0;
    darta=0.409*sin(2.0*pi*doy/365.0-1.37);
    tmp1=-tan(darta)*tan(lat1);
    if (tmp1<= -0.99999)then
        tmp2=pi;
    elseif (tmp1>= 0.99999)then
        tmp2=0.0;
    else
        tmp2=acos(tmp1);
    endif
    omiga=tmp2;
    h=7.639437*omiga;
    end