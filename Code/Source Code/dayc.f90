
    subroutine newdaily_c(lat,ndays,ip,id,iy,tmax,tmin,sw,silt,sand,sfc,tran,ny,num,tmean1,ntree,dl,c,forcn,t_photo,dry_photo,tp_resp,&
    tl_resp,c_soil,n_soil,l_soil,par_tree,parcano,&
    gaparea, npp_c,npp_n,fall_c,fall_n,phenot,&
    cout,wt,gn,fn,fs,fd,fr,ld,&
    ls,ln,ne,nwa,ncell,xn,xs,ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,&
    fj,leafage,nsc0,treeage0,bb1,time,dorm,dist,dd11,d11,cfj,key,dinc1,dinc3,gpp1,ls0,nscstart,woodp,ttc,nscfast,nscd,wt1,npp,ypheno1,leafn)
    integer ntree,i1,i2,it,id,iy,ny,num,fj(ny,3),ip,method,woodp(ny,2)
    real forcn(600,19),t_photo(9),dry_photo(9),tp_resp,tl_resp,&
    c_soil(15),n_soil(15),l_soil(15),par_tree(28,9),&
    parcano(120), gaparea, npp_c,npp_n,fall_c,fall_n,&
    cout(8),dl,c,phenot,grwresp0,tran(ntree),ls0(ntree),kp2(366,1)
    real leaf,root,wood,kl,tlai,ns,am,sl,par,prod,prodn,fineroot,wood1,nscfast(ny,2),lat
    real phen,ll,frl,lln,frln,yxd,tmp1,tmp2,tmp3,tmp4,tmp6,ypheno1,kp1
    real d,h,hr,b,ck,above,below,grwresp,ft,wt,wt1,key(900,3)
    real above1,below1,npp,npp1,npp2,nppn,dleaf,abc,bec,abn,ben,kkk
    real above2,below2,yy,leaf5,k,tmean1,lr,fj1,fj2,fj3,fj4,fj5,fj6,leafn(ndays,3)
    real sn,ss,sd,sr,snc,ssc,sdc,src,gn,gs,gd,gr,fn,fs,fd,fr,ld(ntree),ls(ntree),ln(ntree)
    real xsd,xnd,ncelld,xs(ntree),xn(ntree),ncell(ntree),ann,cr,cs,cet1,cwa,cet0,cwa0
    real nrows0,te0,twa0,ne(ntree),nwa(ntree),ned,nwad,p1,lcell,dcell
    real ttc(4),ff1(ntree,366),ff2(ntree,366),ff3(ntree,366),ff4(ntree,366),&
    ff5(ntree,366),ff6(ntree,366),ff7(ntree,366),ff8(ntree,366),ff9(ntree,366),&
    kp,prod1,nsc0(ntree),ff10(ntree,366),treeage0(ntree),parmax,treeld0,treed!!!!!调试
    real leafage(ntree,2,num),age,pmax1,tlai1,prodtotal(366),llx,frlx,llnx,frlnx
    real nsc2,nsc1,tmp,woodinc,dinc,dinc1(366,50),dinc3(100,50)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!单样地不超过100棵树
    real sw,silt,sand,sfc,tmax,tmin,prod11(86400),tran1(86400),tt(86400),prodk,prodk1
    !real*8 kpmu(864000),dbhf,y00(40,2,5,ntree),y001(1,5,ntree),y0(40,2,5),y01(1,5)
    !real*8 y(40,2,5),y11(1,5),treep(ntree,7),kpmu1(864000,3),vsap
    real rleaf,rabove,rleaf1(86400),rabove1(86400),fjj,nscd(ndays,7)
    real bb1(100),gpp1(ndays,ntree)
    !real ch(12,10000,ntree)
    real dd11(366,ntree),d11(366,ntree),grr(12),cet2,cfj(366,3),vv,numf,cp,nscstart(600)
    !integer numk(366,5,ntree),numk1(4,ntree),knum
    !CHARACTER state (10000,ntree)
    real dorm(ntree),dist(ntree),db,dr,swvol
    real Correct(3,10000,ntree)
    integer time(ntree)
    character(len=3)::dayy
    character(len=64)::filename

    hr=2.0
    do it=1,8
        cout(it)=0.0
    end do
    wt1=0.0
    npp_n=0.0
    ttc=0.0


    parmax=maxval(parcano)

    do it=1,ntree
        if (forcn(it,2) .gt. 0.01 .and. forcn(it,3).gt. 0.01) then
        i1=int(forcn(it,1))
        forcn(it,3)=max(1.0,forcn(it,3))
        i2=int(forcn(it,3))
        d=forcn(it,2)
        h=forcn(it,3)
        hr=2.0
        b=forcn(it,5)
        tlai=forcn(it,7)/forcn(it,13)
        tem=t_photo(i1)
        dry=dry_photo(i1)
        ns=n_soil(15)
        am=par_tree(2,i1)
        sl=par_tree(3,i1)
        kl=par_tree(4,i1)
        par=parcano(i2)
        phen=forcn(it,16)
        leaf=forcn(it,6)
        ft=forcn(it,8)
        ck=par_tree(18,i1)


        !!!!!limitation
        !k=1
        kp1=0.185!0.185
        k=min(1.0,(1-exp(-kp1*forcn(it,9)))/(1-exp(-kp1*nsc0(it))))
        k=max(0.5,k)

        !ls0=9.526
        ann=an(forcn(it,12),i1)
        rr=0.6
        !!!!!!!
        !p1=(xs(it)*(10.0**-3.0)/h)*d
        a0=3.1415*(d**2.0)
        !!!!!!!
        lr=0.08

        leaf5=leaf*max(leaf/max(forcn(it,12),0.001),0.5)
        above=forcn(it,17)+forcn(it,18)!fstem3(d,h,ck)+ftwig3(d,h,b,ck)!
        below=forcn(it,19)!froot3(d,h,hr,ck)!
        wood=above+below

        tmp1=max(above+below,0.00001)

        yyr=forcn(it,11)*15./(abs(forcn(it,11))*15.+tmp1)

        xx=min(max(forcn(it,11)*2.0, tmp1*0.01),tmp1*0.1)/tmp1
        above1=above*xx
        below1=below*xx

        yy=max(min(yyr*0.7+forcn(it,2)*2.0+0.7,3.0),0.01)
        above2=above*yy
        below2=below*yy

        grwresp=0.0!!!!!!!!!!growth respiration
        prod=0.0!!!!!!!!photosynthesis
        prodn=0.0!!!!!!!!!!
        ll=0.0!!!!!!!!!!!!!fall leaf
        lln=0.0!!!!!!!!!!!
        frl=0.0!!!!!!!!!!!!fall fine roots
        frln=0.0!!!!!!!!!!!
        nsc1=forcn(it,9)!!!!!!NSC pool


        !!!!!!!photosynthesis!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        prodtotal=0.0
        prod=0.0

        leafage(it,1,1)=forcn(it,6)
        !if (sn.gt.0.0 .or. ss.gt.0.0 .or. sd.gt.0.0 .or. sr.gt.0.0)then
        if(tmean1.gt.par_tree(22,i1))then
            method=1
            if (id.gt.fj(iy,2).and. id.lt.fj(iy,3))then
                age=90.0
            else
                age=365.0
            endif
            if (i1.eq.1.or.i1.eq.2.or.i1.eq.3.or.i1.eq.4.or.i1.eq.5.or.i1.eq.6.or.i1.eq.7)then
                pmax1=pmax(i1,age)!!!!2.0*
            else
                pmax1=pmax(i1,age)
            endif

            if (lat.gt.50.0)then
                pmax1=0.8*pmax1
            endif

            prod=max(phototree(par,dl,tlai,c,tem,dry,ns,pmax1,sl,&
                kl)*forcn(it,13),0.0)
            prodn=prod/150.0
            n_soil(13)=max(0.0,n_soil(13)-prodn/gaparea*n_soil(13)&
                /max(1.0e-8,n_soil(13)+n_soil(14)))
            n_soil(14)=max(0.0,n_soil(14)-prodn/gaparea*n_soil(14)&
                /max(1.0e-8,n_soil(13)+n_soil(14)))
            n_soil(15)=n_soil(13)+n_soil(14)
        else
            prod=0.0
            prodn=0.0
        endif



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!maintenance respiration
        call resptree(i1,leaf,above1,below1,tp_resp,tl_resp,&
        par_tree,resp_above,resp_below,ns)
        !resp_above=k*resp_above
        !resp_below=k*resp_below

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!falling  
        !call fall(i1,phenot,leaf5,ft,par_tree,ll,frl,lln,frln)   
        if (id<fj(iy,2))then
            call fall(i1,ypheno1,leaf,ft,par_tree,ll,frl,lln,frln)
            if (dry.le.0.2) then !!!!!!!drought condition
                if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then
                    ll=ll+leaf*0.02
                    frl=frl+leaf*0.02
                    lln=ll/par_tree(9,i1)
                    frln=frl/par_tree(9,i1)
                else
                    ll=ll+leaf*0.03
                    frl=frl+leaf*0.03
                    lln=ll/par_tree(9,i1)
                    frln=frl/par_tree(9,i1)
                end if
            endif
        else
            method=1
            if (method.eq.2)then
                do ii=1,id
                    leaf5=leafage(it,1,ii)
                    age=leafage(it,2,ii)
                    if (age>0.0)then!!!!新叶不凋落
                        call fall(i1,ypheno1,leaf5,ft,par_tree,llx,frlx,llnx,frlnx)!!!叶凋落
                        if (dry.le.0.2) then 
                            if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then
                                llx=llx+leaf5*0.02
                                frlx=frlx+leaf5*0.02
                                llnx=llx/par_tree(9,i1)
                                frlnx=frlx/par_tree(9,i1)!!!!!!!!!!!!!!!!!!
                            else
                                llx=llx+leaf5*0.03
                                frlx=frlx+leaf5*0.03
                                llnx=llx/par_tree(9,i1)
                                frlnx=frlx/par_tree(9,i1)!!!!!!!!!!!!!!!!!!!!!!!!!
                            end if
                        end if
                    else
                        llx=0.0
                        frlx=0.0
                        llnx=0.0
                        frlnx=0.0
                    endif
                    leafage(it,1,ii)=leafage(it,1,ii)-llx
                    ll=ll+llx
                    frl=frl+frlx
                    lln=lln+llnx
                    frlnx=frln+frlnx
                enddo
            else
                call fall(i1,ypheno1,leaf,ft,par_tree,llx,frlx,llnx,frlnx)
                ll=ll+llx
                frl=frl+frlx
                lln=lln+llnx
                frlnx=frln+frlnx
            endif
        endif

        !!!!!!!!!!!!!!!!!!!fall wood
        call fallw(i1,above2,below2,par_tree,abc,&
        bec,abn,ben)
        !abc=0.0
        xsd=0.0
        cs=0.0                       


        !!!!!leaf growth!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        method=1
        cr=0.0
        if (method.eq.1)then
            if (id.ge.fj(iy,2) .and. id.le.fj(iy,3))then
                if (i1.eq.5 .or. i1.eq.6 .or. i1.eq.9) then
                    dleaf=max((forcn(it,12)-forcn(it,6))*1.5,0.0)
                    !!bleaf=npp1+forcn(it,9)
                    cleaf=min(forcn(it,12)/15.,forcn(it,12)/15.)
                    dleaf=k*gn*fn*min(dleaf,cleaf)
                else
                    dleaf=max((forcn(it,12)-forcn(it,6))*1.5,0.0)
                    !!bleaf=npp1+forcn(it,9)
                    cleaf=min(forcn(it,12)/30.,forcn(it,12)/15.)!forcn(it,12)/60.
                    dleaf=k*gn*fn*min(dleaf,cleaf)
                endif
                !dleaf=max(0.0,min(forcn(it,9)+prod-resp_above-resp_below-grwresp-cr,dleaf))
            else
                !prod=prod
                !prodn=0.0
                xnd=0.0
                dleaf=0.0
            endif
        endif
        leafn(ip,1)=k
        leafn(ip,2)=gn
        leafn(ip,3)=fn

        grwresp0=dleaf*0.35
        grwresp=grwresp+grwresp0

        !!!!!!fine root growth
        if (id.ge.fj(iy,2) .and. id.le.fj(iy,3))then
            !cr=k*gr*fr*dleaf/3.0
            cr=dleaf/2.0!k*gr*fr*lr*forcn(it,12)/2.0
        else
            cr=0.0
        endif

        grwresp0=cr*0.35
        grwresp=grwresp+grwresp0


        !!!!!!!!!!!!!!!!!!!calculate leaf age          
        if (id.le.fj(iy,2))then
            leafage(it,1,id)=forcn(it,6)
            leafage(it,2,id)=365.0
        else
            leafage(it,1,id)=dleaf
            leafage(it,2,id)=1.0
        endif

        if (id.eq.1)then
            leafage(it,2,1)=1.0
        else
            do ii=2,id
                leafage(it,2,ii-1)=leafage(it,2,ii-1)+1.0
            enddo
        endif

        if (id.eq.1)then
            nscstart(it)=forcn(it,9)
        endif

        tmp=0.0
        woodinc=0.0
        grwresp0=tmp*0.3
        grwresp=grwresp+grwresp0
        npp=prod-resp_above-resp_below
        if (npp.gt.0.0 )  then 
            grwresp=0.25*npp
        else
            grwresp=0.0
        end if
        !npp=npp-grwresp


        if (id.ge.1.and.it.le.366.and.it.eq.4) then
            nscd(ip,1)=forcn(4,6)
            nscd(ip,2)=dleaf!prod
            nscd(ip,3)=ll
            nscd(ip,4)=cr
            nscd(ip,5)=resp_above
            nscd(ip,6)=resp_below
            nscd(ip,7)=grwresp
            npp=prod-resp_above-resp_below-grwresp
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        forcn(it,6)=max(0.0,forcn(it,6)-ll)+dleaf!!!leaf biomass
        forcn(it,8)=max(0.0,forcn(it,8)-frl)+cr!!!fine root biomass
        fjj=forcn(it,17)/(forcn(it,17)+forcn(it,18))
        !forcn(it,17)=max(0.0,(forcn(it,17)-abc*fjj))+(cs+cet1+cwa)*h/(h+hr)*fjj!!!stem biomass
        forcn(it,17)=max(0.0,(forcn(it,17)-abc*fjj))+woodinc*h/(h+hr)*fjj
        !forcn(it,18)=max(0.0,(forcn(it,18)-abc*(1.0-fjj)))+(cs+cet1+cwa)*h/(h+hr)*(1.0-fjj)!!!twig biomas
        forcn(it,18)=max(0.0,(forcn(it,18)-abc*(1.0-fjj)))+woodinc*h/(h+hr)*(1.0-fjj)
        !forcn(it,19)=max(0.0,(forcn(it,19)-bec))+(cs+cet1+cwa)*hr/(h+hr)!!!coarse root biomass
        forcn(it,19)=max(0.0,(forcn(it,19)-bec))+woodinc*hr/(h+hr)

        resp_above=resp_above
        resp_below=resp_below
        forcn(it,9)=forcn(it,9)+prod-resp_above-resp_below-grwresp-dleaf-cr-tmp!!NSC pool
        forcn(it,10)=forcn(it,9)/par_tree(9,i1)
        ff1(it,id)=dleaf
        ff2(it,id)=cs
        ff3(it,id)=cr
        ff4(it,id)=tmp!cet1+cwa
        ff5(it,id)=resp_above
        ff6(it,id)=resp_below
        ff7(it,id)=grwresp
        ff8(it,id)=prod
        ff9(it,id)=forcn(it,3)
        ff10(it,id)=forcn(it,9)

        yxd=forcn(it,9)/par_tree(9,i1)	
        if (forcn(it,10).gt.yxd) then
            yxd0=(forcn(it,10)-yxd)
            forcn(it,10)=forcn(it,10)-yxd0
            abn=yxd0*0.8
            ben=yxd0*0.2
        else	
            forcn(it,10)=forcn(it,10)
            abn=0.0
            ben=0.0
        endif

        if (ll.gt. 1.0e-15) then
            tmp1=0.85-0.018*100.0/par_tree(27,i1)/par_tree(9,i1)
            tmp0=ll/gaparea
            tmp2=tmp1*tmp0
            tmp3=(1.0-tmp1)*tmp0
            l_soil(1)=l_soil(1)+par_tree(27,i1)*0.01*tmp0
            c_soil(1)=c_soil(1)+ tmp2
            c_soil(2)=c_soil(2)+ tmp3
            tmp6=lln/gaparea
            tmp4=tmp6/(1.0+0.2*tmp2/tmp3)
            n_soil(2)=n_soil(2)+tmp4
            n_soil(1)=n_soil(1)+(tmp0-tmp4)
        end if
        if (frl .gt. 1.0e-15) then
            tmp0=frl/gaparea
            tmp2=tmp0!tmp1*tmp0
            tmp3=(1.0-tmp1)*tmp0
            l_soil(3)=l_soil(3)+par_tree(27,i1)*0.01*tmp0
            c_soil(3)=c_soil(3)+ tmp2
            c_soil(4)=c_soil(4)+tmp3
            tmp6=frln/gaparea
            tmp4=tmp6/(1.0+0.2*tmp2/tmp3)
            n_soil(4)=n_soil(4)+tmp4
            n_soil(3)=n_soil(3)+(tmp0-tmp4)
        end if

        c_soil(6)=c_soil(6)+abc/gaparea
        n_soil(6)=n_soil(6)+abn/gaparea
        c_soil(7)=c_soil(7)+bec/gaparea
        n_soil(7)=n_soil(7)+ben/gaparea

        cout(1)=above+cout(1)+forcn(it,6)
        cout(2)=below+cout(2)+forcn(it,8)
        cout(3)=cout(3)+prod
        cout(4)=cout(4)+resp_above+grwresp
        cout(5)=cout(5)+resp_below
        cout(7)=cout(7)+ll+frl+abc+bec
        ttc(1)=ttc(1)+abc 
        ttc(2)=ttc(2)+bec
        ttc(3)=ttc(3)+forcn(it,6)!ttc(3)+ll
        ttc(4)=ttc(4)+frl

        !gpp1(ip,it)=prod
        end if
        wt1=wt1+forcn(it,8)/gaparea
    enddo
    kp2(id,1)=forcn(1,9)
    end