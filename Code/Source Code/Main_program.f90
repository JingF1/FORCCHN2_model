    !!!!!!!!!!!!!!!!!!!!!!!!!Main parogram for call!!!!!!!!!!!!!!!!!!!!
    subroutine forcchn2(fj,yxc2,dayout,yearout,ntrees,ny0,ny,ndays,lat,lon,ele,tmax,tmin,tmean,dayl1,prec,ra,rh,wind,sfc,pwp,vw,sc0,sn0,silt,sand,class,evergr0,deci0,lai0,co2)
    !!!!output: fj: the phenology dates, which included the start time of leaf growth (SOS) and the end time of leaf growth (EOS)
    !!!!output: yxc2: the allocation parameter of each soil pool, which can be used as input instead of the initial soil allocation parameters
    !!!!output: dayout: the daily carbon dynamics, which included above- and belowground biomass, gross primary productivity (GPP), above- and belowground respiration, soil heterotrophic respiration, litter-fall biomass, and soil carbon
    !!!!output: yearout: the yearly carbon dynamics
    !!!!input: the variables can be found in the example run
    integer ny,ny0,ndays,class,iderr,ii3,kdays,nplot
    real ta,dmean,dayl1(ny,366),ds,yb
    real lat,lon,ele,tmax(ny,366),tmin(ny,366), prec(ny,366),ra(ny,366),rh(ny,366),wind(ny,366),sfc,pwp,vw,sc,snn,lai,evergr,deci,co2(ny,366),re_soil,tt(366),tmean(ny,366),tmean1,gpp(ndays,ntrees),nscfast(ny,2),nscd(ndays,7)
    real sc1,sn1,lai1,evergr1,deci1,sc0,sn0,lai0,evergr0,deci0
    real dayout(ndays,8),yearout(ny,8),key(900,3),soc,w_soil1(4)
    real gaparea,par_soil(8),c_soil(15),n_soil(15),l_soil(15),w_soil(4),npp_c,npp_n,fall_c,fall_n,yxc2(11),snc0,forcn(600,19),par_tree(28,9),hold0,wata,parcano(120)
    real tp_resp,tl_resp,t_photo(9),ts_resp,outlai1(999),dry_photo(9),dry_lit,dry_hum,phenfac(366),phenot
    real c_out,n_out,cout(8),hr,year(ny,8),day(ndays,8),outlai(999),fall11(ndays,4),phe1(366),tb,sn11,ypheno1
    real nppmean10,soilresp10,fall10,treeresp10,nepend(12),no3
    integer iy,id,IDUM,nn,ndy,ii,ij,nd1(12),nd2(12),nd(12),lyear
    real nppmean0,soilresp0,fall0,treeresp0,leafwt,leafwt0
    real pp,son,ttc(4),zttc(ndays,4),nva,znva(ny),n_soil11(ny,15),znl11,znl12(ny,366),nav(ndays,15),nim,nim1(ndays),rn1(ndays),no31(ndays),bdh(ny,ntrees), th(ny,ntrees),zbdh(ny,ntrees),everytree(ntrees,5),sw(366),htmax,htmax1,petd1(366)
    real gn,gs,gd,gr,sn,ss,sd,sr,fn,fs,fd,fr,snc,ssc,sdc,src,nscstart(600),ypheno(366),para(1,6),sn00(366,1)
    real xs0,xn0,ne(ntrees),nwa(ntrees),ncell(ntrees),xs(ntrees),xn(ntrees),ff1(ntrees,366),ff2(ntrees,366),ff3(ntrees,366),ff4(ntrees,366),ff5(ntrees,366),ff6(ntrees,366),ff7(ntrees,366),ff8(ntrees,366),ff9(ntrees,366),nsc0(ntrees),ff10(ntrees,366),treeage0(ntrees),tran(ntrees),ht(ntrees,1)
    integer iii,jjj,fj(ny,3),num,dbhnum,itree,woodp(ny,2)
    real ls(ntrees),ls0(ntrees),ln(ntrees),ln0(ntrees),ld0,ldy(ndays,1),ldp(ny,366),ld1(366,ny),lsp,lnp,ld(ny,366,ntrees),ld11(ntrees),pheno12(366,1),tever
    real fj1,fj2,fj3,fj4,fj5,gpp1(ndays,ntrees),npp,npp1(ndays,1)
    real ls01,ln01,leafn(ndays,3),ta1,tb1,ds1
    real bb1(100),dinc1(366,50),dinc2(ny*366,50),dinc3(100,50)
    real dorm(ntrees),dist(ntrees),gmax
    integer time(ntrees)
    data nd1/31,28,31,30,31,30,31,31,30,31,30,31/
    data nd2/31,28,31,30,31,30,31,31,30,31,30,31/
    real,allocatable :: leafage(:,:,:)
    real dd11(366,ntrees),d11(366,ntrees),cfj(366,3)
    
    gaparea=400.0
    ntrees=0
    nplot=1
    call param(par_tree)
    znl12=0.0
    leafwt0=0.0
    nppmean10=0.0
    soilresp10=0.0
    fall10=0.0
    treeresp10=0.0
    nppmean0=0.0
    soilresp0=0.0
    fall0=0.0
    treeresp0=0.0
    n_soil=0.0
    do kk=1,12
        nepend(kk)=0.0
    end do
    do id=1,ndays
        do il=1,8
            dayout(id,il)=0.0
        end do
    end do
    do iy=1,ny
        do il=1,8
            yearout(iy,il)=0.0
        end do
    end do
    do kk=1,ny*12
        outlai1(kk)=0.0
    end do
    gmax=0.15
    do iplot=1,nplot

        !!!!!!!!!!!!!!!!!!!!!!initialization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hr=2.0
        sc1=sc0
        sn1=sn0
        evergr1=evergr0
        deci1=deci0
        !lai1=max(min(lai0,7.0),0.5)

        !!!!!!!!!!!!!!!!vegetation initialization
        call initialer2(trh,sc1,sn1,sfc,pwp,vw,silt,sand,class,evergr1,deci1,lai0,lat,ele,gaparea,par_soil,c_soil,n_soil,l_soil,w_soil,hold0,wata,forcn,ntrees,par_tree,yxc2,tever)

        
        do id=1,ntrees
            nsc0(id)=forcn(id,9)
            treeage0(id)=treeage(forcn(id,2)*2.0,forcn(id,3))
        enddo

        htmax1=0.0
        if (ntrees.eq.1)then
            htmax=forcn(1,3)
        else
            do id=1,ntrees-1
                htmax=max(forcn(id,3),forcn(id+1,3))
                htmax1=max(htmax,htmax1)
            enddo
        endif

        call zeros2(ny,8,year)
        call zeros2(ndays,8,day)

        !!!!!!!!!!!!!predict phenology
        ta1=0.0
        tb1=0.0
        ds1=0.0
        if (fj(1,2).ne.0)then
            do kk=1,fj(1,2)
               ta1=ta1+tmean(1,kk)
            enddo
            ta1=ta1/fj(1,2)
        else
            ta1=0.0
        endif
        if (fj(1,3).ne.0)then
            do kk=1,fj(1,3)
               tb1=tb1+tmean(1,kk)
               ds1=ds1+dayl1(1,kk)
            enddo
            tb1=tb1/fj(1,3)
            ds1=ds1/fj(1,3)
        else
            tb1=0.0
            ds1=0.0
        endif    
        ta=0.0/-1.0692*ta1!!!!the effective temperature threshold of spring
        tb=30.0/8.6060*tb1!!!!the effective temperature threshold of autumn
        ds=12.5/12.0*ds1!!!!the effective photoperiod threshold of autumn
        call pheno11(tmean(1,:),ta,pheno12)
        call phe(phe1,tmean(1,:),dayl1(1,:),tb,ds)
        sn11=-pheno12(fj(1,2),1)!!!!the critical heat theshold of spring
        dmean=pheno12(fj(1,3),1)+sn11!!!!the critical heat theshold of autumn
        yb=phe1(fj(1,3))
        snc0=0.0
        fj=0
        do iy=1,ny
            nn=iy+ny0-1

            if (mod(nn,4).eq.0) then
                ndy=366
            else
                ndy=365
            endif

            do it=1,900
                key(it,1)=0.0
                key(it,2)=0.0
            end do
            do it=1,8
                year(iy,it)=0.0
            end do

            do id=1,ndy
                tt(id)=tmean(iy,id)
            end do
            call pheno(tt,phenfac)
            lyear=0.0
            do id =1,ndy
                lyear=lyear+phenfac(id)
            end do
            lai0=tever!!!!!!!!!!!!

            call pheno11(tmean(iy,:),ta,pheno12)
            call phe(phe1,tmean(iy,:),dayl1(iy,:),tb,ds)
            pheno12=pheno12+sn11
            
 
        do id=1,ndy
            if (pheno12(id,1).ge.0.0)then
                fj(iy,2)=id!!!!!!!!!!!!!!!SOS(unit: doy)
                exit
            endif
        enddo
        do id=1,ndy
            if (phe1(id).le.yb)then
                fj(iy,3)=min(365,id)!!!!!!!!!!!!!!!!!EOS(unit: doy)
            endif
        enddo
        if(fj(iy,2)>=fj(iy,3))then
            goto 30
        endif
        enddo
        do iy=1,ny
           fj(iy,1)=iy+ny0-1
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!begin to predict carbon dynamics and tree growth
        ii=0
        do iy=1,ny
            nn=iy+ny0-1
            if (mod(nn,4).eq.0) then
                ndy=366
            else
                ndy=365
            endif
            do it=1,900
                key(it,1)=0.0
                key(it,2)=0.0
            end do
            do it=1,8
                year(iy,it)=0.0
            end do
            do id=1,ndy
                tt(id)=tmean(iy,id)
            end do
            call pheno(tt,phenfac)
            phenfac(1:fj(iy,2))=0.0
            phenfac(fj(iy,2)+1:fj(iy,3))=1.0
            phenfac(fj(iy,3)+1:366)=0.0

            call pheno(tt,phenfac)
            lyear=0.0
            do id =1,ndy
                lyear=lyear+phenfac(id)
            end do

            call pheno11(tmean(iy,:),ta,pheno12)
            snc0=pheno12(fj(iy,3),1)+sn11
            sn00=pheno12+sn11

            ypheno=0.0
            ypheno(fj(iy,2):fj(iy,3))=1.0
            lai0=tever
            num=ndy
            allocate(leafage(ntrees,2,num))
            leafage=0.0
            sn=0.0
        
            
            !!!!calculate on daily step
            do id=1,ndy
                ii=ii+1
                do ij=1,ntrees
                    ht(ij,1)=forcn(ij,3)
                enddo
                w_soil1=w_soil
                if (iy.eq.1) then
                    !!!soil initialization
                    call newenv(id,lat,ele,lai0,tmax(iy,id),tmin(iy,id),tmean(iy,id),prec(iy,id),ra(iy,id),rh(iy,id),wind(iy,id),w_soil,par_soil,parmax,dl,tday,t_photo,tave,tp_resp,ts_resp,tl_resp,dry_photo,dry_lit,dry_hum,par_tree,runoff,aet_pet,rn,htmax,pwp,petd) !!!!!!!!tp,tl,ts
                else
                    call newenv(id,lat,ele,lai,tmax(iy,id),tmin(iy,id),tmean(iy,id),prec(iy,id),ra(iy,id),rh(iy,id),wind(iy,id),w_soil,par_soil,parmax,dl,tday,t_photo,tave,tp_resp,ts_resp,tl_resp,dry_photo,dry_lit,dry_hum,par_tree,runoff,aet_pet,rn,htmax,pwp,petd)
                endif
                sw(id)=w_soil(1)+w_soil(2)
                sw(id)=min(sw(id),sfc)
                petd1(id)=petd

                !!!!predict phenology of leaf, fine roots, and woods
                call gf(ny,iy,id,tmean(iy,id),tmean(iy,id),sw(id),par_soil,sfc,gn,gs,gd,gr,sn,ss,sd,sr,fn,fs,fd,fr,snc0,ssc,sdc,src,fj)
                call phenoleaf(tmean(iy,id),sn00(id,1),snc0,fn,gn)
                !rn1(ii)=rn
                tmean1=tmean(iy,id)
                
                !!!!predict soil dynamics
                call soil_cn1(tmean1,sw(id),iy,id,sand,runoff,c_soil,n_soil,ph,l_soil,par_soil,tl_resp,dry_lit,dry_hum,re_soil,c_out,n_out,nim,no3)!!!!!!!!!!!!!!!!!!!
                nim1(ii)=nim
                no31(ii)=no3
                day(ii,6)=re_soil

                !!!!allocation of Effective Radiation (layer by layer)
                call newphenology(forcn,ntrees,tday,par_tree,parmax,gaparea,lai,parcano,key,sn,ss,sd,sr,snc0,ssc,sdc,src,phenfac(id))

                !!!!!calculate transpiration
                call tran1(id,lat,ele,ntrees,ht,parcano,lai,tmax(iy,id),tmin(iy,id),tmean(iy,id),prec(iy,id),rh(iy,id),wind(iy,id),tran,htmax,pwp,w_soil1)

                if (mod(id,30).eq.15) then
                    id1=int(id/30)+1+(ny-1)*12
                    outlai(id1)=lai
                endif
                
                !!!!co2 fertilization
                c=(co2(iy,id)-350.0)/(co2(iy,id)+700.0)+1.0
                phenot=phenfac(id)
                ypheno1=ypheno(id)

                !!!!tree daily growth and carbon dynamics
                call newdaily_c(lat,ndays,ii,id,iy,tmax(iy,id),tmin(iy,id),sw(id),silt,sand,sfc,tran,ny,num,tmean1,ntrees,dl,c,forcn,t_photo,dry_photo,tp_resp,tl_resp,c_soil,n_soil,l_soil,par_tree,parcano,gaparea, npp_c,npp_n,fall_c,fall_n,phenot,cout,leafwt,gn,fn,fs,fd,fr,ld11,ls,ln,ne,nwa,ncell,xn,xs,ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,fj,leafage,nsc0,treeage0,bb1,time,dorm,dist,dd11,d11,cfj,key,dinc1,dinc3,gpp1,ls0,nscstart,woodp,ttc,nscfast,nscd,wt1,npp,ypheno1,leafn)
                if (iy.eq.ny) then
                    leafwt0=max(leafwt,leafwt0)
                end if
                npp1(ii,1)=npp
                fall11(ii,1)=ttc(1)/gaparea
                fall11(ii,2)=ttc(2)/gaparea
                fall11(ii,3)=ttc(3)/gaparea
                fall11(ii,4)=ttc(4)/gaparea
                day(ii,1)=cout(1)/gaparea
                day(ii,2)=cout(2)/gaparea
                day(ii,3)=cout(3)/gaparea
                day(ii,4)=cout(4)/gaparea
                day(ii,5)=cout(5)/gaparea
                day(ii,7)=cout(7)/gaparea
                day(ii,8)=0.0
                do ij=1,11
                    day(ii,8)=day(ii,8)+c_soil(ij)
                end do
                year(iy,3)=year(iy,3)+day(ii,3)
                year(iy,4)=year(iy,4)+day(ii,4)
                year(iy,5)=year(iy,5)+day(ii,5)
                year(iy,6)=year(iy,6)+day(ii,6)
                year(iy,7)=year(iy,7)+day(ii,7)
                do ij=1,15
                    nav(ii,ij)=n_soil(ij)
                enddo
                !write(*,*) id
            end do!!!!end the daily step

            do it=1,900
                key(it,1)=key(it,1)/lyear
                key(it,2)=key(it,2)/lyear
            end do

            !!!!!!!!!!!!!!!calculate on yearly step
            call grow2(ntrees,par_tree,forcn,l_soil,c_soil,n_soil,cout,key,gaparea,coutn)
            year(iy,1)=max(cout(1)/gaparea,0.001)
            year(iy,2)=max(cout(2)/gaparea,0.001)
            year(iy,8)=max(cout(8),0.001)
            deallocate(leafage)
        end do
        
        !!!!!!!!!!!!!!!!output results
        do il=1,8
            do iy=1,ny
                yearout(iy,il)=yearout(iy,il)+year(iy,il)/nplot
            end do
            do id=1,ndays
                dayout(id,il)=dayout(id,il)+day(id,il)/nplot
            end do
        end do

        do ii3=1,11
            yxc2(ii3)=c_soil(ii3)/sum(c_soil(1:11))
        end do
30  end do!!!!end the calculation
    return
    end