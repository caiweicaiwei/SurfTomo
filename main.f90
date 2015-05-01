      ! CODE FOR SURFACE WAVE TOMOGRAPHY USING DISPERSION MEASUREMENTS
      ! VERSION:
      !      1.0 
      ! AUTHOR:
      !      HONGJIAN FANG. fanghj@mail.ustc.edu.cn
      ! PURPOSE:
      !      DIRECT INVERT SURFACE WAVE DISPERSION MEASUREMENTS FOR 3-D
      ! STUCTURE WITHOUT THE INTERMEDIATE STEP OF CONSTUCTION THE PHASE
      ! OR GROUP VELOCITY MAPS.
      ! REFERENCE:
      !      Ray tracing based direct inversion of surface dispersion
      !      for 3-D shallow crust structure:methodology and application
      ! HISTORY:
      !       2015/01/31 START TO REORGONIZE THE MESSY CODE   
      !

        program SurfTomo
            use lsmrModule, only:lsmr
            implicit none

! VARIABLE DEFINE

            character inputfile*80
            logical ex
            character dummy*40
            character datafile*80

            integer nx,ny,nz
            real goxd,gozd
            real dvxd,dvzd
            integer nsrc,nrc
            real weight,weight1,weight2,weight3
            real damp,damp1,damp2,damp3
            real minthk
            integer kmax,kmaxRc,kmaxRg,kmaxLc,kmaxLg
            real*8,dimension(:),allocatable:: tRc,tRg,tLc,tLg
            real,dimension(:),allocatable:: depz
        integer itn
        integer nout
        integer localSize
	real mean,std_devs,balances,balanceb
	integer msurf
	integer maxlevel,maxleveld
	real,parameter:: tolr=1e-6
	real,dimension(:),allocatable:: obst,cbst,wt,dtres
	real,dimension(:),allocatable:: pvall,depRp,pvRp
	real sta1_lat,sta1_lon,sta2_lat,sta2_lon
	real dist,dcal
        integer dall
        integer dRp
	integer istep
	real,parameter :: pi=3.1415926535898
	integer checkstat
	integer ii,jj,kk
        real, dimension (:,:), allocatable :: scxf,sczf
        real, dimension (:,:,:), allocatable :: rcxf,rczf
	integer,dimension(:,:),allocatable::wavetype,igrt,nrc1
	integer,dimension(:),allocatable::nsrc1,knum1
	integer,dimension(:,:),allocatable::periods
	real,dimension(:),allocatable::rw
	integer,dimension(:),allocatable::iw,col
        real,dimension(:),allocatable::dv
        real,dimension(:,:,:),allocatable::vsf
        real,dimension(:,:,:),allocatable::vsftrue
	character strf
	integer veltp,wavetp
	real velvalue
	integer knum,knumo,err
	integer istep1,istep2
	integer period
	integer knumi,srcnum,count1
	integer HorizonType,VerticalType
	character line*200
        integer iter,maxiter
        integer iiter,initer
        integer maxnar
        real acond
        real anorm
        real arnorm
        real rnorm
        real xnorm
        character str1
        real atol,btol
        real conlim
        integer istop
        integer itnlim
        integer lenrw,leniw
        integer nar
        integer m,maxvp,n
        integer i,j,k
        real Minvel,MaxVel
        real spfra
        integer sizex,sizey,sizez
        real noiselevel,anomaly
        integer ifsyn
        integer gsini
        real minv,grad,gridsp

! OPEN FILES FIRST TO OUTPUT THE PROCESS
            open(34,file='IterVel.out')
            nout=36
            open(nout,file='lsmr.txt')
            open(66,file='SurfTomo.log')
! OUTPUT PROGRAM INFOMATION            
             write(*,*)
             write(*,*),'                         S U R F  T O M O'
             write(*,*),'PLEASE contact Hongjain Fang &
                 (fanghj@mail.ustc.edu.cn) if you find any bug'
             write(*,*)
             write(66,*)
             write(66,*),'                         S U R F  T O M O'
             write(66,*),'PLEASE contact Hongjain Fang &
                 (fanghj@mail.ustc.edu.cn) if you find any bug'
             write(66,*)

! READ INPUT FILE
            if (iargc() < 1) then
                write(*,*) 'input file [SurfTomo.in(default)]:'
                read(*,'(a)') inputfile
                if (len_trim(inputfile) <=1 ) then
                    inputfile = 'SurfTomo.in'
                else
                    inputfile = inputfile(1:len_trim(inputfile))
                endif
            else
                call getarg(1,inputfile)
            endif
            inquire(file = inputfile, exist = ex)
            if (.not. ex) stop 'unable to open the inputfile'

            open(10,file=inputfile,status='old')
            read(10,'(a30)')dummy
            read(10,'(a30)')dummy
            read(10,'(a30)')dummy
            read(10,*)datafile
            read(10,*) nx,ny,nz
            read(10,*) goxd,gozd
            read(10,*) dvxd,dvzd
            read(10,*) nsrc,nrc
            read(10,*)weight1,weight2,weight3
            read(10,*)damp1,damp2,damp3
!            read(10,*)mode,iflsph
            read(10,*)minthk
            read(10,*)Minvel,Maxvel
	    read(10,*)HorizonType,VerticalType
            read(10,*)maxlevel,maxleveld
            read(10,*)maxiter,initer
            read(10,*)spfra
            read(10,*)kmaxRc
            write(*,*) 'model origin:latitude,longitue'
            write(*,'(2f10.3)') goxd,gozd
            write(*,*) 'grid spacing:latitude,longitue'
            write(*,'(2f10.3)') dvxd,dvzd
            write(*,*) 'model dimension:nx,ny,nz'
            write(*,'(3i5)') nx,ny,nz
            write(66,*) 'model origin:latitude,longitue'
            write(66,'(2f10.3)') goxd,gozd
            write(66,*) 'grid spacing:latitude,longitue'
            write(66,'(2f10.3)') dvxd,dvzd
            write(66,*) 'model dimension:nx,ny,nz'
            write(66,'(3i5)') nx,ny,nz
            if(kmaxRc.gt.0)then
               allocate(tRc(kmaxRc),&
                stat=checkstat)
               if (checkstat > 0) stop 'error allocating RP'
            read(10,*)(tRc(i),i=1,kmaxRc)
            write(*,*)'Rayleigh wave phase velocity used,periods:(s)'
            write(*,'(50f5.1)')(tRc(i),i=1,kmaxRc)
            write(66,*)'Rayleigh wave phase velocity used,periods:(s)'
            write(66,'(50f5.1)')(tRc(i),i=1,kmaxRc)
            endif
            read(10,*)kmaxRg
            if(kmaxRg.gt.0)then
               allocate(tRg(kmaxRg),&
                stat=checkstat)
               if (checkstat > 0) stop 'error allocating RP'
            read(10,*)(tRg(i),i=1,kmaxRg)
            write(*,*)'Rayleigh wave group velocity used,periods:(s)'
            write(*,'(50f5.1)')(tRg(i),i=1,kmaxRg)
            write(66,*)'Rayleigh wave group velocity used,periods:(s)'
            write(66,'(50f5.1)')(tRg(i),i=1,kmaxRg)
            endif
            read(10,*)kmaxLc
            if(kmaxLc.gt.0)then
               allocate(tLc(kmaxLc),&
                stat=checkstat)
               if (checkstat > 0) stop 'error allocating RP'
            read(10,*)(tLc(i),i=1,kmaxLc)
            write(*,*)'Love wave phase velocity used,periods:(s)'
            write(*,'(50f5.1)')(tLc(i),i=1,kmaxLc)
            write(66,*)'Love wave phase velocity used,periods:(s)'
            write(66,'(50f5.1)')(tLc(i),i=1,kmaxLc)
            endif
            read(10,*)kmaxLg
            if(kmaxLg.gt.0)then
               allocate(tLg(kmaxLg),&
                stat=checkstat)
               if (checkstat > 0) stop 'error allocating RP'
            read(10,*)(tLg(i),i=1,kmaxLg)
            write(*,*)'Love wave group velocity used,periods:(s)'
            write(*,'(50f5.1)')(tLg(i),i=1,kmaxLg)
            write(66,*)'Love wave group velocity used,periods:(s)'
            write(66,'(50f5.1)')(tLg(i),i=1,kmaxLg)
            endif
            read(10,*)ifsyn
            read(10,*)sizex,sizey,sizez
            read(10,*)noiselevel,anomaly
            read(10,*)gsini
            read(10,*)minv,grad,gridsp
            close(10)
            kmax=kmaxRc+kmaxRg+kmaxLc+kmaxLg

! READ MEASUREMENTS            
            open(unit=87,file=datafile,status='old')
            allocate(scxf(nsrc,kmax),sczf(nsrc,kmax),&
            rcxf(nrc,nsrc,kmax),rczf(nrc,nsrc,kmax),stat=checkstat)
            if(checkstat > 0)then
            write(6,*)'error with allocate'
            endif
            allocate(periods(nsrc,kmax),wavetype(nsrc,kmax),&
            nrc1(nsrc,kmax),nsrc1(kmax),knum1(kmax),&
            igrt(nsrc,kmax),stat=checkstat)
            if(checkstat > 0)then
            write(6,*)'error with allocate'
            endif
            allocate(obst(nrc*nsrc*kmax),&
            stat=checkstat)
            allocate(pvall(nrc*nsrc*kmax),depRp(nrc*nsrc*kmax),&
            pvRp(nrc*nsrc*kmax),stat=checkstat)
            IF(checkstat > 0)THEN
            write(6,*)'error with allocate'
            ENDIF
            istep=0
            istep2=0
            dall=0
            knumo=12345
            do 
            read(87,'(a)',iostat=err) line
            if(err.eq.0) then
            if(line(1:1).eq.'#') then
            read(line,*) str1,sta1_lat,sta1_lon,period,wavetp,veltp
            if(wavetp.eq.2.and.veltp.eq.0) knum=period
            if(wavetp.eq.2.and.veltp.eq.1) knum=kmaxRc+period
            if(wavetp.eq.1.and.veltp.eq.0) knum=kmaxRg+kmaxRc+period
            if(wavetp.eq.1.and.veltp.eq.1) knum=kmaxLc+kmaxRg+&
               kmaxRc+period
            if(knum.ne.knumo) then
            istep=0
            istep2=istep2+1
            endif
            istep=istep+1
            istep1=0
            sta1_lat=(90.0-sta1_lat)*pi/180.0
            sta1_lon=sta1_lon*pi/180.0
            scxf(istep,knum)=sta1_lat
            sczf(istep,knum)=sta1_lon
            periods(istep,knum)=period
            wavetype(istep,knum)=wavetp
            igrt(istep,knum)=veltp
            nsrc1(knum)=istep
            knum1(istep2)=knum
            knumo=knum
            else
            read(line,*) sta2_lat,sta2_lon,velvalue
            istep1=istep1+1
            dall=dall+1
            sta2_lat=(90.0-sta2_lat)*pi/180.0
            sta2_lon=sta2_lon*pi/180.0
            rcxf(istep1,istep,knum)=sta2_lat
            rczf(istep1,istep,knum)=sta2_lon
            call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist)
            obst(dall)=dist/velvalue
            pvall(dall)=velvalue
            nrc1(istep,knum)=istep1
            endif
            else
            exit
            endif
            enddo
            close(87)
            allocate(depz(nz), stat=checkstat)
           maxnar = dall*nx*ny*nz*spfra!sparsity fraction
           maxvp = (nx-2)*(ny-2)*nz
           allocate(dv(maxvp), stat=checkstat)
           allocate(vsf(nx,ny,nz), stat=checkstat)
           allocate(vsftrue(nx,ny,nz), stat=checkstat)
	
	allocate(rw(maxnar), stat=checkstat)
	if(checkstat > 0)then
	   write(6,*)'error with allocate: real rw'
	endif
	allocate(iw(2*maxnar+1), stat=checkstat)
	if(checkstat > 0)then
	   write(6,*)'error with allocate: integer iw'
	endif
	allocate(col(maxnar), stat=checkstat)
	if(checkstat > 0)then
	   write(6,*)'error with allocate:  integer iw'
	endif
           allocate(cbst(dall+maxvp),wt(dall+maxvp),dtres(dall+maxvp),&
        stat=checkstat)


! ESTIMATE INITIAL MODEL USING RAYLEIGH WAVE PHASE VELOCITY
            istep = 0
            dRp = 0
            do k = 1,kmax
            do j = 1,nsrc1(knum1(k))
            do i = 1,nrc1(j,knum1(k))
            istep = istep + 1
            if(wavetype(j,knum1(k))==2.and.&
               igrt(j,knum1(k))== 0) then
            dRp = dRp + 1
            pvRp(dRp) = pvall(istep)*1.1
            depRp(dRp) = pvall(istep)*tRc(periods(j,knum1(k)))/3.0
            endif
            enddo
            enddo
            enddo
            write(*,'(a,i7)') 'Number of all measurements',istep
            write(*,'(a,i7)') 'Number of Rayleigh wave phase velocity',dRp
            write(66,'(a,i7)') 'Number of measurements',istep
            write(66,'(a,i7)') 'Number of Rayleigh wave phase velocity',dRp


            call EstimateIniVel(vsf,depz,nx,ny,nz,dRp,depRp,pvRp,&
                 gsini,minv,grad,gridsp)
            deallocate(pvall)
            deallocate(pvRp,depRp)

! CHECKERBOARD TEST
            if (ifsyn == 1) then
            write(*,*) 'Checkerboard Resolution Test Begin'
            vsftrue = vsf
            call BuildCheckerboard(vsftrue,nx,ny,nz,&
                 sizex,sizey,sizez,anomaly)
            call synthetic(nx,ny,nz,maxvp,vsftrue,obst,&
                goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                scxf,sczf,rcxf,rczf,nrc1,nsrc1,knum1,kmax,&
                nsrc,nrc,noiselevel)
            endif



! ITERATE UNTILL CONVERGE
            do iter = 1,maxiter
            iw = 0
            rw = 0.0
            col = 0

! COMPUTE SENSITIVITY MATRIX
            call CalSurfG(nx,ny,nz,maxvp,vsf,iw,rw,col,cbst,obst,&
                goxd,gozd,dvxd,dvzd,kmaxRc,kmaxRg,kmaxLc,kmaxLg,&
                tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk,&
                scxf,sczf,rcxf,rczf,nrc1,nsrc1,knum1,kmax,&
                nsrc,nrc,nar,&
                maxlevel,maxleveld,HorizonType,VerticalType)
            
! ADDING REGULARIZATION TERM
        weight = 0.0
        if(iter.le.maxiter/3) then
	weight=weight1
	damp=damp1
	elseif(iter.gt.maxiter/3.and.iter.lt.2*maxiter/3) then
	weight=weight2
	damp=damp2
	else
	weight=weight3
	damp=damp3
	endif

         do i=1,maxvp
         rw(nar+i)=weight
         iw(1+nar+i)=dall+i
         col(nar+i)=i
         cbst(dall+i)=0
         enddo
         nar = nar + maxvp
         iw(1)=nar
         do i=1,nar
         iw(1+nar+i)=col(i)
         enddo
         if (nar > maxnar) stop 'increase sparsity fraction(spfra)'





! CALLING IRLS TO SOLVE THE PROBLEM

          m = dall + maxvp
          n = maxvp
          leniw = 2*nar+1
          lenrw = nar
          dv = 0
          atol = 1e-5
          btol = 1e-5
          conlim = 100
          itnlim = 1000
          istop = 0
          anorm = 0.0
          acond = 0.0
          arnorm = 0.0
          xnorm = 0.0
          

        call LSMR(m, n, leniw, lenrw,iw,rw,cbst, damp,&
           atol, btol, conlim, itnlim, localSize, nout,&
           dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)


        do iiter = 1, initer

        dtres=-cbst
        call aprod(1,m,n,dv,dtres,leniw,lenrw,iw,rw)
        do i=1,m
        if(abs(dtres(i)).lt.tolr) then
        wt(i)= 1.0/sqrt(abs(tolr))
        else
        wt(i)=1.0/sqrt(abs(dtres(i)))
        endif
        enddo
        do i=1,nar
        rw(i)=rw(i)*wt(iw(i+1))
        enddo
        do i=1,m
        dtres(i)=cbst(i)*wt(i)
        enddo

!          m = istep + maxvp
!          n = maxvp
!          leniw = 2*nar+1
!          lenrw = nar
          dv = 0
          atol = 1e-5
          btol = 1e-5
          conlim = 100
          itnlim = 1000
          istop = 0
          anorm = 0.0
          acond = 0.0
          arnorm = 0.0
          xnorm = 0.0
          

        call LSMR(m, n, leniw, lenrw,iw,rw,cbst, damp,&
           atol, btol, conlim, itnlim, localSize, nout,&
           dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)

       do i=1,nar
       rw(i)=rw(i)/wt(iw(i+1))
       enddo

       enddo ! finish inter interations for IRLS
       mean = sum(cbst(1:dall))/dall
       std_devs = sqrt(sum(cbst(1:dall)**2)/dall - mean**2)
       write(*,'(i2,a)'),iter,'th iteration...'
       write(*,'(a,f4.1)'),'weight is:',weight
       write(*,'(a,f5.2,a,f5.2,a)'),'mean and std_devs of residual: ',mean,'s ',std_devs,'s'
       write(66,'(i2,a)'),iter,'th iteration...'
       write(66,'(a,f4.1)'),'weight is:',weight
       write(66,'(a,f5.2,a,f5.2,a)'),'mean and std_devs of residual: ',mean,'s ',std_devs,'s'

	call invwavetrans(nx-2,ny-2,nz,dv,maxlevel,maxleveld,HorizonType,VerticalType)


        do k=1,nz
        do j=1,nx-2
        do i=1,ny-2
        if(dv((k-1)*(nx-2)*(ny-2)+(j-1)*(ny-2)+i).ge.0.5) then
        dv((k-1)*(nx-2)*(ny-2)+(j-1)*(ny-2)+i)=0.5
        endif
        if(dv((k-1)*(nx-2)*(ny-2)+(j-1)*(ny-2)+i).le.-0.5) then
        dv((k-1)*(nx-2)*(ny-2)+(j-1)*(ny-2)+i)=-0.5
        endif
        vsf(i+1,j+1,k)=vsf(i+1,j+1,k)-dv((k-1)*(nx-2)*(ny-2)+(j-1)*(ny-2)+i)
       if(vsf(i+1,j+1,k).lt.Minvel) vsf(i+1,j+1,k)=Minvel
       if(vsf(i+1,j+1,k).gt.Maxvel) vsf(i+1,j+1,k)=Maxvel
!vpf(i+1,j+1,k)=0.9409 + 2.0947*vsf(i+1,j+1,k) - 0.8206*vsf(i+1,j+1,k)**2+ &
!         0.2683*vsf(i+1,j+1,k)**3 - 0.0251*vsf(i+1,j+1,k)**4
!        vpf(i+1,j+1,k)=vsf(i+1,j+1,k)*1.732
        enddo
        enddo
        enddo
!        write(34,*)',OUTPUT P VELOCITY AT ITERATION',iter
!        do k=1,nzf
!        do j=1,nyf
!        write(34,'(100f7.3)') (vpf(i,j,k),i=1,nxf)
!        enddo
!        enddo
        write(34,*)',OUTPUT S VELOCITY AT ITERATION',iter
        do k=1,nz
        do j=1,ny
        write(34,'(100f7.3)') (vsf(i,j,k),i=1,nx)
        enddo
        enddo

            enddo !end iteration

! OUTPUT THE VELOCITY MODEL
                
        write(*,*),'Program finishes successfully'
        write(66,*),'Program finishes successfully'

        if(ifsyn == 1) then
        open(65,file='Vs_model.real')
        open(63,file='Vs_modelSyn.dat')
        do k=1,nz
        do j=1,ny
        write(65,'(100f7.3)') (vsftrue(i,j,k),i=1,nx)
        write(63,'(100f7.3)') (vsf(i,j,k),i=1,nx)
        enddo
        enddo
        close(65)
        close(63)
        write(*,*),'Output True velocity model &
                    to Vs_model.real'
        write(*,*),'Output inverted shear velocity model &
                    to Vs_modelSyn.dat'
        write(66,*),'Output True velocity model &
                    to Vs_model.real'
        write(66,*),'Output inverted shear velocity model &
                    to Vs_modelSyn.dat'
        else
        open(64,file='Vs_modelMeasure.dat')
        do k=1,nz
        do j=1,ny
        write(64,'(100f7.3)') (vsf(i,j,k),i=1,nx)
        enddo
        enddo
        close(64)
        write(*,*),'Output inverted shear velocity model &
                    to Vs_modelMeasure.dat'
        write(66,*),'Output inverted shear velocity model &
                    to Vs_modelMeasure.dat'
        endif

        close(34)
        close(nout) !close lsmr.txt
        close(66) !close surf_tomo.log
        deallocate(obst)
        deallocate(depz)
        deallocate(scxf,sczf)
        deallocate(rcxf,rczf)
        deallocate(wavetype,igrt,nrc1)
        deallocate(nsrc1,knum1,periods)
        deallocate(rw)
        deallocate(iw,col)
        deallocate(cbst,wt,dtres)
        deallocate(dv)
        deallocate(vsf)
        deallocate(vsftrue)
        if(kmaxRc.gt.0) then
	deallocate(tRc)
	endif
	if(kmaxRg.gt.0) then
	deallocate(tRg)
	endif
	if(kmaxLc.gt.0) then
	deallocate(tLc)
	endif
	if(kmaxLg.gt.0) then
	deallocate(tLg)
	endif

            end program            

