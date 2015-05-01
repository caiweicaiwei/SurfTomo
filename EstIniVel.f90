       subroutine EstimateIniVel(vel,depz,nx,ny,nz,nd,dep,pv,&
            gsini,minv,grad,gridsp) 
            implicit none
            integer nd
            real dep(nd)
            real pv(nd)
            integer nx,ny,nz
            real vel(nx,ny,nz)
            real depz(nz)
            integer i
            real v1d(nz)
            integer cn(nz)
            real interval
            real fa,fb,fc
            real ispline
            external ispline
            integer gsini
            real minv,grad,gridsp
           
            vel = 0
            v1d = 0
            cn = 0
            if (gsini == 0) then
            interval = nint(maxval(dep)/nz*10)/10.0
!            write(*,*) 'grid spacing in depth direction:(km)'
!            write(*,'(f5.1)') interval
            do i =1, nd
            if(dep(i)>interval/2) then
            v1d(ceiling((dep(i)-interval/2)/interval)+1)=&
                 v1d(ceiling((dep(i)-interval/2)/interval)+1)+pv(i)
            cn(ceiling((dep(i)-interval/2)/interval)+1)=&
                 cn(ceiling((dep(i)-interval/2)/interval)+1)+1
            endif
            enddo
            do i = 1,nz
            depz(i)=(i-1)*interval
            v1d(i)=v1d(i)/cn(i)
            enddo
            write(*,*) 'grid points in depth direction:(km)'
            write(*,'(50f5.1)') depz

!            call quadfit(nz-2,depz(2:nz-1),v1d(2:nz-1),fa,fb,fc) 
!            do i = 1,nz
!            vel(:,:,i) = fa*depz(i)**2+fb*depz(i)+fc
!            enddo

            call spline(depz(2:nz-1),v1d(2:nz-1),fa,fb,fc,nz-2)
            do i = 1,nz
            vel(:,:,i) = ispline(depz(i),depz(2:nz-1),v1d(2:nz-1),&
                         fa,fb,fc,nz-2)
!            print*, ispline(depz(i),depz(2:nz-1),v1d(2:nz-1),&
!                         fa,fb,fc,nz-2)
            enddo
        else
            do i = 1,nz
            vel(:,:,i) = minv + (i-1)*grad
            depz(i) = (i-1)*gridsp
            enddo
        endif

            end subroutine
