        subroutine BuildCheckerboard(vel,nx,ny,nz,sizex,sizey,sizez,anomaly)
            implicit none
            integer nx,ny,nz
            real vel(nx,ny,nz)
            integer sizex,sizey,sizez
            real anomaly

            integer i,j,k
            integer sz,sy,sx
            real anomalyx,anomalyy,anomalyz

            sz = 0
            anomalyz = anomaly

            do k = 1,nz
            sz = sz + 1
            if(sz == sizez) then
                anomalyz = -anomalyz
                sz = 0
            endif
            
            sy = 0
            anomalyy = anomalyz
            do j = 1,ny
            sy = sy + 1
            if (sy == sizey) then
                anomalyy = -anomalyy
                sy = 0
            endif

            sx = 0
            anomalyx = anomalyy
            do i = 1,nx
            sx = sx + 1
            vel(i,j,k) = vel(i,j,k)*(1+anomalyx)
            if (sx == sizex) then
                anomalyx = -anomalyx
                sx = 0
            endif
            enddo
            enddo
            enddo

            end subroutine
        
