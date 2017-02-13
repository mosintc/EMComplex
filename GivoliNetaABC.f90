! ==== PROBLEM CLASS ====
! Givoli-Neta ABC methods

!initGivoliNetaAuxVariables---------------------------
    subroutine initGivoliNetaAuxVariables(this)
    ! Allocate arrays for Givoli-Neta ABC aux variables
      class(tProblem) :: this;
      integer :: i;
      allocate(this%GNx(0:this%gnNum));
      allocate(this%GNy(0:this%gnNum));
      allocate(this%GNz(0:this%gnNum));
      allocate(this%GNxOld(0:this%gnNum));
      allocate(this%GNyOld(0:this%gnNum));
      allocate(this%GNzOld(0:this%gnNum));
      allocate(this%GNxOld2(0:this%gnNum));
      allocate(this%GNyOld2(0:this%gnNum));
      allocate(this%GNzOld2(0:this%gnNum));
      do i=0,this%gnNum
         call this%GNx(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
         call this%GNy(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
         call this%GNz(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
         call this%GNxOld(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
         call this%GNyOld(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
         call this%GNzOld(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
         call this%GNxOld2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
         call this%GNyOld2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
         call this%GNzOld2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
      end do
      allocate(this%gnA(1:this%gnNum));
      ! Set the angles for Givoli-Neta
      this%gnA(1) = 0.0d0;
      if (this%gnNum>1) then
         do i=2,this%gnNum
            this%gnA(i) = (i-1)*90.d0/this%gnNum;
         enddo
      endif
      allocate(this%GNs(0:this%GNnum));
      allocate(this%GNalp(1:this%GNnum));
      allocate(this%GNbet(0:this%GNnum));
      ! Set the coefficients for Givoli-Neta
      this%GNs(0) = cc;
      do i=1,this%gnNum
         this%GNs(i) = cc / cos(this%gnA(i)/180.0d0*PI);
      enddo
      ! Fill coefficients for Givoli-Neta
      this%GNbet(0) = 1.0d0 / this%GNs(1);
      if (this%gnNum > 1) then
         do i=1,this%gnNum-1
            this%GNalp(i) = 1.0d0/(this%GNs(i)**2) - 1.0d0/(this%GNs(0)**2);
            this%GNbet(i) = 1.0d0/this%GNs(i)+1.0d0/this%GNs(i+1);
         enddo
      endif
    end subroutine initGivoliNetaAuxVariables


!destroyGivoliNetaAuxVariables---------------------------
    subroutine destroyGivoliNetaAuxVariables(this)
    ! Destroy arrays for Givoli-Neta ABC aux variables
      class(tProblem) :: this;
      integer :: i;    
      do i=0,this%gnNum
         call this%GNx(i)%auxmesh_destroy;
         call this%GNy(i)%auxmesh_destroy;
         call this%GNz(i)%auxmesh_destroy;
         call this%GNxOld(i)%auxmesh_destroy;
         call this%GNyOld(i)%auxmesh_destroy;
         call this%GNzOld(i)%auxmesh_destroy;
         call this%GNxOld2(i)%auxmesh_destroy;
         call this%GNyOld2(i)%auxmesh_destroy;
         call this%GNzOld2(i)%auxmesh_destroy;
      end do
      deallocate(this%GNx, this%GNy, this%GNz, this%GNxOld, this%GNyOld, this%GNzOld, this%GNxOld2, this%GNyOld2, this%GNzOld2);
      deallocate(this%gnA);
      deallocate(this%GNs);
      deallocate(this%GNalp);
      deallocate(this%GNbet);
    end subroutine destroyGivoliNetaAuxVariables


!saveGivoliNetaAuxVariables---------------------------
    subroutine saveGivoliNetaAuxVariables(this)
    ! Save arrays for Givoli-Neta ABC aux variables
      class(tProblem) :: this;
      integer :: i, j, k, v;
      ! GNx field filling
      do v=0,this%gnNum
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%GNxOld2(v)%Y(i,j,k)=this%GNxOld(v)%Y(i,j,k);
               this%GNxOld(v)%Y(i,j,k)=this%GNx(v)%Y(i,j,k);
               this%GNxOld2(v)%Z(i,j,k)=this%GNxOld(v)%Z(i,j,k);
               this%GNxOld(v)%Z(i,j,k)=this%GNx(v)%Z(i,j,k);
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! GNy field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%GNyOld2(v)%X(i,j,k)=this%GNyOld(v)%X(i,j,k);
               this%GNyOld(v)%X(i,j,k)=this%GNy(v)%X(i,j,k);
               this%GNyOld2(v)%Z(i,j,k)=this%GNyOld(v)%Z(i,j,k);
               this%GNyOld(v)%Z(i,j,k)=this%GNy(v)%Z(i,j,k);
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! GNz field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%GNzOld2(v)%X(i,j,k)=this%GNzOld(v)%X(i,j,k);
               this%GNzOld(v)%X(i,j,k)=this%GNz(v)%X(i,j,k);
               this%GNzOld2(v)%Y(i,j,k)=this%GNzOld(v)%Y(i,j,k);
               this%GNzOld(v)%Y(i,j,k)=this%GNz(v)%Y(i,j,k);
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      end do
    end subroutine saveGivoliNetaAuxVariables


!fillEBoundaryABCGivoliNeta------------------------------------------------------------------
  subroutine  fillEBoundaryABCGivoliNeta(this,t)
    !Fill E boundary with ABC Givoli-Neta 2nd order Ey_x=0, exact elsewhere
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k, v, it  
    !================ EX =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = this%Source%getEpoint(this%xi05(i), this%yi(0),  this%zi(k), this%Ti(t), i, 0, k, 1);
          this%Ef%X(i, this%Ny, k) = this%Source%getEpoint(this%xi05(i), this%yi(this%Ny),  this%zi(k), this%Ti(t), i, this%Ny, k, 1);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = this%Source%getEpoint(this%xi05(i), this%yi(j),  this%zi(0), this%Ti(t), i, j, 0, 1);
          this%Ef%X(i, j, this%Nz) = this%Source%getEpoint(this%xi05(i), this%yi(j),  this%zi(this%Nz), this%Ti(t), i, j, this%Nz, 1);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !================ EY =================
    ! Fill EY Boundaries
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=1,this%Nx
          this%Ef%Y(i, j, 0) = this%Source%getEpoint(this%xi(i), this%yi05(j),  this%zi(0), this%Ti(t), i, j, 0, 2);
          this%Ef%Y(i, j, this%Nz) = this%Source%getEpoint(this%xi(i), this%yi05(j),  this%zi(this%Nz), this%Ti(t), i, j, this%Nz, 2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
!          this%Ef%Y(0, j, k) = this%Source%getEpoint(this%xi(0), this%yi05(j),  this%zi(k), this%Ti(t), 2);
          this%Ef%Y(this%Nx, j, k) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(j),  this%zi(k), this%Ti(t), this%Nx, j, k, 2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !!!*****************E and Gn0 renew*****************
    call this%GivoliNetaEymain(0);
    call this%GivoliNetaEyaux(0);  
    call this%GivoliNetaEymain(1);
    call this%GivoliNetaEyaux(0);
    !!!****************Gn j renew*****************
    !================ EZ =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)   
    do k=1,this%Nz
       do j=0,this%Ny
          this%Ef%Z(0, j, k) = this%Source%getEpoint(this%xi(0), this%yi(j),  this%zi05(k), this%Ti(t), 0, j, k, 3);
          this%Ef%Z(this%Nx, j, k) = this%Source%getEpoint(this%xi(this%Nx), this%yi(j),  this%zi05(k), this%Ti(t), this%Nx, j, k, 3);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = this%Source%getEpoint(this%xi(i), this%yi(0),  this%zi05(k), this%Ti(t), i, 0, k, 3);
          this%Ef%Z(i, this%Ny, k) = this%Source%getEpoint(this%xi(i), this%yi(this%Ny),  this%zi05(k), this%Ti(t), i, this%Ny, k, 3);
       enddo
    enddo
    !$OMP END PARALLEL DO      
  end subroutine fillEBoundaryABCGivoliNeta


!fillEBoundaryABCGivoliNetaPML------------------------------------------------------------------
  subroutine  fillEBoundaryABCGivoliNetaPML(this,t)
    !Fill E boundary with ABC Givoli-Neta High-order 2th order on Ey_x=0, Unsplit PML elsewhere
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k, v, it  
    !================ EX =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = this%PML%Ef%X(i+this%PML%NumPml, 0+this%PML%NumPml, k+this%PML%NumPml);
          this%Ef%X(i, this%Ny, k) = this%PML%Ef%X(i+this%PML%NumPml, this%Ny+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = this%PML%Ef%X(i+this%PML%NumPml, j+this%PML%NumPml, 0+this%PML%NumPml);
          this%Ef%X(i, j, this%Nz) = this%PML%Ef%X(i+this%PML%NumPml, j+this%PML%NumPml, this%Nz+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !================ EY =================
    ! Fill EY Boundaries
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=1,this%Nx
          this%Ef%Y(i, j, 0) = this%PML%Ef%Y(i+this%PML%NumPml, j+this%PML%NumPml, 0+this%PML%NumPml);
          this%Ef%Y(i, j, this%Nz) = this%PML%Ef%Y(i+this%PML%NumPml, j+this%PML%NumPml, this%Nz+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
!          this%Ef%Y(0, j, k) = this%PML%Ef%Y(0+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
          this%Ef%Y(this%Nx, j, k) = this%PML%Ef%Y(this%Nx+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !!!*****************E and Gn0 renew*****************
    call this%GivoliNetaEymain(0);
    call this%GivoliNetaEyaux(0);  
    call this%GivoliNetaEymain(1);
    call this%GivoliNetaEyaux(0);
    !!!****************Gn j renew*****************
    !================ EZ =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)   
    do k=1,this%Nz
       do j=0,this%Ny
          this%Ef%Z(0, j, k) = this%PML%Ef%Z(0+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
          this%Ef%Z(this%Nx, j, k) = this%PML%Ef%Z(this%Nx+this%PML%NumPml, j+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = this%PML%Ef%Z(i+this%PML%NumPml, 0+this%PML%NumPml, k+this%PML%NumPml);
         this%Ef%Z(i, this%Ny, k) = this%PML%Ef%Z(i+this%PML%NumPml, this%Ny+this%PML%NumPml, k+this%PML%NumPml);
       enddo
    enddo
    !$OMP END PARALLEL DO      
  end subroutine fillEBoundaryABCGivoliNetaPML
  

!fillEBoundaryABCGivoliNetaPure------------------------------------------------------------------
  subroutine  fillEBoundaryABCGivoliNetaPure(this,t)
    !Fill E boundary with ABC Givoli-Neta 2nd order in time on Ey_x=0, big aux problem elsewhere
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k, v, it  
    !================ EX =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = this%buffer%Ef%X(i+this%bufoffsetx, 0+this%bufoffsety, k+this%bufoffsetz);
          this%Ef%X(i, this%Ny, k) = this%buffer%Ef%X(i+this%bufoffsetx, this%Ny+this%bufoffsety, k+this%bufoffsetz);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = this%buffer%Ef%X(i+this%bufoffsetx, j+this%bufoffsety, 0+this%bufoffsetz);
          this%Ef%X(i, j, this%Nz) = this%buffer%Ef%X(i+this%bufoffsetx, j+this%bufoffsety, this%Nz+this%bufoffsetz);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !================ EY =================
    ! Fill EY Boundaries
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=1,this%Nx
          this%Ef%Y(i, j, 0) = this%buffer%Ef%Y(i+this%bufoffsetx, j+this%bufoffsety, 0+this%bufoffsetz);
          this%Ef%Y(i, j, this%Nz) = this%buffer%Ef%Y(i+this%bufoffsetx, j+this%bufoffsety, this%Nz+this%bufoffsetz);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=0,this%Nz
       do j=1,this%Ny
!          this%Ef%Y(0, j, k) = this%buffer%Ef%Y(0+this%bufoffsetx, j+this%bufoffsety, k+this%bufoffsetz);
          this%Ef%Y(this%Nx, j, k) = this%buffer%Ef%Y(this%Nx+this%bufoffsetx, j+this%bufoffsety, k+this%bufoffsetz);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !!!*****************E and Gn0 renew*****************
    call this%GivoliNetaEymain(0);
    call this%GivoliNetaEyaux(0);  
    call this%GivoliNetaEymain(1);
    call this%GivoliNetaEyaux(0);
    !!!****************Gn j renew*****************
    !================ EZ =================
    !$OMP PARALLEL DO SCHEDULE(GUIDED)   
    do k=1,this%Nz
       do j=0,this%Ny
          this%Ef%Z(0, j, k) = this%buffer%Ef%Z(0+this%bufoffsetx, j+this%bufoffsety,k+this%bufoffsetz);
          this%Ef%Z(this%Nx, j, k) = this%buffer%Ef%Z(this%Nx+this%bufoffsetx, j+this%bufoffsety,k+this%bufoffsetz);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = this%buffer%Ef%Z(i+this%bufoffsetx, 0+this%bufoffsety,k+this%bufoffsetz);
         this%Ef%Z(i, this%Ny, k) = this%buffer%Ef%Z(i+this%bufoffsetx, this%Ny+this%bufoffsety,k+this%bufoffsetz);
       enddo
    enddo
    !$OMP END PARALLEL DO      
  end subroutine fillEBoundaryABCGivoliNetaPure
  

!GivoliNetaEYmain------------------------------------------------------------------
  subroutine GivoliNetaEymain(this, mode)
    ! Givoli-Neta 2nd order in time, main update of the field
    class(tProblem) :: this
    integer, intent(in) :: mode
    real*8 :: de, dte, dt2f, dtf, df1, df2
    integer :: i ,j ,k
    !X Faces
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          !EYFace X=0 -----
          select case(this%dxmode)
          case(1)
             de = - (-this%EfOld%Y(0, j, k) + this%EfOld%Y(1, j, k))/(this%hhx);
          case(2)
             de = - (-3.0d0*this%EfOld%Y(0, j, k) + 4.0d0*this%EfOld%Y(1, j, k) - this%EfOld%Y(2, j, k) )/(2.0d0*this%hhx);     ! forward one-sided in X
          end select
          select case(mode)
          case(0)
             dte = this%GNs(1)*(this%gnYold(1)%X(0,j,k) - de);
          case(1)
             dte = this%GNs(1)*(this%gnY(1)%X(0,j,k) - de);
          end select
          select case(this%dtmode)
          case(1)
             this%Ef%Y(0, j, k) = (this%EfOld%Y(0, j, k) + this%ht*dte);               ! in time always backward here
          case(2)
             this%Ef%Y(0, j, k) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, j, k)-this%EfOld2%Y(0, j, k)+2.0d0*this%ht*dte);               ! in time always backward here
          end select
          this%gnY(0)%X(0, j, k) = this%Ef%Y(0, j, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !PLANES
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=2,this%Ny-1
       !{X=0} z=0
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, j, 0) + this%EfOld%Y(1, j, 0))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, j, 0) + 4.0d0*this%EfOld%Y(1, j, 0) - this%EfOld%Y(2, j, 0) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
       case(0)
          dte = this%GNs(1)*(this%gnYold(1)%X(0, j, 0) - de);
       case(1)
          dte = this%GNs(1)*(this%gnY(1)%X(0, j, 0) - de);
       end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, j, 0) = (this%EfOld%Y(0, j, 0) + this%ht*dte);     ! in time always backward here
       case(2)
          this%Ef%Y(0, j, 0) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, j, 0)-this%EfOld2%Y(0, j, 0)+2.0d0*this%ht*dte);     ! in time always backward here
       end select  
       this%gnY(0)%X(0, j, 0) = this%Ef%Y(0, j, 0);
       
       ! {X=0} z=Nz
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, j, this%Nz) + this%EfOld%Y(1, j, this%Nz))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, j, this%Nz) + 4.0d0*this%EfOld%Y(1, j, this%Nz) - this%EfOld%Y(2, j, this%Nz) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
       case(0)
          dte = this%GNs(1)*(this%gnYold(1)%X(0, j, this%Nz) - de);
       case(1)
          dte = this%GNs(1)*(this%gnY(1)%X(0, j, this%Nz) - de);
       end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, j, this%Nz) = (this%EfOld%Y(0, j, this%Nz) + this%ht*dte);               ! in time always backward here
       case(2)
          this%Ef%Y(0, j, this%Nz) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, j, this%Nz)-this%EfOld2%Y(0, j, this%Nz)+2.0d0*this%ht*dte);               ! in time always backward here
       end select  
       this%gnY(0)%X(0, j, this%Nz) = this%Ef%Y(0, j, this%Nz);
    enddo
    !$OMP END PARALLEL DO
    
    !EDGES -------
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       ! {X=0} y=1
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, 1, k) + this%EfOld%Y(1, 1, k))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, 1, k) + 4.0d0*this%EfOld%Y(1, 1, k) - this%EfOld%Y(2, 1, k) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
       case(0)
          dte = this%GNs(1)*(this%gnYold(1)%X(0, 1, k) - de);
       case(1)
          dte = this%GNs(1)*(this%gnY(1)%X(0, 1, k) - de);
       end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, 1, k) = (this%EfOld%Y(0, 1, k) + this%ht*dte);
       case(2)
          this%Ef%Y(0, 1, k) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, 1, k)-this%EfOld2%Y(0, 1, k)+2.0d0*this%ht*dte);               ! in time always backward here
       end select   
       this%gnY(0)%X(0, 1, k) = this%Ef%Y(0, 1, k);
       
       ! {X=0} y=Ny
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, this%Ny, k) + this%EfOld%Y(1, this%Ny, k))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, this%Ny, k) + 4.0d0*this%EfOld%Y(1, this%Ny, k) - this%EfOld%Y(2, this%Ny, k) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
       case(0)
          dte = this%GNs(1)*(this%gnYold(1)%X(0, this%Ny, k) - de);
       case(1)
          dte = this%GNs(1)*(this%gnY(1)%X(0, this%Ny, k) - de);
       end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, this%Ny, k) = (this%EfOld%Y(0, this%Ny, k) + this%ht*dte);  
       case(2)
          this%Ef%Y(0, this%Ny, k) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, this%Ny, k)-this%EfOld2%Y(0, this%Ny, k)+2.0d0*this%ht*dte);                          ! in time always backward here
       end select   
       this%gnY(0)%X(0, this%Ny, k) = this%Ef%Y(0, this%Ny, k);
    enddo
    !$OMP END PARALLEL DO
    
    !CORNERS ----------
    !{X=0} y=1 z=0
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, 1, 0) + this%EfOld%Y(1, 1, 0))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, 1, 0) + 4.0d0*this%EfOld%Y(1, 1, 0) - this%EfOld%Y(2, 1, 0) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
          case(0)
             dte = this%GNs(1)*(this%gnYold(1)%X(0, 1, 0) - de);
          case(1)
             dte = this%GNs(1)*(this%gnY(1)%X(0, 1, 0) - de);
          end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, 1, 0) = (this%EfOld%Y(0, 1, 0) + this%ht*dte);
       case(2)
          this%Ef%Y(0, 1, 0) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, 1, 0)-this%EfOld2%Y(0, 1, 0)+2.0d0*this%ht*dte);               ! in time always backward here
       end select
       this%gnY(0)%X(0, 1, 0) = this%Ef%Y(0, 1, 0);

    !{x=0} y=Ny Z=0
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, this%Ny, 0) + this%EfOld%Y(1, this%Ny, 0))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, this%Ny, 0) + 4.0d0*this%EfOld%Y(1, this%Ny, 0) - this%EfOld%Y(2, this%Ny, 0) )/(2.0d0*this%hhx); ! forward one-sided in X  
       end select
       select case(mode)
          case(0)
             dte = this%GNs(1)*(this%gnYold(1)%X(0, this%Ny, 0) - de);
          case(1)
             dte = this%GNs(1)*(this%gnY(1)%X(0, this%Ny, 0) - de);
          end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, this%Ny, 0) = (this%EfOld%Y(0, this%Ny, 0) + this%ht*dte);
       case(2)
          this%Ef%Y(0, this%Ny, 0) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, this%Ny, 0)-this%EfOld2%Y(0, this%Ny, 0)+2.0d0*this%ht*dte);               ! in time always backward here
       end select
       this%gnY(0)%X(0, this%Ny, 0) = this%Ef%Y(0, this%Ny, 0);
       
    !{X=0} y=1 z=Nz
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, 1, this%Nz) + this%EfOld%Y(1, 1, this%Nz))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, 1, this%Nz) + 4.0d0*this%EfOld%Y(1, 1, this%Nz) - this%EfOld%Y(2, 1, this%Nz) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
          case(0)
             dte = this%GNs(1)*(this%gnYold(1)%X(0, 1, this%Nz) - de);
          case(1)
             dte = this%GNs(1)*(this%gnY(1)%X(0, 1, this%Nz) - de);
          end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, 1, this%Nz) = (this%EfOld%Y(0, 1, this%Nz) + this%ht*dte); 
       case(2)
          this%Ef%Y(0, 1, this%Nz) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, 1, this%Nz)-this%EfOld2%Y(0, 1, this%Nz)+2.0d0*this%ht*dte);               ! in time always backward here
       end select
       this%gnY(0)%X(0, 1, this%Nz) = this%Ef%Y(0, 1, this%Nz);
       
    !{X=0} y=Ny z=Nz
       select case(this%dxmode)
       case(1)
          de = - (-this%EfOld%Y(0, this%Ny, this%Nz) + this%EfOld%Y(1, this%Ny, this%Nz))/(this%hhx);
       case(2)
          de = - (-3.0d0*this%EfOld%Y(0, this%Ny, this%Nz) + 4.0d0*this%EfOld%Y(1, this%Ny, this%Nz) - this%EfOld%Y(2, this%Ny, this%Nz) )/(2.0d0*this%hhx); ! forward one-sided in X
       end select
       select case(mode)
          case(0)
             dte = this%GNs(1)*(this%gnYold(1)%X(0, this%Ny, this%Nz) - de);
          case(1)
             dte = this%GNs(1)*(this%gnY(1)%X(0, this%Ny, this%Nz) - de);
          end select
       select case(this%dtmode)
       case(1)
          this%Ef%Y(0, this%Ny, this%Nz) = (this%EfOld%Y(0, this%Ny, this%Nz) + this%ht*dte);
       case(2)
          this%Ef%Y(0, this%Ny, this%Nz) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, this%Ny, this%Nz)-this%EfOld2%Y(0, this%Ny, this%Nz)+2.0d0*this%ht*dte);            ! in time always backward here
       end select
       this%gnY(0)%X(0, this%Ny, this%Nz) = this%Ef%Y(0, this%Ny, this%Nz);       
   end subroutine GivoliNetaEymain


!GivoliNetaEYaux------------------------------------------------------------------
  subroutine GivoliNetaEYaux(this, mode)
    ! Givoli-Neta 2nd order in time, aux functions update
    class(tProblem) :: this
    integer :: mode;
    real*8 :: de, dte, dt2f, dtf, df1, df2
    integer :: i ,j ,k, v
    if (this%gnNum > 1) then
       do v=1,this%gnNum-1
          !=========================== EY == X BOUNDARIES ============================
          !Faces
          !$OMP PARALLEL DO SCHEDULE(GUIDED)
          do k=1,this%Nz-1
             do j=2,this%Ny-1
                !EYFace X=0 -----
                dt2f = (this%gnY(v-1)%X(0,j,k)-2.0d0*this%gnYold(v-1)%X(0,j,k)+this%gnYold2(v-1)%X(0,j,k))/(this%ht**2); ! central in time around Old
                select case(this%dxauxmode)
                case(1)
                   df1 = (this%gnY(v-1)%X(0,j+1,k)-2.0d0*this%gnY(v-1)%X(0,j,k)+this%gnY(v-1)%X(0,j-1,k))/(this%hhy**2);    ! central in Y
                   df2 = (this%gnY(v-1)%X(0,j,k+1)-2.0d0*this%gnY(v-1)%X(0,j,k)+this%gnY(v-1)%X(0,j,k-1))/(this%hhz**2);    ! central in Z 
                case(2)
                   df1 = (this%gnY(v-1)%X(0,j+1,k)-2.0d0*this%gnY(v-1)%X(0,j,k)+this%gnY(v-1)%X(0,j-1,k))/(this%hhy**2);    ! central in Y
                   df2 = (this%gnY(v-1)%X(0,j,k+1)-2.0d0*this%gnY(v-1)%X(0,j,k)+this%gnY(v-1)%X(0,j,k-1))/(this%hhz**2);    ! central in Z
                end select   
                select case(mode)
                case(0)
                   dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,j,k)+this%GNalp(v)*dt2f+df1+df2);
                case(1)
                   dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,j,k)+this%GNalp(v)*dt2f+df1+df2);
                end select
                select case(this%dtauxmode)
                case(1)
                   this%gnY(v)%X(0,j,k) = (this%gnYold(v)%X(0,j,k) + this%ht*dtf);  ! in time always backward here 
                case(2)
                   this%gnY(v)%X(0,j,k) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,j,k)-this%gnYold2(v)%X(0,j,k)+2.0d0*this%ht*dtf);  ! in time always backward here
                end select   
             enddo
          enddo
          !$OMP END PARALLEL DO
          !EYEdge XZ -------
          !$OMP PARALLEL DO SCHEDULE(GUIDED)
          do j=2,this%Ny-1
             !{X=0} z=0
             dt2f = (this%gnY(v-1)%X(0,j,0)-2.0d0*this%gnYold(v-1)%X(0,j,0)+this%gnYold2(v-1)%X(0,j,0))/(this%ht**2);                                     ! central in time around Old
             select case(this%dxauxmode)
             case(1)
                df1 = (this%gnY(v-1)%X(0,j+1,0)-2.0d0*this%gnY(v-1)%X(0,j,0)+this%gnY(v-1)%X(0,j-1,0))/(this%hhy**2);                                        ! central in Y
                df2 = (this%gnY(v-1)%X(0,j,0)-2.0d0*this%gnY(v-1)%X(0,j,1)+this%gnY(v-1)%X(0,j,2))/(this%hhz**2);                                            ! forward one-sided in Z    
             case(2)
                df1 = (this%gnY(v-1)%X(0,j+1,0)-2.0d0*this%gnY(v-1)%X(0,j,0)+this%gnY(v-1)%X(0,j-1,0))/(this%hhy**2);                                        ! central in Y
                df2 = (2.0d0*this%gnY(v-1)%X(0,j,0)-5.0d0*this%gnY(v-1)%X(0,j,1)+4.0d0*this%gnY(v-1)%X(0,j,2)-this%gnY(v-1)%X(0,j,3))/(this%hhz**2);         ! forward one-sided in Z
             end select  
             select case(mode)
             case(0)
                dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,j,0)+this%GNalp(v)*dt2f+df1+df2);
             case(1)
                dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,j,0)+this%GNalp(v)*dt2f+df1+df2);
             end select
             select case(this%dtauxmode)
             case(1)
                this%gnY(v)%X(0,j,0) = (this%gnYold(v)%X(0,j,0) + this%ht*dtf);  ! in time always backward here
             case(2)
                this%gnY(v)%X(0,j,0) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,j,0)-this%gnYold2(v)%X(0,j,0)+2.0d0*this%ht*dtf);  ! in time always backward here        
             end select   
            
             !{X=0} z=Nz
             dt2f = (this%gnY(v-1)%X(0,j,this%Nz)-2.0d0*this%gnYold(v-1)%X(0,j,this%Nz)+this%gnYold2(v-1)%X(0,j,this%Nz))/(this%ht**2);                                            ! central in time around Old
             select case(this%dxauxmode)
             case(1)
                df1 = (this%gnY(v-1)%X(0,j+1,this%Nz)-2.0d0*this%gnY(v-1)%X(0,j,this%Nz)+this%gnY(v-1)%X(0,j-1,this%Nz))/(this%hhy**2);                                               ! central in Y
                df2 = (this%gnY(v-1)%X(0,j,this%Nz)-2.0d0*this%gnY(v-1)%X(0,j,this%Nz-1)+this%gnY(v-1)%X(0,j,this%Nz-2))/(this%hhz**2);                                               ! backward one-sided Z 
             case(2)   
                df1 = (this%gnY(v-1)%X(0,j+1,this%Nz)-2.0d0*this%gnY(v-1)%X(0,j,this%Nz)+this%gnY(v-1)%X(0,j-1,this%Nz))/(this%hhy**2);                                               ! central in Y
                df2 = (2.0d0*this%gnY(v-1)%X(0,j,this%Nz)-5.0d0*this%gnY(v-1)%X(0,j,this%Nz-1)+4.0d0*this%gnY(v-1)%X(0,j,this%Nz-2)-this%gnY(v-1)%X(0,j,this%Nz-3))/(this%hhz**2);    ! backward one-sided Z
             end select   
             select case(mode)
             case(0)
                dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,j,this%Nz)+this%GNalp(v)*dt2f+df1+df2);
             case(1)
                dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,j,this%Nz)+this%GNalp(v)*dt2f+df1+df2);
             end select
             select case(this%dtauxmode)
             case(1)
                this%gnY(v)%X(0,j,this%Nz) = (this%gnYold(v)%X(0,j,this%Nz) + this%ht*dtf);                                                       ! in time always backward here
             case(2)
                this%gnY(v)%X(0,j,this%Nz) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,j,this%Nz)-this%gnYold2(v)%X(0,j,this%Nz)+2.0d0*this%ht*dtf);  ! in time always backward here
             end select   
          enddo
          !$OMP END PARALLEL DO
          
          !EYEdge XY -------
          !$OMP PARALLEL DO SCHEDULE(GUIDED)
          do k=1,this%Nz-1
             !X=0, Y=1
             dt2f = (this%gnY(v-1)%X(0,1,k)-2.0d0*this%gnYold(v-1)%X(0,1,k)+this%gnYold2(v-1)%X(0,1,k))/(this%ht**2);                                      ! central in time around Old
             select case(this%dxauxmode)
             case(1)
                df1 = (this%gnY(v-1)%X(0,1,k)-2.0d0*this%gnY(v-1)%X(0,2,k)+this%gnY(v-1)%X(0,3,k))/(this%hhy**2);                                             ! forward one-sided in Y
                df2 = (this%gnY(v-1)%X(0,1,k+1)-2.0d0*this%gnY(v-1)%X(0,1,k)+this%gnY(v-1)%X(0,1,k-1))/(this%hhz**2);                                         ! central in Z 
             case(2)
                df1 = (2.0d0*this%gnY(v-1)%X(0,1,k)-5.0d0*this%gnY(v-1)%X(0,2,k)+4.0d0*this%gnY(v-1)%X(0,3,k)-this%gnY(v-1)%X(0,4,k))/(this%hhy**2);          ! forward one-sided in Y
                df2 = (this%gnY(v-1)%X(0,1,k+1)-2.0d0*this%gnY(v-1)%X(0,1,k)+this%gnY(v-1)%X(0,1,k-1))/(this%hhz**2);                                         ! central in Z
             end select  
             select case(mode)
             case(0)
                dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,1,k)+this%GNalp(v)*dt2f+df1+df2);
             case(1)
                dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,1,k)+this%GNalp(v)*dt2f+df1+df2);
             end select
             select case(this%dtauxmode)
             case(1)
                this%gnY(v)%X(0,1,k) = (this%gnYold(v)%X(0,1,k) + this%ht*dtf);  ! in time always backward here
             case(2)   
                this%gnY(v)%X(0,1,k) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,1,k)-this%gnYold2(v)%X(0,1,k)+2.0d0*this%ht*dtf);  ! in time always backward here
             end select   

             !X=0, Y=Ny
             dt2f = (this%gnY(v-1)%X(0,this%Ny,k)-2.0d0*this%gnYold(v-1)%X(0,this%Ny,k)+this%gnYold2(v-1)%X(0,this%Ny,k))/(this%ht**2);                                            ! central in time around Old
             select case(this%dxauxmode)
             case(1)
                df1 = (this%gnY(v-1)%X(0,this%Ny,k)-2.0d0*this%gnY(v-1)%X(0,this%Ny-1,k)+this%gnY(v-1)%X(0,this%Ny-2,k))/(this%hhy**2);                                               ! backward one-sided in Y
                df2 = (this%gnY(v-1)%X(0,this%Ny,k+1)-2.0d0*this%gnY(v-1)%X(0,this%Ny,k)+this%gnY(v-1)%X(0,this%Ny,k-1))/(this%hhz**2);                                               ! central in Z 
             case(2)   
                df1 = (2.0d0*this%gnY(v-1)%X(0,this%Ny,k)-5.0d0*this%gnY(v-1)%X(0,this%Ny-1,k)+4.0d0*this%gnY(v-1)%X(0,this%Ny-2,k)-this%gnY(v-1)%X(0,this%Ny-3,k))/(this%hhy**2);    ! backward one-sided in Y
                df2 = (this%gnY(v-1)%X(0,this%Ny,k+1)-2.0d0*this%gnY(v-1)%X(0,this%Ny,k)+this%gnY(v-1)%X(0,this%Ny,k-1))/(this%hhz**2);                                               ! central in Z
             end select  
             select case(mode)
             case(0)
                dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,this%Ny,k)+this%GNalp(v)*dt2f+df1+df2);
             case(1)
                dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,this%Ny,k)+this%GNalp(v)*dt2f+df1+df2);
             end select
             select case(this%dtauxmode)
             case(1)
                this%gnY(v)%X(0,this%Ny,k) = (this%gnYold(v)%X(0,this%Ny,k) + this%ht*dtf);                                                       ! in time always backward here
             case(2)   
                this%gnY(v)%X(0,this%Ny,k) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,this%Ny,k)-this%gnYold2(v)%X(0,this%Ny,k)+2.0d0*this%ht*dtf);  ! in time always backward here
             end select   
          enddo
          
          !$OMP END PARALLEL DO  
          !EYCorners ----------
          !{X=0} y=1 z=0
          dt2f = (this%gnY(v-1)%X(0,1,0)-2.0d0*this%gnYold(v-1)%X(0,1,0)+this%gnYold2(v-1)%X(0,1,0))/(this%ht**2);                               ! central in time around Old
          select case(this%dxauxmode)
          case(1)
             df1 = (this%gnY(v-1)%X(0,1,0)-2.0d0*this%gnY(v-1)%X(0,2,0)+this%gnY(v-1)%X(0,3,0))/(this%hhy**2);                                      ! forward one-sided in Y
             df2 = (this%gnY(v-1)%X(0,1,0)-2.0d0*this%gnY(v-1)%X(0,1,1)+this%gnY(v-1)%X(0,1,2))/(this%hhz**2);                                      ! forward one-sided in Z 
          case(2)   
             df1 = (2.0d0*this%gnY(v-1)%X(0,1,0)-5.0d0*this%gnY(v-1)%X(0,2,0)+4.0d0*this%gnY(v-1)%X(0,3,0)-this%gnY(v-1)%X(0,4,0))/(this%hhy**2);   ! forward one-sided in Y
             df2 = (2.0d0*this%gnY(v-1)%X(0,1,0)-5.0d0*this%gnY(v-1)%X(0,1,1)+4.0d0*this%gnY(v-1)%X(0,1,2)-this%gnY(v-1)%X(0,1,3))/(this%hhz**2);   ! forward one-sided in Z
          end select    
          select case(mode)
          case(0)
             dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,1,0)+this%GNalp(v)*dt2f+df1+df2);
          case(1)
             dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,1,0)+this%GNalp(v)*dt2f+df1+df2);
          end select
          select case(this%dtauxmode)
          case(1)
             this%gnY(v)%X(0,1,0) = (this%gnYold(v)%X(0,1,0) + this%ht*dtf);                                                                   ! in time always backward here
          case(2)   
             this%gnY(v)%X(0,1,0) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,1,0)-this%gnYold2(v)%X(0,1,0)+2.0d0*this%ht*dtf);                    ! in time always backward here
          end select
          
          !{x=0} y=Ny Z=0
          dt2f = (this%gnY(v-1)%X(0,this%Ny,0)-2.0d0*this%gnYold(v-1)%X(0,this%Ny,0)+this%gnYold2(v-1)%X(0,this%Ny,0))/(this%ht**2);                                           ! central in time around Old
          select case(this%dxauxmode)
          case(1)
             df1 = (this%gnY(v-1)%X(0,this%Ny,0)-2.0d0*this%gnY(v-1)%X(0,this%Ny-1,0)+this%gnY(v-1)%X(0,this%Ny-2,0))/(this%hhy**2);                                              ! backward one-sided in Y
             df2 = (this%gnY(v-1)%X(0,this%Ny,0)-2.0d0*this%gnY(v-1)%X(0,this%Ny,1)+this%gnY(v-1)%X(0,this%Ny,2))/(this%hhz**2);                                                  ! forward one-sided in Z 
          case(2)   
             df1 = (2.0d0*this%gnY(v-1)%X(0,this%Ny,0)-5.0d0*this%gnY(v-1)%X(0,this%Ny-1,0)+4.0d0*this%gnY(v-1)%X(0,this%Ny-2,0)-this%gnY(v-1)%X(0,this%Ny-3,0))/(this%hhy**2);   ! backward one-sided in Y
             df2 = (2.0d0*this%gnY(v-1)%X(0,this%Ny,0)-5.0d0*this%gnY(v-1)%X(0,this%Ny,1)+4.0d0*this%gnY(v-1)%X(0,this%Ny,2)-this%gnY(v-1)%X(0,this%Ny,3))/(this%hhz**2);         ! forward one-sided in Z
          end select    
          select case(mode)
          case(0)
             dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,this%Ny,0)+this%GNalp(v)*dt2f+df1+df2);
          case(1)
             dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,this%Ny,0)+this%GNalp(v)*dt2f+df1+df2);
          end select
          select case(this%dtauxmode)
          case(1)
             this%gnY(v)%X(0,this%Ny,0) = (this%gnYold(v)%X(0,this%Ny,0) + this%ht*dtf);                                                       ! in time always backward here
          case(2)   
             this%gnY(v)%X(0,this%Ny,0) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,this%Ny,0)-this%gnYold2(v)%X(0,this%Ny,0)+2.0d0*this%ht*dtf);  ! in time always backward here
          end select   

          !{X=0} y=1 z=Nz
          dt2f = (this%gnY(v-1)%X(0,1,this%Nz)-2.0d0*this%gnYold(v-1)%X(0,1,this%Nz)+this%gnYold2(v-1)%X(0,1,this%Nz))/(this%ht**2);                                            ! central in time around Old
          select case(this%dxauxmode)
          case(1)
             df1 = (this%gnY(v-1)%X(0,1,this%Nz)-2.0d0*this%gnY(v-1)%X(0,2,this%Nz)+this%gnY(v-1)%X(0,3,this%Nz))/(this%hhy**2);                                                   ! forward one-sided in Y
             df2 = (this%gnY(v-1)%X(0,1,this%Nz)-2.0d0*this%gnY(v-1)%X(0,1,this%Nz-1)+this%gnY(v-1)%X(0,1,this%Nz-2))/(this%hhz**2);                                               ! backward one-sided in Z 
          case(2)   
             df1 = (2.0d0*this%gnY(v-1)%X(0,1,this%Nz)-5.0d0*this%gnY(v-1)%X(0,2,this%Nz)+4.0d0*this%gnY(v-1)%X(0,3,this%Nz)-this%gnY(v-1)%X(0,4,this%Nz))/(this%hhy**2);          ! forward one-sided in Y
             df2 = (2.0d0*this%gnY(v-1)%X(0,1,this%Nz)-5.0d0*this%gnY(v-1)%X(0,1,this%Nz-1)+4.0d0*this%gnY(v-1)%X(0,1,this%Nz-2)-this%gnY(v-1)%X(0,1,this%Nz-3))/(this%hhz**2);    ! backward one-sided in Z
          end select   
          select case(mode)
          case(0)
             dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,1,this%Nz)+this%GNalp(v)*dt2f+df1+df2);
          case(1)
             dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,1,this%Nz)+this%GNalp(v)*dt2f+df1+df2);
          end select
          select case(this%dtauxmode)
          case(1)
             this%gnY(v)%X(0,1,this%Nz) = (this%gnYold(v)%X(0,1,this%Nz) + this%ht*dtf);  ! in time always backward here
          case(2)   
             this%gnY(v)%X(0,1,this%Nz) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,1,this%Nz)-this%gnYold2(v)%X(0,1,this%Nz)+2.0d0*this%ht*dtf);  ! in time always backward here
          end select
          
          !{X=0} y=Ny z=Nz
          dt2f = (this%gnY(v-1)%X(0,this%Ny,this%Nz)-2.0d0*this%gnYold(v-1)%X(0,this%Ny,this%Nz)+this%gnYold2(v-1)%X(0,this%Ny,this%Nz))/(this%ht**2);                                                ! central in time around Old
          select case(this%dxauxmode)
          case(1)
             df1 = (this%gnY(v-1)%X(0,this%Ny,this%Nz)-2.0d0*this%gnY(v-1)%X(0,this%Ny-1,this%Nz)+this%gnY(v-1)%X(0,this%Ny-2,this%Nz))/(this%hhy**2);                                                   ! backward one-sided in Y
             df2 = (this%gnY(v-1)%X(0,this%Ny,this%Nz)-2.0d0*this%gnY(v-1)%X(0,this%Ny,this%Nz-1)+this%gnY(v-1)%X(0,this%Ny,this%Nz-2))/(this%hhz**2);                                                   ! backward one-sided in Z 
          case(2)   
             df1 = (2.0d0*this%gnY(v-1)%X(0,this%Ny,this%Nz)-5.0d0*this%gnY(v-1)%X(0,this%Ny-1,this%Nz)+4.0d0*this%gnY(v-1)%X(0,this%Ny-2,this%Nz)-this%gnY(v-1)%X(0,this%Ny-3,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
             df2 = (2.0d0*this%gnY(v-1)%X(0,this%Ny,this%Nz)-5.0d0*this%gnY(v-1)%X(0,this%Ny,this%Nz-1)+4.0d0*this%gnY(v-1)%X(0,this%Ny,this%Nz-2)-this%gnY(v-1)%X(0,this%Ny,this%Nz-3))/(this%hhz**2);  ! backward one-sided in Z
          end select   
          select case(mode)
          case(0)
             dtf = 1.0d0/this%GNbet(v)*(this%gnYOld(v+1)%X(0,this%Ny,this%Nz)+this%GNalp(v)*dt2f+df1+df2);
          case(1)
             dtf = 1.0d0/this%GNbet(v)*(this%gnY(v+1)%X(0,this%Ny,this%Nz)+this%GNalp(v)*dt2f+df1+df2);
          end select
          select case(this%dtauxmode)
          case(1)
             this%gnY(v)%X(0,this%Ny,this%Nz) = (this%gnYold(v)%X(0,this%Ny,this%Nz) + this%ht*dtf);  ! in time always backward here
          case(2)   
             this%gnY(v)%X(0,this%Ny,this%Nz) = 1.0d0/3.0d0*(4.0d0*this%gnYold(v)%X(0,this%Ny,this%Nz)-this%gnYold2(v)%X(0,this%Ny,this%Nz)+2.0d0*this%ht*dtf);  ! in time always backward here
          end select   
       enddo
    endif
  end subroutine GivoliNetaEYaux
