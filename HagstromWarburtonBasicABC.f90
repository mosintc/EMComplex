!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains part of implementation of tProblem class

! ==== PROBLEM CLASS ====
! Hagstrom-Warburton Basic ABC methods

!initHWAuxVariables---------------------------
    subroutine initHWAuxVariables(this)
    ! Allocate arrays for Hagstrom-Warburton ABC aux variables
      class(tProblem) :: this;
      integer :: i, j, info;
      real*8, allocatable :: ipiv(:), Ematr(:,:);
      allocate(this%hwX(0:this%hwNum));
      allocate(this%hwY(0:this%hwNum));
      allocate(this%hwZ(0:this%hwNum));
      allocate(this%hwXold(0:this%hwNum));
      allocate(this%hwYold(0:this%hwNum));
      allocate(this%hwZold(0:this%hwNum));
      allocate(this%hwXold2(0:this%hwNum));
      allocate(this%hwYold2(0:this%hwNum));
      allocate(this%hwZold2(0:this%hwNum));

      do i=0,this%hwNum
         call this%hwX(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
         call this%hwY(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
         call this%hwZ(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
         call this%hwXold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
         call this%hwYold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
         call this%hwZold(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
         call this%hwXold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 0, this%problem_id);
         call this%hwYold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 1, this%problem_id);
         call this%hwZold2(i)%auxmesh_init(this%Nx, this%Ny, this%Nz, 2, this%problem_id);
      end do

      ! Set angles for HW
      allocate(this%hwAng(0:this%hwNum-1));
      this%hwAng(0) = 0.0d0;
      if (this%hwNum>1) then
         do i=2,this%hwNum
            this%hwAng(i-1) = 0.0d0!(i-1)*90.0d0/this%hwNum;
         enddo
      endif

      ! Set basic coefficients for HW
      allocate(this%hwAC(0:this%hwNum-1));
      do i=0,this%hwNum-1
         this%hwAC(i) = cos(this%hwAng(i)/180.0d0*PI);
      enddo

      !Fill coefficients for HW
      allocate(this%hwL(0:this%hwNum-1, 0:this%hwNum-1));
      allocate(this%hwA(0:this%hwNum-1, 0:this%hwNum-1));
      allocate(this%hwB(0:this%hwNum-1, 0:this%hwNum-1));
      allocate(this%hwC(0:this%hwNum-1, 0:this%hwNum-1));
      allocate(this%hwD(0:this%hwNum-1, 0:this%hwNum-1));
      allocate(Ematr(0:this%hwNum-1, 0:this%hwNum-1));
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            this%hwL(i,j)=0.0d0;
            this%hwA(i,j)=0.0d0;
            this%hwB(i,j)=0.0d0;
            this%hwC(i,j)=0.0d0;
            this%hwD(i,j)=0.0d0;
            Ematr(i,j)=0.0d0;
            if (i==j) then
               Ematr(i,j)=1.0d0;
            endif
         enddo
      enddo
      
      if (this%hwNum>1) then
         this%hwL(1,0) = 2.0d0*this%hwAC(1)*(1.0d0-this%hwAC(0)**2);
         this%hwL(1,1) = 1.0d0+this%hwAC(1)**2+2.0d0*this%hwAC(0)*this%hwAC(1);
         this%hwD(1,0) = 2.0d0*this%hwAC(1)*(cc**2)*(this%ht**2);
         this%hwD(1,1) = 1.0d0*(cc**2)*(this%ht**2);        
         if (this%hwNum>2) then
            this%hwL(1,2) =  (1.0d0-this%hwAC(1)**2);
            this%hwD(1,2) = 1.0d0*(cc**2)*(this%ht**2);        
            do i=2,this%hwNum-1
               this%hwL(i,i-1) = this%hwAC(i)*(1.0d0-this%hwAC(i-1)**2);
               this%hwL(i,i) = this%hwAC(i)*(1.0d0+this%hwAC(i-1)**2)+this%hwAC(i-1)*(1.0d0+this%hwAC(i)**2);
               this%hwD(i,i-1) = this%hwAC(i)*(cc**2)*(this%ht**2);
               this%hwD(i,i) = (this%hwAC(i)+this%hwAC(i-1))*(cc**2)*(this%ht**2);
               if (i/=this%hwNum-1) then
                  this%hwL(i,i+1) = this%hwAC(i-1)*(1.0d0-this%hwAC(i)**2);
                  this%hwD(i,i+1) = this%hwAC(i-1)*(cc**2)*(this%ht**2);
               endif
            enddo
         endif     
      endif

      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            this%hwA(i,j)=this%hwL(i,j);
            this%hwB(i,j)=2.0d0*this%hwL(i,j);
            this%hwC(i,j)=this%hwL(i,j);
         enddo
      enddo

      

      select case(this%dtmode)
      case(1)
          this%hwA(0,0)=1.0d0*this%hwAC(0);
          this%hwB(0,0)=1.0d0*this%hwAC(0);
          this%hwC(0,0)=0.0d0;
          if (this%hwNum>1) then
             this%hwA(0,1)=-1.0d0;
             this%hwB(0,1)=-1.0d0;
             this%hwC(0,1)=0.0d0;
          endif
      case(2)
          this%hwA(0,0)=3.0d0*this%hwAC(0);
          this%hwB(0,0)=4.0d0*this%hwAC(0);
          this%hwC(0,0)=this%hwAC(0);
          if (this%hwNum>1) then
             this%hwA(0,1)=-3.0d0;
             this%hwB(0,1)=-4.0d0;
             this%hwC(0,1)=-1.0d0;
          endif
      end select

      write(600,*) 'A array';
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            WRITE(600, 521, advance="no") this%hwA(i,j)
         enddo
         write(600,*) '';
      enddo
      write(600,*) 'L array';
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            WRITE(600, 521, advance="no") this%hwL(i,j)
         enddo
         write(600,*) '';
      enddo     
      write(600,*) 'B array';
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            WRITE(600, 521, advance="no") this%hwB(i,j)
         enddo
         write(600,*) '';
      enddo
      write(600,*) 'C array';
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            WRITE(600, 521, advance="no") this%hwC(i,j)
         enddo
         write(600,*) '';
      enddo
      write(600,*) 'D array';
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            WRITE(600, 521, advance="no") this%hwD(i,j)
         enddo
         write(600,*) '';
      enddo
      
      allocate(ipiv(0:this%hwNum));
      call dgesv(this%hwNum,this%hwNum,this%hwA,this%hwNum,ipiv,Ematr,this%hwNum,info);

      do i=0,this%hwNum-1 
         do j=0,this%hwNum-1
            this%hwA(i,j)=Ematr(i,j);
         enddo
      enddo

      write(600,*) 'A^-1 array';
      do i=0,this%hwNum-1
         do j=0,this%hwNum-1
            WRITE(600, 521, advance="no") this%hwA(i,j)
         enddo
         write(600,*) '';
      enddo

      521   format (E11.4, ' ');   
      
      allocate(this%hwFnew(0:this%hwNum-1), this%hwF(0:this%hwNum-1), this%hwFold(0:this%hwNum-1), this%hwGF(0:this%hwNum-1));
      do i=0,this%hwNum-1
         this%hwFnew(i)=0.0d0;
         this%hwF(i)=0.0d0;
         this%hwFold(i)=0.0d0;
         this%hwGF(i)=0.0d0;
      enddo

      deallocate(ipiv);
      deallocate(Ematr);

    end subroutine initHWAuxVariables


!destroyHWAuxVariables---------------------------
    subroutine destroyHWAuxVariables(this)
    ! Destroy arrays for Givoli-Neta ABC aux variables
      class(tProblem) :: this;
      integer :: i;    
      do i=0,this%hwNum
         call this%hwX(i)%auxmesh_destroy;
         call this%hwY(i)%auxmesh_destroy;
         call this%hwZ(i)%auxmesh_destroy;
         call this%hwXOld(i)%auxmesh_destroy;
         call this%hwYOld(i)%auxmesh_destroy;
         call this%hwZOld(i)%auxmesh_destroy;
         call this%hwXOld2(i)%auxmesh_destroy;
         call this%hwYOld2(i)%auxmesh_destroy;
         call this%hwZOld2(i)%auxmesh_destroy;
      end do
      deallocate(this%hwX, this%hwY, this%hwZ, this%hwXOld, this%hwYOld, this%hwZOld, this%hwXOld2, this%hwYOld2, this%hwZOld2);
      deallocate(this%hwAng);
      deallocate(this%hwAC);
      deallocate(this%hwL);
      deallocate(this%hwA);
      deallocate(this%hwB);
      deallocate(this%hwC);
      deallocate(this%hwD);
      deallocate(this%hwFnew, this%hwF, this%hwFold, this%hwGF);
    end subroutine destroyHWAuxVariables


!saveHWAuxVariables---------------------------
    subroutine saveHWAuxVariables(this)
    ! Save arrays for HW ABC aux variables
      class(tProblem) :: this;
      integer :: i, j, k, v;
      do v=0,this%hwNum
      ! HWx field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%HWxOld2(v)%Y(i,j,k)=this%HWxOld(v)%Y(i,j,k);
               this%HWxOld(v)%Y(i,j,k)=this%HWx(v)%Y(i,j,k);
               this%HWxOld2(v)%Z(i,j,k)=this%HWxOld(v)%Z(i,j,k);
               this%HWxOld(v)%Z(i,j,k)=this%HWx(v)%Z(i,j,k);
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! HWy field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%HWyOld2(v)%X(i,j,k)=this%HWyOld(v)%X(i,j,k);
               this%HWyOld(v)%X(i,j,k)=this%HWy(v)%X(i,j,k);
               this%HWyOld2(v)%Z(i,j,k)=this%HWyOld(v)%Z(i,j,k);
               this%HWyOld(v)%Z(i,j,k)=this%HWy(v)%Z(i,j,k);
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      ! HWz field filling
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%HWzOld2(v)%X(i,j,k)=this%HWzOld(v)%X(i,j,k);
               this%HWzOld(v)%X(i,j,k)=this%HWz(v)%X(i,j,k);
               this%HWzOld2(v)%Y(i,j,k)=this%HWzOld(v)%Y(i,j,k);
               this%HWzOld(v)%Y(i,j,k)=this%HWz(v)%Y(i,j,k);
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      enddo
    end subroutine saveHWAuxVariables   


!fillEBoundaryABCHW------------------------------------------------------------------
  subroutine  fillEBoundaryABCHW(this,t)
    !Fill E boundary with ABC Hagstrom-Warburton Ey_x=0, exact elsewhere
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
    !      this%Ef%Y(0, j, k) = this%Source%getEpoint(this%xi(0), this%yi05(j),  this%zi(k), this%Ti(t), 2);
          this%Ef%Y(this%Nx, j, k) = this%Source%getEpoint(this%xi(this%Nx), this%yi05(j),  this%zi(k), this%Ti(t), this%Nx, j, k, 2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !!!*****************Ey_x=0 update*****************
     call this%HWEyUpdate;
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
  end subroutine fillEBoundaryABCHW

  
!fillEBoundaryABCHWPML------------------------------------------------------------------
  subroutine  fillEBoundaryABCHWPML(this,t)
    !Fill E boundary with ABC Hagstrom-Warburton High-order on Ey_x=0, Unsplit PML elsewhere
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k, v, it
    write(*,*) 'using pml';
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
    !!!*****************E renew*****************
    call this%HWEyUpdate;
    !!!****************E renew*****************
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
  end subroutine fillEBoundaryABCHWPML  

  
!HWEyUpdate------------------------------------------------------------------
  subroutine HWEyUpdate(this)
    !Fill Ey_x=0 boundary with ABC Hagstrom-Warburton
    class(tProblem) :: this
    integer :: j ,k, v, w
    real*8 :: de, df1, df2, tp;
    character*8 :: trans;
    real*8, allocatable :: temp(:);
    allocate(temp(0:this%hwNum-1));
    !=========================== EY == X BOUNDARIES ============================
    !EYFace X=0 -----
    do k=0,this%Nz
       do j=1,this%Ny
          select case(this%dxmode)
          case(1)
             de = (-this%EfOld%Y(0, j, k) + this%EfOld%Y(1, j, k))/(this%hhx);
          case(2)
             de = (-3.0d0*this%EfOld%Y(0, j, k) + 4.0d0*this%EfOld%Y(1, j, k) - this%EfOld%Y(2, j, k) )/(2.0d0*this%hhx);
          end select
          do v=0,this%hwNum-1 ! Vectors formation
             this%hwF(v) = this%hwYOld(v)%X(0,j,k);
             this%hwFold(v) = this%hwYold2(v)%X(0,j,k);
             select case(this%dxauxmode)
             case(1) !===1st order===
               if (j==1.AND.k==0) then                             ! j=1, k=0
                  df1 = (this%hwYOld(v)%X(0,1,0)-2.0d0*this%hwYOld(v)%X(0,2,0)+this%hwYOld(v)%X(0,3,0))/(this%hhy**2);   ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(0,1,0)-2.0d0*this%hwYOld(v)%X(0,1,1)+this%hwYOld(v)%X(0,1,2))/(this%hhz**2);   ! forward one-sided in Z
               elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                  df1 = (this%hwYOld(v)%X(0,this%Ny,0)-2.0d0*this%hwYOld(v)%X(0,this%Ny-1,0)+this%hwYOld(v)%X(0,this%Ny-2,0))/(this%hhy**2);     ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(0,this%Ny,0)-2.0d0*this%hwYOld(v)%X(0,this%Ny,1)+this%hwYOld(v)%X(0,this%Ny,2))/(this%hhz**2);         ! forward one-sided in Z 
               elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                  df1 = (this%hwYOld(v)%X(0,1,this%Nz)-2.0d0*this%hwYOld(v)%X(0,2,this%Nz)+this%hwYOld(v)%X(0,3,this%Nz))/(this%hhy**2);         ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(0,1,this%Nz)-2.0d0*this%hwYOld(v)%X(0,1,this%Nz-1)+this%hwYOld(v)%X(0,1,this%Nz-2))/(this%hhz**2);     ! backward one-sided in Z 
               elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(0,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(0,this%Ny-1,this%Nz)+this%hwYOld(v)%X(0,this%Ny-2,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(0,this%Ny,this%Nz)-2.0d0*this%hwYOld(v)%X(0,this%Ny,this%Nz-1)+this%hwYOld(v)%X(0,this%Ny,this%Nz-2))/(this%hhz**2);  ! backward one-sided in Z 
               elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                  df1 = (this%hwYOld(v)%X(0,j+1,0)-2.0d0*this%hwYOld(v)%X(0,j,0)+this%hwYOld(v)%X(0,j-1,0))/(this%hhy**2);                                      ! central in Y
                  df2 = (this%hwYOld(v)%X(0,j,0)-2.0d0*this%hwYOld(v)%X(0,j,1)+this%hwYOld(v)%X(0,j,2))/(this%hhz**2);                                          ! forward one-sided in Z      
               elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                  df1 = (this%hwYOld(v)%X(0,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(0,j,this%Nz)+this%hwYOld(v)%X(0,j-1,this%Nz))/(this%hhy**2);                    ! central in Y
                  df2 = (this%hwYOld(v)%X(0,j,this%Nz)-2.0d0*this%hwYOld(v)%X(0,j,this%Nz-1)+this%hwYOld(v)%X(0,j,this%Nz-2))/(this%hhz**2);                    ! backward one-sided Z 
               elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(0,1,k)-2.0d0*this%hwYOld(v)%X(0,2,k)+this%hwYOld(v)%X(0,3,k))/(this%hhy**2);                                           ! forward one-sided in Y
                  df2 = (this%hwYOld(v)%X(0,1,k+1)-2.0d0*this%hwYOld(v)%X(0,1,k)+this%hwYOld(v)%X(0,1,k-1))/(this%hhz**2);                                       ! central in Z 
               elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                  df1 = (this%hwYOld(v)%X(0,this%Ny,k)-2.0d0*this%hwYOld(v)%X(0,this%Ny-1,k)+this%hwYOld(v)%X(0,this%Ny-2,k))/(this%hhy**2);    ! backward one-sided in Y
                  df2 = (this%hwYOld(v)%X(0,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(0,this%Ny,k)+this%hwYOld(v)%X(0,this%Ny,k-1))/(this%hhz**2);    ! central in Z 
               else
                  df1 = (this%hwYOld(v)%X(0,j+1,k)-2.0d0*this%hwYOld(v)%X(0,j,k)+this%hwYOld(v)%X(0,j-1,k))/(this%hhy**2);    ! central in Y
                  df2 = (this%hwYOld(v)%X(0,j,k+1)-2.0d0*this%hwYOld(v)%X(0,j,k)+this%hwYOld(v)%X(0,j,k-1))/(this%hhz**2);    ! central in Z 
               endif
             case(2) !===2nd order===
                if (j==1.AND.k==0) then                             ! j=1, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(0,1,0)-5.0d0*this%hwYOld(v)%X(0,2,0)+4.0d0*this%hwYOld(v)%X(0,3,0)-this%hwYOld(v)%X(0,4,0))/(this%hhy**2);   ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(0,1,0)-5.0d0*this%hwYOld(v)%X(0,1,1)+4.0d0*this%hwYOld(v)%X(0,1,2)-this%hwYOld(v)%X(0,1,3))/(this%hhz**2);   ! forward one-sided in Z
                elseif (j==this%Ny.AND.k==0) then                   ! j=Ny, k=0
                   df1 = (2.0d0*this%hwYOld(v)%X(0,this%Ny,0)-5.0d0*this%hwYOld(v)%X(0,this%Ny-1,0)+4.0d0*this%hwYOld(v)%X(0,this%Ny-2,0)-this%hwYOld(v)%X(0,this%Ny-3,0))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(0,this%Ny,0)-5.0d0*this%hwYOld(v)%X(0,this%Ny,1)+4.0d0*this%hwYOld(v)%X(0,this%Ny,2)-this%hwYOld(v)%X(0,this%Ny,3))/(this%hhz**2);          ! forward one-sided in Z 
                elseif (j==1.AND.k==this%Nz) then                   ! j=1, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(0,1,this%Nz)-5.0d0*this%hwYOld(v)%X(0,2,this%Nz)+4.0d0*this%hwYOld(v)%X(0,3,this%Nz)-this%hwYOld(v)%X(0,4,this%Nz))/(this%hhy**2);          ! forward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(0,1,this%Nz)-5.0d0*this%hwYOld(v)%X(0,1,this%Nz-1)+4.0d0*this%hwYOld(v)%X(0,1,this%Nz-2)-this%hwYOld(v)%X(0,1,this%Nz-3))/(this%hhz**2);    ! backward one-sided in Z 
                elseif (j==this%Ny.AND.k==this%Nz) then             ! j=Ny, k=Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(0,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(0,this%Ny-1,this%Nz)+4.0d0*this%hwYOld(v)%X(0,this%Ny-2,this%Nz)-this%hwYOld(v)%X(0,this%Ny-3,this%Nz))/(this%hhy**2);  ! backward one-sided in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(0,this%Ny,this%Nz)-5.0d0*this%hwYOld(v)%X(0,this%Ny,this%Nz-1)+4.0d0*this%hwYOld(v)%X(0,this%Ny,this%Nz-2)-this%hwYOld(v)%X(0,this%Ny,this%Nz-3))/(this%hhz**2);  ! backward one-sided in Z 
                elseif (j>1.AND.j<this%Ny.AND.k==0) then            ! j>1, j<Ny, k=0
                   df1 = (this%hwYOld(v)%X(0,j+1,0)-2.0d0*this%hwYOld(v)%X(0,j,0)+this%hwYOld(v)%X(0,j-1,0))/(this%hhy**2);                                                                  ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(0,j,0)-5.0d0*this%hwYOld(v)%X(0,j,1)+4.0d0*this%hwYOld(v)%X(0,j,2)-this%hwYOld(v)%X(0,j,3))/(this%hhz**2);                                  ! forward one-sided in Z      
                elseif (j>1.AND.j<this%Ny.AND.k==this%Nz) then      ! j>1, j<Ny, k=Nz
                   df1 = (this%hwYOld(v)%X(0,j+1,this%Nz)-2.0d0*this%hwYOld(v)%X(0,j,this%Nz)+this%hwYOld(v)%X(0,j-1,this%Nz))/(this%hhy**2);                                                ! central in Y
                   df2 = (2.0d0*this%hwYOld(v)%X(0,j,this%Nz)-5.0d0*this%hwYOld(v)%X(0,j,this%Nz-1)+4.0d0*this%hwYOld(v)%X(0,j,this%Nz-2)-this%hwYOld(v)%X(0,j,this%Nz-3))/(this%hhz**2);    ! backward one-sided Z 
                elseif (j==1.AND.k>0.AND.k<this%Nz) then            ! j=1, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(0,1,k)-5.0d0*this%hwYOld(v)%X(0,2,k)+4.0d0*this%hwYOld(v)%X(0,3,k)-this%hwYOld(v)%X(0,4,k))/(this%hhy**2);                                  ! forward one-sided in Y
                   df2 = (this%hwYOld(v)%X(0,1,k+1)-2.0d0*this%hwYOld(v)%X(0,1,k)+this%hwYOld(v)%X(0,1,k-1))/(this%hhz**2);                                                                  ! central in Z 
                elseif (j==this%Ny.AND.k>0.AND.k<this%Nz) then      ! j=Ny, k>0, k<Nz
                   df1 = (2.0d0*this%hwYOld(v)%X(0,this%Ny,k)-5.0d0*this%hwYOld(v)%X(0,this%Ny-1,k)+4.0d0*this%hwYOld(v)%X(0,this%Ny-2,k)-this%hwYOld(v)%X(0,this%Ny-3,k))/(this%hhy**2);    ! backward one-sided in Y
                   df2 = (this%hwYOld(v)%X(0,this%Ny,k+1)-2.0d0*this%hwYOld(v)%X(0,this%Ny,k)+this%hwYOld(v)%X(0,this%Ny,k-1))/(this%hhz**2);                                                ! central in Z 
                else
                   df1 = (this%hwYOld(v)%X(0,j+1,k)-2.0d0*this%hwYOld(v)%X(0,j,k)+this%hwYOld(v)%X(0,j-1,k))/(this%hhy**2);    ! central in Y
                   df2 = (this%hwYOld(v)%X(0,j,k+1)-2.0d0*this%hwYOld(v)%X(0,j,k)+this%hwYOld(v)%X(0,j,k-1))/(this%hhz**2);    ! central in Z 
                endif
             end select
             this%hwGF(v) = df1+df2;
          enddo

          do v=0,this%hwNum-1
             temp(v) = 0.0d0;
          enddo
          ! B*Y^(n)
          do v=0,this%hwNum-1
             tp = 0.0d0;
             do w=0,this%hwNum-1
                tp = tp + this%hwB(v,w)*this%hwF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! C*Y^(n-1)
          do v=0,this%hwNum-1
             tp = 0.0d0;
             do w=0,this%hwNum-1
                tp = tp - this%hwC(v,w)*this%hwFold(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          ! D*GY^(n)
          do v=0,this%hwNum-1
             tp = 0.0d0;
             do w=0,this%hwNum-1
               tp = tp + this%hwD(v,w)*this%hwGF(w);
             enddo
             temp(v) = temp(v)+tp;
          enddo
          !add dx_u
          select case(this%dtmode)
          case(1)
             temp(0) = temp(0)+1.0d0*cc*this%ht*de;
          case(2)
             temp(0) = temp(0)+2.0d0*cc*this%ht*de;
          end select
          ! A^(-1)*()
          do v=0,this%hwNum-1
             tp = 0.0d0;
             do w=0,this%hwNum-1
                tp=tp + this%hwA(v,w)*temp(w);
             enddo
             this%hwY(v)%X(0,j,k) = tp;
          enddo
          this%Ef%Y(0, j, k) = this%hwY(0)%X(0,j,k);
       enddo
    enddo
!    Example of the First order ABC
!    do k=0,this%Nz
!       do j=1,this%Ny
!          de = (-3.0d0*this%EfOld%Y(0, j, k) + 4.0d0*this%EfOld%Y(1, j, k) - this%EfOld%Y(2, j, k) )/(2.0d0*this%hhx); 
!          this%Ef%Y(0, j, k) = 1.0d0/3.0d0*(4.0d0*this%EfOld%Y(0, j, k)-this%EfOld2%Y(0, j, k) + 2.0d0*this%ht*de);
!          !de = (-this%EfOld%Y(0, j, k) + this%EfOld%Y(1, j, k))/(this%hhx); 
!          !this%Ef%Y(0, j, k) = (this%EfOld%Y(0, j, k) + this%ht*de);
!       enddo
!    enddo
  end subroutine HWEyUpdate 
