! ==== PROBLEM CLASS ====
! Fill E boundaries with simple ABCs


!fillEBoundaryABCSommerfeld------------------------------------------------------------------
  subroutine  fillEBoundaryABCSommerfeld(this,t)
    !Fill E boundary with Sommerfeld 0th order ABC
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    real*8 :: gx, gy, gz, gm, dg;
    real*8 :: tm1, tm2;
    gx = (1.0d0-cc*this%ht/this%hhx)/(1.0d0+cc*this%ht/this%hhx);
    gy = (1.0d0-cc*this%ht/this%hhy)/(1.0d0+cc*this%ht/this%hhy);
    gz = (1.0d0-cc*this%ht/this%hhz)/(1.0d0+cc*this%ht/this%hhz);
    dg = (1.0d0-cc*this%ht/(sqrt(2.0d0)*this%hhx))/(1.0d0+cc*this%ht/(sqrt(2.0d0)*this%hhx));
    !===== EX =====
    ! Faces
    do k=1,this%Nz-1
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = gy*this%Ef%X(i, 0, k) - gy*this%Ef%X(i, 1, k) + this%Efold%X(i, 1, k);
          this%Ef%X(i, this%Ny, k) = gy*this%Ef%X(i, this%Ny, k) - gy*this%Ef%X(i, this%Ny-1, k) + this%Efold%X(i, this%Ny-1, k);
       enddo
    enddo
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = gz*this%Ef%X(i, j, 0) - gz*this%Ef%X(i, j, 1) + this%Efold%X(i, j, 1);
          this%Ef%X(i, j, this%Nz) = gz*this%Ef%X(i, j, this%Nz) - gz*this%Ef%X(i, j, this%Nz-1) + this%Efold%X(i, j, this%Nz-1);
       enddo
    enddo
    ! Edges
    do i=1,this%Nx
       tm1 = gy*this%Ef%X(i, 0, 0) - gy*this%Ef%X(i, 1, 0) + this%Efold%X(i, 1, 0);
       tm2 = gz*this%Ef%X(i, 0, 0) - gz*this%Ef%X(i, 0, 1) + this%Efold%X(i, 0, 1);
       this%Ef%X(i, 0, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, 0, 0) = dg*this%Ef%X(i, 0, 0) - dg*this%Ef%X(i, 1, 1) + this%Efold%X(i, 1, 1);
       tm1 = gy*this%Ef%X(i, this%Ny, 0) - gy*this%Ef%X(i, this%Ny-1, 0) + this%Efold%X(i, this%Ny-1, 0);
       tm2 = gz*this%Ef%X(i, this%Ny, 0) - gz*this%Ef%X(i, this%Ny, 1) + this%Efold%X(i, this%Ny, 1);
       this%Ef%X(i, this%Ny, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, this%Ny, 0) = dg*this%Ef%X(i, this%Ny, 0) - dg*this%Ef%X(i, this%Ny-1, 1) + this%Efold%X(i, this%Ny-1, 1);
       tm1 = gy*this%Ef%X(i, 0, this%Nz) - gy*this%Ef%X(i, 1, this%Nz) + this%Efold%X(i, 1, this%Nz);
       tm2 = gz*this%Ef%X(i, 0, this%Nz) - gz*this%Ef%X(i, 0, this%Nz-1) + this%Efold%X(i, 0, this%Nz-1);
       this%Ef%X(i, 0, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, 0, this%Nz) = dg*this%Ef%X(i, 0, this%Nz) - dg*this%Ef%X(i, 1, this%Nz-1) + this%Efold%X(i, 1, this%Nz-1);
       tm1 = gy*this%Ef%X(i, this%Ny, this%Nz) - gy*this%Ef%X(i, this%Ny-1, this%Nz) + this%Efold%X(i, this%Ny-1, this%Nz);
       tm2 = gz*this%Ef%X(i, this%Ny, this%Nz) - gz*this%Ef%X(i, this%Ny, this%Nz-1) + this%Efold%X(i, this%Ny, this%Nz-1);
       this%Ef%X(i, this%Ny, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, this%Ny, this%Nz) = dg*this%Ef%X(i, this%Ny, this%Nz) - dg*this%Ef%X(i, this%Ny-1, this%Nz-1) + this%Efold%X(i, this%Ny-1, this%Nz-1);
    enddo
    
    !===== EY =====
    ! Faces
    do j=1,this%Ny
       do i=1,this%Nx-1
          this%Ef%Y(i, j, 0) = gz*this%Ef%Y(i, j, 0)-gz*this%Ef%Y(i, j, 1)+this%Efold%Y(i, j, 1);
          this%Ef%Y(i, j, this%Nz) = gz*this%Ef%Y(i, j, this%Nz)-gz*this%Ef%Y(i, j, this%Nz-1)+this%Efold%Y(i, j, this%Nz-1);
       enddo
    enddo     
    do k=1,this%Nz-1
       do j=1,this%Ny
          this%Ef%Y(0, j, k) = gx*this%Ef%Y(0, j, k)-gx*this%Ef%Y(1, j, k)+this%Efold%Y(1, j, k);
          this%Ef%Y(this%Nx, j, k) = gx*this%Ef%Y(this%Nx, j, k)-gx*this%Ef%Y(this%Nx-1, j, k)+this%Efold%Y(this%Nx-1, j, k);
       enddo
    enddo
    ! Edges
    do j=1,this%Ny
       tm1 = gx*this%Ef%Y(0, j, 0) - gx*this%Ef%Y(1, j, 0) + this%Efold%Y(1, j, 0);
       tm2 = gz*this%Ef%Y(0, j, 0) - gz*this%Ef%Y(0, j, 1) + this%Efold%Y(0, j, 1);
       this%Ef%Y(0, j, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(0, j, 0) = dg*this%Ef%Y(0, j, 0) - dg*this%Ef%Y(1, j, 1) + this%Efold%Y(1, j, 1);
       tm1 = gx*this%Ef%Y(this%Nx, j, 0) - gx*this%Ef%Y(this%Nx-1, j, 0) + this%Efold%Y(this%Nx-1, j, 0);
       tm2 = gz*this%Ef%Y(this%Nx, j, 0) - gz*this%Ef%Y(this%Nx, j, 1) + this%Efold%Y(this%Nx, j, 1);
       this%Ef%Y(this%Nx, j, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(this%Nx, j, 0) = dg*this%Ef%Y(this%Nx, j, 0) - dg*this%Ef%Y(this%Nx-1, j, 1) + this%Efold%Y(this%Nx-1, j, 1);
       tm1 = gx*this%Ef%Y(0, j, this%Nz) - gx*this%Ef%Y(1, j, this%Nz) + this%Efold%Y(1, j, this%Nz);
       tm2 = gz*this%Ef%Y(0, j, this%Nz) - gz*this%Ef%Y(0, j, this%Nz-1) + this%Efold%Y(0, j, this%Nz-1);
       this%Ef%Y(0, j, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(0, j, this%Nz) = dg*this%Ef%Y(0, j, this%Nz) - dg*this%Ef%Y(1, j, this%Nz-1) + this%Efold%Y(1, j, this%Nz-1);   
       tm1 = gx*this%Ef%Y(this%Nx, j, this%Nz) - gx*this%Ef%Y(this%Nx-1, j, this%Nz) + this%Efold%Y(this%Nx-1, j, this%Nz);
       tm2 = gz*this%Ef%Y(this%Nx, j, this%Nz) - gz*this%Ef%Y(this%Nx, j, this%Nz-1) + this%Efold%Y(this%Nx, j, this%Nz-1);
       this%Ef%Y(this%Nx, j, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(this%Nx, j, this%Nz) = dg*this%Ef%Y(this%Nx, j, this%Nz) - dg*this%Ef%Y(this%Nx-1, j, this%Nz-1) + this%Efold%Y(this%Nx-1, j, this%Nz-1);
    enddo
    
    !===== EZ =====
    !Faces
    do k=1,this%Nz
       do j=1,this%Ny-1
          this%Ef%Z(0, j, k) = gx*this%Ef%Z(0, j, k)-gx*this%Ef%Z(1, j, k)+this%Efold%Z(1, j, k);
          this%Ef%Z(this%Nx, j, k) = gx*this%Ef%Z(this%Nx, j, k)-gx*this%Ef%Z(this%Nx-1, j, k)+this%Efold%Z(this%Nx-1, j, k);
       enddo
    enddo
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = gy*this%Ef%Z(i, 0, k)-gy*this%Ef%Z(i, 1, k)+this%Efold%Z(i, 1, k);
          this%Ef%Z(i, this%Ny, k) = gy*this%Ef%Z(i, this%Ny, k)-gy*this%Ef%Z(i, this%Ny-1, k)+this%Efold%Z(i, this%Ny-1, k);
       enddo
    enddo
    ! Edges
    do k=1,this%Nz
       tm1 = gx*this%Ef%Z(0, 0, k) - gx*this%Ef%Z(1, 0, k) + this%Efold%Z(1, 0, k);
       tm2 = gy*this%Ef%Z(0, 0, k) - gy*this%Ef%Z(0, 1, k) + this%Efold%Z(0, 1, k);
       this%Ef%Z(0, 0, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(0, 0, k) = gx*this%Ef%Z(0, 0, k) - gx*this%Ef%Z(1, 1, k) + this%Efold%Z(1, 1, k);
       tm1 = gx*this%Ef%Z(this%Nx, 0, k) - gx*this%Ef%Z(this%Nx-1, 0, k) + this%Efold%Z(this%Nx-1, 0, k);
       tm2 = gy*this%Ef%Z(this%Nx, 0, k) - gy*this%Ef%Z(this%Nx, 1, k) + this%Efold%Z(this%Nx, 1, k);
       this%Ef%Z(this%Nx, 0, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(this%Nx, 0, k) = gx*this%Ef%Z(this%Nx, 0, k) - gx*this%Ef%Z(this%Nx-1, 1, k) + this%Efold%Z(this%Nx-1, 1, k);
       tm1 = gx*this%Ef%Z(0, this%Ny, k) - gx*this%Ef%Z(1, this%Ny, k) + this%Efold%Z(1, this%Ny, k);
       tm2 = gy*this%Ef%Z(0, this%Ny, k) - gy*this%Ef%Z(0, this%Ny-1, k) + this%Efold%Z(0, this%Ny-1, k);
       this%Ef%Z(0, this%Ny, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(0, this%Ny, k) = gx*this%Ef%Z(0, this%Ny, k) - gx*this%Ef%Z(1, this%Ny-1, k) + this%Efold%Z(1, this%Ny-1, k);
       tm1 = gx*this%Ef%Z(this%Nx, this%Ny, k) - gx*this%Ef%Z(this%Nx-1, this%Ny, k) + this%Efold%Z(this%Nx-1, this%Ny, k);
       tm2 = gy*this%Ef%Z(this%Nx, this%Ny, k) - gy*this%Ef%Z(this%Nx, this%Ny-1, k) + this%Efold%Z(this%Nx, this%Ny-1, k);
       this%Ef%Z(this%Nx, this%Ny, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(this%Nx, this%Ny, k) = gx*this%Ef%Z(this%Nx, this%Ny, k) - gx*this%Ef%Z(this%Nx-1, this%Ny-1, k) + this%Efold%Z(this%Nx-1, this%Ny-1, k);
    enddo   
  end subroutine fillEBoundaryABCSommerfeld


!fillEBoundaryABCHigdon------------------------------------------------------------------
  subroutine  fillEBoundaryABCHigdon(this,t)
    !Fill E boundary with Higdon 2 angles ABC
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    real*8 :: gx1, gx2, gy1, gy2, gz1, gz2;
    real*8 :: tm1, tm2;
    gx1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx));
    gx2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx));
    gy1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy));
    gy2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy));
    gz1 = (1.0d0-cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhz));
    gz2 = (1.0d0-cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx))/(1.0d0+cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhz));
    !===== EX =====
    ! Faces
     !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = -(gy1+gy2)*this%Ef%X(i, 1, k)-gy1*gy2*this%Ef%X(i, 2, k)  &
          +(gy1+gy2)*this%Efold%X(i, 0, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%X(i, 1, k)+(gy1+gy2)*this%Efold%X(i, 2, k)  &
          -gy1*gy2*this%Efold2%X(i, 0, k)-(gy1+gy2)*this%Efold2%X(i, 1, k)-this%Efold2%X(i, 2, k);

          this%Ef%X(i, this%Ny, k) = -(gy1+gy2)*this%Ef%X(i, this%Ny-1, k)-gy1*gy2*this%Ef%X(i, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%X(i, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%X(i, this%Ny-1, k)+(gy1+gy2)*this%Efold%X(i, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%X(i, this%Ny, k)-(gy1+gy2)*this%Efold2%X(i, this%Ny-1, k)-this%Efold2%X(i, this%Ny-2, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = -(gz1+gz2)*this%Ef%X(i, j, 1)-gz1*gz2*this%Ef%X(i, j, 2)  &
          +(gz1+gz2)*this%Efold%X(i, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%X(i, j, 1)+(gz1+gz2)*this%Efold%X(i, j, 2)  &
          -gz1*gz2*this%Efold2%X(i, j, 0)-(gz1+gz2)*this%Efold2%X(i, j, 1)-this%Efold2%X(i, j, 2);
          
          this%Ef%X(i, j, this%Nz) = -(gz1+gz2)*this%Ef%X(i, j, this%Nz-1)-gz1*gz2*this%Ef%X(i, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%X(i, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%X(i, j, this%Nz-1)+(gz1+gz2)*this%Efold%X(i, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%X(i, j, this%Nz)-(gz1+gz2)*this%Efold2%X(i, j, this%Nz-1)-this%Efold2%X(i, j, this%Nz-2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do i=1,this%Nx
       tm1 = -(gy1+gy2)*this%Ef%X(i, 1, 0)-gy1*gy2*this%Ef%X(i, 2, 0)  &
          +(gy1+gy2)*this%Efold%X(i, 0, 0)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%X(i, 1, 0)+(gy1+gy2)*this%Efold%X(i, 2, 0)  &
          -gy1*gy2*this%Efold2%X(i, 0, 0)-(gy1+gy2)*this%Efold2%X(i, 1, 0)-this%Efold2%X(i, 2, 0);       
       tm2 = -(gz1+gz2)*this%Ef%X(i, 0, 1)-gz1*gz2*this%Ef%X(i, 0, 2)  &
          +(gz1+gz2)*this%Efold%X(i, 0, 0)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%X(i, 0, 1)+(gz1+gz2)*this%Efold%X(i, 0, 2)  &
          -gz1*gz2*this%Efold2%X(i, 0, 0)-(gz1+gz2)*this%Efold2%X(i, 0, 1)-this%Efold2%X(i, 0, 2);
       this%Ef%X(i, 0, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, 0, 0) = dg*this%Ef%X(i, 0, 0) - dg*this%Ef%X(i, 1, 1) + this%Efold%X(i, 1, 1);
       
       tm1 = -(gy1+gy2)*this%Ef%X(i, this%Ny-1, 0)-gy1*gy2*this%Ef%X(i, this%Ny-2, 0)  &
          +(gy1+gy2)*this%Efold%X(i, this%Ny, 0)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%X(i, this%Ny-1, 0)+(gy1+gy2)*this%Efold%X(i, this%Ny-2, 0)  &
          -gy1*gy2*this%Efold2%X(i, this%Ny, 0)-(gy1+gy2)*this%Efold2%X(i, this%Ny-1, 0)-this%Efold2%X(i, this%Ny-2, 0);
       tm2 = -(gz1+gz2)*this%Ef%X(i, this%Ny, 1)-gz1*gz2*this%Ef%X(i, this%Ny, 2)  &
          +(gz1+gz2)*this%Efold%X(i, this%Ny, 0)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%X(i, this%Ny, 1)+(gz1+gz2)*this%Efold%X(i, this%Ny, 2)  &
          -gz1*gz2*this%Efold2%X(i, this%Ny, 0)-(gz1+gz2)*this%Efold2%X(i, this%Ny, 1)-this%Efold2%X(i, this%Ny, 2);
       this%Ef%X(i, this%Ny, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, this%Ny, 0) = dg*this%Ef%X(i, this%Ny, 0) - dg*this%Ef%X(i, this%Ny-1, 1) + this%Efold%X(i, this%Ny-1, 1);

       tm1 = -(gy1+gy2)*this%Ef%X(i, 1, this%Nz)-gy1*gy2*this%Ef%X(i, 2, this%Nz)  &
          +(gy1+gy2)*this%Efold%X(i, 0, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%X(i, 1, this%Nz)+(gy1+gy2)*this%Efold%X(i, 2, this%Nz)  &
          -gy1*gy2*this%Efold2%X(i, 0, this%Nz)-(gy1+gy2)*this%Efold2%X(i, 1, this%Nz)-this%Efold2%X(i, 2, this%Nz);    
       tm2 = -(gz1+gz2)*this%Ef%X(i, 0, this%Nz-1)-gz1*gz2*this%Ef%X(i, 0, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%X(i, 0, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%X(i, 0, this%Nz-1)+(gz1+gz2)*this%Efold%X(i, 0, this%Nz-2)  &
          -gz1*gz2*this%Efold2%X(i, 0, this%Nz)-(gz1+gz2)*this%Efold2%X(i, 0, this%Nz-1)-this%Efold2%X(i, 0, this%Nz-2);
       this%Ef%X(i, 0, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, 0, this%Nz) = dg*this%Ef%X(i, 0, this%Nz) - dg*this%Ef%X(i, 1, this%Nz-1) + this%Efold%X(i, 1, this%Nz-1);

       tm1 = -(gy1+gy2)*this%Ef%X(i, this%Ny-1, this%Nz)-gy1*gy2*this%Ef%X(i, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%Efold%X(i, this%Ny, this%Nz)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%X(i, this%Ny-1, this%Nz)+(gy1+gy2)*this%Efold%X(i, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%Efold2%X(i, this%Ny, this%Nz)-(gy1+gy2)*this%Efold2%X(i, this%Ny-1, this%Nz)-this%Efold2%X(i, this%Ny-2, this%Nz);
       tm2 = -(gz1+gz2)*this%Ef%X(i, this%Ny, this%Nz-1)-gz1*gz2*this%Ef%X(i, this%Ny, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%X(i, this%Ny, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%X(i, this%Ny, this%Nz-1)+(gz1+gz2)*this%Efold%X(i, this%Ny, this%Nz-2)  &
          -gz1*gz2*this%Efold2%X(i, this%Ny, this%Nz)-(gz1+gz2)*this%Efold2%X(i, this%Ny, this%Nz-1)-this%Efold2%X(i, this%Ny, this%Nz-2);
       this%Ef%X(i, this%Ny, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, this%Ny, this%Nz) = dg*this%Ef%X(i, this%Ny, this%Nz) - dg*this%Ef%X(i, this%Ny-1, this%Nz-1) + this%Efold%X(i, this%Ny-1, this%Nz-1);
    enddo
    !$OMP END PARALLEL DO
    !===== EY =====
    ! Faces
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=1,this%Nx-1
          this%Ef%Y(i, j, 0) = -(gz1+gz2)*this%Ef%Y(i, j, 1)-gz1*gz2*this%Ef%Y(i, j, 2)  &
          +(gz1+gz2)*this%Efold%Y(i, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%Y(i, j, 1)+(gz1+gz2)*this%Efold%Y(i, j, 2)  &
          -gz1*gz2*this%Efold2%Y(i, j, 0)-(gz1+gz2)*this%Efold2%Y(i, j, 1)-this%Efold2%Y(i, j, 2);

          this%Ef%Y(i, j, this%Nz) = -(gz1+gz2)*this%Ef%Y(i, j, this%Nz-1)-gz1*gz2*this%Ef%Y(i, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%Y(i, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%Y(i, j, this%Nz-1)+(gz1+gz2)*this%Efold%Y(i, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%Y(i, j, this%Nz)-(gz1+gz2)*this%Efold2%Y(i, j, this%Nz-1)-this%Efold2%Y(i, j, this%Nz-2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=1,this%Ny
          this%Ef%Y(0, j, k) = -(gx1+gx2)*this%Ef%Y(1, j, k)-gx1*gx2*this%Ef%Y(2, j, k)  &
          +(gx1+gx2)*this%Efold%Y(0, j, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Y(1, j, k)+(gx1+gx2)*this%Efold%Y(2, j, k)  &
          -gx1*gx2*this%Efold2%Y(0, j, k)-(gx1+gx2)*this%Efold2%Y(1, j, k)-this%Efold2%Y(2, j, k);
 
          this%Ef%Y(this%Nx, j, k) = -(gx1+gx2)*this%Ef%Y(this%Nx-1, j, k)-gx1*gx2*this%Ef%Y(this%Nx-2, j, k)  &
          +(gx1+gx2)*this%Efold%Y(this%Nx, j, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Y(this%Nx-1, j, k)+(gx1+gx2)*this%Efold%Y(this%Nx-2, j, k)  &
          -gx1*gx2*this%Efold2%Y(this%Nx, j, k)-(gx1+gx2)*this%Efold2%Y(this%Nx-1, j, k)-this%Efold2%Y(this%Nx-2, j, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       tm1 = -(gx1+gx2)*this%Ef%Y(1, j, 0)-gx1*gx2*this%Ef%Y(2, j, 0)  &
          +(gx1+gx2)*this%Efold%Y(0, j, 0)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Y(1, j, 0)+(gx1+gx2)*this%Efold%Y(2, j, 0)  &
          -gx1*gx2*this%Efold2%Y(0, j, 0)-(gx1+gx2)*this%Efold2%Y(1, j, 0)-this%Efold2%Y(2, j, 0);
       tm2 = -(gz1+gz2)*this%Ef%Y(0, j, 1)-gz1*gz2*this%Ef%Y(0, j, 2)  &
          +(gz1+gz2)*this%Efold%Y(0, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%Y(0, j, 1)+(gz1+gz2)*this%Efold%Y(0, j, 2)  &
          -gz1*gz2*this%Efold2%Y(0, j, 0)-(gz1+gz2)*this%Efold2%Y(0, j, 1)-this%Efold2%Y(0, j, 2);
       this%Ef%Y(0, j, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(0, j, 0) = dg*this%Ef%Y(0, j, 0) - dg*this%Ef%Y(1, j, 1) + this%Efold%Y(1, j, 1);
       
       tm1 = -(gx1+gx2)*this%Ef%Y(this%Nx-1, j, 0)-gx1*gx2*this%Ef%Y(this%Nx-2, j, 0)  &
          +(gx1+gx2)*this%Efold%Y(this%Nx, j, 0)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Y(this%Nx-1, j, 0)+(gx1+gx2)*this%Efold%Y(this%Nx-2, j, 0)  &
          -gx1*gx2*this%Efold2%Y(this%Nx, j, 0)-(gx1+gx2)*this%Efold2%Y(this%Nx-1, j, 0)-this%Efold2%Y(this%Nx-2, j, 0);
       tm2 = -(gz1+gz2)*this%Ef%Y(this%Nx, j, 1)-gz1*gz2*this%Ef%Y(this%Nx, j, 2)  &
          +(gz1+gz2)*this%Efold%Y(this%Nx, j, 0)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%Y(this%Nx, j, 1)+(gz1+gz2)*this%Efold%Y(this%Nx, j, 2)  &
          -gz1*gz2*this%Efold2%Y(this%Nx, j, 0)-(gz1+gz2)*this%Efold2%Y(this%Nx, j, 1)-this%Efold2%Y(this%Nx, j, 2);
       this%Ef%Y(this%Nx, j, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(this%Nx, j, 0) = dg*this%Ef%Y(this%Nx, j, 0) - dg*this%Ef%Y(this%Nx-1, j, 1) + this%Efold%Y(this%Nx-1, j, 1);

       tm1 = -(gx1+gx2)*this%Ef%Y(1, j, this%Nz)-gx1*gx2*this%Ef%Y(2, j, this%Nz)  &
          +(gx1+gx2)*this%Efold%Y(0, j, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Y(1, j, this%Nz)+(gx1+gx2)*this%Efold%Y(2, j, this%Nz)  &
          -gx1*gx2*this%Efold2%Y(0, j, this%Nz)-(gx1+gx2)*this%Efold2%Y(1, j, this%Nz)-this%Efold2%Y(2, j, this%Nz);
       tm2 = -(gz1+gz2)*this%Ef%Y(0, j, this%Nz-1)-gz1*gz2*this%Ef%Y(0, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%Y(0, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%Y(0, j, this%Nz-1)+(gz1+gz2)*this%Efold%Y(0, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%Y(0, j, this%Nz)-(gz1+gz2)*this%Efold2%Y(0, j, this%Nz-1)-this%Efold2%Y(0, j, this%Nz-2);     
       this%Ef%Y(0, j, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(0, j, this%Nz) = dg*this%Ef%Y(0, j, this%Nz) - dg*this%Ef%Y(1, j, this%Nz-1) + this%Efold%Y(1, j, this%Nz-1);   

       tm1 = -(gx1+gx2)*this%Ef%Y(this%Nx-1, j, this%Nz)-gx1*gx2*this%Ef%Y(this%Nx-2, j, this%Nz)  &
          +(gx1+gx2)*this%Efold%Y(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Y(this%Nx-1, j, this%Nz)+(gx1+gx2)*this%Efold%Y(this%Nx-2, j, this%Nz)  &
          -gx1*gx2*this%Efold2%Y(this%Nx, j, this%Nz)-(gx1+gx2)*this%Efold2%Y(this%Nx-1, j, this%Nz)-this%Efold2%Y(this%Nx-2, j, this%Nz);
       tm2 = -(gz1+gz2)*this%Ef%Y(this%Nx, j, this%Nz-1)-gz1*gz2*this%Ef%Y(this%Nx, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%Y(this%Nx, j, this%Nz)+(2.0d0+2.0d0*gz1*gz2)*this%Efold%Y(this%Nx, j, this%Nz-1)+(gz1+gz2)*this%Efold%Y(this%Nx, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%Y(this%Nx, j, this%Nz)-(gz1+gz2)*this%Efold2%Y(this%Nx, j, this%Nz-1)-this%Efold2%Y(this%Nx, j, this%Nz-2);
       this%Ef%Y(this%Nx, j, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%Y(this%Nx, j, this%Nz) = dg*this%Ef%Y(this%Nx, j, this%Nz) - dg*this%Ef%Y(this%Nx-1, j, this%Nz-1) + this%Efold%Y(this%Nx-1, j, this%Nz-1);
    enddo
    !$OMP END PARALLEL DO
    
    !===== EZ =====
    ! Faces
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=1,this%Ny-1
          this%Ef%Z(0, j, k) = -(gx1+gx2)*this%Ef%Z(1, j, k)-gx1*gx2*this%Ef%Z(2, j, k)  &
          +(gx1+gx2)*this%Efold%Z(0, j, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Z(1, j, k)+(gx1+gx2)*this%Efold%Z(2, j, k)  &
          -gx1*gx2*this%Efold2%Z(0, j, k)-(gx1+gx2)*this%Efold2%Z(1, j, k)-this%Efold2%Z(2, j, k);
 
          this%Ef%Z(this%Nx, j, k) = -(gx1+gx2)*this%Ef%Z(this%Nx-1, j, k)-gx1*gx2*this%Ef%Z(this%Nx-2, j, k)  &
          +(gx1+gx2)*this%Efold%Z(this%Nx, j, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Z(this%Nx-1, j, k)+(gx1+gx2)*this%Efold%Z(this%Nx-2, j, k)  &
          -gx1*gx2*this%Efold2%Z(this%Nx, j, k)-(gx1+gx2)*this%Efold2%Z(this%Nx-1, j, k)-this%Efold2%Z(this%Nx-2, j, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = -(gy1+gy2)*this%Ef%Z(i, 1, k)-gy1*gy2*this%Ef%Z(i, 2, k)  &
          +(gy1+gy2)*this%Efold%Z(i, 0, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%Z(i, 1, k)+(gy1+gy2)*this%Efold%Z(i, 2, k)  &
          -gy1*gy2*this%Efold2%Z(i, 0, k)-(gy1+gy2)*this%Efold2%Z(i, 1, k)-this%Efold2%Z(i, 2, k);

          this%Ef%Z(i, this%Ny, k) = -(gy1+gy2)*this%Ef%Z(i, this%Ny-1, k)-gy1*gy2*this%Ef%Z(i, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%Z(i, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%Z(i, this%Ny-1, k)+(gy1+gy2)*this%Efold%Z(i, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%Z(i, this%Ny, k)-(gy1+gy2)*this%Efold2%Z(i, this%Ny-1, k)-this%Efold2%Z(i, this%Ny-2, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       tm1 = -(gx1+gx2)*this%Ef%Z(1, 0, k)-gx1*gx2*this%Ef%Z(2, 0, k)  &
          +(gx1+gx2)*this%Efold%Z(0, 0, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Z(1, 0, k)+(gx1+gx2)*this%Efold%Z(2, 0, k)  &
          -gx1*gx2*this%Efold2%Z(0, 0, k)-(gx1+gx2)*this%Efold2%Z(1, 0, k)-this%Efold2%Z(2, 0, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(0, 1, k)-gy1*gy2*this%Ef%Z(0, 2, k)  &
          +(gy1+gy2)*this%Efold%Z(0, 0, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%Z(0, 1, k)+(gy1+gy2)*this%Efold%Z(0, 2, k)  &
          -gy1*gy2*this%Efold2%Z(0, 0, k)-(gy1+gy2)*this%Efold2%Z(0, 1, k)-this%Efold2%Z(0, 2, k);
       this%Ef%Z(0, 0, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(0, 0, k) = gx*this%Ef%Z(0, 0, k) - gx*this%Ef%Z(1, 1, k) + this%Efold%Z(1, 1, k);

       tm1 = -(gx1+gx2)*this%Ef%Z(this%Nx-1, 0, k)-gx1*gx2*this%Ef%Z(this%Nx-2, 0, k)  &
          +(gx1+gx2)*this%Efold%Z(this%Nx, 0, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Z(this%Nx-1, 0, k)+(gx1+gx2)*this%Efold%Z(this%Nx-2, 0, k)  &
          -gx1*gx2*this%Efold2%Z(this%Nx, 0, k)-(gx1+gx2)*this%Efold2%Z(this%Nx-1, 0, k)-this%Efold2%Z(this%Nx-2, 0, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(this%Nx, 1, k)-gy1*gy2*this%Ef%Z(this%Nx, 2, k)  &
          +(gy1+gy2)*this%Efold%Z(this%Nx, 0, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%Z(this%Nx, 1, k)+(gy1+gy2)*this%Efold%Z(this%Nx, 2, k)  &
          -gy1*gy2*this%Efold2%Z(this%Nx, 0, k)-(gy1+gy2)*this%Efold2%Z(this%Nx, 1, k)-this%Efold2%Z(this%Nx, 2, k);
       this%Ef%Z(this%Nx, 0, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(this%Nx, 0, k) = gx*this%Ef%Z(this%Nx, 0, k) - gx*this%Ef%Z(this%Nx-1, 1, k) + this%Efold%Z(this%Nx-1, 1, k);

       tm1 = -(gx1+gx2)*this%Ef%Z(1, this%Ny, k)-gx1*gx2*this%Ef%Z(2, this%Ny, k)  &
          +(gx1+gx2)*this%Efold%Z(0, this%Ny, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Z(1, this%Ny, k)+(gx1+gx2)*this%Efold%Z(2, this%Ny, k)  &
          -gx1*gx2*this%Efold2%Z(0, this%Ny, k)-(gx1+gx2)*this%Efold2%Z(1, this%Ny, k)-this%Efold2%Z(2, this%Ny, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(0, this%Ny-1, k)-gy1*gy2*this%Ef%Z(0, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%Z(0, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%Z(0, this%Ny-1, k)+(gy1+gy2)*this%Efold%Z(0, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%Z(0, this%Ny, k)-(gy1+gy2)*this%Efold2%Z(0, this%Ny-1, k)-this%Efold2%Z(0, this%Ny-2, k);
       this%Ef%Z(0, this%Ny, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(0, this%Ny, k) = gx*this%Ef%Z(0, this%Ny, k) - gx*this%Ef%Z(1, this%Ny-1, k) + this%Efold%Z(1, this%Ny-1, k);

       tm1 = -(gx1+gx2)*this%Ef%Z(this%Nx-1, this%Ny, k)-gx1*gx2*this%Ef%Z(this%Nx-2, this%Ny, k)  &
          +(gx1+gx2)*this%Efold%Z(this%Nx, this%Ny, k)+(2.0d0+2.0d0*gx1*gx2)*this%Efold%Z(this%Nx-1, this%Ny, k)+(gx1+gx2)*this%Efold%Z(this%Nx-2, this%Ny, k)  &
          -gx1*gx2*this%Efold2%Z(this%Nx, this%Ny, k)-(gx1+gx2)*this%Efold2%Z(this%Nx-1, this%Ny, k)-this%Efold2%Z(this%Nx-2, this%Ny, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(this%Nx, this%Ny-1, k)-gy1*gy2*this%Ef%Z(this%Nx, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%Z(this%Nx, this%Ny, k)+(2.0d0+2.0d0*gy1*gy2)*this%Efold%Z(this%Nx, this%Ny-1, k)+(gy1+gy2)*this%Efold%Z(this%Nx, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%Z(this%Nx, this%Ny, k)-(gy1+gy2)*this%Efold2%Z(this%Nx, this%Ny-1, k)-this%Efold2%Z(this%Nx, this%Ny-2, k);
       this%Ef%Z(this%Nx, this%Ny, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(this%Nx, this%Ny, k) = gx*this%Ef%Z(this%Nx, this%Ny, k) - gx*this%Ef%Z(this%Nx-1, this%Ny-1, k) + this%Efold%Z(this%Nx-1, this%Ny-1, k);
    enddo
    !$OMP END PARALLEL DO
  end subroutine fillEBoundaryABCHigdon


!fillEBoundaryABCBetzMittra------------------------------------------------------------------
  subroutine  fillEBoundaryABCBetzMittra(this,t)
    !Fill E boundary with Betz-Mittra ABC
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    real*8 :: alp, ro1x, ro2x, betax, gx1, gx2, ro1y, ro2y, betay, gy1, gy2, ro1z, ro2z, betaz, gz1, gz2;
    real*8 :: tm1, tm2;
    alp = 0.1d0/this%hhx;
    ro1x = cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhx);
    ro2x = cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhx);
    betax = (1.0d0+ro1x)/(1.0d0+ro1x*(1.0d0+alp*this%hhx));
    gx1 = (1.0d0-ro1x)/(1.0d0+ro1x*(1.0d0+alp*this%hhx));
    gx2 = (1.0d0-ro2x)/(1.0d0+ro2x*(1.0d0+alp*this%hhx));

    ro1y = cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhy);
    ro2y = cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhy);
    betay = (1.0d0+ro1y)/(1.0d0+ro1y*(1.0d0+alp*this%hhy));
    gy1 = (1.0d0-ro1y)/(1.0d0+ro1y*(1.0d0+alp*this%hhy));
    gy2 = (1.0d0-ro2y)/(1.0d0+ro2y*(1.0d0+alp*this%hhy));

    ro1z = cc*this%ht/(cos(this%hig1/180.0d0*PI)*this%hhz);
    ro2z = cc*this%ht/(cos(this%hig2/180.0d0*PI)*this%hhz);
    betaz = (1.0d0+ro1z)/(1.0d0+ro1z*(1.0d0+alp*this%hhz));
    gz1 = (1.0d0-ro1z)/(1.0d0+ro1z*(1.0d0+alp*this%hhz));
    gz2 = (1.0d0-ro2z)/(1.0d0+ro2z*(1.0d0+alp*this%hhz));

    !===== EX =====
    ! Faces
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do i=1,this%Nx
          this%Ef%X(i, 0, k) = -(gy1+gy2)*this%Ef%X(i, 1, k)-gy1*gy2*this%Ef%X(i, 2, k)  &
          +(gy1+gy2)*this%Efold%X(i, 0, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%X(i, 1, k)+(gy1+gy2*betay)*this%Efold%X(i, 2, k)  &
          -gy1*gy2*this%Efold2%X(i, 0, k)-(gy1+gy2*betay)*this%Efold2%X(i, 1, k)-betay*this%Efold2%X(i, 2, k);

          this%Ef%X(i, this%Ny, k) = -(gy1+gy2)*this%Ef%X(i, this%Ny-1, k)-gy1*gy2*this%Ef%X(i, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%X(i, this%Ny, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%X(i, this%Ny-1, k)+(gy1+gy2*betay)*this%Efold%X(i, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%X(i, this%Ny, k)-(gy1+gy2*betay)*this%Efold2%X(i, this%Ny-1, k)-betay*this%Efold2%X(i, this%Ny-2, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=1,this%Nx
          this%Ef%X(i, j, 0) = -(gz1+gz2)*this%Ef%X(i, j, 1)-gz1*gz2*this%Ef%X(i, j, 2)  &
          +(gz1+gz2)*this%Efold%X(i, j, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%X(i, j, 1)+(gz1+gz2*betaz)*this%Efold%X(i, j, 2)  &
          -gz1*gz2*this%Efold2%X(i, j, 0)-(gz1+gz2*betaz)*this%Efold2%X(i, j, 1)-betaz*this%Efold2%X(i, j, 2);
          
          this%Ef%X(i, j, this%Nz) = -(gz1+gz2)*this%Ef%X(i, j, this%Nz-1)-gz1*gz2*this%Ef%X(i, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%X(i, j, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%X(i, j, this%Nz-1)+(gz1+gz2*betaz)*this%Efold%X(i, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%X(i, j, this%Nz)-(gz1+gz2*betaz)*this%Efold2%X(i, j, this%Nz-1)-betaz*this%Efold2%X(i, j, this%Nz-2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do i=1,this%Nx
       tm1 = -(gy1+gy2)*this%Ef%X(i, 1, 0)-gy1*gy2*this%Ef%X(i, 2, 0)  &
          +(gy1+gy2)*this%Efold%X(i, 0, 0)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%X(i, 1, 0)+(gy1+gy2*betay)*this%Efold%X(i, 2, 0)  &
          -gy1*gy2*this%Efold2%X(i, 0, 0)-(gy1+gy2*betay)*this%Efold2%X(i, 1, 0)-betay*this%Efold2%X(i, 2, 0);
       
       tm2 = -(gz1+gz2)*this%Ef%X(i, 0, 1)-gz1*gz2*this%Ef%X(i, 0, 2)  &
          +(gz1+gz2)*this%Efold%X(i, 0, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%X(i, 0, 1)+(gz1+gz2*betay)*this%Efold%X(i, 0, 2)  &
          -gz1*gz2*this%Efold2%X(i, 0, 0)-(gz1+gz2*betaz)*this%Efold2%X(i, 0, 1)-betaz*this%Efold2%X(i, 0, 2);
       this%Ef%X(i, 0, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, 0, 0) = dg*this%Ef%X(i, 0, 0) - dg*this%Ef%X(i, 1, 1) + this%Efold%X(i, 1, 1);
       
       tm1 = -(gy1+gy2)*this%Ef%X(i, this%Ny-1, 0)-gy1*gy2*this%Ef%X(i, this%Ny-2, 0)  &
          +(gy1+gy2)*this%Efold%X(i, this%Ny, 0)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%X(i, this%Ny-1, 0)+(gy1+gy2*betay)*this%Efold%X(i, this%Ny-2, 0)  &
          -gy1*gy2*this%Efold2%X(i, this%Ny, 0)-(gy1+gy2*betay)*this%Efold2%X(i, this%Ny-1, 0)-betay*this%Efold2%X(i, this%Ny-2, 0);
       
       tm2 = -(gz1+gz2)*this%Ef%X(i, this%Ny, 1)-gz1*gz2*this%Ef%X(i, this%Ny, 2)  &
          +(gz1+gz2)*this%Efold%X(i, this%Ny, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%X(i, this%Ny, 1)+(gz1+gz2*betaz)*this%Efold%X(i, this%Ny, 2)  &
          -gz1*gz2*this%Efold2%X(i, this%Ny, 0)-(gz1+gz2*betaz)*this%Efold2%X(i, this%Ny, 1)-betaz*this%Efold2%X(i, this%Ny, 2);
       this%Ef%X(i, this%Ny, 0) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, this%Ny, 0) = dg*this%Ef%X(i, this%Ny, 0) - dg*this%Ef%X(i, this%Ny-1, 1) + this%Efold%X(i, this%Ny-1, 1);

       tm1 = -(gy1+gy2)*this%Ef%X(i, 1, this%Nz)-gy1*gy2*this%Ef%X(i, 2, this%Nz)  &
          +(gy1+gy2)*this%Efold%X(i, 0, this%Nz)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%X(i, 1, this%Nz)+(gy1+gy2*betay)*this%Efold%X(i, 2, this%Nz)  &
          -gy1*gy2*this%Efold2%X(i, 0, this%Nz)-(gy1+gy2*betay)*this%Efold2%X(i, 1, this%Nz)-betay*this%Efold2%X(i, 2, this%Nz);
       
       tm2 = -(gz1+gz2)*this%Ef%X(i, 0, this%Nz-1)-gz1*gz2*this%Ef%X(i, 0, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%X(i, 0, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%X(i, 0, this%Nz-1)+(gz1+gz2*betaz)*this%Efold%X(i, 0, this%Nz-2)  &
          -gz1*gz2*this%Efold2%X(i, 0, this%Nz)-(gz1+gz2*betaz)*this%Efold2%X(i, 0, this%Nz-1)-betaz*this%Efold2%X(i, 0, this%Nz-2);
       this%Ef%X(i, 0, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, 0, this%Nz) = dg*this%Ef%X(i, 0, this%Nz) - dg*this%Ef%X(i, 1, this%Nz-1) + this%Efold%X(i, 1, this%Nz-1);

       tm1 = -(gy1+gy2)*this%Ef%X(i, this%Ny-1, this%Nz)-gy1*gy2*this%Ef%X(i, this%Ny-2, this%Nz)  &
          +(gy1+gy2)*this%Efold%X(i, this%Ny, this%Nz)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%X(i, this%Ny-1, this%Nz)+(gy1+gy2*betay)*this%Efold%X(i, this%Ny-2, this%Nz)  &
          -gy1*gy2*this%Efold2%X(i, this%Ny, this%Nz)-(gy1+gy2*betay)*this%Efold2%X(i, this%Ny-1, this%Nz)-betay*this%Efold2%X(i, this%Ny-2, this%Nz);
       
       tm2 = -(gz1+gz2)*this%Ef%X(i, this%Ny, this%Nz-1)-gz1*gz2*this%Ef%X(i, this%Ny, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%X(i, this%Ny, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%X(i, this%Ny, this%Nz-1)+(gz1+gz2*betaz)*this%Efold%X(i, this%Ny, this%Nz-2)  &
          -gz1*gz2*this%Efold2%X(i, this%Ny, this%Nz)-(gz1+gz2*betaz)*this%Efold2%X(i, this%Ny, this%Nz-1)-betaz*this%Efold2%X(i, this%Ny, this%Nz-2);
       this%Ef%X(i, this%Ny, this%Nz) = (tm1 + tm2)/2.0d0;
       !this%Ef%X(i, this%Ny, this%Nz) = dg*this%Ef%X(i, this%Ny, this%Nz) - dg*this%Ef%X(i, this%Ny-1, this%Nz-1) + this%Efold%X(i, this%Ny-1, this%Nz-1);
    enddo
    !$OMP END PARALLEL DO
    !===== EY =====
    ! Faces
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       do i=1,this%Nx-1
          this%Ef%Y(i, j, 0) = -(gz1+gz2)*this%Ef%Y(i, j, 1)-gz1*gz2*this%Ef%Y(i, j, 2)  &
          +(gz1+gz2)*this%Efold%Y(i, j, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%Y(i, j, 1)+(gz1+gz2*betaz)*this%Efold%Y(i, j, 2)  &
          -gz1*gz2*this%Efold2%Y(i, j, 0)-(gz1+gz2*betaz)*this%Efold2%Y(i, j, 1)-betaz*this%Efold2%Y(i, j, 2);

          this%Ef%Y(i, j, this%Nz) = -(gz1+gz2)*this%Ef%Y(i, j, this%Nz-1)-gz1*gz2*this%Ef%Y(i, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%Y(i, j, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%Y(i, j, this%Nz-1)+(gz1+gz2*betaz)*this%Efold%Y(i, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%Y(i, j, this%Nz)-(gz1+gz2*betaz)*this%Efold2%Y(i, j, this%Nz-1)-betaz*this%Efold2%Y(i, j, this%Nz-2);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=1,this%Ny
          this%Ef%Y(0, j, k) = -(gx1+gx2)*this%Ef%Y(1, j, k)-gx1*gx2*this%Ef%Y(2, j, k)  &
          +(gx1+gx2)*this%Efold%Y(0, j, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Y(1, j, k)+(gx1+gx2*betax)*this%Efold%Y(2, j, k)  &
          -gx1*gx2*this%Efold2%Y(0, j, k)-(gx1+gx2*betax)*this%Efold2%Y(1, j, k)-betax*this%Efold2%Y(2, j, k);
 
          this%Ef%Y(this%Nx, j, k) = -(gx1+gx2)*this%Ef%Y(this%Nx-1, j, k)-gx1*gx2*this%Ef%Y(this%Nx-2, j, k)  &
          +(gx1+gx2)*this%Efold%Y(this%Nx, j, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Y(this%Nx-1, j, k)+(gx1+gx2*betax)*this%Efold%Y(this%Nx-2, j, k)  &
          -gx1*gx2*this%Efold2%Y(this%Nx, j, k)-(gx1+gx2*betax)*this%Efold2%Y(this%Nx-1, j, k)-betax*this%Efold2%Y(this%Nx-2, j, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny
       tm1 = -(gx1+gx2)*this%Ef%Y(1, j, 0)-gx1*gx2*this%Ef%Y(2, j, 0)  &
          +(gx1+gx2)*this%Efold%Y(0, j, 0)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Y(1, j, 0)+(gx1+gx2*betax)*this%Efold%Y(2, j, 0)  &
          -gx1*gx2*this%Efold2%Y(0, j, 0)-(gx1+gx2*betax)*this%Efold2%Y(1, j, 0)-betax*this%Efold2%Y(2, j, 0);
       tm2 = -(gz1+gz2)*this%Ef%Y(0, j, 1)-gz1*gz2*this%Ef%Y(0, j, 2)  &
          +(gz1+gz2)*this%Efold%Y(0, j, 0)+((betax+1.0d0)+2.0d0*gz1*gz2)*this%Efold%Y(0, j, 1)+(gz1+gz2*betaz)*this%Efold%Y(0, j, 2)  &
          -gz1*gz2*this%Efold2%Y(0, j, 0)-(gz1+gz2*betaz)*this%Efold2%Y(0, j, 1)-betaz*this%Efold2%Y(0, j, 2);
       this%Ef%Y(0, j, 0) = (tm1 + tm2)/2;
       !this%Ef%Y(0, j, 0) = dg*this%Ef%Y(0, j, 0) - dg*this%Ef%Y(1, j, 1) + this%Efold%Y(1, j, 1);
       
       tm1 = -(gx1+gx2)*this%Ef%Y(this%Nx-1, j, 0)-gx1*gx2*this%Ef%Y(this%Nx-2, j, 0)  &
          +(gx1+gx2)*this%Efold%Y(this%Nx, j, 0)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Y(this%Nx-1, j, 0)+(gx1+gx2*betax)*this%Efold%Y(this%Nx-2, j, 0)  &
          -gx1*gx2*this%Efold2%Y(this%Nx, j, 0)-(gx1+gx2*betax)*this%Efold2%Y(this%Nx-1, j, 0)-betax*this%Efold2%Y(this%Nx-2, j, 0);
       tm2 = -(gz1+gz2)*this%Ef%Y(this%Nx, j, 1)-gz1*gz2*this%Ef%Y(this%Nx, j, 2)  &
          +(gz1+gz2)*this%Efold%Y(this%Nx, j, 0)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%Y(this%Nx, j, 1)+(gz1+gz2*betaz)*this%Efold%Y(this%Nx, j, 2)  &
          -gz1*gz2*this%Efold2%Y(this%Nx, j, 0)-(gz1+gz2*betaz)*this%Efold2%Y(this%Nx, j, 1)-betaz*this%Efold2%Y(this%Nx, j, 2);
       this%Ef%Y(this%Nx, j, 0) = (tm1 + tm2)/2;
       !this%Ef%Y(this%Nx, j, 0) = dg*this%Ef%Y(this%Nx, j, 0) - dg*this%Ef%Y(this%Nx-1, j, 1) + this%Efold%Y(this%Nx-1, j, 1);

       tm1 = -(gx1+gx2)*this%Ef%Y(1, j, this%Nz)-gx1*gx2*this%Ef%Y(2, j, this%Nz)  &
          +(gx1+gx2)*this%Efold%Y(0, j, this%Nz)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Y(1, j, this%Nz)+(gx1+gx2*betax)*this%Efold%Y(2, j, this%Nz)  &
          -gx1*gx2*this%Efold2%Y(0, j, this%Nz)-(gx1+gx2*betax)*this%Efold2%Y(1, j, this%Nz)-betax*this%Efold2%Y(2, j, this%Nz);
       tm2 = -(gz1+gz2)*this%Ef%Y(0, j, this%Nz-1)-gz1*gz2*this%Ef%Y(0, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%Y(0, j, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%Y(0, j, this%Nz-1)+(gz1+gz2*betaz)*this%Efold%Y(0, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%Y(0, j, this%Nz)-(gz1+gz2*betaz)*this%Efold2%Y(0, j, this%Nz-1)-betaz*this%Efold2%Y(0, j, this%Nz-2);     
       this%Ef%Y(0, j, this%Nz) = (tm1 + tm2)/2;
       !this%Ef%Y(0, j, this%Nz) = dg*this%Ef%Y(0, j, this%Nz) - dg*this%Ef%Y(1, j, this%Nz-1) + this%Efold%Y(1, j, this%Nz-1);   

       tm1 = -(gx1+gx2)*this%Ef%Y(this%Nx-1, j, this%Nz)-gx1*gx2*this%Ef%Y(this%Nx-2, j, this%Nz)  &
          +(gx1+gx2)*this%Efold%Y(this%Nx, j, this%Nz)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Y(this%Nx-1, j, this%Nz)+(gx1+gx2*betax)*this%Efold%Y(this%Nx-2, j, this%Nz)  &
          -gx1*gx2*this%Efold2%Y(this%Nx, j, this%Nz)-(gx1+gx2*betax)*this%Efold2%Y(this%Nx-1, j, this%Nz)-betax*this%Efold2%Y(this%Nx-2, j, this%Nz);
       tm2 = -(gz1+gz2)*this%Ef%Y(this%Nx, j, this%Nz-1)-gz1*gz2*this%Ef%Y(this%Nx, j, this%Nz-2)  &
          +(gz1+gz2)*this%Efold%Y(this%Nx, j, this%Nz)+((betaz+1.0d0)+2.0d0*gz1*gz2)*this%Efold%Y(this%Nx, j, this%Nz-1)+(gz1+gz2*betaz)*this%Efold%Y(this%Nx, j, this%Nz-2)  &
          -gz1*gz2*this%Efold2%Y(this%Nx, j, this%Nz)-(gz1+gz2*betaz)*this%Efold2%Y(this%Nx, j, this%Nz-1)-betaz*this%Efold2%Y(this%Nx, j, this%Nz-2);
       this%Ef%Y(this%Nx, j, this%Nz) = (tm1 + tm2)/2;
       !this%Ef%Y(this%Nx, j, this%Nz) = dg*this%Ef%Y(this%Nx, j, this%Nz) - dg*this%Ef%Y(this%Nx-1, j, this%Nz-1) + this%Efold%Y(this%Nx-1, j, this%Nz-1);
    enddo
    !$OMP END PARALLEL DO
    
    !===== EZ =====
    ! Faces
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do j=1,this%Ny-1
          this%Ef%Z(0, j, k) = -(gx1+gx2)*this%Ef%Z(1, j, k)-gx1*gx2*this%Ef%Z(2, j, k)  &
          +(gx1+gx2)*this%Efold%Z(0, j, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Z(1, j, k)+(gx1+gx2*betax)*this%Efold%Z(2, j, k)  &
          -gx1*gx2*this%Efold2%Z(0, j, k)-(gx1+gx2*betax)*this%Efold2%Z(1, j, k)-betax*this%Efold2%Z(2, j, k);
 
          this%Ef%Z(this%Nx, j, k) = -(gx1+gx2)*this%Ef%Z(this%Nx-1, j, k)-gx1*gx2*this%Ef%Z(this%Nx-2, j, k)  &
          +(gx1+gx2)*this%Efold%Z(this%Nx, j, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Z(this%Nx-1, j, k)+(gx1+gx2*betax)*this%Efold%Z(this%Nx-2, j, k)  &
          -gx1*gx2*this%Efold2%Z(this%Nx, j, k)-(gx1+gx2*betax)*this%Efold2%Z(this%Nx-1, j, k)-betax*this%Efold2%Z(this%Nx-2, j, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = -(gy1+gy2)*this%Ef%Z(i, 1, k)-gy1*gy2*this%Ef%Z(i, 2, k)  &
          +(gy1+gy2)*this%Efold%Z(i, 0, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%Z(i, 1, k)+(gy1+gy2*betay)*this%Efold%Z(i, 2, k)  &
          -gy1*gy2*this%Efold2%Z(i, 0, k)-(gy1+gy2*betay)*this%Efold2%Z(i, 1, k)-betay*this%Efold2%Z(i, 2, k);

          this%Ef%Z(i, this%Ny, k) = -(gy1+gy2)*this%Ef%Z(i, this%Ny-1, k)-gy1*gy2*this%Ef%Z(i, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%Z(i, this%Ny, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%Z(i, this%Ny-1, k)+(gy1+gy2*betay)*this%Efold%Z(i, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%Z(i, this%Ny, k)-(gy1+gy2*betay)*this%Efold2%Z(i, this%Ny-1, k)-betay*this%Efold2%Z(i, this%Ny-2, k);
       enddo
    enddo
    !$OMP END PARALLEL DO
    ! Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz
       tm1 = -(gx1+gx2)*this%Ef%Z(1, 0, k)-gx1*gx2*this%Ef%Z(2, 0, k)  &
          +(gx1+gx2)*this%Efold%Z(0, 0, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Z(1, 0, k)+(gx1+gx2*betax)*this%Efold%Z(2, 0, k)  &
          -gx1*gx2*this%Efold2%Z(0, 0, k)-(gx1+gx2*betax)*this%Efold2%Z(1, 0, k)-betax*this%Efold2%Z(2, 0, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(0, 1, k)-gy1*gy2*this%Ef%Z(0, 2, k)  &
          +(gy1+gy2)*this%Efold%Z(0, 0, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%Z(0, 1, k)+(gy1+gy2*betay)*this%Efold%Z(0, 2, k)  &
          -gy1*gy2*this%Efold2%Z(0, 0, k)-(gy1+gy2*betay)*this%Efold2%Z(0, 1, k)-betay*this%Efold2%Z(0, 2, k);
       this%Ef%Z(0, 0, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(0, 0, k) = gx*this%Ef%Z(0, 0, k) - gx*this%Ef%Z(1, 1, k) + this%Efold%Z(1, 1, k);

       tm1 = -(gx1+gx2)*this%Ef%Z(this%Nx-1, 0, k)-gx1*gx2*this%Ef%Z(this%Nx-2, 0, k)  &
          +(gx1+gx2)*this%Efold%Z(this%Nx, 0, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Z(this%Nx-1, 0, k)+(gx1+gx2*betax)*this%Efold%Z(this%Nx-2, 0, k)  &
          -gx1*gx2*this%Efold2%Z(this%Nx, 0, k)-(gx1+gx2*betax)*this%Efold2%Z(this%Nx-1, 0, k)-betax*this%Efold2%Z(this%Nx-2, 0, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(this%Nx, 1, k)-gy1*gy2*this%Ef%Z(this%Nx, 2, k)  &
          +(gy1+gy2)*this%Efold%Z(this%Nx, 0, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%Z(this%Nx, 1, k)+(gy1+gy2*betay)*this%Efold%Z(this%Nx, 2, k)  &
          -gy1*gy2*this%Efold2%Z(this%Nx, 0, k)-(gy1+gy2*betay)*this%Efold2%Z(this%Nx, 1, k)-betay*this%Efold2%Z(this%Nx, 2, k);
       this%Ef%Z(this%Nx, 0, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(this%Nx, 0, k) = gx*this%Ef%Z(this%Nx, 0, k) - gx*this%Ef%Z(this%Nx-1, 1, k) + this%Efold%Z(this%Nx-1, 1, k);

       tm1 = -(gx1+gx2)*this%Ef%Z(1, this%Ny, k)-gx1*gx2*this%Ef%Z(2, this%Ny, k)  &
          +(gx1+gx2)*this%Efold%Z(0, this%Ny, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Z(1, this%Ny, k)+(gx1+gx2*betax)*this%Efold%Z(2, this%Ny, k)  &
          -gx1*gx2*this%Efold2%Z(0, this%Ny, k)-(gx1+gx2*betax)*this%Efold2%Z(1, this%Ny, k)-betax*this%Efold2%Z(2, this%Ny, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(0, this%Ny-1, k)-gy1*gy2*this%Ef%Z(0, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%Z(0, this%Ny, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%Z(0, this%Ny-1, k)+(gy1+gy2*betay)*this%Efold%Z(0, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%Z(0, this%Ny, k)-(gy1+gy2*betay)*this%Efold2%Z(0, this%Ny-1, k)-betay*this%Efold2%Z(0, this%Ny-2, k);
       this%Ef%Z(0, this%Ny, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(0, this%Ny, k) = gx*this%Ef%Z(0, this%Ny, k) - gx*this%Ef%Z(1, this%Ny-1, k) + this%Efold%Z(1, this%Ny-1, k);

       tm1 = -(gx1+gx2)*this%Ef%Z(this%Nx-1, this%Ny, k)-gx1*gx2*this%Ef%Z(this%Nx-2, this%Ny, k)  &
          +(gx1+gx2)*this%Efold%Z(this%Nx, this%Ny, k)+((betax+1.0d0)+2.0d0*gx1*gx2)*this%Efold%Z(this%Nx-1, this%Ny, k)+(gx1+gx2*betax)*this%Efold%Z(this%Nx-2, this%Ny, k)  &
          -gx1*gx2*this%Efold2%Z(this%Nx, this%Ny, k)-(gx1+gx2*betax)*this%Efold2%Z(this%Nx-1, this%Ny, k)-betax*this%Efold2%Z(this%Nx-2, this%Ny, k);
       tm2 = -(gy1+gy2)*this%Ef%Z(this%Nx, this%Ny-1, k)-gy1*gy2*this%Ef%Z(this%Nx, this%Ny-2, k)  &
          +(gy1+gy2)*this%Efold%Z(this%Nx, this%Ny, k)+((betay+1.0d0)+2.0d0*gy1*gy2)*this%Efold%Z(this%Nx, this%Ny-1, k)+(gy1+gy2*betay)*this%Efold%Z(this%Nx, this%Ny-2, k)  &
          -gy1*gy2*this%Efold2%Z(this%Nx, this%Ny, k)-(gy1+gy2*betay)*this%Efold2%Z(this%Nx, this%Ny-1, k)-betay*this%Efold2%Z(this%Nx, this%Ny-2, k);
       this%Ef%Z(this%Nx, this%Ny, k) = (tm1 + tm2)/2.0d0;
       !this%Ef%Z(this%Nx, this%Ny, k) = gx*this%Ef%Z(this%Nx, this%Ny, k) - gx*this%Ef%Z(this%Nx-1, this%Ny-1, k) + this%Efold%Z(this%Nx-1, this%Ny-1, k);
    enddo
    !$OMP END PARALLEL DO
  end subroutine fillEBoundaryABCBetzMittra


!fillEBoundaryABCMur------------------------------------------------------------------
  subroutine  fillEBoundaryABCMur(this,t)
    !Fill E boundary with Mur ABC
    class(tProblem) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    real*8 :: gq1, gq2, gq3, gq4, gx, gy, gz;
    real*8 :: tm1, tm2, tm3;
    gx = (1.0d0-cc*this%ht/this%hhx)/(1.0d0+cc*this%ht/this%hhx);
    gy = (1.0d0-cc*this%ht/this%hhy)/(1.0d0+cc*this%ht/this%hhy);
    gz = (1.0d0-cc*this%ht/this%hhz)/(1.0d0+cc*this%ht/this%hhz);
    !================ EX ================= 
    !====== Faces
    !=On Y
    gq1 = (cc*this%ht-this%hhy)/(cc*this%ht+this%hhy);
    gq2 = 2.0d0*this%hhy/(cc*this%ht+this%hhy);
    gq3 = this%hhy*((cc*this%ht)**2) / (2.0d0*(this%hhx**2)*(cc*this%ht+this%hhy));
    gq4 = this%hhy*((cc*this%ht)**2) / (2.0d0*(this%hhz**2)*(cc*this%ht+this%hhy));
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do i=2,this%Nx-1
          this%Ef%X(i, 0, k) = - this%Efold2%X(i, 1, k) &
                               + gq1*(this%Ef%X(i, 1, k)+this%Efold2%X(i, 0, k)) &
                               + gq2*(this%Efold%X(i, 0, k)+this%Efold%X(i, 1, k)) &
                               + gq3*(this%Efold%X(i+1, 0, k)-2.0d0*this%Efold%X(i, 0, k)+this%Efold%X(i-1, 0, k)+this%Efold%X(i+1, 1, k)-2.0d0*this%Efold%X(i, 1, k)+this%Efold%X(i-1, 1, k)) &
                               + gq4*(this%Efold%X(i, 0, k+1)-2.0d0*this%Efold%X(i, 0, k)+this%Efold%X(i, 0, k-1)+this%Efold%X(i, 1, k+1)-2.0d0*this%Efold%X(i, 1, k)+this%Efold%X(i, 1, k-1))
          this%Ef%X(i, this%Ny, k) = - this%Efold2%X(i, this%Ny-1, k) &
                                     + gq1*(this%Ef%X(i, this%Ny-1, k)+this%Efold2%X(i, this%Ny, k)) &
                                     + gq2*(this%Efold%X(i, this%Ny, k)+this%Efold%X(i, this%Ny-1, k)) &
                                     + gq3*(this%Efold%X(i+1, this%Ny, k)-2.0d0*this%Efold%X(i, this%Ny, k)+this%Efold%X(i-1, this%Ny, k)+this%Efold%X(i+1, this%Ny-1, k)-2.0d0*this%Efold%X(i, this%Ny-1, k)+this%Efold%X(i-1, this%Ny-1, k)) &
                                     + gq4*(this%Efold%X(i, this%Ny, k+1)-2.0d0*this%Efold%X(i, this%Ny, k)+this%Efold%X(i, this%Ny, k-1)+this%Efold%X(i, this%Ny-1, k+1)-2.0d0*this%Efold%X(i, this%Ny-1, k)+this%Efold%X(i, this%Ny-1, k-1))
       enddo
    enddo
    !$OMP END PARALLEL DO
    !=On Z
    gq1 = (cc*this%ht-this%hhz)/(cc*this%ht+this%hhz);
    gq2 = 2.0d0*this%hhz/(cc*this%ht+this%hhz);
    gq3 = this%hhz*((cc*this%ht)**2) / (2.0d0*(this%hhx**2)*(cc*this%ht+this%hhz));
    gq4 = this%hhz*((cc*this%ht)**2) / (2.0d0*(this%hhy**2)*(cc*this%ht+this%hhz));
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1
       do i=2,this%Nx-1
          this%Ef%X(i, j, 0) = - this%Efold2%X(i, j, 1) &
                               + gq1*(this%Ef%X(i, j, 1)+this%Efold2%X(i, j, 0)) &
                               + gq2*(this%Efold%X(i, j, 0)+this%Efold%X(i, j, 1)) &
                               + gq3*(this%Efold%X(i+1, j, 0)-2.0d0*this%Efold%X(i, j, 0)+this%Efold%X(i-1, j, 0)+this%Efold%X(i+1, j, 1)-2.0d0*this%Efold%X(i, j, 1)+this%Efold%X(i-1, j, 1)) &
                               + gq4*(this%Efold%X(i, j+1, 0)-2.0d0*this%Efold%X(i, j, 0)+this%Efold%X(i, j-1, 0)+this%Efold%X(i, j+1, 1)-2.0d0*this%Efold%X(i, j, 1)+this%Efold%X(i, j-1, 1))
          this%Ef%X(i, j, this%Nz) = - this%Efold2%X(i, j, this%Nz-1) &
                                     + gq1*(this%Ef%X(i, j, this%Nz-1)+this%Efold2%X(i, j, this%Nz)) &
                                     + gq2*(this%Efold%X(i, j, this%Nz)+this%Efold%X(i, j, this%Nz-1)) &
                                     + gq3*(this%Efold%X(i+1, j, this%Nz)-2.0d0*this%Efold%X(i, j, this%Nz)+this%Efold%X(i-1, j, this%Nz)+this%Efold%X(i+1, j, this%Nz-1)-2.0d0*this%Efold%X(i, j, this%Nz-1)+this%Efold%X(i-1, j, this%Nz-1)) &
                                     + gq4*(this%Efold%X(i, j+1, this%Nz)-2.0d0*this%Efold%X(i, j, this%Nz)+this%Efold%X(i, j-1, this%Nz)+this%Efold%X(i, j+1, this%Nz-1)-2.0d0*this%Efold%X(i, j, this%Nz-1)+this%Efold%X(i, j-1, this%Nz-1))
       enddo
    enddo
    !$OMP END PARALLEL DO
    !====== Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do i=2,this%Nx-1 !YZ edges
       tm1 = gy*this%Ef%X(i, 0, 0) - gy*this%Ef%X(i, 1, 0) + this%Efold%X(i, 1, 0);
       tm2 = gz*this%Ef%X(i, 0, 0) - gz*this%Ef%X(i, 0, 1) + this%Efold%X(i, 0, 1);
       this%Ef%X(i, 0, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%X(i, this%Ny, 0) - gy*this%Ef%X(i, this%Ny-1, 0) + this%Efold%X(i, this%Ny-1, 0);
       tm2 = gz*this%Ef%X(i, this%Ny, 0) - gz*this%Ef%X(i, this%Ny, 1) + this%Efold%X(i, this%Ny, 1);
       this%Ef%X(i, this%Ny, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%X(i, 0, this%Nz) - gy*this%Ef%X(i, 1, this%Nz) + this%Efold%X(i, 1, this%Nz);
       tm2 = gz*this%Ef%X(i, 0, this%Nz) - gz*this%Ef%X(i, 0, this%Nz-1) + this%Efold%X(i, 0, this%Nz-1);
       this%Ef%X(i, 0, this%Nz) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%X(i, this%Ny, this%Nz) - gy*this%Ef%X(i, this%Ny-1, this%Nz) + this%Efold%X(i, this%Ny-1, this%Nz);
       tm2 = gz*this%Ef%X(i, this%Ny, this%Nz) - gz*this%Ef%X(i, this%Ny, this%Nz-1) + this%Efold%X(i, this%Ny, this%Nz-1);
       this%Ef%X(i, this%Ny, this%Nz) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1 !XZ edges
       tm1 = gx*this%Ef%X(1, j, 0) - gx*this%Ef%X(2, j, 0) + this%Efold%X(2, j, 0);
       tm2 = gz*this%Ef%X(1, j, 0) - gz*this%Ef%X(1, j, 1) + this%Efold%X(1, j, 1);
       this%Ef%X(1, j, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%X(this%Nx, j, 0) - gx*this%Ef%X(this%Nx-1, j, 0) + this%Efold%X(this%Nx-1, j, 0);
       tm2 = gz*this%Ef%X(this%Nx, j, 0) - gz*this%Ef%X(this%Nx, j, 1) + this%Efold%X(this%Nx, j, 1);
       this%Ef%X(this%Nx, j, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%X(1, j, this%Nz) - gx*this%Ef%X(2, j, this%Nz) + this%Efold%X(2, j, this%Nz);
       tm2 = gz*this%Ef%X(1, j, this%Nz) - gz*this%Ef%X(1, j, this%Nz-1) + this%Efold%X(1, j, this%Nz-1);
       this%Ef%X(1, j, this%Nz) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%X(this%Nx, j, this%Nz) - gx*this%Ef%X(this%Nx-1, j, this%Nz) + this%Efold%X(this%Nx-1, j, this%Nz);
       tm2 = gz*this%Ef%X(this%Nx, j, this%Nz) - gz*this%Ef%X(this%Nx, j, this%Nz-1) + this%Efold%X(this%Nx, j, this%Nz-1);
       this%Ef%X(this%Nx, j, this%Nz) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1 !XY edges
       tm1 = gx*this%Ef%X(1, 0, k) - gx*this%Ef%X(2, 0, k) + this%Efold%X(2, 0, k);
       tm2 = gy*this%Ef%X(1, 0, k) - gy*this%Ef%X(1, 1, k) + this%Efold%X(1, 1, k);
       this%Ef%X(1, 0, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%X(this%Nx, 0, k) - gx*this%Ef%X(this%Nx-1, 0, k) + this%Efold%X(this%Nx-1, 0, k);
       tm2 = gy*this%Ef%X(this%Nx, 0, k) - gy*this%Ef%X(this%Nx, 1, k) + this%Efold%X(this%Nx, 1, k);
       this%Ef%X(this%Nx, 0, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%X(1, this%Ny, k) - gx*this%Ef%X(2, this%Ny, k) + this%Efold%X(2, this%Ny, k);
       tm2 = gy*this%Ef%X(1, this%Ny, k) - gy*this%Ef%X(1, this%Ny-1, k) + this%Efold%X(1, this%Ny-1, k);
       this%Ef%X(1, this%Ny, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%X(this%Nx, this%Ny, k) - gx*this%Ef%X(this%Nx-1, this%Ny, k) + this%Efold%X(this%Nx-1, this%Ny, k);
       tm2 = gy*this%Ef%X(this%Nx, this%Ny, k) - gy*this%Ef%X(this%Nx, this%Ny-1, k) + this%Efold%X(this%Nx, this%Ny-1, k);
       this%Ef%X(this%Nx, this%Ny, k) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !==== Corners
       tm1 = gx*this%Ef%X(1, 0, 0) - gx*this%Ef%X(2, 0, 0) + this%Efold%X(2, 0, 0);
       tm2 = gy*this%Ef%X(1, 0, 0) - gy*this%Ef%X(1, 1, 0) + this%Efold%X(1, 1, 0);
       tm3 = gz*this%Ef%X(1, 0, 0) - gz*this%Ef%X(1, 0, 1) + this%Efold%X(1, 0, 1);
       this%Ef%X(1, 0, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(1, this%Ny, 0) - gx*this%Ef%X(2, this%Ny, 0) + this%Efold%X(2, this%Ny, 0);
       tm2 = gy*this%Ef%X(1, this%Ny, 0) - gy*this%Ef%X(1, this%Ny-1, 0) + this%Efold%X(1, this%Ny-1, 0);
       tm3 = gz*this%Ef%X(1, this%Ny, 0) - gz*this%Ef%X(1, this%Ny, 1) + this%Efold%X(1, this%Ny, 1);
       this%Ef%X(1, this%Ny, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(1, 0, this%Nz) - gx*this%Ef%X(2, 0, this%Nz) + this%Efold%X(2, 0, this%Nz);
       tm2 = gy*this%Ef%X(1, 0, this%Nz) - gy*this%Ef%X(1, 1, this%Nz) + this%Efold%X(1, 1, this%Nz);
       tm3 = gz*this%Ef%X(1, 0, this%Nz) - gz*this%Ef%X(1, 0, this%Nz-1) + this%Efold%X(1, 0, this%Nz-1);
       this%Ef%X(1, 0, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(1, this%Ny, this%Nz) - gx*this%Ef%X(2, this%Ny, this%Nz) + this%Efold%X(2, this%Ny, this%Nz);
       tm2 = gy*this%Ef%X(1, this%Ny, this%Nz) - gy*this%Ef%X(1, this%Ny-1, this%Nz) + this%Efold%X(1, this%Ny-1, this%Nz);
       tm3 = gz*this%Ef%X(1, this%Ny, this%Nz) - gz*this%Ef%X(1, this%Ny, this%Nz-1) + this%Efold%X(1, this%Ny, this%Nz-1);
       this%Ef%X(1, this%Ny, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(this%Nx, 0, 0) - gx*this%Ef%X(this%Nx-1, 0, 0) + this%Efold%X(this%Nx-1, 0, 0);
       tm2 = gy*this%Ef%X(this%Nx, 0, 0) - gy*this%Ef%X(this%Nx, 1, 0) + this%Efold%X(this%Nx, 1, 0);
       tm3 = gz*this%Ef%X(this%Nx, 0, 0) - gz*this%Ef%X(this%Nx, 0, 1) + this%Efold%X(this%Nx, 0, 1);
       this%Ef%X(this%Nx, 0, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(this%Nx, this%Ny, 0) - gx*this%Ef%X(this%Nx-1, this%Ny, 0) + this%Efold%X(this%Nx-1, this%Ny, 0);
       tm2 = gy*this%Ef%X(this%Nx, this%Ny, 0) - gy*this%Ef%X(this%Nx, this%Ny-1, 0) + this%Efold%X(this%Nx, this%Ny-1, 0);
       tm3 = gz*this%Ef%X(this%Nx, this%Ny, 0) - gz*this%Ef%X(this%Nx, this%Ny, 1) + this%Efold%X(this%Nx, this%Ny, 1);
       this%Ef%X(this%Nx, this%Ny, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(this%Nx, 0, this%Nz) - gx*this%Ef%X(this%Nx-1, 0, this%Nz) + this%Efold%X(this%Nx-1, 0, this%Nz);
       tm2 = gy*this%Ef%X(this%Nx, 0, this%Nz) - gy*this%Ef%X(this%Nx, 1, this%Nz) + this%Efold%X(this%Nx, 1, this%Nz);
       tm3 = gz*this%Ef%X(this%Nx, 0, this%Nz) - gz*this%Ef%X(this%Nx, 0, this%Nz-1) + this%Efold%X(this%Nx, 0, this%Nz-1);
       this%Ef%X(this%Nx, 0, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%X(this%Nx, this%Ny, this%Nz) - gx*this%Ef%X(this%Nx-1, this%Ny, this%Nz) + this%Efold%X(this%Nx-1, this%Ny, this%Nz);
       tm2 = gy*this%Ef%X(this%Nx, this%Ny, this%Nz) - gy*this%Ef%X(this%Nx, this%Ny-1, this%Nz) + this%Efold%X(this%Nx, this%Ny-1, this%Nz);
       tm3 = gz*this%Ef%X(this%Nx, this%Ny, this%Nz) - gz*this%Ef%X(this%Nx, this%Ny, this%Nz-1) + this%Efold%X(this%Nx, this%Ny, this%Nz-1);
       this%Ef%X(this%Nx, this%Ny, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;
    
    !================ EY ================= 
    !==== Faces
    !=On X
    gq1 = (cc*this%ht-this%hhx)/(cc*this%ht+this%hhx);
    gq2 = 2.0d0*this%hhx/(cc*this%ht+this%hhx);
    gq3 = this%hhx*((cc*this%ht)**2) / (2.0d0*(this%hhy**2)*(cc*this%ht+this%hhx));
    gq4 = this%hhx*((cc*this%ht)**2) / (2.0d0*(this%hhz**2)*(cc*this%ht+this%hhx));
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1
       do j=2,this%Ny-1
          this%Ef%Y(0, j, k) = - this%Efold2%Y(1, j, k) &
                               + gq1*(this%Ef%Y(1, j, k)+this%Efold2%Y(0, j, k)) &
                               + gq2*(this%Efold%Y(0, j, k)+this%Efold%Y(1, j, k)) &
                               + gq3*(this%Efold%Y(0, j+1, k)-2.0d0*this%Efold%Y(0, j, k)+this%Efold%Y(0, j-1, k)+this%Efold%Y(1, j+1, k)-2.0d0*this%Efold%Y(1, j, k)+this%Efold%Y(1, j-1, k)) &
                               + gq4*(this%Efold%Y(0, j, k+1)-2.0d0*this%Efold%Y(0, j, k)+this%Efold%Y(0, j, k-1)+this%Efold%Y(1, j, k+1)-2.0d0*this%Efold%Y(1, j, k)+this%Efold%Y(1, j, k-1))

          this%Ef%Y(this%Nx, j, k) = - this%Efold2%Y(this%Nx-1, j, k) &
                                     + gq1*(this%Ef%Y(this%Nx-1, j, k)+this%Efold2%Y(this%Nx, j, k)) &
                                     + gq2*(this%Efold%Y(this%Nx, j, k)+this%Efold%Y(this%Nx-1, j, k)) &
                                     + gq3*(this%Efold%Y(this%Nx, j+1, k)-2.0d0*this%Efold%Y(this%Nx, j, k)+this%Efold%Y(this%Nx, j-1, k)+this%Efold%Y(this%Nx-1, j+1, k)-2.0d0*this%Efold%Y(this%Nx-1, j, k)+this%Efold%Y(this%Nx-1, j-1, k)) &
                                     + gq4*(this%Efold%Y(this%Nx, j, k+1)-2.0d0*this%Efold%Y(this%Nx, j, k)+this%Efold%Y(this%Nx, j, k-1)+this%Efold%Y(this%Nx-1, j, k+1)-2.0d0*this%Efold%Y(this%Nx-1, j, k)+this%Efold%Y(this%Nx-1, j, k-1));
       enddo
    enddo
    !$OMP END PARALLEL DO
    !=On Z
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=2,this%Ny-1
       do i=1,this%Nx-1
          this%Ef%Y(i, j, 0) = - this%Efold2%Y(i, j, 1) &
                               + gq1*(this%Ef%Y(i, j, 1)+this%Efold2%Y(i, j, 0)) &
                               + gq2*(this%Efold%Y(i, j, 0)+this%Efold%Y(i, j, 1)) &
                               + gq3*(this%Efold%Y(i+1, j, 0)-2.0d0*this%Efold%Y(i, j, 0)+this%Efold%Y(i-1, j, 0)+this%Efold%Y(i+1, j, 1)-2.0d0*this%Efold%Y(i, j, 1)+this%Efold%Y(i-1, j, 1)) &
                               + gq4*(this%Efold%Y(i, j+1, 0)-2.0d0*this%Efold%Y(i, j, 0)+this%Efold%Y(i, j-1, 0)+this%Efold%Y(i, j+1, 1)-2.0d0*this%Efold%Y(i, j, 1)+this%Efold%Y(i, j-1, 1))
          this%Ef%Y(i, j, this%Nz) = - this%Efold2%Y(i, j, this%Nz-1) &
                                     + gq1*(this%Ef%Y(i, j, this%Nz-1)+this%Efold2%Y(i, j, this%Nz)) &
                                     + gq2*(this%Efold%Y(i, j, this%Nz)+this%Efold%Y(i, j, this%Nz-1)) &
                                     + gq3*(this%Efold%Y(i+1, j, this%Nz)-2.0d0*this%Efold%Y(i, j, this%Nz)+this%Efold%Y(i-1, j, this%Nz)+this%Efold%Y(i+1, j, this%Nz-1)-2.0d0*this%Efold%Y(i, j, this%Nz-1)+this%Efold%Y(i-1, j, this%Nz-1)) &
                                     + gq4*(this%Efold%Y(i, j+1, this%Nz)-2.0d0*this%Efold%Y(i, j, this%Nz)+this%Efold%Y(i, j-1, this%Nz)+this%Efold%Y(i, j+1, this%Nz-1)-2.0d0*this%Efold%Y(i, j, this%Nz-1)+this%Efold%Y(i, j-1, this%Nz-1))
       enddo
    enddo
    !$OMP END PARALLEL DO
    !====== Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do i=1,this%Nx-1 !YZ edges
       tm1 = gy*this%Ef%Y(i, 1, 0) - gy*this%Ef%Y(i, 2, 0) + this%Efold%Y(i, 2, 0);
       tm2 = gz*this%Ef%Y(i, 1, 0) - gz*this%Ef%Y(i, 1, 1) + this%Efold%Y(i, 1, 1);
       this%Ef%Y(i, 1, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%Y(i, this%Ny, 0) - gy*this%Ef%Y(i, this%Ny-1, 0) + this%Efold%Y(i, this%Ny-1, 0);
       tm2 = gz*this%Ef%Y(i, this%Ny, 0) - gz*this%Ef%Y(i, this%Ny, 1) + this%Efold%Y(i, this%Ny, 1);
       this%Ef%Y(i, this%Ny, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%Y(i, 1, this%Nz) - gy*this%Ef%Y(i, 2, this%Nz) + this%Efold%Y(i, 2, this%Nz);
       tm2 = gz*this%Ef%Y(i, 1, this%Nz) - gz*this%Ef%Y(i, 1, this%Nz-1) + this%Efold%Y(i, 1, this%Nz-1);
       this%Ef%Y(i, 1, this%Nz) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%Y(i, this%Ny, this%Nz) - gy*this%Ef%Y(i, this%Ny-1, this%Nz) + this%Efold%Y(i, this%Ny-1, this%Nz);
       tm2 = gz*this%Ef%Y(i, this%Ny, this%Nz) - gz*this%Ef%Y(i, this%Ny, this%Nz-1) + this%Efold%Y(i, this%Ny, this%Nz-1);
       this%Ef%Y(i, this%Ny, this%Nz) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=2,this%Ny-1 !XZ edges
       tm1 = gx*this%Ef%Y(0, j, 0) - gx*this%Ef%Y(1, j, 0) + this%Efold%Y(1, j, 0);
       tm2 = gz*this%Ef%Y(0, j, 0) - gz*this%Ef%Y(0, j, 1) + this%Efold%Y(0, j, 1);
       this%Ef%Y(0, j, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, j, 0) - gx*this%Ef%Y(this%Nx-1, j, 0) + this%Efold%Y(this%Nx-1, j, 0);
       tm2 = gz*this%Ef%Y(this%Nx, j, 0) - gz*this%Ef%Y(this%Nx, j, 1) + this%Efold%Y(this%Nx, j, 1);
       this%Ef%Y(this%Nx, j, 0) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Y(0, j, this%Nz) - gx*this%Ef%Y(1, j, this%Nz) + this%Efold%Y(1, j, this%Nz);
       tm2 = gz*this%Ef%Y(0, j, this%Nz) - gz*this%Ef%Y(0, j, this%Nz-1) + this%Efold%Y(0, j, this%Nz-1);
       this%Ef%Y(0, j, this%Nz) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, j, this%Nz) - gx*this%Ef%Y(this%Nx-1, j, this%Nz) + this%Efold%Y(this%Nx-1, j, this%Nz);
       tm2 = gz*this%Ef%Y(this%Nx, j, this%Nz) - gz*this%Ef%Y(this%Nx, j, this%Nz-1) + this%Efold%Y(this%Nx, j, this%Nz-1);
       this%Ef%Y(this%Nx, j, this%Nz) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=1,this%Nz-1 !XY edges
       tm1 = gx*this%Ef%Y(0, 1, k) - gx*this%Ef%Y(1, 1, k) + this%Efold%Y(1, 1, k);
       tm2 = gy*this%Ef%Y(0, 1, k) - gy*this%Ef%Y(0, 2, k) + this%Efold%Y(0, 2, k);
       this%Ef%Y(0, 1, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, 1, k) - gx*this%Ef%Y(this%Nx-1, 1, k) + this%Efold%Y(this%Nx-1, 1, k);
       tm2 = gy*this%Ef%Y(this%Nx, 1, k) - gy*this%Ef%Y(this%Nx, 2, k) + this%Efold%Y(this%Nx, 2, k);
       this%Ef%Y(this%Nx, 1, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Y(0, this%Ny, k) - gx*this%Ef%Y(1, this%Ny, k) + this%Efold%Y(1, this%Ny, k);
       tm2 = gy*this%Ef%Y(0, this%Ny, k) - gy*this%Ef%Y(0, this%Ny-1, k) + this%Efold%Y(0, this%Ny-1, k);
       this%Ef%Y(0, this%Ny, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, this%Ny, k) - gx*this%Ef%Y(this%Nx-1, this%Ny, k) + this%Efold%Y(this%Nx-1, this%Ny, k);
       tm2 = gy*this%Ef%Y(this%Nx, this%Ny, k) - gy*this%Ef%Y(this%Nx, this%Ny-1, k) + this%Efold%Y(this%Nx, this%Ny-1, k);
       this%Ef%Y(this%Nx, this%Ny, k) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !==== Corners
       tm1 = gx*this%Ef%Y(0, 1, 0) - gx*this%Ef%Y(1, 1, 0) + this%Efold%Y(1, 1, 0);
       tm2 = gy*this%Ef%Y(0, 1, 0) - gy*this%Ef%Y(0, 2, 0) + this%Efold%Y(0, 2, 0);
       tm3 = gz*this%Ef%Y(0, 1, 0) - gz*this%Ef%Y(0, 1, 1) + this%Efold%Y(0, 1, 1);
       this%Ef%Y(0, 1, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, 1, 0) - gx*this%Ef%Y(this%Nx-1, 1, 0) + this%Efold%Y(this%Nx-1, 1, 0);
       tm2 = gy*this%Ef%Y(this%Nx, 1, 0) - gy*this%Ef%Y(this%Nx, 2, 0) + this%Efold%Y(this%Nx, 2, 0);
       tm3 = gz*this%Ef%Y(this%Nx, 1, 0) - gz*this%Ef%Y(this%Nx, 1, 1) + this%Efold%Y(this%Nx, 1, 1);
       this%Ef%Y(this%Nx, 1, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Y(0, 1, this%Nz) - gx*this%Ef%Y(1, 1, this%Nz) + this%Efold%Y(1, 1, this%Nz);
       tm2 = gy*this%Ef%Y(0, 1, this%Nz) - gy*this%Ef%Y(0, 2, this%Nz) + this%Efold%Y(0, 2, this%Nz);
       tm3 = gz*this%Ef%Y(0, 1, this%Nz) - gz*this%Ef%Y(0, 1, this%Nz-1) + this%Efold%Y(0, 1, this%Nz-1);
       this%Ef%Y(0, 1, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, 1, this%Nz) - gx*this%Ef%Y(this%Nx-1, 1, this%Nz) + this%Efold%Y(this%Nx-1, 1, this%Nz);
       tm2 = gy*this%Ef%Y(this%Nx, 1, this%Nz) - gy*this%Ef%Y(this%Nx, 2, this%Nz) + this%Efold%Y(this%Nx, 2, this%Nz);
       tm3 = gz*this%Ef%Y(this%Nx, 1, this%Nz) - gz*this%Ef%Y(this%Nx, 1, this%Nz-1) + this%Efold%Y(this%Nx, 1, this%Nz-1);
       this%Ef%Y(this%Nx, 1, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;
       
       tm1 = gx*this%Ef%Y(0, this%Ny, 0) - gx*this%Ef%Y(1, this%Ny, 0) + this%Efold%Y(1, this%Ny, 0);
       tm2 = gy*this%Ef%Y(0, this%Ny, 0) - gy*this%Ef%Y(0, this%Ny-1, 0) + this%Efold%Y(0, this%Ny-1, 0);
       tm3 = gz*this%Ef%Y(0, this%Ny, 0) - gz*this%Ef%Y(0, this%Ny, 1) + this%Efold%Y(0, this%Ny, 1);
       this%Ef%Y(0, this%Ny, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, this%Ny, 0) - gx*this%Ef%Y(this%Nx-1, this%Ny, 0) + this%Efold%Y(this%Nx-1, this%Ny, 0);
       tm2 = gy*this%Ef%Y(this%Nx, this%Ny, 0) - gy*this%Ef%Y(this%Nx, this%Ny-1, 0) + this%Efold%Y(this%Nx, this%Ny-1, 0);
       tm3 = gz*this%Ef%Y(this%Nx, this%Ny, 0) - gz*this%Ef%Y(this%Nx, this%Ny, 1) + this%Efold%Y(this%Nx, this%Ny, 1);
       this%Ef%Y(this%Nx, this%Ny, 0) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Y(0, this%Ny, this%Nz) - gx*this%Ef%Y(1, this%Ny, this%Nz) + this%Efold%Y(1, this%Ny, this%Nz);
       tm2 = gy*this%Ef%Y(0, this%Ny, this%Nz) - gy*this%Ef%Y(0, this%Ny-1, this%Nz) + this%Efold%Y(0, this%Ny-1, this%Nz);
       tm3 = gz*this%Ef%Y(0, this%Ny, this%Nz) - gz*this%Ef%Y(0, this%Ny, this%Nz-1) + this%Efold%Y(0, this%Ny, this%Nz-1);
       this%Ef%Y(0, this%Ny, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Y(this%Nx, this%Ny, this%Nz) - gx*this%Ef%Y(this%Nx-1, this%Ny, this%Nz) + this%Efold%Y(this%Nx-1, this%Ny, this%Nz);
       tm2 = gy*this%Ef%Y(this%Nx, this%Ny, this%Nz) - gy*this%Ef%Y(this%Nx, this%Ny-1, this%Nz) + this%Efold%Y(this%Nx, this%Ny-1, this%Nz);
       tm3 = gz*this%Ef%Y(this%Nx, this%Ny, this%Nz) - gz*this%Ef%Y(this%Nx, this%Ny, this%Nz-1) + this%Efold%Y(this%Nx, this%Ny, this%Nz-1);
       this%Ef%Y(this%Nx, this%Ny, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;
       
    !================ EZ =================
    !==== Faces
    !=On X
    gq1 = (cc*this%ht-this%hhx)/(cc*this%ht+this%hhx);
    gq2 = 2.0d0*this%hhx/(cc*this%ht+this%hhx);
    gq3 = this%hhx*((cc*this%ht)**2) / (2.0d0*(this%hhy**2)*(cc*this%ht+this%hhx));
    gq4 = this%hhx*((cc*this%ht)**2) / (2.0d0*(this%hhz**2)*(cc*this%ht+this%hhx));
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=2,this%Nz-1
       do j=1,this%Ny-1
          this%Ef%Z(0, j, k) = - this%Efold2%Z(1, j, k) &
                               + gq1*(this%Ef%Z(1, j, k)+this%Efold2%Z(0, j, k)) &
                               + gq2*(this%Efold%Z(0, j, k)+this%Efold%Z(1, j, k)) &
                               + gq3*(this%Efold%Z(0, j+1, k)-2.0d0*this%Efold%Z(0, j, k)+this%Efold%Z(0, j-1, k)+this%Efold%Z(1, j+1, k)-2.0d0*this%Efold%Z(1, j, k)+this%Efold%Z(1, j-1, k)) &
                               + gq4*(this%Efold%Z(0, j, k+1)-2.0d0*this%Efold%Z(0, j, k)+this%Efold%Z(0, j, k-1)+this%Efold%Z(1, j, k+1)-2.0d0*this%Efold%Z(1, j, k)+this%Efold%Z(1, j, k-1))

          this%Ef%Z(this%Nx, j, k) = - this%Efold2%Z(this%Nx-1, j, k) &
                                     + gq1*(this%Ef%Z(this%Nx-1, j, k)+this%Efold2%Z(this%Nx, j, k)) &
                                     + gq2*(this%Efold%Z(this%Nx, j, k)+this%Efold%Z(this%Nx-1, j, k)) &
                                     + gq3*(this%Efold%Z(this%Nx, j+1, k)-2.0d0*this%Efold%Z(this%Nx, j, k)+this%Efold%Z(this%Nx, j-1, k)+this%Efold%Z(this%Nx-1, j+1, k)-2.0d0*this%Efold%Z(this%Nx-1, j, k)+this%Efold%Z(this%Nx-1, j-1, k)) &
                                     + gq4*(this%Efold%Z(this%Nx, j, k+1)-2.0d0*this%Efold%Z(this%Nx, j, k)+this%Efold%Z(this%Nx, j, k-1)+this%Efold%Z(this%Nx-1, j, k+1)-2.0d0*this%Efold%Z(this%Nx-1, j, k)+this%Efold%Z(this%Nx-1, j, k-1));
        enddo
     enddo
     !$OMP END PARALLEL DO
    !=On Y
    gq1 = (cc*this%ht-this%hhy)/(cc*this%ht+this%hhy);
    gq2 = 2.0d0*this%hhy/(cc*this%ht+this%hhy);
    gq3 = this%hhy*((cc*this%ht)**2) / (2.0d0*(this%hhx**2)*(cc*this%ht+this%hhy));
    gq4 = this%hhy*((cc*this%ht)**2) / (2.0d0*(this%hhz**2)*(cc*this%ht+this%hhy));
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=2,this%Nz-1
       do i=1,this%Nx-1
          this%Ef%Z(i, 0, k) = - this%Efold2%Z(i, 1, k) &
                               + gq1*(this%Ef%Z(i, 1, k)+this%Efold2%Z(i, 0, k)) &
                               + gq2*(this%Efold%Z(i, 0, k)+this%Efold%Z(i, 1, k)) &
                               + gq3*(this%Efold%Z(i+1, 0, k)-2.0d0*this%Efold%Z(i, 0, k)+this%Efold%Z(i-1, 0, k)+this%Efold%Z(i+1, 1, k)-2.0d0*this%Efold%Z(i, 1, k)+this%Efold%Z(i-1, 1, k)) &
                               + gq4*(this%Efold%Z(i, 0, k+1)-2.0d0*this%Efold%Z(i, 0, k)+this%Efold%Z(i, 0, k-1)+this%Efold%Z(i, 1, k+1)-2.0d0*this%Efold%Z(i, 1, k)+this%Efold%Z(i, 1, k-1))
          this%Ef%Z(i, this%Ny, k) = - this%Efold2%Z(i, this%Ny-1, k) &
                                     + gq1*(this%Ef%Z(i, this%Ny-1, k)+this%Efold2%Z(i, this%Ny, k)) &
                                     + gq2*(this%Efold%Z(i, this%Ny, k)+this%Efold%Z(i, this%Ny-1, k)) &
                                     + gq3*(this%Efold%Z(i+1, this%Ny, k)-2.0d0*this%Efold%Z(i, this%Ny, k)+this%Efold%Z(i-1, this%Ny, k)+this%Efold%Z(i+1, this%Ny-1, k)-2.0d0*this%Efold%Z(i, this%Ny-1, k)+this%Efold%Z(i-1, this%Ny-1, k)) &
                                     + gq4*(this%Efold%Z(i, this%Ny, k+1)-2.0d0*this%Efold%Z(i, this%Ny, k)+this%Efold%Z(i, this%Ny, k-1)+this%Efold%Z(i, this%Ny-1, k+1)-2.0d0*this%Efold%Z(i, this%Ny-1, k)+this%Efold%Z(i, this%Ny-1, k-1));
       enddo
    enddo
    !$OMP END PARALLEL DO
    !==== Edges
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do k=2,this%Nz-1 !XY edges
       tm1 = gx*this%Ef%Z(0, 0, k) - gx*this%Ef%Z(1, 0, k) + this%Efold%Z(1, 0, k);
       tm2 = gy*this%Ef%Z(0, 0, k) - gy*this%Ef%Z(0, 1, k) + this%Efold%Z(0, 1, k);
       this%Ef%Z(0, 0, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Z(this%Nx, 0, k) - gx*this%Ef%Z(this%Nx-1, 0, k) + this%Efold%Z(this%Nx-1, 0, k);
       tm2 = gy*this%Ef%Z(this%Nx, 0, k) - gy*this%Ef%Z(this%Nx, 1, k) + this%Efold%Z(this%Nx, 1, k);
       this%Ef%Z(this%Nx, 0, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Z(0, this%Ny, k) - gx*this%Ef%Z(1, this%Ny, k) + this%Efold%Z(1, this%Ny, k);
       tm2 = gy*this%Ef%Z(0, this%Ny, k) - gy*this%Ef%Z(0, this%Ny-1, k) + this%Efold%Z(0, this%Ny-1, k);
       this%Ef%Z(0, this%Ny, k) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Z(this%Nx, this%Ny, k) - gx*this%Ef%Z(this%Nx-1, this%Ny, k) + this%Efold%Z(this%Nx-1, this%Ny, k);
       tm2 = gy*this%Ef%Z(this%Nx, this%Ny, k) - gy*this%Ef%Z(this%Nx, this%Ny-1, k) + this%Efold%Z(this%Nx, this%Ny-1, k);
       this%Ef%Z(this%Nx, this%Ny, k) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do j=1,this%Ny-1 !XZ edges
       tm1 = gx*this%Ef%Z(0, j, 1) - gx*this%Ef%Z(1, j, 1) + this%Efold%Z(1, j, 1);
       tm2 = gz*this%Ef%Z(0, j, 1) - gz*this%Ef%Z(0, j, 2) + this%Efold%Z(0, j, 2);
       this%Ef%Z(0, j, 1) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Z(0, j, this%Nz) - gx*this%Ef%Z(1, j, this%Nz) + this%Efold%Z(1, j, this%Nz);
       tm2 = gz*this%Ef%Z(0, j, this%Nz) - gz*this%Ef%Z(0, j, this%Nz-1) + this%Efold%Z(0, j, this%Nz-1);
       this%Ef%Z(0, j, this%Nz) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Z(this%Nx, j, 1) - gx*this%Ef%Z(this%Nx-1, j, 1) + this%Efold%Z(this%Nx-1, j, 1);
       tm2 = gz*this%Ef%Z(this%Nx, j, 1) - gz*this%Ef%Z(this%Nx, j, 2) + this%Efold%Z(this%Nx, j, 2);
       this%Ef%Z(this%Nx, j, 1) = (tm1 + tm2)/2.0d0;

       tm1 = gx*this%Ef%Z(this%Nx, j, this%Nz) - gx*this%Ef%Z(this%Nx-1, j, this%Nz) + this%Efold%Z(this%Nx-1, j, this%Nz);
       tm2 = gz*this%Ef%Z(this%Nx, j, this%Nz) - gz*this%Ef%Z(this%Nx, j, this%Nz-1) + this%Efold%Z(this%Nx, j, this%Nz-1);
       this%Ef%Z(this%Nx, j, this%Nz) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SCHEDULE(GUIDED)
    do i=1,this%Nx-1 !YZ edges
       tm1 = gy*this%Ef%Z(i, 0, 1) - gy*this%Ef%Z(i, 1, 1) + this%Efold%Z(i, 1, 1);
       tm2 = gz*this%Ef%Z(i, 0, 1) - gz*this%Ef%Z(i, 0, 2) + this%Efold%Z(i, 0, 2);
       this%Ef%Z(i, 0, 1) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%Z(i, 0, this%Nz) - gy*this%Ef%Z(i, 1, this%Nz) + this%Efold%Z(i, 1, this%Nz);
       tm2 = gz*this%Ef%Z(i, 0, this%Nz) - gz*this%Ef%Z(i, 0, this%Nz-1) + this%Efold%Z(i, 0, this%Nz-1);
       this%Ef%Z(i, 0, this%Nz) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%Z(i, this%Ny, 1) - gy*this%Ef%Z(i, this%Ny-1, 1) + this%Efold%Z(i, this%Ny-1, 1);
       tm2 = gz*this%Ef%Z(i, this%Ny, 1) - gz*this%Ef%Z(i, this%Ny, 2) + this%Efold%Z(i, this%Ny, 2);
       this%Ef%Z(i, this%Ny, 1) = (tm1 + tm2)/2.0d0;

       tm1 = gy*this%Ef%Z(i, this%Ny, this%Nz) - gy*this%Ef%Z(i, this%Ny-1, this%Nz) + this%Efold%Z(i, this%Ny-1, this%Nz);
       tm2 = gz*this%Ef%Z(i, this%Ny, this%Nz) - gz*this%Ef%Z(i, this%Ny, this%Nz-1) + this%Efold%Z(i, this%Ny, this%Nz-1);
       this%Ef%Z(i, this%Ny,  this%Nz) = (tm1 + tm2)/2.0d0;
    enddo
    !$OMP END PARALLEL DO
    !==== Corners
       tm1 = gx*this%Ef%Z(0, 0, 1) - gx*this%Ef%Z(1, 0, 1) + this%Efold%Z(1, 0, 1);
       tm2 = gy*this%Ef%Z(0, 0, 1) - gy*this%Ef%Z(0, 1, 1) + this%Efold%Z(0, 1, 1);
       tm3 = gz*this%Ef%Z(0, 0, 1) - gz*this%Ef%Z(0, 0, 2) + this%Efold%Z(0, 0, 2);
       this%Ef%Z(0, 0, 1) = (tm1 + tm2 + tm3)/3.0d0;
       
       tm1 = gx*this%Ef%Z(0, this%Ny, 1) - gx*this%Ef%Z(1, this%Ny, 1) + this%Efold%Z(1, this%Ny, 1);
       tm2 = gy*this%Ef%Z(0, this%Ny, 1) - gy*this%Ef%Z(0, this%Ny-1, 1) + this%Efold%Z(0, this%Ny-1, 1);
       tm3 = gz*this%Ef%Z(0, this%Ny, 1) - gz*this%Ef%Z(0, this%Ny, 2) + this%Efold%Z(0, this%Ny, 2);
       this%Ef%Z(0, this%Ny, 1) = (tm1 + tm2 + tm3)/3.0d0;
       
       tm1 = gx*this%Ef%Z(this%Nx, 0, 1) - gx*this%Ef%Z(this%Nx-1, 0, 1) + this%Efold%Z(this%Nx-1, 0, 1);
       tm2 = gy*this%Ef%Z(this%Nx, 0, 1) - gy*this%Ef%Z(this%Nx, 1, 1) + this%Efold%Z(this%Nx, 1, 1);
       tm3 = gz*this%Ef%Z(this%Nx, 0, 1) - gz*this%Ef%Z(this%Nx, 0, 2) + this%Efold%Z(this%Nx, 0, 2);
       this%Ef%Z(this%Nx, 0, 1) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Z(this%Nx, this%Ny, 1) - gx*this%Ef%Z(this%Nx-1, this%Ny, 1) + this%Efold%Z(this%Nx-1, this%Ny, 1);
       tm2 = gy*this%Ef%Z(this%Nx, this%Ny, 1) - gy*this%Ef%Z(this%Nx, this%Ny-1, 1) + this%Efold%Z(this%Nx, this%Ny-1, 1);
       tm3 = gz*this%Ef%Z(this%Nx, this%Ny, 1) - gz*this%Ef%Z(this%Nx, this%Ny, 2) + this%Efold%Z(this%Nx, this%Ny, 2);
       this%Ef%Z(this%Nx, this%Ny, 1) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Z(0, 0, this%Nz) - gx*this%Ef%Z(1, 0, this%Nz) + this%Efold%Z(1, 0, this%Nz);
       tm2 = gy*this%Ef%Z(0, 0, this%Nz) - gy*this%Ef%Z(0, 1, this%Nz) + this%Efold%Z(0, 1, this%Nz);
       tm3 = gz*this%Ef%Z(0, 0, this%Nz) - gz*this%Ef%Z(0, 0, this%Nz-1) + this%Efold%Z(0, 0, this%Nz-1);
       this%Ef%Z(0, 0, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;
       
       tm1 = gx*this%Ef%Z(0, this%Ny, this%Nz) - gx*this%Ef%Z(1, this%Ny, this%Nz) + this%Efold%Z(1, this%Ny, this%Nz);
       tm2 = gy*this%Ef%Z(0, this%Ny, this%Nz) - gy*this%Ef%Z(0, this%Ny-1, this%Nz) + this%Efold%Z(0, this%Ny-1, this%Nz);
       tm3 = gz*this%Ef%Z(0, this%Ny, this%Nz) - gz*this%Ef%Z(0, this%Ny, this%Nz-1) + this%Efold%Z(0, this%Ny, this%Nz-1);
       this%Ef%Z(0, this%Ny, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;
       
       tm1 = gx*this%Ef%Z(this%Nx, 0, this%Nz) - gx*this%Ef%Z(this%Nx-1, 0, this%Nz) + this%Efold%Z(this%Nx-1, 0, this%Nz);
       tm2 = gy*this%Ef%Z(this%Nx, 0, this%Nz) - gy*this%Ef%Z(this%Nx, 1, this%Nz) + this%Efold%Z(this%Nx, 1, this%Nz);
       tm3 = gz*this%Ef%Z(this%Nx, 0, this%Nz) - gz*this%Ef%Z(this%Nx, 0, this%Nz-1) + this%Efold%Z(this%Nx, 0, this%Nz-1);
       this%Ef%Z(this%Nx, 0, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;

       tm1 = gx*this%Ef%Z(this%Nx, this%Ny, this%Nz) - gx*this%Ef%Z(this%Nx-1, this%Ny, this%Nz) + this%Efold%Z(this%Nx-1, this%Ny, this%Nz);
       tm2 = gy*this%Ef%Z(this%Nx, this%Ny, this%Nz) - gy*this%Ef%Z(this%Nx, this%Ny-1, this%Nz) + this%Efold%Z(this%Nx, this%Ny-1, this%Nz);
       tm3 = gz*this%Ef%Z(this%Nx, this%Ny, this%Nz) - gz*this%Ef%Z(this%Nx, this%Ny, this%Nz-1) + this%Efold%Z(this%Nx, this%Ny, this%Nz-1);
       this%Ef%Z(this%Nx, this%Ny, this%Nz) = (tm1 + tm2 + tm3)/3.0d0;
  end subroutine fillEBoundaryABCMur
