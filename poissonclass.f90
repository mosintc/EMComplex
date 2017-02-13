module poissonclass
  ! Poisson problem class
  use commonvars
  use meshclass
  use parallel
  !use fftsg3d
  
  implicit none;
  
  type, public :: tPoisson
     ! Buffer
     type(TMesh) :: Ef;
     type(TMesh) :: Eps, Ec1, Ec0;
     ! Neumann problem parameters
     integer :: npx, npy, npz;
     integer :: Nin;
     real*8, allocatable :: rhsE(:,:,:) !for E- Poisson
     real*8, allocatable ::  phiE(:,:,:)
     real*8, allocatable ::  lambdax(:), lambday(:), lambdaz(:);! electric
     ! Problem parameters
     integer :: Nx, Ny, Nz;
     real*8 :: hhx, hhy, hhz, ht;
   contains
      procedure poisson_init;
      procedure poisson_destroy;
      procedure Esolver;
      procedure gradientEcalc;
      procedure Ecpmst1_formation;
      procedure Ecmpst1_to_Ecmpst0;
      procedure divCalc;
   end type tPoisson

contains

!poisson_init------------------------------------------------------------------   
    subroutine poisson_Init(this, Nx, Ny, Nz, hhx, hhy, hhz, ht, pid, npx, npy, npz)
    !Initialize the Poisson problem
      class(tPoisson) :: this
      integer :: Nx, Ny, Nz;
      real*8 :: hhx, hhy, hhz, ht;
      integer :: pid;
      integer :: npx, npy, npz;
      integer mx, my, mz, i, j, k
      this%npx = npx;
      this%npy = npy;
      this%npz = npz;
      this%Nx = Nx;
      this%Ny = Ny;
      this%Nz = Nz;
      this%hhx = hhx;
      this%hhy = hhy;
      this%hhz = hhz;
      this%ht = ht;
      this%Nin = (this%Nx-this%npx)/2;
      call this%Ef%mesh_init(this%Nx,this%Ny,this%Nz,0,1000+pid);
      call this%Eps%mesh_init(this%npx,this%npy,this%npz,0,1100+pid);
      call this%Ec0%mesh_init(this%Nx,this%Ny,this%Nz,0,1200+pid);
      call this%Ec1%mesh_init(this%Nx,this%Ny,this%Nz,0,1300+pid);
      allocate(this%rhsE(1:this%npx+1,1:this%npy+1,1:this%npz+1));
      allocate(this%phiE(0:this%npx,0:this%npy,0:this%npz));
      do i=1,this%npx+1
         do j=1,this%npy+1
            do k=1,this%npz+1
               this%rhsE(i,j,k)=0.0d0;
               this%phiE(i-1,j-1,k-1)=0.0d0;
            enddo
         enddo
      enddo
      allocate(this%lambdax(0:this%npx));
      allocate(this%lambday(0:this%npy));
      allocate(this%lambdaz(0:this%npz));
      !write(*,*) 'Lambdax-----'
      do mx=0,this%npx ! 
         this%lambdax(mx)=-((2.0D0*dsin(0.5D0*Pi*mx/(this%npx+1)))/this%hhx)**2
      !   write(*,*) this%lambdax(mx)
      end do
      !write(*,*) 'Lambday-----'
      do my=0,this%npy !
         this%lambday(my)=-((2.0D0*dsin(0.5D0*Pi*my/(this%npy+1)))/this%hhy)**2
      !   write(*,*) this%lambdax(my)
      end do
      !write(*,*) 'Lambdaz-----'
      do mz=0,this%npz !
         this%lambdaz(mz)=-((2.0D0*dsin(0.5D0*Pi*mz/(this%npz+1)))/this%hhz)**2
      !   write(*,*) this%lambdax(mz)
      end do
    end subroutine poisson_Init;

    
!poisson_destroy------------------------------------------------------------------      
    subroutine poisson_Destroy(this)
    !Destroy the Poisson problem
      class(tPoisson) :: this
      call this%Ef%mesh_destroy;
      call this%Eps%mesh_destroy;
      call this%Ec0%mesh_destroy;
      call this%Ec1%mesh_destroy;
      deallocate(this%rhsE, this%phiE);
    end subroutine poisson_Destroy;
    

!Esolver------------------------------------------------------------------      
    subroutine Esolver(this,t)
    ! Calculate E currents with Poisson problem
      class(tPoisson) :: this
      integer :: t;
      integer i,k,j,mx,my,mz,j1,j2,j3
      integer nmax,nmaxsqrt,n1,n2,n3
      parameter (nmax = 2048) !128
      parameter (nmaxsqrt = 256)
      integer ips(0 : nmaxsqrt + 1)
      real*8  tps(0:8*nmax-1),wps(0:nmax*3/2-1);
      real*8 Deltapsix,Deltapsiy,Deltapsiz
      real*8 lambxyzsum
      real*8 rhsEk(0:this%npx,0:this%npy,0:this%npz)
      real*8 psiE(0:this%npx+2,0:this%npy+2,0:this%npz+2)
      
      ips(0) = 0;
      psiE(:,:,:)=0.0d0;
      psiE(this%npx+2,1:this%npy+1,1:this%npz+1)=-this%Ef%X(this%npx+1+this%Nin,0+this%Nin:this%npy+this%Nin,0+this%Nin:this%npz+this%Nin)*this%hhx;
      psiE(0,1:this%npy+1,1:this%npz+1)= this%Ef%X(0+this%Nin,0+this%Nin:this%npy+this%Nin,0+this%Nin:this%npz+this%Nin)*this%hhx;
      psiE(1:this%npx+1,this%npy+2,1:this%npz+1)= -this%Ef%Y(0+this%Nin:this%npx+this%Nin,this%npy+1+this%Nin,0+this%Nin:this%npz+this%Nin)*this%hhy;
      psiE(1:this%npx+1,0,1:this%npz+1)= this%Ef%Y(0+this%Nin:this%npx+this%Nin,0+this%Nin,0+this%Nin:this%npz+this%Nin)*this%hhy;
      psiE(1:this%npx+1,1:this%npy+1,this%npz+2)=-this%Ef%Z(0+this%Nin:this%npx+this%Nin,0+this%Nin:this%npy+this%Nin,this%npz+1+this%Nin)*this%hhz;
      psiE(1:this%npx+1,1:this%npy+1,0)= this%Ef%Z(0+this%Nin:this%npx+this%Nin,0+this%Nin:this%npy+this%Nin,0+this%Nin)*this%hhz;

      do i=1,this%npx+1 !1,n+1 originally
         do j=1,this%npy+1
            do k=1,this%npz+1
               Deltapsix=(psiE(i-1,j,k)+psiE(i+1,j,k)-2.0d0*psiE(i,j,k))/(this%hhx**2);
               Deltapsiy=(psiE(i,j-1,k)+psiE(i,j+1,k)-2.0d0*psiE(i,j,k))/(this%hhy**2);
               Deltapsiz=(psiE(i,j,k-1)+psiE(i,j,k+1)-2.0d0*psiE(i,j,k))/(this%hhz**2);
               this%rhsE(i,j,k)=-(Deltapsix+Deltapsiy+Deltapsiz)
            enddo
         enddo  
      end do

      rhsEk(0:this%npx,0:this%npy,0:this%npz)=this%rhsE(1:this%npx+1,1:this%npy+1,1:this%npz+1);

      call ddct3d(this%npx+1,this%npy+1,this%npz+1,this%npy+1,this%npz+1,-1,rhsEk(0:this%npx,0:this%npy,0:this%npz),tps,ips,wps);

        
      do mx=0,this%npx !1,n+1 originally
         do my=0,this%npy
            do mz=0,this%npz
               if ((mx==0).AND.(my==0).AND.(mz==0)) then
                  rhsEk(0,0,0)=10.0d0
               else
                  lambxyzsum=this%lambdax(mx)+this%lambday(my)+this%lambdaz(mz)
                  rhsEk(mx,my,mz)=rhsEk(mx,my,mz)/lambxyzsum
               endif
            enddo
         enddo
      enddo

      do j3 = 0, this%npz ! - 1
         do j2 = 0, this%npy ! - 1
            rhsEk(0, j2, j3) = rhsEk(0, j2, j3)*0.5d0
         end do
         do j1 = 0, this%npx !- 1
            rhsEk(j1,0,j3)=rhsEk(j1,0,j3)*0.5d0
         end do
      end do
      do j2 = 0, this%npy ! - 1
         do j1 = 0, this%npx !- 1
            rhsEk(j1,j2,0) = rhsEk(j1,j2,0)*0.5d0
         end do
      end do

      call ddct3d(this%npx+1,this%npy+1,this%npx+1,this%npy+1,this%npz+1,1,rhsEk(0:this%npx,0:this%npy,0:this%npz),tps,ips,wps);

      do j3 = 0, this%npz
         do j2 = 0, this%npy
            do j1 = 0, this%npx
               rhsEk(j1,j2,j3) = rhsEk(j1,j2,j3)*(8.0d0/dfloat(this%npx+1)/dfloat(this%npy+1)/dfloat(this%npz+1))
            end do
         end do
      end do

      this%phiE(0:this%npx,0:this%npy,0:this%npz)=rhsEk(0:this%npx,0:this%npy,0:this%npz);

    end subroutine Esolver;


!gradientEcalc----------------------------------------------------------------
    subroutine gradientEcalc(this,t)
    ! Calculate E field components as a gradient of potential
      class(tPoisson) :: this
      integer :: t;
      integer :: ix, iy, iz;
      do ix=1,this%npx
         this%Eps%X(ix,0:this%npy,0:this%npz)=-(this%phiE(ix,0:this%npy,0:this%npz)-this%phiE(ix-1,0:this%npy,0:this%npz))/this%hhx;
      enddo
      do iy=1,this%npy
         this%Eps%Y(0:this%npx,iy,0:this%npz)=-(this%phiE(0:this%npx,iy,0:this%npz)-this%phiE(0:this%npx,iy-1,0:this%npz))/this%hhy;
      enddo
      do iz=1,this%npz
         this%Eps%Z(0:this%npx,0:this%npy,iz)=-(this%phiE(0:this%npx,0:this%npy,iz)-this%phiE(0:this%npx,0:this%npy,iz-1))/this%hhz;
      enddo      
    end subroutine gradientEcalc


!Ecpmst1_formation------------------------------------------------------------    
    subroutine  Ecpmst1_formation(this,t)
    ! Composite field formation
      class(tPoisson) :: this;
      integer t;
      integer :: ix, iy, iz;
      this%Ec1%X(1:this%Nx,0:this%Ny,0:this%Nz)=this%Ef%X(1:this%Nx,0:this%Ny,0:this%Nz)
      this%Ec1%Y(0:this%Nx,1:this%Ny,0:this%Nz)=this%Ef%Y(0:this%Nx,1:this%Ny,0:this%Nz)
      this%Ec1%Z(0:this%Nx,0:this%Ny,1:this%Nz)=this%Ef%Z(0:this%Nx,0:this%Ny,1:this%Nz)
      this%Ec1%X(1+this%Nin:this%npx+this%Nin,0+this%Nin:this%npy+this%Nin,0+this%Nin:this%npz+this%Nin)=this%Eps%X(1:this%npx,0:this%npy,0:this%npz);
      this%Ec1%Y(0+this%Nin:this%npx+this%Nin,1+this%Nin:this%npy+this%Nin,0+this%Nin:this%npz+this%Nin)=this%Eps%Y(0:this%npx,1:this%npy,0:this%npz);
      this%Ec1%Z(0+this%Nin:this%npx+this%Nin,0+this%Nin:this%npy+this%Nin,1+this%Nin:this%npz+this%Nin)=this%Eps%Z(0:this%npx,0:this%npy,1:this%npz);
      call this%divCalc(t);
    end subroutine  Ecpmst1_formation


!DivCalc------------------------------------------------------------    
    subroutine divCalc(this,t)
    ! Calculate divergence of the composite field
      class(tPoisson) :: this;
      integer t;
      integer :: ix, iy, iz;
      real*8, allocatable :: diE(:,:,:)
      real*8 ::  maxdiv
      integer :: mx, my, mz;
      allocate(diE(1:this%Nx-1,1:this%Ny-1,1:this%Nz-1));
      maxdiv = 0;
      diE(:,:,:) = 0.0d0;
      do ix=1,this%Nx-1
         do iy=1,this%Ny-1
            do iz=1,this%Nz-1
               diE(ix,iy,iz)=(this%Ec1%X(ix+1,iy,iz)-this%Ec1%X(ix,iy,iz))/this%hhx+(this%Ec1%Y(ix,iy+1,iz)-this%Ec1%Y(ix,iy,iz))/this%hhy+(this%Ec1%Z(ix,iy,iz+1)-this%Ec1%Z(ix,iy,iz))/this%hhz;
               if (diE(ix,iy,iz)>maxdiv) then
                  maxdiv = diE(ix,iy,iz);
                  mx = ix;
                  my = iy;
                  mz = iz;
               endif   
            enddo
         enddo
      enddo
      ! write(*,*) 'DivE Composite:', maxval(abs(diE(:,:,:)));
      if (file_compositedivergence==1) then
         if (file_writelocations==1) then
            write(unit=300,fmt=10) maxval(abs(diE(:,:,:))), maxloc(abs(diE(:,:,:)));
         else
            write(unit=300,fmt=20) maxval(abs(diE(:,:,:)));
         endif
         10 format(E22.15,',',I2,',',I2,',',I2);
         20 format(E22.15);
      endif
      deallocate(diE);
    end subroutine divCalc

    
!Ecpmst1_formation------------------------------------------------------------   
    subroutine  Ecmpst1_to_Ecmpst0(this)
    ! Copy Ec1 to Ec0
      class(tPoisson) :: this;
      integer :: i,j,k;
      do k=0,this%Nz
         do j=0,this%Ny
            do i=1,this%Nx
               this%Ec0%X(i,j,k)=this%Ec1%X(i,j,k);
            enddo
         enddo
      enddo
      do k=0,this%Nz
         do j=1,this%Ny
            do i=0,this%Nx
               this%Ec0%Y(i,j,k)=this%Ec1%Y(i,j,k);
            enddo
         enddo
      enddo
      do k=1,this%Nz
         do j=0,this%Ny
            do i=0,this%Nx
               this%Ec0%Z(i,j,k)=this%Ec1%Z(i,j,k);
            enddo
         enddo
      enddo
     end subroutine  Ecmpst1_to_Ecmpst0
    
end module poissonclass
