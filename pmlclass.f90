module pmlclass
! PML class that holds all the preferences of the PML 
  use commonvars
  use meshclass
  use parallel
   implicit none;
   type, public :: tPML
     integer :: pml_id; 
     integer :: EboundarytypePML;
     integer :: HboundarytypePML;
      ! 0 - zeros on the bound
      ! 1 - analytical solution as a boundary condition
      ! 2 - scheme
     integer :: NumPml;                ! Thickness of the PML in points
     integer :: Nx, Ny, Nz;            ! Number of points in PML at all
     real*8  :: hhx, hhy, hhz, ht;     ! differentials
     real*8 :: Apml;                   ! Dissipation factor of the MPL
     real*8, allocatable :: xi(:), yi(:), zi(:);
     real*8, allocatable :: xi05(:), yi05(:), zi05(:);
     real*8, allocatable :: sgmxi(:), sgmyi(:), sgmzi(:);
     real*8, allocatable :: sgmxi05(:), sgmyi05(:), sgmzi05(:);
     type(tMesh) Ef, Hf;
     type(tMesh) Epml, Hpml;
     type(tMesh) Ean, Han;
     type(tMesh) JfE, JfH;
     real*8 :: err, abserr;
     type(tFieldError) :: Ex_err, Ey_err, Ez_err, Hx_err, Hy_err, Hz_err;
     character(len=2) :: err_f;
    contains
      procedure pml_init;
      procedure pml_destroy;
      procedure pml_initphyspoints;
      procedure pml_initsigmaarrays;
      procedure pml_dostep
      ! Update functions
      procedure updateEInteriorPML
      procedure updateHInteriorPML
      procedure updateEInteriorPMLeverywhere;
      procedure updateHInteriorPMLeverywhere;
      procedure fillEBoundaryPML;
      procedure fillHBoundaryPML;
      ! Point functions
      procedure updateEpointPML;
      procedure updateHpointPML;
      procedure fillEBoundaryPointPML;
      procedure fillHBoundaryPointPML;
      ! Sigma functions
      procedure sgmpml_x;
      procedure sgmpml_y;
      procedure sgmpml_z;
  end type tPML

contains

!PML_init------------------------------------------------------------------   
  subroutine PML_init(this, Px, Py, Pz, nhhx, nhhy, nhhz, nht, thickness, PMLparam, new_id)
  ! Initialize PML  
    class(tPML) :: this
    integer, intent(in) :: Px, Py, Pz;
    real*8, intent(in) :: nhhx, nhhy, nhhz, nht;
    integer, intent(in) :: thickness;
    real*8, intent(in) :: PMLParam;
    integer, intent(in) :: new_id;
    this%pml_id = new_id;
    this%EboundarytypePML = 0;
    this%HboundarytypePML = 0;
    this%Nx = Px + 2*thickness;
    this%Ny = Py + 2*thickness;
    this%Nz = Pz + 2*thickness;
    this%hhx = nhhx;
    this%hhy = nhhy;
    this%hhz = nhhz;
    this%ht = nht;
    this%NumPML = thickness;
    this%Apml = PMLparam;
    call this%Ef%mesh_init(this%Nx,this%Ny,this%Nz,0,this%pml_id);
    call this%Hf%mesh_init(this%Nx,this%Ny,this%Nz,1,this%pml_id);
    call this%Epml%mesh_init(this%Nx,this%Ny,this%Nz,0,this%pml_id);
    call this%Hpml%mesh_init(this%Nx,this%Ny,this%Nz,1,this%pml_id);
    call this%JfE%mesh_init(this%Nx,this%Ny,this%Nz,0,this%pml_id);
    call this%JfH%mesh_init(this%Nx,this%Ny,this%Nz,1,this%pml_id);
    call this%Ean%mesh_init(this%Nx,this%Ny,this%Nz,0,this%pml_id);
    call this%Han%mesh_init(this%Nx,this%Ny,this%Nz,1,this%pml_id);
    allocate(this%xi(0:this%Nx));
    allocate(this%yi(0:this%Ny));
    allocate(this%zi(0:this%Nz));
    allocate(this%xi05(0:this%Nx+1));
    allocate(this%yi05(0:this%Ny+1));
    allocate(this%zi05(0:this%Nz+1));
    allocate(this%sgmxi(0:this%Nx+1));
    allocate(this%sgmyi(0:this%Ny+1));
    allocate(this%sgmzi(0:this%Nz+1));
    allocate(this%sgmxi05(0:this%Nx+1));
    allocate(this%sgmyi05(0:this%Ny+1));
    allocate(this%sgmzi05(0:this%Nz+1));
    call this%pml_initPhysPoints;
    call this%pml_initSigmaArrays;
  end subroutine PML_init

  
!PML_destroy------------------------------------------------------------------     
    subroutine PML_Destroy(this)
    !Destroy the problem structures and free the memory  
      class(tPML) :: this
      call this%Ef%mesh_destroy;
      call this%Hf%mesh_destroy;
      call this%Epml%mesh_destroy;
      call this%Hpml%mesh_destroy;
      call this%JfE%mesh_destroy;
      call this%JfH%mesh_destroy;
      call this%Ean%mesh_destroy;
      call this%Han%mesh_destroy;
      deallocate(this%xi,this%yi,this%zi, this%xi05, this%yi05, this%zi05);
      deallocate(this%sgmxi,this%sgmyi,this%sgmzi, this%sgmxi05, this%sgmyi05, this%sgmzi05);
    end subroutine PML_Destroy;
    

!PML_initPhysPoints-----------------------------------------------------------      
    subroutine PML_initPhysPoints(this)
    !Initialize arrays with physical points
      class(tPML) :: this
      integer :: k
      real*8 :: fx;
      do k=0,this%Nx
         this%Xi(k)=-this%Nx*this%hhx*cntx+this%hhx*(k);
         this%Xi05(k)=-this%Nx*this%hhx*cntx+this%hhx*(k-0.5d0);    
      enddo
      this%Xi05(this%Nx+1)=-this%Nx*this%hhx*cntx+this%hhx*(this%Nx+1.0d0-0.5d0);  
      do k=0,this%Ny
         this%Yi(k)=-this%Ny*this%hhy*cnty+this%hhy*(k);
         this%Yi05(k)=-this%Ny*this%hhy*cnty+this%hhy*(k-0.5d0);
      enddo
      this%Yi05(this%Ny+1)=-this%Ny*this%hhy*cnty+this%hhy*(this%Ny+1.0d0-0.5d0); 
      do k=0,this%Nz
         this%Zi(k)=-this%Nz*this%hhz*cntz+this%hhz*(k);
         this%Zi05(k)=-this%Nz*this%hhz*cntz+this%hhz*(k-0.5d0);
      enddo
      this%Zi05(this%Nz+1)=-this%Nz*this%hhz*cntz+this%hhz*(this%Nz+1.0d0-0.5d0);
    end subroutine PML_initPhysPoints;

    
!PML_initSigmaArrays-----------------------------------------------------------      
    subroutine PML_initSigmaArrays(this)
    !Initialize arrays with physical points
      class(tPML) :: this
      integer :: k
      do k=0,this%Nx
         this%sgmXi(k)=this%sgmpml_x(this%Xi(k));
         this%sgmXi05(k)=this%sgmpml_x(this%Xi05(k));       
      enddo
      this%sgmXi05(this%Nx+1)=this%sgmpml_x(this%Xi05(this%Nx+1));  
      do k=0,this%Ny
         this%sgmYi(k)=this%sgmpml_y(this%Yi(k));
         this%sgmYi05(k)=this%sgmpml_y(this%Yi05(k));
      enddo
      this%sgmYi05(this%Ny+1)=this%sgmpml_y(this%Yi05(this%Ny+1)); 
      do k=0,this%Nz
         this%sgmZi(k)=this%sgmpml_z(this%Zi(k));
         this%sgmZi05(k)=this%sgmpml_z(this%Zi05(k));
      enddo
      this%sgmZi05(this%Nz+1)=this%sgmpml_z(this%Zi05(this%Nz+1)); 
    end subroutine PML_initSigmaArrays;

    
!PML_DoStep------------------------------------------------------------------     
    subroutine PML_DoStep(this, t)
    !Do update step in PML 
      class(tPML) :: this
      integer, intent(in) :: t;
       call cpu_time(tstarts(15))
       !$ tstarts(15) = OMP_get_wtime();
          call this%updateEInteriorPML;        ! Calculate E field in PML everywhere except outer boundary
       call cpu_time(tends(15))
       !$ tends(15) = OMP_get_wtime();

       call cpu_time(tstarts(16))
       !$ tstarts(16) = OMP_get_wtime();
          call this%updateHInteriorPML;        ! Calculate H field in PML everywhere except outer boundary
       call cpu_time(tends(16))
       !$ tends(16) = OMP_get_wtime();
       
      !call this%fillEBoundaryPML(t);       ! Fill E boundary in PML with current Eboundary condition
      !call this%fillHBoundaryPML(t);       ! Fill H boundary in PML with current Hboundary condition
    end subroutine PML_DoStep;    

    
!SGMPML_X------------------------------------------------------------------   
  real*8 function sgmpml_x(this, x)
  ! Get sigmaX function in the point
    class(tPML) :: this;
    real*8, intent(in) :: x;
    if (x > this%xi(this%Nx-this%NumPML)) then
       sgmpml_x=this%Apml*(x-this%xi(this%Nx-this%NumPML))**2
    elseif (x < this%xi(this%NumPml)) then
       sgmpml_x=this%Apml*(x-this%xi(this%Numpml))**2
    else
       sgmpml_x=0.0d0
    endif
  end function sgmpml_x

!SGMPML_Y------------------------------------------------------------------    
  real*8 function sgmpml_y(this, y)
  ! Get sigmaY function in the point
    class(tPML) :: this;
    real*8, intent(in) :: y;
    if (y > this%yi(this%Ny-this%NumPML)) then
       sgmpml_y=this%Apml*(y-this%yi(this%Ny-this%NumPML))**2
    elseif (y < this%yi(this%NumPml)) then
       sgmpml_y=this%Apml*(y-this%yi(this%Numpml))**2
    else
       sgmpml_y=0.0d0
    endif
  end function sgmpml_y

!SGMPML_Z------------------------------------------------------------------  
  real*8 function sgmpml_z(this, z)
  ! Get sigmaZ function in the point
    class(tPML) :: this;
    real*8, intent(in) :: z;
    if (z > this%zi(this%Nz-this%NumPML)) then
       sgmpml_z=this%Apml*(z-this%zi(this%Nz-this%NumPML))**2
    elseif (z < this%zi(this%NumPml)) then
       sgmpml_z=this%Apml*(z-this%zi(this%Numpml))**2
    else
       sgmpml_z=0.0d0
    endif
  end function sgmpml_z


!updateEpointPML------------------------------------------------------------------
    subroutine updateEpointPML(this, Px, Py, Pz, component)
    ! Update E point with scheme
      class(tPML) :: this
      integer, intent(in) :: Px, Py, Pz, component;
      integer :: inc;
      real*8 :: sg_z, sg_y, sg_x, a11, a1, a2, a3;
      real*8 :: DH, fldpml0;
      inc = 0;
      select case(component)
      case(1) ! EX update
      !  if (checkbounds==1) then
      !     ! Check all the necessary points are correct
      !     inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, component);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py+1, Pz, 3);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 3);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz+1, 2);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 2);
      !  endif
      !  if (inc==0) then
           sg_x =  this%sgmXi05(Px);! this%sgmpml_x(this%xi05(Px));
           sg_y =  this%sgmYi(Py);! this%sgmpml_y(this%yi(Py));
           sg_z =  this%sgmZi(Pz);! this%sgmpml_z(this%zi(Pz));
           a11 = (2.0d0 + cc*this%ht*sg_y);
           a1 = (2.0d0 - cc*this%ht*sg_y) / a11;
           a2 = (2.0d0 + cc*this%ht*sg_x) / a11;
           a3 = (cc*this%ht*sg_x - 2.0d0) / a11;
           DH = cc*this%ht*( (this%Hf%Z(Px,Py+1,Pz)-this%Hf%Z(Px,Py,Pz))/this%hhy - (this%Hf%Y(Px,Py,Pz+1)-this%Hf%Y(Px,Py,Pz))/this%hhz );
           ! update EpmlX field
           fldpml0=this%Epml%X(Px,Py,Pz);
           this%Epml%X(Px,Py,Pz)=(DH-this%ht*this%JfE%X(Px,Py,Pz)*4.0d0*Pi/cc+this%Epml%X(Px,Py,Pz)*(1.0d0-0.5d0*sg_z*this%ht*cc))/(1.0d0+0.5d0*sg_z*this%ht*cc);
           ! real EX field
           this%Ef%X(Px,Py,Pz)=this%Ef%X(Px,Py,Pz)*a1+this%Epml%X(Px,Py,Pz)*a2+fldpml0*a3;
      !     if (checkupdated==1) then
      !        this%Ef%chX(Px,Py,Pz) = this%Ef%chX(Px,Py,Pz)+1;   ! this string is important for update checking
      !     endif   
      !  else
      !      fatalerr = 301;
      !      write(*,*) 'Alert! EX PML point update bounds error.'
      !  endif 
     case(2) ! EY update
      !  if (checkbounds==1) then
      !      ! Check all the necessary points are correct
      !      inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, component);
      !      inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz+1, 1);
      !      inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 1);
      !      inc = inc + this%Hf%mesh_checkpoint(Px+1, Py, Pz, 3);
      !      inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 3);
      !  endif
      !  if (inc==0) then
           sg_x = this%sgmXi(Px);  !this%sgmpml_x(this%xi(Px));
           sg_y = this%sgmYi05(Py);!this%sgmpml_y(this%yi05(Py));
           sg_z = this%sgmZi(Pz);  !this%sgmpml_z(this%zi(Pz));
           a11 = (2.0d0 + cc*this%ht*sg_z);
           a1 = (2.0d0 - cc*this%ht*sg_z) / a11;
           a2 = (2.0d0 + cc*this%ht*sg_y) / a11;
           a3 = (cc*this%ht*sg_y - 2.0d0) / a11;
           DH = cc*this%ht*( (this%Hf%X(Px,Py,Pz+1)-this%Hf%X(Px,Py,Pz))/this%hhz - (this%Hf%Z(Px+1,Py,Pz)-this%Hf%Z(Px,Py,Pz))/this%hhx );
           ! update EpmlY field
           fldpml0=this%Epml%Y(Px,Py,Pz);
           this%Epml%Y(Px,Py,Pz)=(DH-this%ht*this%JfE%Y(Px,Py,Pz)*4.0d0*Pi/cc+this%Epml%Y(Px,Py,Pz)*(1.0d0-0.5d0*sg_x*this%ht*cc))/(1.0d0+0.5d0*sg_x*this%ht*cc);
           ! real EY field
           this%Ef%Y(Px,Py,Pz)=this%Ef%Y(Px,Py,Pz)*a1+this%Epml%Y(Px,Py,Pz)*a2+fldpml0*a3;
      !     if (checkupdated==1) then
      !        this%Ef%chY(Px,Py,Pz) = this%Ef%chY(Px,Py,Pz)+1;   ! this string is important for update checking
      !     endif   
      !  else
      !    fatalerr = 301;
      !    write(*,*) 'Alert! EY PML point update bounds error.'
      !  endif 
      case(3) ! EZ update
      !  if (checkbounds==1) then
      !   ! Check all the necessary points are correct
      !     inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, component);
      !     inc = inc + this%Hf%mesh_checkpoint(Px+1, Py, Pz, 2);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 2);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py+1, Pz, 1);
      !     inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, 1);
      !  endif  
      !  if (inc==0) then
           sg_x = this%sgmXi(Px); !this%sgmpml_x(this%xi(Px));
           sg_y = this%sgmYi(Py); !this%sgmpml_y(this%yi(Py));
           sg_z = this%sgmZi05(Pz); !this%sgmpml_z(this%zi05(Pz));
           a11 = (2.0d0 + cc*this%ht*sg_x);
           a1 = (2.0d0 - cc*this%ht*sg_x) / a11;
           a2 = (2.0d0 + cc*this%ht*sg_z) / a11;
           a3 = (cc*this%ht*sg_z - 2.0d0) / a11;
           DH = cc*this%ht*( (this%Hf%Y(Px+1,Py,Pz)-this%Hf%Y(Px,Py,Pz))/this%hhx - (this%Hf%X(Px,Py+1,Pz)-this%Hf%X(Px,Py,Pz))/this%hhy );
           ! update EpmlZ field
           fldpml0=this%Epml%Z(Px,Py,Pz);
           this%Epml%Z(Px,Py,Pz)=(DH-this%ht*this%JfE%Z(Px,Py,Pz)*4.0d0*Pi/cc+this%Epml%Z(Px,Py,Pz)*(1.0d0-0.5d0*sg_y*this%ht*cc))/(1.0d0+0.5d0*sg_y*this%ht*cc);
           ! real EZ field
           this%Ef%Z(Px,Py,Pz)=this%Ef%Z(Px,Py,Pz)*a1+this%Epml%Z(Px,Py,Pz)*a2+fldpml0*a3;
      !     if (checkupdated==1) then
      !        this%Ef%chZ(Px,Py,Pz) = this%Ef%chZ(Px,Py,Pz)+1;   ! this string is important for update checking
      !     endif   
      !  else
      !     fatalerr = 301;
      !     write(*,*) 'Alert! EZ PML point update bounds error.'
      !  endif
      case default ! Incorrect
         fatalerr = 303;
         write(*,*) 'Alert! Incorrect use of Component in updateEpoint function!', component  
      end select
    end subroutine updateEpointPML;

    
!updateHpointPML------------------------------------------------------------------
    subroutine updateHpointPML(this, Px, Py, Pz, component)
    ! Update H point with scheme
      class(tPML) :: this
      integer, intent(in) :: Px, Py, Pz, component;
      integer :: inc;
      real*8 :: sg_z, sg_y, sg_x, a0, a1, a2, a3, a4, a5, a6, a11;
      real*8 :: DH, fldpml0;
      inc = 0;
      select case(component)
      case(1) ! HX update
        ! if (checkbounds==1) then
        !    ! Check all the necessary points are correct
        !    inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, component);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 3);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py-1, Pz, 3);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 2);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz-1, 2);
        ! endif
        ! if (inc==0) then
            sg_x = this%sgmXi(Px); !this%sgmpml_x(this%xi(Px));
            sg_y = this%sgmYi05(Py); !this%sgmpml_y(this%yi05(Py));
            sg_z = this%sgmZi05(Pz); !this%sgmpml_z(this%zi05(Pz));
            a11=(2.0d0+cc*this%ht*sg_z)
            a1=(2.0d0-cc*this%ht*sg_z)/a11
            a2=2.0d0*cc*this%ht/a11
            a0=2.0d0+cc*this%ht*sg_y
            a4=(2.0d0-cc*this%ht*sg_y)/a0
            ! update HpmlX field
            DH = ( (this%Ef%Z(Px,Py,Pz)-this%Ef%Z(Px,Py-1,Pz))/this%hhy - (this%Ef%Y(Px,Py,Pz)-this%Ef%Y(Px,Py,Pz-1))/this%hhz );
            fldpml0=this%Hpml%X(Px,Py,Pz);
            this%Hpml%X(Px,Py,Pz)=this%Hpml%X(Px,Py,Pz)*a1-(DH+this%JfH%X(Px,Py,Pz)*4.0d0*Pi/cc)*a2;
            ! real HX field
            a5=(2.0+cc*this%ht*sg_x)/a0
            a6=(-2.0+cc*this%ht*sg_x)/a0
            this%Hf%X(Px,Py,Pz)=this%Hf%X(Px,Py,Pz)*a4+this%Hpml%X(Px,Py,Pz)*a5+fldpml0*a6;
        !    if (checkupdated==1) then
        !       this%Hf%chX(Px,Py,Pz) = this%Hf%chX(Px,Py,Pz)+1; ! this string is important for update checking
        !    endif   
        ! else
        !    fatalerr = 302;
        !    write(*,*) 'Alert! HX PML point update bounds error.'
        ! endif 
      case(2) ! HY update
        ! if (checkbounds==1) then
        !    ! Check all the necessary points are correct
        !    inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, component);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 1);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz-1, 1);
        !    inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 3);
        !    inc = inc + this%Ef%mesh_checkpoint(Px-1, Py, Pz, 3);
        ! endif  
        ! if (inc==0) then
            sg_x = this%sgmXi05(Px); !this%sgmpml_x(this%xi05(Px));
            sg_y = this%sgmYi(Py);   !this%sgmpml_y(this%yi(Py));
            sg_z = this%sgmZi05(Pz); !this%sgmpml_z(this%zi05(Pz));
            a11=(2.0d0+cc*this%ht*sg_z)
            a4=(2.0d0-cc*this%ht*sg_z)/a11
            a5=(2.0d0+cc*this%ht*sg_y)/a11
            a6=(-2.0d0+cc*this%ht*sg_y)/a11
            a0=2.0d0+cc*this%ht*sg_x
            a1=(2.0d0-cc*this%ht*sg_x)/a0
            a2=2.0d0*cc*this%ht/a0
            ! update HpmlY field
            DH = ( (this%Ef%X(Px,Py,Pz)-this%Ef%X(Px,Py,Pz-1))/this%hhz - (this%Ef%Z(Px,Py,Pz)-this%Ef%Z(Px-1,Py,Pz))/this%hhx );
            fldpml0=this%Hpml%Y(Px,Py,Pz);
            this%Hpml%Y(Px,Py,Pz)=-(DH+this%JfH%Y(Px,Py,Pz)*4.0d0*Pi/cc)*a2+this%Hpml%Y(Px,Py,Pz)*a1;
            ! real HY field
            this%Hf%Y(Px,Py,Pz)=this%Hf%Y(Px,Py,Pz)*a4+this%Hpml%Y(Px,Py,Pz)*a5+fldpml0*a6;
        !    if (checkupdated==1) then
        !       this%Hf%chY(Px,Py,Pz) = this%Hf%chY(Px,Py,Pz)+1; ! this string is important for update checking
        !    endif   
        ! else
        !    fatalerr = 302;
        !    write(*,*) 'Alert! HY PML point update bounds error.'
        ! endif 
      case(3) ! HZ update
        !if (checkbounds==1) then 
        !   ! Check all the necessary points are correct
        !   inc = inc + this%Hf%mesh_checkpoint(Px, Py, Pz, component);
        !   inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 2);
        !   inc = inc + this%Ef%mesh_checkpoint(Px-1, Py, Pz, 2);
        !   inc = inc + this%Ef%mesh_checkpoint(Px, Py, Pz, 1);
        !   inc = inc + this%Ef%mesh_checkpoint(Px, Py-1, Pz, 1);
        !endif   
        !if (inc==0) then
           sg_x = this%sgmXi05(Px); !this%sgmpml_x(this%xi05(Px));
           sg_y = this%sgmYi05(Py); !this%sgmpml_y(this%yi05(Py));
           sg_z = this%sgmZi(Pz);   !this%sgmpml_z(this%zi(Pz));
           a1=(2.0d0-cc*this%ht*sg_y)/(2.0d0+cc*this%ht*sg_y)
           a2=2.0d0*cc*this%ht/(2.0d0+cc*this%ht*sg_y)
           a11=(2.0d0+cc*this%ht*sg_x)
           a4=(2.0d0-cc*this%ht*sg_x)/a11
           a5=(2.0d0+cc*this%ht*sg_z)/a11
           a6=(-2.0d0+cc*this%ht*sg_z)/a11
           ! update HpmlZ field
           DH = ( (this%Ef%Y(Px,Py,Pz)-this%Ef%Y(Px-1,Py,Pz))/this%hhx - (this%Ef%X(Px,Py,Pz)-this%Ef%X(Px,Py-1,Pz))/this%hhy );
           fldpml0=this%Hpml%Z(Px,Py,Pz);
           this%Hpml%Z(Px,Py,Pz)=-(DH+this%JfH%Z(Px,Py,Pz)*4.0d0*Pi/cc)*a2+this%Hpml%Z(Px,Py,Pz)*a1;
           ! real HZ field
           this%Hf%Z(Px,Py,Pz)=this%Hf%Z(Px,Py,Pz)*a4+this%Hpml%Z(Px,Py,Pz)*a5+fldpml0*a6;
        !   if (checkupdated==1) then
        !      this%Hf%chZ(Px,Py,Pz) = this%Hf%chZ(Px,Py,Pz)+1; ! this string is important for update checking
        !   endif   
        ! else
        !   fatalerr = 302;
        !   write(*,*) 'Alert! HZ PML point update bounds error.'
        !endif
      case default ! Incorrect
         fatalerr = 304;
         write(*,*) 'Alert! Incorrect use of Component in updateEpoint function!', component 
      end select
    end subroutine updateHpointPML;

!updateEInteriorPML------------------------------------------------------------------
    subroutine updateEInteriorPML(this)
    ! Update all the E fields in the interior
      class(tPML) :: this
      integer :: i, j ,k
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny-1
            do i=1,this%Nx
              if (.NOT.((i>=this%NumPML+1).AND.(i<=this%Nx-this%NumPML).AND.(j>=this%NumPML+1).AND.(j<=this%Ny-this%NumPML-1).AND.(k>=this%NumPML+1).AND.(k<=this%Nz-this%NumPML-1))) then
                 call this%updateEpointPML(i, j, k, 1); ! Ex update
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny
            do i=1,this%Nx-1
               if (.NOT.((i>=this%NumPML+1).AND.(i<=this%Nx-this%NumPML-1).AND.(j>=this%NumPML+1).AND.(j<=this%Ny-this%NumPML).AND.(k>=this%NumPML+1).AND.(k<=this%Nz-this%NumPML-1))) then
                  call this%updateEpointPML(i, j, k, 2); ! Ey update
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny-1
            do i=1,this%Nx-1
               if (.NOT.((i>=this%NumPML+1).AND.(i<=this%Nx-this%NumPML-1).AND.(j>=this%NumPML+1).AND.(j<=this%Ny-this%NumPML-1).AND.(k>=this%NumPML+1).AND.(k<=this%Nz-this%NumPML))) then
                  call this%updateEpointPML(i, j, k, 3); ! Ez update
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine updateEInteriorPML;


!updateHInterior------------------------------------------------------------------
    subroutine updateHInteriorPML(this)
    ! Update all the H fields in the interior
      class(tPML) :: this
      integer :: i, j ,k
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny
            do i=1,this%Nx-1
              if (.NOT.((i>=this%NumPML+1).AND.(i<=this%Nx-this%NumPML-1).AND.(j>=this%NumPML+2).AND.(j<=this%Ny-this%NumPML-1).AND.(k>=this%NumPML+2).AND.(k<=this%Nz-this%NumPML-1))) then
                  call this%updateHpointPML(i, j, k, 1); ! Hx update
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny-1
            do i=1,this%Nx
               if (.NOT.((i>=this%NumPML+2).AND.(i<=this%Nx-this%NumPML-1).AND.(j>=this%NumPML+1).AND.(j<=this%Ny-this%NumPML-1).AND.(k>=this%NumPML+2).AND.(k<=this%Nz-this%NumPML-1))) then
                  call this%updateHpointPML(i, j, k, 2); ! Hy update
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny
            do i=1,this%Nx
               if (.NOT.((i>=this%NumPML+2).AND.(i<=this%Nx-this%NumPML-1).AND.(j>=this%NumPML+2).AND.(j<=this%Ny-this%NumPML-1).AND.(k>=this%NumPML+1).AND.(k<=this%Nz-this%NumPML-1))) then
                  call this%updateHpointPML(i, j, k, 3); ! Hz update
               endif
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine updateHInteriorPML;


!updateEInteriorPMLEverywhere------------------------------------------------------------------
    subroutine updateEInteriorPMLEverywhere(this)
    ! Update all the E fields in the interior everywhere
      class(tPML) :: this
      integer :: i, j ,k
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny-1
            do i=1,this%Nx
               call this%updateEpointPML(i, j, k, 1); ! Ex update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny
            do i=1,this%Nx-1
               call this%updateEpointPML(i, j, k, 2); ! Ey update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny-1
            do i=1,this%Nx-1
               call this%updateEpointPML(i, j, k, 3); ! Ez update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine updateEInteriorPMLEverywhere;


!updateHInterior------------------------------------------------------------------
    subroutine updateHInteriorPMLEverywhere(this)
    ! Update all the H fields in the interior everywhere
      class(tPML) :: this
      integer :: i, j ,k
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny
            do i=1,this%Nx-1
               call this%updateHpointPML(i, j, k, 1); ! Hx update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz
         do j=1,this%Ny-1
            do i=1,this%Nx
               call this%updateHpointPML(i, j, k, 2); ! Hy update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do k=1,this%Nz-1
         do j=1,this%Ny
            do i=1,this%Nx
               call this%updateHpointPML(i, j, k, 3); ! Hz update
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine updateHInteriorPMLEverywhere;
    

!fillEBoundaryPML------------------------------------------------------------------
  subroutine  fillEBoundaryPML(this,t)
    !Fill E boundary with some boundary condition
    class(tPML) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    ! Fill EX Boundaries
    !$OMP PARALLEL DO SHARED(balance)
    do k=0,this%Nz
       do i=1,this%Nx
          call this%fillEBoundaryPointPML(t, i, 0, k, 1);
          call this%fillEBoundaryPointPML(t, i, this%Ny, k, 1);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(balance)
    do j=1,this%Ny-1
       do i=1,this%Nx
          call this%fillEBoundaryPointPML(t, i, j, 0, 1); 
          call this%fillEBoundaryPointPML(t, i, j, this%Nz, 1);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Fill EY Boundaries
    !$OMP PARALLEL DO SHARED(balance)
    do j=1,this%Ny
       do i=0,this%Nx
          call this%fillEBoundaryPointPML(t, i, j, 0, 2);
          call this%fillEBoundaryPointPML(t, i, j, this%Nz, 2);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(balance)
    do k=1,this%Nz-1
       do j=1,this%Ny
          call this%fillEBoundaryPointPML(t, 0, j, k, 2);
          call this%fillEBoundaryPointPML(t, this%Nx, j, k, 2);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Fill EZ Boundaries
    !$OMP PARALLEL DO SHARED(balance)
    do k=1,this%Nz
       do j=0,this%Ny
          call this%fillEBoundaryPointPML(t, 0, j, k, 3);
          call this%fillEBoundaryPointPML(t, this%Nx, j, k, 3);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(balance)
    do k=1,this%Nz
       do i=1,this%Nx-1
          call this%fillEBoundaryPointPML(t, i, 0, k, 3);
          call this%fillEBoundaryPointPML(t, i, this%Ny, k, 3);
          !$ call AddParBalance;
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine fillEBoundaryPML


!fillHBoundaryPML------------------------------------------------------------------
  subroutine  fillHBoundaryPML(this, t)
    ! Update H boundary with scheme
    class(tPML) :: this
    integer, intent(in) :: t
    integer :: i ,j ,k
    ! HX Boundaries Update
    !$OMP PARALLEL DO SHARED(balance)
    do j=1,this%Ny
       do i=0,this%Nx
          call this%fillHBoundaryPointPML(t, i, j, 1, 1);
          call this%fillHBoundaryPointPML(t, i, j, this%Nz, 1);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(balance)
    do k=2,this%Nz-1
       do i=0,this%Nx
          call this%fillHBoundaryPointPML(t, i, 1, k, 1);
          call this%fillHBoundaryPointPML(t, i, this%Ny, k, 1);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
  !  !$OMP PARALLEL DO SHARED(balance)
  !  do k=2,this%Nz-1
  !     do j=2,this%Ny-1
  !        call this%fillHBoundaryPointPML(t, 0, j, k, 1);
  !        call this%fillHBoundaryPointPML(t, this%Nx, j, k, 1);
  !        !$ call AddParBalance;
  !     enddo
  !  enddo
  !  !$OMP END PARALLEL DO
    
    ! HY Boundaries Update
    !$OMP PARALLEL DO SHARED(balance)
    do k=1,this%Nz
       do j=0,this%Ny
          call this%fillHBoundaryPointPML(t, 1, j, k, 2);
          call this%fillHBoundaryPointPML(t, this%Nx, j, k, 2);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(balance)
    do j=0,this%Ny
       do i=2,this%Nx-1
          call this%fillHBoundaryPointPML(t, i, j, 1, 2);
          call this%fillHBoundaryPointPML(t, i, j, this%Nz, 2);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
   ! !$OMP PARALLEL DO SHARED(balance)
   ! do k=2,this%Nz-1
   !    do i=2,this%Nx-1
   !       call this%fillHBoundaryPointPML(t, i, 0, k, 2);
   !       call this%fillHBoundaryPointPML(t, i, this%Ny, k, 2);
   !       !$ call AddParBalance;
   !    enddo
   ! enddo
   ! !$OMP END PARALLEL DO
    
    ! HZ Boundaries Update
    !$OMP PARALLEL DO SHARED(balance)
    do k=0,this%Nz
       do i=1,this%Nx
          call this%fillHBoundaryPointPML(t, i, 1, k, 3);
          call this%fillHBoundaryPointPML(t, i, this%Ny, k, 3);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO SHARED(balance)
    do k=0,this%Nz
       do j=2,this%Ny-1
          call this%fillHBoundaryPointPML(t, 1, j, k, 3);
          call this%fillHBoundaryPointPML(t, this%Nx, j, k, 3);
          !$ call AddParBalance;
       enddo
    enddo
    !$OMP END PARALLEL DO
  !  !$OMP PARALLEL DO SHARED(balance)
  !  do j=2,this%Ny-1
  !     do i=2,this%Nx-1
  !        call this%fillHBoundaryPointPML(t, i, j, 0, 3);
  !        call this%fillHBoundaryPointPML(t, i, j, this%Nz, 3);
  !        !$ call AddParBalance;
  !     enddo
  !  enddo
  !  !$OMP END PARALLEL DO
  end subroutine fillHBoundaryPML


!fillEBoundaryPointPML------------------------------------------------------------------
  subroutine  fillEBoundaryPointPML(this, t, Px, Py, Pz, component)
    ! Fill E point with some boundary condition
    class(tPML) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    select case (this%EboundarytypePML)
    case (0) ! Zeros
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz) = 0;
       case(2)
          this%Ef%Y(Px, Py, Pz) = 0;
       case(3)
          this%Ef%Z(Px, Py, Pz) = 0;
       case default ! illegal type
          fatalerr = 307;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPointPML function!', component           
       end select
    case (1) ! Analytical solution   
       select case (component)
       case(1)
          this%Ef%X(Px, Py, Pz) = this%Ean%X(Px, Py, Pz);
       case(2)
          this%Ef%Y(Px, Py, Pz) = this%Ean%Y(Px, Py, Pz);
       case(3)
          this%Ef%Z(Px, Py, Pz) = this%Ean%Z(Px, Py, Pz);
       case default ! illegal type
          fatalerr = 207;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPoint function!', component    
       end select   
    case (2) ! Scheme
       select case (component)
       case(1)
          call this%updateEpointPML(Px, Py, Pz, 1);
       case(2)
          call this%updateEpointPML(Px, Py, Pz, 2);
       case(3)
          call this%updateEpointPML(Px, Py, Pz, 3);
       case default ! illegal type
          fatalerr = 307;
          write(*,*) 'Alert! Incorrect use of Component in fillEboundaryPointPML function!', component    
       end select
    case default ! illegal type
       fatalerr = 305;
       write(*,*) 'Alert! Illegal boundary condition for E in fillEboundaryPML function!', this%EboundarytypePML     
    end select
    if (checkupdated==1) then
       if (this%EboundarytypePML /= 2) then
          select case (component)
          case(1)
             this%Ef%chX(Px, Py, Pz) = this%Ef%chX(Px, Py, Pz)+1;
          case(2)
             this%Ef%chY(Px, Py, Pz) = this%Ef%chY(Px, Py, Pz)+1;
          case(3)
             this%Ef%chZ(Px, Py, Pz) = this%Ef%chZ(Px, Py, Pz)+1;
          end select
       endif
    endif  
  end subroutine fillEBoundaryPointPML


!fillHBoundaryPointPML------------------------------------------------------------------
  subroutine  fillHBoundaryPointPML(this, t, Px, Py, Pz, component)
    ! Fill H point with some boundary condition
    class(tPML) :: this
    integer, intent(in) :: t;
    integer, intent(in) :: Px, Py, Pz, component;
    select case (this%HboundarytypePML)
    case (0) ! fill with zero
       select case (component)
       case(1)
          this%Hf%X(Px, Py, Pz) = 0;
       case(2)
          this%Hf%Y(Px, Py, Pz) = 0;
       case(3)
          this%Hf%Z(Px, Py, Pz) = 0;
       case default ! illegal type
          fatalerr = 308;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPointPML function!', component        
       end select
    case (1) ! fill from the analytical solution
       select case (component)
       case(1)
          this%Hf%X(Px, Py, Pz) = this%Han%X(Px, Py, Pz);
       case(2)
          this%Hf%Y(Px, Py, Pz) = this%Han%Y(Px, Py, Pz);
       case(3)
          this%Hf%Z(Px, Py, Pz) = this%Han%Z(Px, Py, Pz);
       case default ! illegal type
          fatalerr = 208;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPoint function!', component       
       end select   
    case (2) ! fill with scheme
       select case (component)
       case(1)
          call this%updateHpointPML(Px, Py, Pz, 1);
       case(2)
          call this%updateHpointPML(Px, Py, Pz, 2);
       case(3)
          call this%updateHpointPML(Px, Py, Pz, 3);
       case default ! illegal type
          fatalerr = 308;
          write(*,*) 'Alert! Incorrect use of Component in fillHboundaryPointPML function!', component    
       end select
    case default ! illegal type
       fatalerr = 306;
       write(*,*) 'Alert! Illegal boundary condition for H in fillHboundaryPML function!', this%HboundarytypePML
    end select
    if (checkupdated==1) then
       if (this%HboundarytypePML /= 2) then  
          select case (component)
          case(1)
             this%Hf%chX(Px, Py, Pz) = this%Hf%chX(Px, Py, Pz)+1;
          case(2)
             this%Hf%chY(Px, Py, Pz) = this%Hf%chY(Px, Py, Pz)+1;
          case(3)
             this%Hf%chZ(Px, Py, Pz) = this%Hf%chZ(Px, Py, Pz)+1;
          end select
       endif
    endif
  end subroutine fillHBoundaryPointPML  
  
end module pmlclass
