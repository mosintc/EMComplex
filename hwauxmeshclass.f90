module auxmeshclass
! Class that contains arrays with field components  
  use commonvars
   implicit none;
  type, public :: tAuxMesh
     integer :: amtype;
     integer :: holder;
     real*8, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:); ! actual arrays
  contains
      procedure :: auxmesh_init
      procedure :: auxmesh_destroy
      procedure :: auxmesh_clear
  end type tAuxMesh    
   
contains
  
!AUXMESH_INIT------------------------------------------------------------------   
    subroutine auxmesh_init(this, Px, Py, Pz, amtype, msholder)
    !Initialization of the auxiliary mesh  
      class(tAuxMesh) :: this
      integer, intent(in) :: Px, Py, Pz;
      integer, intent(in) :: amtype;
      integer, intent(in) :: msholder;
      this%amtype = amtype;
      this%holder = msholder;
      select case(this%amtype)
      case (0) 
         allocate(this%Y(1:Px, 0:Py, 0:Pz), this%Z(1:Px, 0:Py, 0:Pz));        
      case (1) 
         allocate(this%X(0:Px, 1:Py, 0:Pz), this%Z(0:Px, 1:Py, 0:Pz));
      case (2) 
         allocate(this%X(0:Px, 0:Py, 1:Pz), this%Y(0:Px, 0:Py, 1:Pz));
      end select
      call this%auxmesh_clear;
    end subroutine auxmesh_init;

!AUXMESH_DESTROY------------------------------------------------------------------   
    subroutine auxmesh_destroy(this) 
      !Destroy Mesh structure and deallocate the memor
      class(tAuxMesh) :: this
      select case(this%amtype)
      case (0) 
         deallocate(this%Y, this%Z);        
      case (1) 
         deallocate(this%X, this%Z);
      case (2) 
         deallocate(this%X, this%Y);
      end select
    end subroutine auxmesh_destroy;


!AuxMesh_Clear------------------------------------------------------------------      
    subroutine auxmesh_clear(this) 
    !Clear check arrays
      class(tAuxMesh) :: this
      select case(this%amtype)
      case (0)
         this%Y(:,:,:)=0.0d0;
         this%Z(:,:,:)=0.0d0;        
      case (1) 
         this%X(:,:,:)=0.0d0;
         this%Z(:,:,:)=0.0d0;
      case (2)
         this%X(:,:,:)=0.0d0;
         this%Y(:,:,:)=0.0d0;
      end select
    end subroutine auxmesh_clear;
       
end module auxmeshclass
