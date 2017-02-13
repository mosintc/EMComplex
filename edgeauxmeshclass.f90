module EdgeAuxMeshclass
! Class that contains arrays for auxiliary variables on the edges
  use commonvars
   implicit none;
  type, public :: tEdgeAuxMesh
     integer :: amtype;
     integer :: holder;
     real*8, allocatable :: fbottom(:), fleft(:), fup(:), fright(:); ! actual arrays
  contains
      procedure :: EdgeAuxMesh_init
      procedure :: EdgeAuxMesh_destroy
      procedure :: EdgeAuxMesh_clear
  end type tEdgeAuxMesh    
   
contains
  
!EdgeAUXMESH_INIT------------------------------------------------------------------   
    subroutine EdgeAuxMesh_init(this, Vx, Vy, amtype, msholder)
    !Initialization of the edge auxiliary meshes  
      class(tEdgeAuxMesh) :: this
      integer, intent(in) :: Vx, Vy;
      integer, intent(in) :: amtype;
      integer, intent(in) :: msholder;
      this%amtype = amtype;
      this%holder = msholder;
      select case(this%amtype)
      case (0) 
         allocate(this%fbottom(0:Vx), this%fleft(0:Vy), this%fup(0:Vx), this%fright(0:Vy));        
      case (1) 
         allocate(this%fbottom(1:Vx), this%fleft(0:Vy), this%fup(1:Vx), this%fright(0:Vy));  
      end select
      call this%EdgeAuxMesh_clear;
    end subroutine EdgeAuxMesh_init;

!EdgeAUXMESH_DESTROY------------------------------------------------------------------   
    subroutine Edgeauxmesh_destroy(this) 
      !Destroy Mesh structure and deallocate the memor
      class(tEdgeAuxMesh) :: this
      deallocate(this%fbottom, this%fleft, this%fup, this%fright);        
    end subroutine Edgeauxmesh_destroy;


!EdgeAuxMesh_Clear------------------------------------------------------------------      
    subroutine Edgeauxmesh_clear(this) 
    !Clear check arrays
      class(tEdgeAuxMesh) :: this
      this%fbottom(:)=0.0d0;
      this%fleft(:)=0.0d0;
      this%fup(:)=0.0d0;
      this%fright(:)=0.0d0;
    end subroutine Edgeauxmesh_clear;
       
end module EdgeAuxMeshclass
