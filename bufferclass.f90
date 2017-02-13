module bufferclass
  ! Transport to transfer data from problems to problems
  use commonvars
  use meshclass
  implicit none;
  
  type, public :: tBuffer
     type(TMesh) :: Ef, Hf, Je, Jh;
     type(TMesh) :: Ec0, Ec1;
     type(TMesh) :: mainMu;
     type(tMesh) :: Estatic;
     type(tMesh) :: Hstatic;
   contains
      procedure buffer_init;
      procedure buffer_destroy;
      procedure clearFields;
   end type tBuffer

contains

!buffer_init------------------------------------------------------------------   
    subroutine buffer_Init(this, Px, Py, Pz)
    !Initialize the buffer
      class(tBuffer) :: this
      integer, intent(in) :: Px, Py, Pz;
      call this%Ef%mesh_init(Px, Py, Pz, 0, -1);
      call this%Hf%mesh_init(Px, Py, Pz, 1, -1);
      call this%Je%mesh_init(Px, Py, Pz, 0, -1);
      call this%Jh%mesh_init(Px, Py, Pz, 1, -1);
      call this%Estatic%mesh_init(Px, Py, Pz, 0, -1);
      call this%Hstatic%mesh_init(Px, Py, Pz, 1, -1);
      call this%Ec0%mesh_init(Px, Py, Pz, 0, -1);
      call this%Ec1%mesh_init(Px, Py, Pz, 0, -1);
      call this%mainMu%mesh_init(Px, Py, Pz, 0, -1);
      this%mainMU%X(1:Px,0:Py,0:Pz)=1;
      this%mainMU%Y(0:Px,1:Py,0:Pz)=1;
      this%mainMU%Z(0:Px,0:Py,1:Pz)=1;
    end subroutine buffer_Init;

!buffer_destroy------------------------------------------------------------------      
    subroutine buffer_Destroy(this)
    !Destroy the buffer structures and free the memory  
      class(tBuffer) :: this
      call this%Ef%mesh_destroy;
      call this%Hf%mesh_destroy;
      call this%Je%mesh_destroy;
      call this%Jh%mesh_destroy;
      call this%Estatic%mesh_destroy;
      call this%Hstatic%mesh_destroy;
      call this%Ec0%mesh_destroy;
      call this%Ec1%mesh_destroy;
      call this%mainMu%mesh_destroy;
    end subroutine buffer_Destroy;
    
!buffer_clear------------------------------------------------------------------      
    subroutine clearFields(this)
    !Clear all the data in the buffer 
      class(tBuffer) :: this
      call this%Ef%mesh_clear;
      call this%Hf%mesh_clear;
      call this%Je%mesh_clear;
      call this%Jh%mesh_clear;
    end subroutine clearFields;
   
end module bufferclass
