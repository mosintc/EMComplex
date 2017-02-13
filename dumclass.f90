module dumclass
  ! DUM class
  implicit none;
  
  type, public :: tDUM
     integer :: pt(0:2);
   contains
      procedure init;
   end type tDUM
   
contains

!init------------------------------------------------------------------   
    subroutine init(this)
      class(tDUM) :: this;
      integer :: i;
      do i=0,2
         this%pt(i)=0;
      enddo
    end subroutine init;

end module dumclass
  
