module sourceclass
! Source class that describes analytical solution and sources  
  use commonvars
  use newdipole
  use sdipole
   implicit none;
  type, public :: tSource
     integer :: sourcetype;
     integer :: problem_type
     type(tNewDipole) :: newdip;
    contains
      procedure getEpoint;
      procedure getHpoint;
      procedure getJpoint;
  end type tSource

contains

!getEpoint------------------------------------------------------------------  
  real*8 function getEpoint(this, x, y, z, t, xp, yp, zp, type)
  ! Get E field value from the analytic solution  
    class(tSource) :: this
    real*8, intent(in) :: x, y, z, t
    integer, intent(in) :: xp, yp, zp;
    integer, intent(in) :: type
    real*8 :: res;
    res = 0;
    select case (this%sourcetype)
    case (0)
        res = this%newdip%getNewDipoleEpoint(x, y, z, xp, yp, zp, t, type, this%problem_type);
       ! res = getDipoleEpoint(x,y,z,t,type);
    case (1)
       res = getSDipoleEpoint(x,y,z,t,type);
    end select
    getEpoint = res;
  end function getEpoint

!getHpoint------------------------------------------------------------------    
  real*8 function getHpoint(this, x, y, z, t, xp, yp, zp, type)
  ! Get H field value from the analytic solution  
    class(tSource) :: this
    real*8, intent(in) :: x, y, z, t
    integer, intent(in) :: xp, yp, zp;
    integer, intent(in) :: type
    real*8 :: res;
    res = 0;
    select case (this%sourcetype)
    case (0)
       res = this%newdip%getNewDipoleHpoint(x, y, z, xp, yp, zp, t, type, this%problem_type);
       ! res = getDipoleHpoint(x,y,z,t,type);
    case (1)
       res = getSDipoleHpoint(x,y,z,t,type);
    end select   
    getHpoint = res;
  end function getHpoint

!getJpoint------------------------------------------------------------------    
  real*8 function getJpoint(this, x, y, z, t, xp, yp, zp, type)
  ! Get J source value from the analytic solution    
    class(tSource) :: this
    real*8, intent(in) :: x, y, z, t
    integer, intent(in) :: xp, yp, zp;
    integer, intent(in) :: type
    real*8 :: res;
    res = 0;
    select case (this%sourcetype)
    case (0)
       res = this%newdip%getNewDipoleJpoint(x, y, z, xp, yp, zp, t, type, this%problem_type);
       ! res = getDipoleJpoint(x,y,z,t,type);
    case (1)
       res = getSDipoleJpoint(x,y,z,t,type);
    end select    
    getJpoint = res;
  end function getJpoint
  
end module sourceclass
