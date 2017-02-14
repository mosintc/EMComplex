!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains implementation of tTimers class

module timersclass
  ! The set of timers to control the problems
  use commonvars
  implicit none;
  
  type, public :: tTimers
     integer :: Nprobs;
     integer, allocatable :: isused(:);
     integer, allocatable :: timer(:);
   contains
      procedure timers_init;
      procedure timers_destroy;
      procedure getfreetimer;
      procedure closetimer;
      procedure howmanyisused;
   end type tTimers

contains

!timers_init------------------------------------------------------------------   
    subroutine timers_init(this, Nprobs)
    !Initialize the timers arrays
      class(tTimers) :: this
      integer, intent(in) :: Nprobs;
      integer :: i;
      this%NProbs = Nprobs-1;
      allocate(this%isused(0:Nprobs-1));
      allocate(this%timer(0:Nprobs-1));
      do i=0,Nprobs-1
         this%isused(i) = 0;
         this%timer(i) = 0;
      enddo
    end subroutine timers_init;

!timers_destroy------------------------------------------------------------------      
    subroutine timers_destroy(this)
    !Destroy the timers arrays
      class(tTimers) :: this
      deallocate(this%isused, this%timer);
    end subroutine timers_destroy;

!getfreetimer------------------------------------------------------------------      
    integer function getfreetimer(this)
    !Find the first available timer
      class(tTimers) :: this
      integer :: i, res, ch;
      i = 0;
      do while (i<=this%Nprobs.AND.this%isused(i)==1)
         i = i+1;
      end do
      
      if (i>this%Nprobs) then
         getfreetimer = -1;
      else
         write(*,*) 'We choosed', i;
         getfreetimer = i;
      endif 
    end function getfreetimer;


!howmanyisused------------------------------------------------------------------      
    integer function howmanyisused(this)
    !Set the timer to the free parameters
      class(tTimers) :: this
      integer :: i, num;
      num = 0;
      do i=0,this%Nprobs
         if (this%isused(i) == 1) then
            num = num+1;
         endif
      enddo
      howmanyisused = num;
    end function howmanyisused; 
    

!closetimer------------------------------------------------------------------      
    subroutine closetimer(this, N)
    !Set the timer to the free parameters
      class(tTimers) :: this
      integer, intent(in) :: N;
      this%isused(N) = 0;
      this%timer(N) = 0;
    end subroutine closetimer;    
    
end module timersclass
