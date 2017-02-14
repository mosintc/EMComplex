!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains implementation of tMesh class

module meshclass
! Class that contains arrays with field components  
  use commonvars
   implicit none;
  type, public :: tMesh
     integer :: shifttype;
     integer :: holder;
     integer, allocatable :: bounds(:, :, :);
     real*8, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:); ! actual arrays
     integer, allocatable :: chx(:,:,:), chy(:,:,:), chz(:,:,:);
  contains
      procedure :: mesh_init
      procedure :: mesh_destroy
      procedure :: mesh_clear
      procedure :: mesh_clearChA
      procedure :: mesh_checkarrays
      procedure :: mesh_checkpoint
  end type tMesh    
   
contains
!MESH_INIT------------------------------------------------------------------   
    subroutine mesh_init(this, Px, Py, Pz, stype, msholder)
    !Initialization of the mesh  
      class(tMesh) :: this
      integer, intent(in) :: Px, Py, Pz;
      integer, intent(in) :: stype;
      integer, intent(in) :: msholder;
      this%shifttype = stype;
      this%holder = msholder;
      allocate(this%bounds(1:3, 1:3, 1:2));
      select case(this%shifttype)
      case (0) ! Points on the centers of edges
         allocate(this%X(1:Px, 0:Py, 0:Pz), this%Y(0:Px, 1:Py, 0:Pz), this%Z(0:Px, 0:Py, 1:Pz));
         allocate(this%chX(1:Px, 0:Py, 0:Pz), this%chY(0:Px, 1:Py, 0:Pz), this%chZ(0:Px, 0:Py, 1:Pz));
         this%bounds(1,1,1)=1;
         this%bounds(1,1,2)=Px;
         this%bounds(1,2,1)=0;
         this%bounds(1,2,2)=Py;
         this%bounds(1,3,1)=0;
         this%bounds(1,3,2)=Pz;
         
         this%bounds(2,1,1)=0;
         this%bounds(2,1,2)=Px;
         this%bounds(2,2,1)=1;
         this%bounds(2,2,2)=Py;
         this%bounds(2,3,1)=0;
         this%bounds(2,3,2)=Pz;
         
         this%bounds(3,1,1)=0;
         this%bounds(3,1,2)=Px;
         this%bounds(3,2,1)=0;
         this%bounds(3,2,2)=Py;
         this%bounds(3,3,1)=1;
         this%bounds(3,3,2)=Pz;
      case (1) ! Points in the centers of faces
         allocate(this%X(0:Px, 1:Py, 1:Pz), this%Y(1:Px, 0:Py, 1:Pz), this%Z(1:Px, 1:Py, 0:Pz));
         allocate(this%chX(0:Px, 1:Py, 1:Pz), this%chY(1:Px, 0:Py, 1:Pz), this%chZ(1:Px, 1:Py, 0:Pz));
         this%bounds(1,1,1)=0;
         this%bounds(1,1,2)=Px;
         this%bounds(1,2,1)=1;
         this%bounds(1,2,2)=Py;
         this%bounds(1,3,1)=1;
         this%bounds(1,3,2)=Pz;
         
         this%bounds(2,1,1)=1;
         this%bounds(2,1,2)=Px;
         this%bounds(2,2,1)=0;
         this%bounds(2,2,2)=Py;
         this%bounds(2,3,1)=1;
         this%bounds(2,3,2)=Pz;
         
         this%bounds(3,1,1)=1;
         this%bounds(3,1,2)=Px;
         this%bounds(3,2,1)=1;
         this%bounds(3,2,2)=Py;
         this%bounds(3,3,1)=0;
         this%bounds(3,3,2)=Pz;
      end select
      call this%mesh_clear;
      call this%mesh_clearChA;
    end subroutine mesh_init;

!MESH_DESTROY------------------------------------------------------------------   
    subroutine mesh_destroy(this) 
      !Destroy Mesh structure and deallocate the memor
      class(tMesh) :: this
      deallocate(this%x, this%y, this%z, this%chx, this%chy, this%chz, this%bounds);
    end subroutine mesh_destroy;


!Mesh_CheckPointBounds------------------------------------------------------------------  
    integer function mesh_checkpoint(this, Px, Py, Pz, component) 
    !Clear check arrays
      class(tMesh) :: this
      integer :: Px, Py, Pz, component, res
      res = 0;
      if ((Px<this%bounds(component,1,1)).OR.(Px>this%bounds(component,1,2))) then
         res = 1;
      endif
      if ((Py<this%bounds(component,2,1)).OR.(Py>this%bounds(component,2,2))) then
         res = 1;
      endif
      if ((Pz<this%bounds(component,3,1)).OR.(Pz>this%bounds(component,3,2))) then
         res = 1;
      endif
      mesh_checkpoint = res;
    end function mesh_checkpoint;


!Mesh_Clear------------------------------------------------------------------      
    subroutine mesh_clear(this) 
    !Clear check arrays
      class(tMesh) :: this
      integer :: i,j,k;
      this%X(:,:,:)=0.0d0;   
      this%Y(:,:,:)=0.0d0; 
      this%Z(:,:,:)=0.0d0;
    end subroutine mesh_clear;
    

!Mesh_ClearCheckArrays------------------------------------------------------------------      
    subroutine mesh_clearChA(this) 
    !Clear check arrays
      class(tMesh) :: this
      integer :: i,j,k;
       do k=this%bounds(1,3,1),this%bounds(1,3,2)
         do j=this%bounds(1,2,1),this%bounds(1,2,2)
            do i=this%bounds(1,1,1),this%bounds(1,1,2)
               this%chX(i,j,k)=0;   
            enddo
         enddo
      enddo
      do k=this%bounds(2,3,1),this%bounds(2,3,2)
         do j=this%bounds(2,2,1),this%bounds(2,2,2)
            do i=this%bounds(2,1,1),this%bounds(2,1,2)
               this%chY(i,j,k)=0; 
            enddo
         enddo
      enddo
      do k=this%bounds(3,3,1),this%bounds(3,3,2)
         do j=this%bounds(3,2,1),this%bounds(3,2,2)
            do i=this%bounds(3,1,1),this%bounds(3,1,2)
               this%chZ(i,j,k)=0; 
            enddo
         enddo
      enddo
    end subroutine mesh_clearChA;

!Mesh_CheckArrays------------------------------------------------------------------         
    subroutine mesh_checkarrays(this) 
    !Clear check the update results of the mesh
      class(tMesh) :: this
      integer :: i,j,k;
      character :: fld;
      select case (this%shifttype)
      case (0)
         fld='E';
      case (1)
         fld='H';
      end select
      do k=this%bounds(1,3,1),this%bounds(1,3,2)
         do j=this%bounds(1,2,1),this%bounds(1,2,2)
            do i=this%bounds(1,1,1),this%bounds(1,1,2)
               if (this%chX(i,j,k)==0) then
                  write(*,fmt=41) this%holder, fld, i, j, k
               endif
               if (this%chX(i,j,k)>1) then
                 write(*,fmt=44) this%holder, fld, this%chX(i,j,k), i, j, k
               endif        
            enddo
         enddo
      enddo
      do k=this%bounds(2,3,1),this%bounds(2,3,2)
         do j=this%bounds(2,2,1),this%bounds(2,2,2)
            do i=this%bounds(2,1,1),this%bounds(2,1,2)
               if (this%chY(i,j,k)==0) then
                  write(*,fmt=42) this%holder, fld, i, j, k
               endif
               if (this%chY(i,j,k)>1) then
                  write(*,fmt=45) this%holder, fld, this%chY(i,j,k), i, j, k
               endif
            enddo
         enddo
      enddo
      do k=this%bounds(3,3,1),this%bounds(3,3,2)
         do j=this%bounds(3,2,1),this%bounds(3,2,2)
            do i=this%bounds(3,1,1),this%bounds(3,1,2)
               if (this%chZ(i,j,k)==0) then
                  write(*,fmt=43) this%holder, fld, i, j, k
               endif
               if (this%chZ(i,j,k)>1) then
                  write(*,fmt=46) this%holder, fld, this%chZ(i,j,k), i, j, k
               endif
            enddo
         enddo
      enddo
41    format ('Problem:',I2,2x,A,'X component is not updated!: (',I3,',',I3,',',I3,')');
42    format ('Problem:',I2,2x,A,'Y component is not updated!: (',I3,',',I3,',',I3,')');
43    format ('Problem:',I2,2x,A,'Z component is not updated!: (',I3,',',I3,',',I3,')');
44    format ('Problem:',I2,2x,A,'X component is updated ',I3,' times! (',I3,',',I3,',',I3,')');      
45    format ('Problem:',I2,2x,A,'Y component is updated ',I3,' times! (',I3,',',I3,',',I3,')');    
46    format ('Problem:',I2,2x,A,'Z component is updated ',I3,' times! (',I3,',',I3,',',I3,')');         
    end subroutine mesh_checkarrays;    
    
end module meshclass
