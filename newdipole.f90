!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains implementation of divergence free source

module newdipole
! Computate the analytic solution of the dipole
  use commonvars
  use meshclass
  
  implicit none;
  type, public :: tNewDipole
   contains
     procedure getNewDipoleEpoint;
     procedure getNewDipoleHpoint;
     procedure getNewDipoleJpoint;
  end type tNewDipole
contains

!getF------------------------------------------
    real*8 function getF(t)
    ! Calculate the F      
      real*8, intent(in) :: t
      getF=a*sin(v*t)
    end function getF

!getDF------------------------------------------
    real*8 function getDF(t,type)
    ! Calculate the derivative of F  
      real*8, intent(in) :: t
      integer, intent(in) :: type
      getDF = 0;
      if (type==91) then            ! by t 
         getDF = a*cos(v*t)*v;
      else if (type == 92) then      ! by 2t 
         getDF = -a*sin(v*t)*(v**2);
      else if (type == 93) then      ! by 3t 
         getDF = -a*cos(v*t)*(v**3);
      else
         write(*,*) 'Alert: getDF get an unexpected type:', type 
      end if
    end function getDF

!getNu------------------------------------------
    real*8 function getNu(t)
    ! Calculate the Nu 
      real*8, intent(in) :: t
      getNu = 0.0d0;
      if ((q*t>0.0d0).and.(q*t<1.0d0)) then           
         getNu = -3432.0d0*(t*q)**(15)+25740.0d0*(t*q)**(14)-83160.0d0*(t*q)**(13)+150150.0d0*(t*q)**(12)-163800.0d0*(t*q)**(11)+108108.0d0*(t*q)**(10)-40040.0d0*(t*q)**(9)+6435.0d0*(t*q)**(8)
      else if (q*t<= 0.0d0 ) then     
         getNu = 0.0d0;
      else if (q*t >= 1.0d0) then   
         getNu = 1.0d0;
      end if
    end function getNu
    
!getDNu------------------------------------------
    real*8 function getDNu(t, type)
    ! Calculate the derivative of Nu  
      real*8, intent(in) :: t
      integer, intent(in) :: type
      getDNu = 0.0d0;
      if (type==91) then            ! by t 
        if ((q*t>0.0d0).and.(q*t<1.0d0)) then           
           getDNu = -51480.0d0*(q**15)*(t**14)+360360.0d0*(q**14)*(t**13)-1081080.0d0*(q**13)*(t**12)+1801800.0d0*(q**12)*(t**11)-1801800.0d0*(q**11)*(t**10)+1081080.0d0*(q**10)*(t**9)-360360.0d0*(q**9)*(t**8)+51480.0d0*(q**8)*(t**7) 
        else      
           getDNu = 0.0d0
        end if
      else if (type == 92) then      ! by 2t 
        if ((q*t>0.0d0).and.(q*t<1.0d0)) then           
           getDNu = -720720.0d0*(q**15)*(t**13)+4684680.0d0*(q**14)*(t**12)-12972960.0d0*(q**13)*(t**11)+19819800.0d0*(q**12)*(t**10)-18018000.0d0*(q**11)*(t**9)+9729720.0d0*(q**10)*(t**8)-2882880.0d0*(q**9)*(t**7)+360360.0d0*(q**8)*(t**6) 
        else      
           getDNu = 0.0d0
        end if
      else if (type == 93) then      ! by 3t 
         if ((q*t>0.0d0).and.(q*t<1.0d0)) then           
           getDNu = -9369360.0d0*(q**15)*(t**12)+56216160.0d0*(q**14)*(t**11)-142702560.0d0*(q**13)*(t**10)+198198000.0d0*(q**12)*(t**9)-162162000.0d0*(q**11)*(t**8)+77837760.0d0*(q**10)*(t**7)-20180160.0d0*(q**9)*(t**6)+2162160.0d0*(q**8)*(t**5)
         else      
           getDNu = 0.0d0
          end if
      else      
         write(*,*) 'Alert: getDNu get an unexpected type:', type 
      end if
    end function getDNu

!getM------------------------------------------
   real*8 function getM(x,y,z)
    ! Calculate M function  
      real*8, intent(in) :: x, y, z
      real*8 :: R
      R=sqrt(x**2+y**2+z**2)
      getM=0.0d0
      if (R*p<1.0d0) then
         getM = -3432.0d0*(R*p)**(15)+25740.0d0*(R*p)**(14)-83160.0d0*(R*p)**(13)+150150.0d0*(R*p)**(12)-163800.0d0*(R*p)**(11)+108108.0d0*(R*p)**(10)-40040.0d0*(R*p)**(9)+6435.0d0*(R*p)**(8);
      else if (R*p>=1.0d0) then
         getM = 1.0d0
      end if
    end function getM

!getDM------------------------------------------
    real*8 function getDM(x,y,z,type)
    ! Calculate derivatives of M function  
      real*8, intent(in) :: x, y, z
      integer, intent(in) :: type
      real*8 :: R, detm
      R=sqrt(x**2+y**2+z**2)  
      getDM=0.0d0
      if (R*p<1.0d0) then           
         if (type==1) then             !by X
           detm = -51480.0d0*(p**15)*R**(13)+360360.0d0*(p**14)*R**(12)-1081080.0d0*(p**13)*R**(11)+1801800.0d0*(p**12)*R**(10)-1801800.0d0*(p**11)*R**(9)+1081080.0d0*(p**10)*R**(8)-360360.0d0*(p**9)*R**(7)+51480.0d0*(p**8)*(R**6) 
           getDM = x*detm
         elseif (type==2) then          !by Y
           detm = -51480.0d0*(p**15)*R**(13)+360360.0d0*(p**14)*R**(12)-1081080.0d0*(p**13)*R**(11)+1801800.0d0*(p**12)*R**(10)-1801800.0d0*(p**11)*R**(9)+1081080.0d0*(p**10)*R**(8)-360360.0d0*(p**9)*R**(7)+51480.0d0*(p**8)*(R**6)
           getDM = y*detm
         elseif (type==3) then          !by Z
           detm = -51480.0d0*(p**15)*R**(13)+360360.0d0*(p**14)*R**(12)-1081080.0d0*(p**13)*R**(11)+1801800.0d0*(p**12)*R**(10)-1801800.0d0*(p**11)*R**(9)+1081080.0d0*(p**10)*R**(8)-360360.0d0*(p**9)*R**(7)+51480.0d0*(p**8)*(R**6) 
           getDM = z*detm
         elseif (type==31) then          !by 2zx
           detm = -669240.0d0*(p**15)*(R**11)+4324320.0d0*(p**14)*(R**10)-11891880.0d0*(p**13)*(R**9)+18018000.0d0*(p**12)*(R**8)-16216200.0d0*(p**11)*(R**7)+8648640.0d0*(p**10)*(R**6)-2522520.0d0*(p**9)*(R**5)+308880.0d0*(p**8)*(R**4)
           getDM = x*z*detm
         elseif (type==32) then          !by 2zy
           detm = -669240.0d0*(p**15)*(R**11)+4324320.0d0*(p**14)*(R**10)-11891880.0d0*(p**13)*(R**9)+18018000.0d0*(p**12)*(R**8)-16216200.0d0*(p**11)*(R**7)+8648640.0d0*(p**10)*(R**6)-2522520.0d0*(p**9)*(R**5)+308880.0d0*(p**8)*(R**4)
           getDM = z*y*detm  
         elseif (type==21) then         !by 2yx
           detm = -669240.0d0*(p**15)*(R**11)+4324320.0d0*(p**14)*(R**10)-11891880.0d0*(p**13)*(R**9)+18018000.0d0*(p**12)*(R**8)-16216200.0d0*(p**11)*(R**7)+8648640.0d0*(p**10)*(R**6)-2522520.0d0*(p**9)*(R**5)+308880.0d0*(p**8)*(R**4)
           getDM = x*y*detm
         elseif (type==11) then         !by 2x
           getDM = -669240.0d0*(R**11)*(p**15)*(x**2)-51480.0d0*(R**13)*(p**15)+4324320.0d0*(R**10)*(p**14)*(x**2)+360360.0d0*(R**12)*(p**14)-11891880.0d0*(R**9)*(p**13)*(x**2)-1081080.0d0*(R**11)*(p**13)+18018000.0d0*(R**8)*(p**12)*(x**2)+1801800.0d0*(R**10)*(p**12)-16216200.0d0*(R**7)*(p**11)*(x**2)-1801800.0d0*(R**9)*(p**11)+8648640.0d0*(R**6)*(p**10)*(x**2)+1081080.0d0*(R**8)*(p**10)-2522520.0d0*(R**5)*(p**9)*(x**2)-360360.0d0*(R**7)*(p**9)+308880.0d0*(R**4)*(p**8)*(x**2)+51480.0d0*(R**6)*(p**8)  
         elseif (type==22) then         !by 2y
           getDM = -669240.0d0*(R**11)*(p**15)*(y**2)-51480.0d0*(R**13)*(p**15)+4324320.0d0*(R**10)*(p**14)*(y**2)+360360.0d0*(R**12)*(p**14)-11891880.0d0*(R**9)*(p**13)*(y**2)-1081080.0d0*(R**11)*(p**13)+18018000.0d0*(R**8)*(p**12)*(y**2)+1801800.0d0*(R**10)*(p**12)-16216200.0d0*(R**7)*(p**11)*(y**2)-1801800.0d0*(R**9)*(p**11)+8648640.0d0*(R**6)*(p**10)*(y**2)+1081080.0d0*(R**8)*(p**10)-2522520.0d0*(R**5)*(p**9)*(y**2)-360360.0d0*(R**7)*(p**9)+308880.0d0*(R**4)*(p**8)*(y**2)+51480.0d0*(R**6)*(p**8)
         elseif (type==33) then         !by 2z
           getDM = -669240.0d0*(R**11)*(p**15)*(z**2)-51480.0d0*(R**13)*(p**15)+4324320.0d0*(R**10)*(p**14)*(z**2)+360360.0d0*(R**12)*(p**14)-11891880.0d0*(R**9)*(p**13)*(z**2)-1081080.0d0*(R**11)*(p**13)+18018000.0d0*(R**8)*(p**12)*(z**2)+1801800.0d0*(R**10)*(p**12)-16216200.0d0*(R**7)*(p**11)*(z**2)-1801800.0d0*(R**9)*(p**11)+8648640.0d0*(R**6)*(p**10)*(z**2)+1081080.0d0*(R**8)*(p**10)-2522520.0d0*(R**5)*(p**9)*(z**2)-360360.0d0*(R**7)*(p**9)+308880.0d0*(R**4)*(p**8)*(z**2)+51480.0d0*(R**6)*(p**8)
          elseif (type==222) then       !by 3y
           getDM = -7361640.0d0*(R**9)*(p**15)*(y**3)-2007720.0d0*(R**11)*(p**15)*y+43243200.0d0*(R**8)*(p**14)*(y**3)+12972960.0d0*(R**10)*(p**14)*y-107026920.0d0*(R**7)*(p**13)*(y**3)-35675640.0d0*(R**9)*(p**13)*y+144144000.0d0*(R**6)*(p**12)*(y**3)+54054000.0d0*(R**8)*(p**12)*y-113513400.0d0*(R**5)*(p**11)*(y**3)-48648600.0d0*(R**7)*(p**11)*y+51891840.0d0*(R**4)*(p**10)*(y**3)+25945920.0d0*(R**6)*(p**10)*y-12612600.0d0*(R**3)*(p**9)*(y**3)-7567560.0d0*(R**5)*(p**9)*y+1235520.0d0*(R**2)*(p**8)*(y**3)+926640.0d0*(R**4)*(p**8)*y
          elseif (type==333) then       !by 3z
           getDM = -7361640.0d0*(R**9)*(p**15)*(z**3)-2007720.0d0*(R**11)*(p**15)*z+43243200.0d0*(R**8)*(p**14)*(z**3)+12972960.0d0*(R**10)*(p**14)*z-107026920.0d0*(R**7)*(p**13)*(z**3)-35675640.0d0*(R**9)*(p**13)*z+144144000.0d0*(R**6)*(p**12)*(z**3)+54054000.0d0*(R**8)*(p**12)*z-113513400.0d0*(R**5)*(p**11)*(z**3)-48648600.0d0*(R**7)*(p**11)*z+51891840.0d0*(R**4)*(p**10)*(z**3)+25945920.0d0*(R**6)*(p**10)*z-12612600.0d0*(R**3)*(p**9)*(z**3)-7567560.0d0*(R**5)*(p**9)*z+1235520.0d0*(R**2)*(p**8)*(z**3)+926640.0d0*(R**4)*(p**8)*z
          elseif (type==322) then       !by zyy
           getDM = -7361640.0d0*(R**9)*(p**15)*z*(y**2)-669240.0d0*(R**11)*(p**15)*z+43243200.0d0*(R**8)*(p**14)*z*(y**2)+4324320.0d0*(R**10)*(p**14)*z-107026920.0d0*(R**7)*(p**13)*z*(y**2)-11891880.0d0*(R**9)*(p**13)*z+144144000.0d0*(R**6)*(p**12)*z*(y**2)+18018000.0d0*(R**8)*(p**12)*z-113513400.0d0*(R**5)*(p**11)*z*(y**2)-16216200.0d0*(R**7)*(p**11)*z+51891840.0d0*(R**4)*(p**10)*z*(y**2)+8648640.0d0*(R**6)*(p**10)*z-12612600.0d0*(R**3)*(p**9)*z*(y**2)-2522520.0d0*(R**5)*(p**9)*z+1235520.0d0*(R**2)*(p**8)*z*(y**2)+308880.0d0*(R**4)*(p**8)*z
          elseif (type==311) then       !by zxx
           getDM = -7361640.0d0*(R**9)*(p**15)*z*(x**2)-669240.0d0*(R**11)*(p**15)*z+43243200.0d0*(R**8)*(p**14)*z*(x**2)+4324320.0d0*(R**10)*(p**14)*z-107026920.0d0*(R**7)*(p**13)*z*(x**2)-11891880.0d0*(R**9)*(p**13)*z+144144000.0d0*(R**6)*(p**12)*z*(x**2)+18018000.0d0*(R**8)*(p**12)*z-113513400.0d0*(R**5)*(p**11)*z*(x**2)-16216200.0d0*(R**7)*(p**11)*z+51891840.0d0*(R**4)*(p**10)*z*(x**2)+8648640.0d0*(R**6)*(p**10)*z-12612600.0d0*(R**3)*(p**9)*z*(x**2)-2522520.0d0*(R**5)*(p**9)*z+1235520.0d0*(R**2)*(p**8)*z*(x**2)+308880.0d0*(R**4)*(p**8)*z
          elseif (type==211) then       !by yxx
             getDM = -7361640.0d0*(R**9)*(p**15)*y*(x**2)-669240.0d0*(R**11)*(p**15)*y+43243200.0d0*(R**8)*(p**14)*y*(x**2)+4324320.0d0*(R**10)*(p**14)*y-107026920.0d0*(R**7)*(p**13)*y*(x**2)-11891880.0d0*(R**9)*(p**13)*y+144144000.0d0*(R**6)*(p**12)*y*(x**2)+18018000.0d0*(R**8)*(p**12)*y-113513400.0d0*(R**5)*(p**11)*y*(x**2)-16216200.0d0*(R**7)*(p**11)*y+51891840.0d0*(R**4)*(p**10)*y*(x**2)+8648640.0d0*(R**6)*(p**10)*y-12612600.0d0*(R**3)*(p**9)*y*(x**2)-2522520.0d0*(R**5)*(p**9)*y+1235520.0d0*(R**2)*(p**8)*y*(x**2)+308880.0d0*(R**4)*(p**8)*y
          elseif (type==332) then       !by zzy
             getDM = -7361640.0d0*(R**9)*(p**15)*y*(z**2)-669240.0d0*(R**11)*(p**15)*y+43243200.0d0*(R**8)*(p**14)*y*(z**2)+4324320.0d0*(R**10)*(p**14)*y-107026920.0d0*(R**7)*(p**13)*y*(z**2)-11891880.0d0*(R**9)*(p**13)*y+144144000.0d0*(R**6)*(p**12)*y*(z**2)+18018000.0d0*(R**8)*(p**12)*y-113513400.0d0*(R**5)*(p**11)*y*(z**2)-16216200.0d0*(R**7)*(p**11)*y+51891840.0d0*(R**4)*(p**10)*y*(z**2)+8648640.0d0*(R**6)*(p**10)*y-12612600.0d0*(R**3)*(p**9)*y*(z**2)-2522520.0d0*(R**5)*(p**9)*y+1235520.0d0*(R**2)*(p**8)*y*(z**2)+308880.0d0*(R**4)*(p**8)*y
          else
             write(*,*) 'Alert: getDM get an unexpected type:', type
          end if
      else
         getDM = 0.0d0
      end if
    end function getDM

!getHi------------------------------------------
    real*8 function getHi(t)
    ! Calculate Hi function  
      real*8, intent(in) :: t
      getHi = getF(t)*getNu(t)
    end function getHi

!getDHi------------------------------------------
    real*8 function getDHi(t, type)
    ! Calculate the derivative of Nu  
      real*8, intent(in) :: t
      integer, intent(in) :: type
      getDHi = 0;
      if (type==91) then            ! by t 
         getDHi = getDF(t,91)*getNu(t)+getF(t)*getDNu(t,91)
      else if (type == 92) then      ! by 2t 
         getDHi = getDNu(t,92)*getF(t)+2*getDNu(t,91)*getDF(t,91)+getNu(t)*getDF(t,92)
      else if (type == 93) then      ! by 3t 
         getDHi = getDNu(t,93)*getF(t)+3*getDNu(t,92)*getDF(t,91)+3*getDNu(t,91)*getDF(t,92)+getNu(t)*getDF(t,93)
      else
         write(*,*) 'Alert: getDHi get an unexpected type:', type
      end if
    end function getDHi

!getW------------------------------------------
    real*8 function getW(x,y,z,t)
    ! Calculate W function  
      real*8, intent(in) :: x,y,z,t
      real*8 :: R, ct
      R = sqrt(x**2 + y**2 + z**2)
      ct = t-R/cc
      if (ct >= 0.0d0) then
         getW = -1.0d0/R*(getHi(ct)/(R**2)+1.0d0/(cc*R)*getDHi(ct,91))  
      else
         getW = 0.0d0
      end if
    end function getW

!getDW------------------------------------------
    real*8 function getDW(x,y,z,t,type)
    ! Calculate derivatives of M function  
      real*8, intent(in) :: x, y, z, t
      integer, intent(in) :: type
      real*8 :: R, detw, ct
      R=sqrt(x**2+y**2+z**2)
      ct = t-R/cc
      getDW=0
      if (type==91) then             !by t
          if (ct > 0.0d0) then
             getDW = -1.0d0/R*(getDHi(ct,91)/(R**2)+getDHi(ct,92)/(R*cc));
          end if    
      elseif (type==92) then          !by 2t
         if (ct > 0.0d0) then
           getDW = -1.0d0/R*(getDHi(ct,92)/(R**2)+getDHi(ct,93)/(R*cc));
         end if 
      elseif (type==1) then          !by x
         if (ct > 0) then
           detw = 1.0d0/(R**3)*(getHi(ct)/(R**2)+getDHi(ct,91)/(R*cc))+1.0d0/R*(2.0d0*getDHi(ct,91)/((R**3)*cc)+2.0d0*getHi(ct)/(R**4)+getDHi(ct,92)/((cc**2)*(R**2)));
           getDW = x*detw
         end if 
      elseif (type==2) then          !by y
         if (ct > 0.0d0) then
           detw = 1.0d0/(R**3)*(getHi(ct)/(R**2)+getDHi(ct,91)/(R*cc))+1.0d0/R*(2.0d0*getDHi(ct,91)/((R**3)*cc)+2.0d0*getHi(ct)/(R**4)+getDHi(ct,92)/((cc**2)*(R**2)));
           getDW = y*detw
         end if
      elseif (type==3) then           !by z
         if (ct > 0.0d0) then
           detw = 1.0d0/(R**3)*(getHi(ct)/(R**2)+getDHi(ct,91)/(R*cc))+1.0d0/R*(2.0d0*getDHi(ct,91)/((R**3)*cc)+2.0d0*getHi(ct)/(R**4)+getDHi(ct,92)/((cc**2)*(R**2)));
           getDW = z*detw
         end if
      elseif (type==31) then         !by 2zx
         if (ct > 0.0d0) then
           detw = -3.0d0/(R**5)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))+2.0d0/(R**3)*(-2.0d0*getDHi(ct,91)/((R**3)*cc)-2.0d0*getHi(ct)/(R**4)-getDHi(ct,92)/((cc**2)*(R**2)))-1.0d0/R*(4.0d0*getDHi(ct,92)/((R**4)*(cc**2))+8.0d0*getDHi(ct,91)/((R**5)*cc)+8.0d0*getHi(ct)/(R**6)+getDHi(ct,93)/((cc**3)*(R**3)))
           getDW = z*x*detw
        end if
      elseif (type==32) then         !by 2zy
         if (ct > 0.0d0) then
           detw = -3.0d0/(R**5)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))+2.0d0/(R**3)*(-2.0d0*getDHi(ct,91)/((R**3)*cc)-2.0d0*getHi(ct)/(R**4)-getDHi(ct,92)/((cc**2)*(R**2)))-1.0d0/R*(4*getDHi(ct,92)/((R**4)*(cc**2))+8.0d0*getDHi(ct,91)/((R**5)*cc)+8.0d0*getHi(ct)/(R**6)+getDHi(ct,93)/((cc**3)*(R**3)))
           getDW = z*y*detw
        end if  
      elseif (type==21) then         !by 2yx
         if (ct > 0.0d0) then
           detw = -3.0d0/(R**5)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))+2.0d0/(R**3)*(-2.0d0*getDHi(ct,91)/((R**3)*cc)-2.0d0*getHi(ct)/(R**4)-getDHi(ct,92)/((cc**2)*(R**2)))-1.0d0/R*(4*getDHi(ct,92)/((R**4)*(cc**2))+8.0d0*getDHi(ct,91)/((R**5)*cc)+8.0d0*getHi(ct)/(R**6)+getDHi(ct,93)/((cc**3)*(R**3)))
           getDW = y*x*detw
        end if
      elseif (type==11) then         !by 2x
         if (ct > 0.0d0) then
           getDW =  -3.0d0*(x**2)/(R**5)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))-2.0d0*(x**2)/(R**3)*(2.0d0*getDHi(ct,91)/((R**3)*cc)+2.0d0*getHi(ct)/(R**4)+getDHi(ct,92)/((cc**2)*(R**2)))+1.0d0/(R**3)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))-1.0d0/R*(4.0d0*getDHi(ct,92)*(x**2)/((R**4)*(cc**2))+8.0d0*getDHi(ct,91)*(x**2)/((R**5)*cc)-2.0d0*getDHi(ct,91)/((R**3)*cc)+8.0d0*getHi(ct)*(x**2)/(R**6)-2.0d0*getHi(ct)/(R**4)+getDHi(ct,93)*(x**2)/((cc**3)*(R**3))-getDHi(ct,92)/((cc**2)*(R**2))) 
         end if  
      elseif (type==22) then         !by 2y
         if (ct > 0.0d0) then
           getDW =  -3.0d0*(y**2)/(R**5)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))-2.0d0*(y**2)/(R**3)*(2.0d0*getDHi(ct,91)/((R**3)*cc)+2.0d0*getHi(ct)/(R**4)+getDHi(ct,92)/((cc**2)*(R**2)))+1.0d0/(R**3)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))-1.0d0/R*(4.0d0*getDHi(ct,92)*(y**2)/((R**4)*(cc**2))+8.0d0*getDHi(ct,91)*(y**2)/((R**5)*cc)-2.0d0*getDHi(ct,91)/((R**3)*cc)+8.0d0*getHi(ct)*(y**2)/(R**6)-2.0d0*getHi(ct)/(R**4)+getDHi(ct,93)*(y**2)/((cc**3)*(R**3))-getDHi(ct,92)/((cc**2)*(R**2))) 
         end if
      elseif (type==33) then         !by 2z
         if (ct > 0.0d0) then
           getDW =  -3.0d0*(z**2)/(R**5)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))-2.0d0*(z**2)/(R**3)*(2.0d0*getDHi(ct,91)/((R**3)*cc)+2.0d0*getHi(ct)/(R**4)+getDHi(ct,92)/((cc**2)*(R**2)))+1.0d0/(R**3)*(getHi(ct)/(R**2)+getDHi(ct,91)/(cc*R))-1.0d0/R*(4.0d0*getDHi(ct,92)*(z**2)/((R**4)*(cc**2))+8.0d0*getDHi(ct,91)*(z**2)/((R**5)*cc)-2.0d0*getDHi(ct,91)/((R**3)*cc)+8.0d0*getHi(ct)*(z**2)/(R**6)-2.0d0*getHi(ct)/(R**4)+getDHi(ct,93)*(z**2)/((cc**3)*(R**3))-getDHi(ct,92)/((cc**2)*(R**2))) 
         end if
       elseif (type==39) then         !by zt
         if (ct > 0.0d0) then
           detw = 1.0d0/(R**3)*(getDHi(ct,91)/(R**2)+getDHi(ct,92)/(R*cc))+1.0d0/R*(2.0d0*getDHi(ct,92)/((R**3)*cc)+2.0d0*getDHi(ct,91)/(R**4)+getDHi(ct,93)/((cc**2)*(R**2)));
           getDW = z*detw 
         end if
       elseif (type==29) then         !by yt
         if (ct > 0.0d0) then
           detw = 1.0d0/(R**3)*(getDHi(ct,91)/(R**2)+getDHi(ct,92)/(R*cc))+1.0d0/R*(2.0d0*getDHi(ct,92)/((R**3)*cc)+2.0d0*getDHi(ct,91)/(R**4)+getDHi(ct,93)/((cc**2)*(R**2)));
           getDW = y*detw 
        end if
      else
        write(*,*) 'Alert: getDW get an unexpected type:', type
      end if
    end function getDW  

!getNewDipoleEpoint------------------------------------------
    real*8 function getNewDipoleEpoint(this, x, y, z, xp, yp, zp, t, type, problem_type)
      ! Calculate Electric pole at the time-point
      class(tNewDipole) :: this;
      real*8, intent(in) :: x, y, z;
      integer, intent(in) :: xp, yp, zp;
      real*8, intent(in) :: t;
      integer, intent(in) :: type;
      integer, intent(in) :: problem_type;
      real*8 :: CgetW, getNuCT, getDHiCT91, getDHiCT92, getDHiCT93, getDNu91CT, getDNu92CT, getDNu93CT, R, ct,  detw;
      real*8 :: getDWCT91, getDWCT39, getDWCT29;
      real*8 :: ct5, ct6, ct7, ct8, ct9, ct10, ct11, ct12, ct13, ct14, ct15, cosvct, sinvct;
      if (type/=1) then
         R = sqrt(x**2 + y**2 + z**2)
         ct = t-R/cc
         ct5 = ct**5;
         ct6 = ct**6;
         ct7 = ct**7;
         ct8 = ct**8;
         ct9 = ct**9;
         ct10 = ct**10;
         ct11 = ct**11;
         ct12 = ct**12;
         ct13 = ct**13;
         ct14 = ct**14;
         ct15 = ct**15;
         
         getNuCT = 0.0d0;
         if ((q*ct>0.0d0).and.(q*ct<1.0d0)) then           
            getNuCT = cNu(0)*(ct15)+cNu(1)*(ct14)+cNu(2)*(ct13)+cNu(3)*(ct12)+cNu(4)*(ct11)+cNu(5)*(ct10)+cNu(6)*(ct9)+cNu(7)*(ct8);
         else if (q*ct<= 0.0d0 ) then     
            getNuCT = 0.0d0;
         else if (q*ct >= 1.0d0) then   
            getNuCT = 1.0d0;
         end if
       
         if ((q*ct>0.0d0).and.(q*ct<1.0d0)) then           
            getDNu91CT = cDNu(0)*(ct14)+cDNu(1)*(ct13)+cDNu(2)*(ct12)+cDNu(3)*(ct11)+cDNu(4)*(ct10)+cDNu(5)*(ct9)+cDNu(6)*(ct8)+cDNu(7)*(ct7)
            getDNu92CT = cDNu(8)*(ct13)+cDNu(9)*(ct12)+cDNu(10)*(ct11)+cDNu(11)*(ct10)+cDNu(12)*(ct9)+cDNu(13)*(ct8)+cDNu(14)*(ct7)+cDNu(15)*(ct6)
            getDNu93CT = cDNu(16)*(ct12)+cDNu(17)*(ct11)+cDNu(18)*(ct10)+cDNu(19)*(ct9)+cDNu(20)*(ct8)+cDNu(21)*(ct7)+cDNu(22)*(ct6)+cDNu(23)*(ct5)
         else      
            getDNu91CT = 0.0d0
            getDNu92CT = 0.0d0
            getDNu93CT = 0.0d0
         end if
         cosvct = cos(v*ct);
         sinvct = sin(v*ct);
         getDHiCT91 = a*cosvct*v*getNuCT+a*sinvct*getDNu91CT;
         getDHiCT92 = getDNu92CT*a*sinvct+2*getDNu91CT*a*cosvct*v-getNuCT*a*sinvct*(v**2);
         getDHiCT93 = getDNu93CT*a*sinvct+3*getDNu92CT*a*cosvct*v-3*getDNu91CT*a*sinvct*(v**2)-getNuCT*a*cosvct*(v**3);

         getDWCT91 = 0.0d0;
         getDWCT39 = 0.0d0;
         getDWCT29 = 0.0d0;
         if (ct > 0.0d0) then
            getDWCT91 = -1.0d0/R*(getDHiCT91/(R**2)+getDHiCT92/(R*cc));
            detw = 1.0d0/(R**3)*(getDHiCT91/(R**2)+getDHiCT92/(R*cc))+1.0d0/R*(2.0d0*getDHiCT92/((R**3)*cc)+2.0d0*getDHiCT91/(R**4)+getDHiCT93/((cc**2)*(R**2)));
            getDWCT39 = z*detw 
            getDWCT29 = y*detw 
         end if   
      endif
      getNewDipoleEpoint = 0.0d0;

      if (problem_type == 0) then
         if (type==1) then            ! EX
            getNewDipoleEpoint = 0.0d0
         else if (type == 2) then      ! EY
            getNewDipoleEpoint = -x/cc*(MEJy(xp,yp,zp)*getDWCT39+DMEJy(2,xp,yp,zp)*getDWCT91);
         else if (type == 3) then      ! EZ
            getNewDipoleEpoint = x/cc*(DMEJz(1,xp,yp,zp)*getDWCT91+MEJz(xp,yp,zp)*getDWCT29);
         else
            write(*,*) 'Alert: getEpoint get an unexpected type:', type
         end if
      else if (problem_type == 1) then
         if (type==1) then            ! EX
            getNewDipoleEpoint = 0.0d0
         else if (type == 2) then      ! EY
            getNewDipoleEpoint = -x/cc*(auxMEJy(xp,yp,zp)*getDWCT39+auxDMEJy(2,xp,yp,zp)*getDWCT91);
         else if (type == 3) then      ! EZ
            getNewDipoleEpoint = x/cc*(auxDMEJz(1,xp,yp,zp)*getDWCT91+auxMEJz(xp,yp,zp)*getDWCT29);
         else
            write(*,*) 'Alert: getEpoint get an unexpected type:', type
         end if
      endif   
    end function getNewDipoleEpoint

    
!Hpoint------------------------------------------
    real*8 function getNewDipoleHpoint(this, x, y, z, xp, yp, zp, t, type, problem_type)
      ! Calculate Magnetic pole at the time-point
      class(tNewDipole) :: this;
      real*8, intent(in) :: x, y, z;
      integer, intent(in) :: xp, yp, zp;
      real*8, intent(in) :: t
      integer, intent(in) :: type
      integer, intent(in) :: problem_type;
      real*8 :: CgetW, getNuCT, getDHiCT91, getDHiCT92, getDHiCT93, getDNu91CT, getDNu92CT, getDNu93CT, R, ct,  detw;
      real*8 :: getDWCT1, getDWCT2, getDWCT3, getDWCT31, getDWCT21, getDWCT11, getDWCT22, getDWCT33;
      real*8 :: ct5, ct6, ct7, ct8, ct9, ct10, ct11, ct12, ct13, ct14, ct15, cosvct, sinvct;
      R = sqrt(x**2 + y**2 + z**2)
      ct = t-R/cc
      getNuCT = 0.0d0;
      if ((q*ct>0.0d0).and.(q*ct<1.0d0)) then           
         getNuCT = cNu(0)*(ct**15)+cNu(1)*(ct**14)+cNu(2)*(ct**13)+cNu(3)*(ct**12)+cNu(4)*(ct**11)+cNu(5)*(ct**10)+cNu(6)*(ct**9)+cNu(7)*(ct**8);
      else if (q*ct<= 0.0d0 ) then     
         getNuCT = 0.0d0;
      else if (q*ct >= 1.0d0) then   
         getNuCT = 1.0d0;
      end if
       
      if ((q*ct>0.0d0).and.(q*ct<1.0d0)) then           
         getDNu91CT = cDNu(0)*(ct**14)+cDNu(1)*(ct**13)+cDNu(2)*(ct**12)+cDNu(3)*(ct**11)+cDNu(4)*(ct**10)+cDNu(5)*(ct**9)+cDNu(6)*(ct**8)+cDNu(7)*(ct**7)
         getDNu92CT = cDNu(8)*(ct**13)+cDNu(9)*(ct**12)+cDNu(10)*(ct**11)+cDNu(11)*(ct**10)+cDNu(12)*(ct**9)+cDNu(13)*(ct**8)+cDNu(14)*(ct**7)+cDNu(15)*(ct**6) 
         getDNu93CT = cDNu(16)*(ct**12)+cDNu(17)*(ct**11)+cDNu(18)*(ct**10)+cDNu(19)*(ct**9)+cDNu(20)*(ct**8)+cDNu(21)*(ct**7)+cDNu(22)*(ct**6)+cDNu(23)*(ct**5)
      else      
         getDNu91CT = 0.0d0
         getDNu92CT = 0.0d0
         getDNu93CT = 0.0d0
      end if

      getDHiCT91 = a*cos(v*ct)*v*getNuCT+a*sin(v*ct)*getDNu91CT;
      getDHiCT92 = getDNu92CT*a*sin(v*ct)+2*getDNu91CT*a*cos(v*ct)*v-getNuCT*a*sin(v*ct)*(v**2);
      getDHiCT93 = getDNu93CT*a*sin(v*ct)+3*getDNu92CT*a*cos(v*ct)*v-3*getDNu91CT*a*sin(v*ct)*(v**2)-getNuCT*a*cos(v*ct)*(v**3);

      if (ct >= 0.0d0) then
         CgetW = -1.0d0/R*(a*sin(v*ct)*getNuCT/(R**2)+1.0d0/(cc*R)*getDHiCT91);
      else
         CgetW = 0.0d0
      end if

      getDWCT1 = 0.0d0;
      getDWCT2 = 0.0d0;
      getDWCT3 = 0.0d0;
      getDWCT31 = 0.0d0;
      getDWCT21 = 0.0d0;
      getDWCT22 = 0.0d0;
      getDWCT33 = 0.0d0;
      if (ct > 0.0d0) then
         detw = 1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(R*cc))+1.0d0/R*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)));
         getDWCT1 = x*detw;
         getDWCT2 = y*detw;
         getDWCT3 = z*detw
             detw = -3.0d0/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))+2.0d0/(R**3)*(-2.0d0*getDHiCT91/((R**3)*cc)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)-getDHiCT92/((cc**2)*(R**2)))-1.0d0/R*(4.0d0*getDHiCT92/((R**4)*(cc**2))+8.0d0*getDHiCT91/((R**5)*cc)+8.0d0*a*sin(v*ct)*getNuCT/(R**6)+getDHiCT93/((cc**3)*(R**3)))
         getDWCT31 = z*x*detw;
         getDWCT21 = y*x*detw
         getDWCT22 =  -3.0d0*(y**2)/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-2.0d0*(y**2)/(R**3)*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)))+1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-1.0d0/R*(4.0d0*getDHiCT92*(y**2)/((R**4)*(cc**2))+8.0d0*getDHiCT91*(y**2)/((R**5)*cc)-2.0d0*getDHiCT91/((R**3)*cc)+8.0d0*a*sin(v*ct)*getNuCT*(y**2)/(R**6)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT93*(y**2)/((cc**3)*(R**3))-getDHiCT92/((cc**2)*(R**2))) 
         getDWCT33 =  -3.0d0*(z**2)/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-2.0d0*(z**2)/(R**3)*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)))+1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-1.0d0/R*(4.0d0*getDHiCT92*(z**2)/((R**4)*(cc**2))+8.0d0*getDHiCT91*(z**2)/((R**5)*cc)-2.0d0*getDHiCT91/((R**3)*cc)+8.0d0*a*sin(v*ct)*getNuCT*(z**2)/(R**6)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT93*(z**2)/((cc**3)*(R**3))-getDHiCT92/((cc**2)*(R**2))) 
      end if   

      if (problem_type == 0) then
         getNewDipoleHpoint = 0.0d0
         if (type==1) then             ! HX
            getNewDipoleHpoint = -x*(getDWCT22*MHx(xp,yp,zp)+2.0d0*getDWCT2*DMHx(1,xp,yp,zp)+CgetW*DMHx(7,xp,yp,zp)+getDWCT33*MHx(xp,yp,zp)+2.0d0*getDWCT3*DMHx(2,xp,yp,zp)+CgetW*DMHx(8,xp,yp,zp));
         else if (type == 2) then      ! HY
            getNewDipoleHpoint = getDWCT2*MHy(xp,yp,zp)+x*getDWCT21*MHy(xp,yp,zp)+x*getDWCT2*DMHy(0,xp,yp,zp)+CgetW*DMHy(1,xp,yp,zp)+x*getDWCT1*DMHy(1,xp,yp,zp)+x*CgetW*DMHy(5,xp,yp,zp);
         else if (type == 3) then      ! HZ
            getNewDipoleHpoint = getDWCT3*MHz(xp,yp,zp)+x*getDWCT31*MHz(xp,yp,zp)+x*getDWCT3*DMHz(0,xp,yp,zp)+CgetW*DMHz(2,xp,yp,zp)+x*getDWCT1*DMHz(2,xp,yp,zp)+x*CgetW*DMHz(3,xp,yp,zp);
         else
            write(*,*) 'Alert: getHpoint get an unexpected type!', type   
         end if
      else if (problem_type == 1) then
         getNewDipoleHpoint = 0.0d0
         if (type==1) then             ! HX
            getNewDipoleHpoint = -x*(getDWCT22*auxMHx(xp,yp,zp)+2.0d0*getDWCT2*auxDMHx(1,xp,yp,zp)+CgetW*auxDMHx(7,xp,yp,zp)+getDWCT33*auxMHx(xp,yp,zp)+2.0d0*getDWCT3*auxDMHx(2,xp,yp,zp)+CgetW*auxDMHx(8,xp,yp,zp));
         else if (type == 2) then      ! HY
            getNewDipoleHpoint = getDWCT2*auxMHy(xp,yp,zp)+x*getDWCT21*auxMHy(xp,yp,zp)+x*getDWCT2*auxDMHy(0,xp,yp,zp)+CgetW*auxDMHy(1,xp,yp,zp)+x*getDWCT1*auxDMHy(1,xp,yp,zp)+x*CgetW*auxDMHy(5,xp,yp,zp);
         else if (type == 3) then      ! HZ
            getNewDipoleHpoint = getDWCT3*auxMHz(xp,yp,zp)+x*getDWCT31*auxMHz(xp,yp,zp)+x*getDWCT3*auxDMHz(0,xp,yp,zp)+CgetW*auxDMHz(2,xp,yp,zp)+x*getDWCT1*auxDMHz(2,xp,yp,zp)+x*CgetW*auxDMHz(3,xp,yp,zp);
         else
            write(*,*) 'Alert: getHpoint get an unexpected type!', type   
         end if
      endif
    end function getNewDipoleHpoint


!Jpoint------------------------------------------
    real*8 function getNewDipoleJpoint(this, x, y, z, xp, yp, zp, t, type, problem_type)
      ! Calculate Currents at the time-point
      class(tNewDipole) :: this;
      real*8, intent(in) :: x, y, z;
      integer, intent(in) :: xp, yp, zp;
      real*8, intent(in) :: t
      integer, intent(in) :: type
      integer, intent(in) :: problem_type;
      real*8 :: res;
      real*8 :: CgetW, getNuCT, getDHiCT91, getDHiCT92, getDHiCT93, getDNu91CT, getDNu92CT, getDNu93CT, R, ct,  detw;
      real*8 :: getDWCT91, getDWCT92, getDWCT1, getDWCT2, getDWCT3, getDWCT31, getDWCT32, getDWCT21, getDWCT11, getDWCT22, getDWCT33, getDWCT39, getDWCT29;
      real*8 :: ct5, ct6, ct7, ct8, ct9, ct10, ct11, ct12, ct13, ct14, ct15, cosvct, sinvct;
      if (type/=1) then         
         R = sqrt(x**2 + y**2 + z**2)
         ct = t-R/cc
         ct5 = ct**5;
         ct6 = ct**6;
         ct7 = ct**7;
         ct8 = ct**8;
         ct9 = ct**9;
         ct10 = ct**10;
         ct11 = ct**11;
         ct12 = ct**12;
         ct13 = ct**13;
         ct14 = ct**14;
         ct15 = ct**15;
         getNuCT = 0.0d0;
         if ((q*ct>0.0d0).and.(q*ct<1.0d0)) then           
            getNuCT = cNu(0)*(ct15)+cNu(1)*(ct14)+cNu(2)*(ct13)+cNu(3)*(ct12)+cNu(4)*(ct11)+cNu(5)*(ct10)+cNu(6)*(ct9)+cNu(7)*(ct8);
         else if (q*ct<= 0.0d0 ) then     
            getNuCT = 0.0d0;
         else if (q*ct >= 1.0d0) then   
            getNuCT = 1.0d0;
         end if
       
         getDNu91CT = 0.0d0;
         getDNu92CT = 0.0d0;
         getDNu93CT = 0.0d0;
          if ((q*ct>0.0d0).and.(q*ct<1.0d0)) then           
             getDNu91CT = cDNu(0)*(ct14)+cDNu(1)*(ct13)+cDNu(2)*(ct12)+cDNu(3)*(ct11)+cDNu(4)*(ct10)+cDNu(5)*(ct9)+cDNu(6)*(ct8)+cDNu(7)*(ct7)
             getDNu92CT = cDNu(8)*(ct13)+cDNu(9)*(ct12)+cDNu(10)*(ct11)+cDNu(11)*(ct10)+cDNu(12)*(ct9)+cDNu(13)*(ct8)+cDNu(14)*(ct7)+cDNu(15)*(ct6)
             getDNu93CT = cDNu(16)*(ct12)+cDNu(17)*(ct11)+cDNu(18)*(ct10)+cDNu(19)*(ct9)+cDNu(20)*(ct8)+cDNu(21)*(ct7)+cDNu(22)*(ct6)+cDNu(23)*(ct5)
          else      
             getDNu91CT = 0.0d0
             getDNu92CT = 0.0d0
             getDNu93CT = 0.0d0
          end if

          cosvct = cos(v*ct);
          sinvct = sin(v*ct);
          getDHiCT91 = a*cos(v*ct)*v*getNuCT+a*sin(v*ct)*getDNu91CT;
          getDHiCT92 = getDNu92CT*a*sin(v*ct)+2*getDNu91CT*a*cos(v*ct)*v-getNuCT*a*sin(v*ct)*(v**2);
          getDHiCT93 = getDNu93CT*a*sin(v*ct)+3*getDNu92CT*a*cos(v*ct)*v-3*getDNu91CT*a*sin(v*ct)*(v**2)-getNuCT*a*cos(v*ct)*(v**3);

          if (ct >= 0.0d0) then
             CgetW = -1.0d0/R*(a*sin(v*ct)*getNuCT/(R**2)+1.0d0/(cc*R)*getDHiCT91);
          else
             CgetW = 0.0d0
          end if

          
          getDWCT91 = 0.0d0;
          getDWCT92 = 0.0d0;
          getDWCT1 = 0.0d0;
          getDWCT2 = 0.0d0;
          getDWCT3 = 0.0d0;
          getDWCT31 = 0.0d0;
          getDWCT32 = 0.0d0;
          getDWCT21 = 0.0d0;
          getDWCT11 = 0.0d0;
          getDWCT22 = 0.0d0;
          getDWCT33 = 0.0d0;
          getDWCT39 = 0.0d0;
          getDWCT29 = 0.0d0;
          
          if (ct > 0.0d0) then
             getDWCT91 = -1.0d0/R*(getDHiCT91/(R**2)+getDHiCT92/(R*cc));
             getDWCT92 = -1.0d0/R*(getDHiCT92/(R**2)+getDHiCT93/(R*cc));
             detw = 1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(R*cc))+1.0d0/R*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)));
             getDWCT1 = x*detw;
             getDWCT2 = y*detw;
             getDWCT3 = z*detw
             detw = -3.0d0/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))+2.0d0/(R**3)*(-2.0d0*getDHiCT91/((R**3)*cc)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)-getDHiCT92/((cc**2)*(R**2)))-1.0d0/R*(4.0d0*getDHiCT92/((R**4)*(cc**2))+8.0d0*getDHiCT91/((R**5)*cc)+8.0d0*a*sin(v*ct)*getNuCT/(R**6)+getDHiCT93/((cc**3)*(R**3)))
             getDWCT31 = z*x*detw;
             getDWCT32 = z*y*detw
             getDWCT21 = y*x*detw
             getDWCT11 =  -3.0d0*(x**2)/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-2.0d0*(x**2)/(R**3)*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)))+1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-1.0d0/R*(4.0d0*getDHiCT92*(x**2)/((R**4)*(cc**2))+8.0d0*getDHiCT91*(x**2)/((R**5)*cc)-2.0d0*getDHiCT91/((R**3)*cc)+8.0d0*a*sin(v*ct)*getNuCT*(x**2)/(R**6)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT93*(x**2)/((cc**3)*(R**3))-getDHiCT92/((cc**2)*(R**2))) 
             getDWCT22 =  -3.0d0*(y**2)/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-2.0d0*(y**2)/(R**3)*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)))+1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-1.0d0/R*(4.0d0*getDHiCT92*(y**2)/((R**4)*(cc**2))+8.0d0*getDHiCT91*(y**2)/((R**5)*cc)-2.0d0*getDHiCT91/((R**3)*cc)+8.0d0*a*sin(v*ct)*getNuCT*(y**2)/(R**6)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT93*(y**2)/((cc**3)*(R**3))-getDHiCT92/((cc**2)*(R**2))) 
             getDWCT33 =  -3.0d0*(z**2)/(R**5)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-2.0d0*(z**2)/(R**3)*(2.0d0*getDHiCT91/((R**3)*cc)+2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT92/((cc**2)*(R**2)))+1.0d0/(R**3)*(a*sin(v*ct)*getNuCT/(R**2)+getDHiCT91/(cc*R))-1.0d0/R*(4.0d0*getDHiCT92*(z**2)/((R**4)*(cc**2))+8.0d0*getDHiCT91*(z**2)/((R**5)*cc)-2.0d0*getDHiCT91/((R**3)*cc)+8.0d0*a*sin(v*ct)*getNuCT*(z**2)/(R**6)-2.0d0*a*sin(v*ct)*getNuCT/(R**4)+getDHiCT93*(z**2)/((cc**3)*(R**3))-getDHiCT92/((cc**2)*(R**2))) 
             detw = 1.0d0/(R**3)*(getDHiCT91/(R**2)+getDHiCT92/(R*cc))+1.0d0/R*(2.0d0*getDHiCT92/((R**3)*cc)+2.0d0*getDHiCT91/(R**4)+getDHiCT93/((cc**2)*(R**2)));
             getDWCT39 = z*detw 
             getDWCT29 = y*detw 
        end if   
     endif

     if (problem_type == 0) then
        getNewDipoleJpoint = 0.0d0
        if (type==1) then             ! JX
           getNewDipoleJpoint = 0.0d0
        else if (type == 2) then      ! JY         
           getNewDipoleJpoint = -x/(cc**2)*DMEJy(2,xp,yp,zp)*getDWCT92+x*DMEJy(11,xp,yp,zp)*CgetW+x*DMEJy(7,xp,yp,zp)*getDWCT3+2.0d0*x*DMEJy(4,xp,yp,zp)*getDWCT2+2.0d0*x*DMEJy(1,xp,yp,zp)*getDWCT32+x*DMEJy(2,xp,yp,zp)*getDWCT22+x*DMEJy(10,xp,yp,zp)*CgetW+3.0d0*x*DMEJy(8,xp,yp,zp)*getDWCT3+3.0d0*x*DMEJy(2,xp,yp,zp)*getDWCT33+2.0d0*DMEJy(3,xp,yp,zp)*CgetW+2.0d0*DMEJy(2,xp,yp,zp)*getDWCT1+x*DMEJy(12,xp,yp,zp)*CgetW+2.0d0*x*DMEJy(3,xp,yp,zp)*getDWCT1+x*DMEJy(2,xp,yp,zp)*getDWCT11+2.0d0*DMEJy(0,xp,yp,zp)*getDWCT3+x*DMEJy(6,xp,yp,zp)*getDWCT3+2.0d0*x*DMEJy(0,xp,yp,zp)*getDWCT31;
        else if (type == 3) then      ! JZ
           getNewDipoleJpoint = x/(cc**2)*DMEJz(1,xp,yp,zp)*getDWCT92-2.0d0*DMEJz(5,xp,yp,zp)*CgetW-2.0d0*DMEJz(1,xp,yp,zp)*getDWCT1-x*DMEJz(13,xp,yp,zp)*CgetW-2.0d0*x*DMEJz(5,xp,yp,zp)*getDWCT1-x*DMEJz(1,xp,yp,zp)*getDWCT11-2.0d0*DMEJz(0,xp,yp,zp)*getDWCT2-x*DMEJz(6,xp,yp,zp)*getDWCT2-2.0d0*x*DMEJz(0,xp,yp,zp)*getDWCT21-x*DMEJz(9,xp,yp,zp)*CgetW-3.0d0*x*DMEJz(7,xp,yp,zp)*getDWCT2-3.0d0*x*DMEJz(1,xp,yp,zp)*getDWCT22-x*DMEJz(14,xp,yp,zp)*CgetW-x*DMEJz(8,xp,yp,zp)*getDWCT2-2.0d0*x*DMEJz(4,xp,yp,zp)*getDWCT3-2.0d0*x*DMEJz(2,xp,yp,zp)*getDWCT32-x*DMEJz(1,xp,yp,zp)*getDWCT33;
        else
           write(*,*) 'Alert: getJpoint get an unexpected type!', type  
        end if
     else if (problem_type == 1) then
        getNewDipoleJpoint = 0.0d0
        if (type==1) then             ! JX
           getNewDipoleJpoint = 0.0d0
        else if (type == 2) then      ! JY         
           getNewDipoleJpoint = -x/(cc**2)*auxDMEJy(2,xp,yp,zp)*getDWCT92+x*auxDMEJy(11,xp,yp,zp)*CgetW+x*auxDMEJy(7,xp,yp,zp)*getDWCT3+2.0d0*x*auxDMEJy(4,xp,yp,zp)*getDWCT2+2.0d0*x*auxDMEJy(1,xp,yp,zp)*getDWCT32+x*auxDMEJy(2,xp,yp,zp)*getDWCT22+x*auxDMEJy(10,xp,yp,zp)*CgetW+3.0d0*x*auxDMEJy(8,xp,yp,zp)*getDWCT3+3.0d0*x*auxDMEJy(2,xp,yp,zp)*getDWCT33+2.0d0*auxDMEJy(3,xp,yp,zp)*CgetW+2.0d0*auxDMEJy(2,xp,yp,zp)*getDWCT1+x*auxDMEJy(12,xp,yp,zp)*CgetW+2.0d0*x*auxDMEJy(3,xp,yp,zp)*getDWCT1+x*auxDMEJy(2,xp,yp,zp)*getDWCT11+2.0d0*auxDMEJy(0,xp,yp,zp)*getDWCT3+x*auxDMEJy(6,xp,yp,zp)*getDWCT3+2.0d0*x*auxDMEJy(0,xp,yp,zp)*getDWCT31;
        else if (type == 3) then      ! JZ
           getNewDipoleJpoint = x/(cc**2)*auxDMEJz(1,xp,yp,zp)*getDWCT92-2.0d0*auxDMEJz(5,xp,yp,zp)*CgetW-2.0d0*auxDMEJz(1,xp,yp,zp)*getDWCT1-x*auxDMEJz(13,xp,yp,zp)*CgetW-2.0d0*x*auxDMEJz(5,xp,yp,zp)*getDWCT1-x*auxDMEJz(1,xp,yp,zp)*getDWCT11-2.0d0*auxDMEJz(0,xp,yp,zp)*getDWCT2-x*auxDMEJz(6,xp,yp,zp)*getDWCT2-2.0d0*x*auxDMEJz(0,xp,yp,zp)*getDWCT21-x*auxDMEJz(9,xp,yp,zp)*CgetW-3.0d0*x*auxDMEJz(7,xp,yp,zp)*getDWCT2-3.0d0*x*auxDMEJz(1,xp,yp,zp)*getDWCT22-x*auxDMEJz(14,xp,yp,zp)*CgetW-x*auxDMEJz(8,xp,yp,zp)*getDWCT2-2.0d0*x*auxDMEJz(4,xp,yp,zp)*getDWCT3-2.0d0*x*auxDMEJz(2,xp,yp,zp)*getDWCT32-x*auxDMEJz(1,xp,yp,zp)*getDWCT33;
        else
           write(*,*) 'Alert: getJpoint get an unexpected type!', type  
        end if
     endif
     
      getNewDipoleJpoint = -getNewDipoleJpoint*cc/4.0d0/Pi;
    end function getNewDipoleJpoint
    
end module newdipole
