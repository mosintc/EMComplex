!-----------------------------------------------------------------------
!  Copyright 2017 Mikhail Osintcev
!  This file is part of the EMtool developed at NCSU
!-----------------------------------------------------------------------
! This module contains implementation of divergence free source (old and slow implementation)

module dipole
! Computate the analytic solution of the dipole
  use commonvars
   implicit none;
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

!Epoint------------------------------------------
    real*8 function getDipoleEpoint(x,y,z,t,type)
    ! Calculate Electric pole in the time-point  
      real*8, intent(in) :: x,y,z,t
      integer, intent(in) :: type
      getDipoleEpoint = 0.0d0;
      if (type==1) then            ! EX
         getDipoleEpoint = 0.0d0
      else if (type == 2) then      ! EY
         getDipoleEpoint = -x/cc*(getM(x,y,z)*getDW(x,y,z,t,39)+getDM(x,y,z,3)*getDW(x,y,z,t,91))    
      else if (type == 3) then      ! EZ
         getDipoleEpoint = x/cc*(getDM(x,y,z,2)*getDW(x,y,z,t,91)+getM(x,y,z)*getDW(x,y,z,t,29))
      else
         write(*,*) 'Alert: getEpoint get an unexpected type:', type
      end if
    end function getDipoleEpoint

!Hpoint------------------------------------------
    real*8 function getDipoleHpoint(x,y,z,t,type)
    ! Calculate Magnetic pole in the time-point  
      real*8, intent(in) :: x,y,z,t
      integer, intent(in) :: type
      getDipoleHpoint = 0.0d0
      if (type==1) then             ! HX
         getDipoleHpoint = -x*(getDW(x,y,z,t,22)*getM(x,y,z)+2*getDW(x,y,z,t,2)*getDM(x,y,z,2)+getW(x,y,z,t)*getDM(x,y,z,22)+getDW(x,y,z,t,33)*getM(x,y,z)+2*getDW(x,y,z,t,3)*getDM(x,y,z,3)+getW(x,y,z,t)*getDM(x,y,z,33));
      else if (type == 2) then      ! HY
         getDipoleHpoint = getDW(x,y,z,t,2)*getM(x,y,z)+x*getDW(x,y,z,t,21)*getM(x,y,z)+x*getDW(x,y,z,t,2)*getDM(x,y,z,1)+getW(x,y,z,t)*getDM(x,y,z,2)+x*getDW(x,y,z,t,1)*getDM(x,y,z,2)+x*getW(x,y,z,t)*getDM(x,y,z,21);
      else if (type == 3) then      ! HZ
         getDipoleHpoint = getDW(x,y,z,t,3)*getM(x,y,z)+x*getDW(x,y,z,t,31)*getM(x,y,z)+x*getDW(x,y,z,t,3)*getDM(x,y,z,1)+getW(x,y,z,t)*getDM(x,y,z,3)+x*getDW(x,y,z,t,1)*getDM(x,y,z,3)+x*getW(x,y,z,t)*getDM(x,y,z,31);
      else
         write(*,*) 'Alert: getHpoint get an unexpected type!', type   
      end if
    end function getDipoleHpoint

!Jpoint------------------------------------------
    real*8 function getDipoleJpoint(x,y,z,t,type)
    ! Calculate Currents in the time-point  
      real*8, intent(in) :: x,y,z,t
      integer, intent(in) :: type
      real*8 :: res
      getDipoleJpoint = 0.0d0
      if (type==1) then             ! JX
         getDipoleJpoint = 0.0d0
      else if (type == 2) then      ! JY         
         getDipoleJpoint = -x/(cc**2)*getDM(x,y,z,3)*getDW(x,y,z,t,92)+x*getDM(x,y,z,322)*getW(x,y,z,t)+x*getDM(x,y,z,22)*getDW(x,y,z,t,3)+2.0d0*x*getDM(x,y,z,32)*getDW(x,y,z,t,2)+2.0d0*x*getDM(x,y,z,2)*getDW(x,y,z,t,32)+x*getDM(x,y,z,3)*getDW(x,y,z,t,22)+x*getDM(x,y,z,333)*getW(x,y,z,t)+3.0d0*x*getDM(x,y,z,33)*getDW(x,y,z,t,3)+3.0d0*x*getDM(x,y,z,3)*getDW(x,y,z,t,33)+2.0d0*getDM(x,y,z,31)*getW(x,y,z,t)+2.0d0*getDM(x,y,z,3)*getDW(x,y,z,t,1)+x*getDM(x,y,z,311)*getW(x,y,z,t)+2.0d0*x*getDM(x,y,z,31)*getDW(x,y,z,t,1)+x*getDM(x,y,z,3)*getDW(x,y,z,t,11)+2.0d0*getDM(x,y,z,1)*getDW(x,y,z,t,3)+x*getDM(x,y,z,11)*getDW(x,y,z,t,3)+2.0d0*x*getDM(x,y,z,1)*getDW(x,y,z,t,31);
      else if (type == 3) then      ! JZ
         getDipoleJpoint = x/(cc**2)*getDM(x,y,z,2)*getDW(x,y,z,t,92)-2.0d0*getDM(x,y,z,21)*getW(x,y,z,t)-2.0d0*getDM(x,y,z,2)*getDW(x,y,z,t,1)-x*getDM(x,y,z,211)*getW(x,y,z,t)-2.0d0*x*getDM(x,y,z,21)*getDW(x,y,z,t,1)-x*getDM(x,y,z,2)*getDW(x,y,z,t,11)-2.0d0*getDM(x,y,z,1)*getDW(x,y,z,t,2)-x*getDM(x,y,z,11)*getDW(x,y,z,t,2)-2.0d0*x*getDM(x,y,z,1)*getDW(x,y,z,t,21)-x*getDM(x,y,z,222)*getW(x,y,z,t)-3.0d0*x*getDM(x,y,z,22)*getDW(x,y,z,t,2)-3.0d0*x*getDM(x,y,z,2)*getDW(x,y,z,t,22)-x*getDM(x,y,z,332)*getW(x,y,z,t)-x*getDM(x,y,z,33)*getDW(x,y,z,t,2)-2.0d0*x*getDM(x,y,z,32)*getDW(x,y,z,t,3)-2.0d0*x*getDM(x,y,z,3)*getDW(x,y,z,t,32)-x*getDM(x,y,z,2)*getDW(x,y,z,t,33)
      else
        write(*,*) 'Alert: getJpoint get an unexpected type!', type  
      end if
      getDipoleJpoint = -getDipoleJpoint*cc/4.0d0/Pi;
    end function getDipoleJpoint
    
end module dipole
