module useful_screening
   use replay_types, only: rk,ok,invalid
   use replay_measure, only: measurement_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: axis_cost_bounds
contains
   pure subroutine axis_cost_bounds(m,axis,coarse,mid,lower,upper,status)
      type(measurement_t),intent(in) :: m
      integer,intent(in) :: axis
      logical,intent(in) :: coarse
      real(rk),intent(out) :: mid,lower,upper
      integer,intent(out) :: status
      real(rk) :: windows(9),v,warm,calibration
      integer :: n,i,j
      status=invalid;mid=0;lower=0;upper=0
      n=m%completed_windows
      if(n<1.or.n>9.or.n/=m%requested_windows) return
      select case(axis)
      case(1)
         windows(:n)=m%cpu_windows(:n);warm=m%warm_cpu;calibration=m%calibration_cpu
      case(2)
         if(any(m%calls(:n)<=0)) return
         windows(:n)=m%wall_seconds(:n)/real(m%calls(:n),rk)
         warm=m%warm_wall;calibration=m%calibration_wall
      case default
         return
      end select
      if(.not.all(ieee_is_finite(windows(:n)))) return
      if(minval(windows(:n))<=0) return
      do i=2,n
         v=windows(i);j=i-1
         do while(j>=1)
            if(windows(j)<=v) exit
            windows(j+1)=windows(j);j=j-1
         end do
         windows(j+1)=v
      end do
      mid=windows((n+1)/2)
      if(mod(n,2)==0) mid=(windows(n/2)+windows(n/2+1))/2
      lower=windows(1);upper=windows(n)
      if(coarse) then
         if(.not.m%warm_recorded.or..not.m%calibration_recorded) return
         if(.not.all(ieee_is_finite([warm,calibration]))) return
         lower=min(lower,warm,calibration);upper=max(upper,warm,calibration)
         if(lower<=0) return
      end if
      status=ok
   end subroutine
end module
