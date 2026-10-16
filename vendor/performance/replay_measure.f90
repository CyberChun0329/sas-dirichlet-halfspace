module replay_measure
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use replay_types, only: rk,ok,invalid,numerical_failure
   implicit none
   private
   public :: measurement_t,batch_iface,measure_batch
   type measurement_t
      real(rk) :: median_cpu=0,low_cpu=0,high_cpu=0,median_wall=0
      real(rk) :: cpu_windows(9)=0,window_seconds(9)=0,wall_seconds(9)=0
      integer(int64) :: calls(9)=0
      integer :: window_count=0,requested_windows=0,completed_windows=0
      real(rk) :: minimum_window=0,calibration_cpu=0,calibration_wall=0,warm_cpu=0,warm_wall=0
      logical :: warm_recorded=.false.,calibration_recorded=.false.
      character(48) :: limit_reason=''
   end type
   abstract interface
      subroutine batch_iface(values,status)
         import rk
         real(rk),intent(out) :: values(:)
         integer,intent(out) :: status
      end subroutine
   end interface
contains
   ! Each call recomputes the batch. Preparation/total semantics belong to fn.
   ! max_windows is the requested count (1..9); defaults are 3 windows and 0.001 s.
   subroutine measure_batch(fn,values,measurement,status,min_window,max_windows)
      procedure(batch_iface) :: fn
      real(rk),intent(out) :: values(:)
      type(measurement_t),intent(out) :: measurement
      integer,intent(out) :: status
      real(rk),optional,intent(in) :: min_window
      integer,optional,intent(in) :: max_windows
      integer(int64),parameter :: call_limit=10000000_int64
      real(rk),parameter :: time_limit=30.0_rk
      real(rk) :: baseline(size(values)),c0,c1,elapsed,estimate,wall,target
      real(rk) :: cw(9),ww(9)
      integer(int64) :: t0,t1,rate,clockmax,chunk,j,count
      integer :: nw,w,st
      measurement=measurement_t();status=invalid;nw=9;target=.004_rk
      if(present(min_window)) target=min_window
      if(present(max_windows)) nw=max_windows
      measurement%requested_windows=nw;measurement%minimum_window=target
      if(size(values)==0 .or. nw<1 .or. nw>9) return
      if(.not.ieee_is_finite(target) .or. target<=0 .or. target>time_limit) return
      call system_clock(t0,rate,clockmax)
      if(rate<=0) then
         measurement%limit_reason='SYSTEM_CLOCK_UNAVAILABLE';return
      end if
      ! Warm code once, then calibrate a separate actual call.
      call system_clock(t0);call cpu_time(c0);call fn(values,st);call cpu_time(c1);call system_clock(t1)
      measurement%warm_recorded=.true.;measurement%warm_cpu=c1-c0
      measurement%warm_wall=clock_elapsed(t0,t1,rate,clockmax)
      if(st/=ok .or. .not.all(ieee_is_finite(values))) then
         status=numerical_failure;measurement%limit_reason='WARM_RESPONSE_FAILURE';return
      end if
      if(c0<0.or.c1<c0) then
         measurement%limit_reason='CPU_TIME_UNAVAILABLE';return
      end if
      if(measurement%warm_cpu>time_limit.or.measurement%warm_wall>time_limit) then
         status=numerical_failure;measurement%limit_reason='WARM_TIME_LIMIT';return
      end if
      baseline=values
      call system_clock(t0);call cpu_time(c0);call fn(values,st);call cpu_time(c1);call system_clock(t1)
      measurement%calibration_recorded=.true.;measurement%calibration_cpu=c1-c0
      measurement%calibration_wall=clock_elapsed(t0,t1,rate,clockmax)
      if(st/=ok .or. .not.all(ieee_is_finite(values)) .or. any(values/=baseline)) then
         status=numerical_failure;measurement%limit_reason='CALIBRATION_RESPONSE_FAILURE';return
      end if
      if(c0<0 .or. c1<c0) then
         measurement%limit_reason='CPU_TIME_UNAVAILABLE';return
      end if
      if(measurement%calibration_cpu>time_limit.or.measurement%calibration_wall>time_limit) then
         status=numerical_failure;measurement%limit_reason='CALIBRATION_TIME_LIMIT';return
      end if
      estimate=max(c1-c0,1e-7_rk)
      do w=1,nw
         measurement%window_count=w
         count=0
         chunk=max(1_int64,int(min(real(call_limit,rk),target/estimate),int64))
         call system_clock(t0);call cpu_time(c0)
         do
            chunk=min(chunk,call_limit-count)
            if(chunk<=0) then
               call record_window(measurement,w,count,c0,t0,rate,clockmax,c1,elapsed,wall)
               status=numerical_failure;measurement%limit_reason='BATCH_LIMIT';return
            end if
            do j=1,chunk
               call fn(values,st)
               count=count+1
               if(st/=ok) then
                  call record_window(measurement,w,count,c0,t0,rate,clockmax,c1,elapsed,wall)
                  status=numerical_failure;measurement%limit_reason='CALLBACK_FAILURE';return
               end if
            end do
            call record_window(measurement,w,count,c0,t0,rate,clockmax,c1,elapsed,wall)
            if(c1<c0 .or. elapsed<0 .or. .not.ieee_is_finite(elapsed)) then
               status=numerical_failure;measurement%limit_reason='INVALID_CPU_WINDOW';return
            end if
            if(elapsed>time_limit .or. wall>time_limit) then
               status=numerical_failure;measurement%limit_reason='TIME_LIMIT';return
            end if
            if(elapsed>=target) exit
            if(elapsed>0) then
               estimate=elapsed/real(count,rk)
               chunk=max(1_int64,int(min(real(call_limit-count,rk),(target-elapsed)/estimate+1),int64))
            else
               chunk=min(call_limit-count,2*chunk)
            end if
         end do
         ! Values are actual final responses, never a cached answer used by fn.
         if(.not.all(ieee_is_finite(values)) .or. any(values/=baseline)) then
            status=numerical_failure;measurement%limit_reason='WINDOW_RESPONSE_CHANGED';return
         end if
         measurement%completed_windows=w
         cw(w)=measurement%cpu_windows(w);ww(w)=wall/real(count,rk)
         estimate=cw(w)
      end do
      call sort(cw(:nw));call sort(ww(:nw))
      measurement%median_cpu=median(cw(:nw));measurement%median_wall=median(ww(:nw))
      measurement%low_cpu=cw(1);measurement%high_cpu=cw(nw);status=ok
   end subroutine

   ! Snapshot actual attempted calls and elapsed times before any censored-window return.
   subroutine record_window(m,w,count,c0,t0,rate,clockmax,c1,elapsed,wall)
      type(measurement_t),intent(inout) :: m
      integer,intent(in) :: w
      integer(int64),intent(in) :: count,t0,rate,clockmax
      real(rk),intent(in) :: c0
      real(rk),intent(out) :: c1,elapsed,wall
      integer(int64) :: t1
      call cpu_time(c1);call system_clock(t1)
      elapsed=c1-c0;wall=clock_elapsed(t0,t1,rate,clockmax)
      m%calls(w)=count;m%window_seconds(w)=elapsed;m%wall_seconds(w)=wall
      if(count>0) m%cpu_windows(w)=elapsed/real(count,rk)
   end subroutine

   pure real(rk) function clock_elapsed(a,b,rate,clockmax) result(v)
      integer(int64),intent(in) :: a,b,rate,clockmax
      if(b>=a) then
         v=real(b-a,rk)/real(rate,rk)
      else
         v=(real(clockmax-a,rk)+real(b,rk)+1)/real(rate,rk)
      end if
   end function

   pure subroutine sort(a)
      real(rk),intent(inout) :: a(:)
      real(rk) :: v
      integer :: i,j
      do i=2,size(a)
         v=a(i);j=i-1
         do while(j>=1)
            if(a(j)<=v) exit
            a(j+1)=a(j);j=j-1
         end do
         a(j+1)=v
      end do
   end subroutine
   pure real(rk) function median(a) result(v)
      real(rk),intent(in) :: a(:)
      integer :: n
      n=size(a);v=a((n+1)/2)
      if(mod(n,2)==0) v=(a(n/2)+a(n/2+1))/2
   end function
end module
