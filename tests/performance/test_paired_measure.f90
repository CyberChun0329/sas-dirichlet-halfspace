! Synthetic callback accounting only: no scientific input bank or resource tuning.
program test_paired_measure
   use, intrinsic :: iso_fortran_env, only: int64
   use replay_types, only: rk,ok,invalid
   use unified_policy, only: policy_t
   use replay_measure, only: measurement_t
   use paired_measure, only: measure_paired,pair_trace_t
   implicit none
   type(policy_t) :: p
   type(measurement_t) :: m(2)
   type(pair_trace_t) :: trace(2)
   real(rk) :: a(1),b(1)
   integer :: status(2),left_calls,right_calls,scope,scope_calls(3,2),j
   left_calls=0;right_calls=0;scope_calls=0;p%formal_allowed=.true.
   do scope=1,3
      call measure_paired(left,right,a,b,m,status,p,.true.,0,trace)
      call require(all(status==ok),'synthetic formal-shaped pair completed')
      call require(all(m%completed_windows==9),'nine actual windows per side')
      do j=1,2
         call require(all(m(j)%calls>0_int64),'actual calls retained')
         call require(all(m(j)%window_seconds>=p%formal_min_window),'each CPU window resolved')
         call require(trace(j)%rate>0.and.all(trace(j)%window_end>=trace(j)%window_start),'clock trace retained')
      end do
      call require(trace(1)%window_start(1)<trace(2)%window_start(1),'first window left then right')
      call require(trace(2)%window_start(2)<trace(1)%window_start(2),'second window right then left')
   end do
   call require(all(scope_calls>0),'independent PREP QUERY TOTAL callbacks exercised')
   scope=1;left_calls=0;right_calls=0
   call measure_paired(failing,right,a,b,m,status,p,.false.,0,trace)
   call require(status(1)/=ok.and.status(2)==ok,'one side failure preserves peer evidence')
   call require(m(1)%completed_windows==0.and.m(1)%calls(1)==1,'callback failure attempt retained')
   call require(m(2)%completed_windows==1,'peer coarse scope completed')
   p%call_limit=1_int64
   call measure_paired(tiny,tiny,a,b,m,status,p,.false.,0,trace)
   call require(all(status/=ok),'call budget censors unresolved windows')
   call require(all(m%completed_windows==0),'censored window not marked complete')
   call require(all(m(1)%calls(:1)==1_int64).and.all(m(2)%calls(:1)==1_int64),'censored actual calls retained')
   ! Warm/calibration are tiny, so their estimate schedules a multi-call chunk.
   ! Only its first callback is slow: a chunk-boundary check would call again.
   p=policy_t();p%callback_limit_seconds=.01_rk
   left_calls=0;right_calls=0
   call measure_paired(slow_after_calibration,right,a,b,m,status,p,.false.,0,trace)
   call require(status(1)/=ok.and.status(2)==ok,'returned callback time limit preserves peer evidence')
   call require(m(1)%warm_recorded.and.m(1)%calibration_recorded,'fast warm and calibration retained')
   call require(m(1)%warm_cpu<p%callback_limit_seconds.and. &
      m(1)%calibration_cpu<p%coarse_min_window/2,'fast calibration schedules more than one chunk call')
   call require(left_calls==3.and.m(1)%calls(1)==1_int64,'first over-limit callback stops chunk remainder')
   call require(m(1)%completed_windows==0.and.m(1)%window_seconds(1)>p%callback_limit_seconds, &
      'over-limit CPU window is retained and censored')
   call require(m(1)%wall_seconds(1)>0.and.trace(1)%window_end(1)>trace(1)%window_start(1), &
      'over-limit wall time and clock endpoints retained')
   call require(m(1)%limit_reason=='CALLBACK_OR_WINDOW_TIME_LIMIT','explicit returned-callback limit reason')
   call require(m(2)%completed_windows==1.and.m(2)%calls(1)>0_int64,'peer successful window retained')
   print '(a)','PASS paired measurement synthetic success, rotating order, independent scopes, failure and resource-limit censor'
contains
   subroutine require(test,label)
      logical,intent(in) :: test
      character(*),intent(in) :: label
      if(.not.test) then
         print *,trim(label);error stop 'paired test failed'
      end if
   end subroutine
   subroutine work(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      real(rk) :: total
      integer :: k
      total=0
      do k=1,1024
         total=total+sin(real(k,rk))
      end do
      values=total+scope;status=ok
   end subroutine
   subroutine left(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      left_calls=left_calls+1;scope_calls(scope,1)=scope_calls(scope,1)+1
      call work(values,status)
   end subroutine
   subroutine right(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      right_calls=right_calls+1;scope_calls(scope,2)=scope_calls(scope,2)+1
      call work(values,status)
   end subroutine
   subroutine failing(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      left_calls=left_calls+1;values=1;status=ok
      if(left_calls>=3) status=invalid
   end subroutine
   subroutine tiny(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      values=1;status=ok
   end subroutine
   subroutine slow_after_calibration(values,status)
      real(rk),intent(out) :: values(:)
      integer,intent(out) :: status
      real(rk) :: started,now
      left_calls=left_calls+1;values=1;status=ok
      if(left_calls/=3) return
      call cpu_time(started)
      do
         call cpu_time(now)
         if(now-started>=2*p%callback_limit_seconds) exit
      end do
   end subroutine
end program
