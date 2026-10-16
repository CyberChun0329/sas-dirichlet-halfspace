program test_replay_measure
   use, intrinsic :: iso_fortran_env, only: int64
   use replay_types, only: rk,ok,numerical_failure
   use replay_measure
   implicit none
   type(measurement_t) :: query,total,failed
   real(rk) :: values(3),reference(3),sorted(9),tmp
   integer(int64) :: actual_calls=0,prep_calls=0
   integer :: status,i,j,mode
   real(rk),volatile :: work_sink=0
   mode=1
   call measure_batch(batch,values,query,status)
   call require(status==ok,'query measurement succeeds')
   reference=values
   call check_measurement(query)
   call require(actual_calls==2+sum(query%calls),'warm/calibration plus actual query calls')
   actual_calls=0;mode=2
   call measure_batch(batch,values,total,status)
   call require(status==ok,'direct total measurement succeeds')
   call require(all(values==reference),'total preserves actual query values')
   call check_measurement(total)
   call require(actual_calls==2+sum(total%calls).and.prep_calls==actual_calls,'prep reruns for every total call')
   ! The third and later responses change, so calibration succeeds and the first window fails.
   actual_calls=0;mode=3
   call measure_batch(batch,values,failed,status)
   call require(status==numerical_failure.and.failed%limit_reason=='WINDOW_RESPONSE_CHANGED', &
      'changed timed-window response rejected')
   call require(failed%requested_windows==9.and.failed%window_count==1.and.failed%completed_windows==0, &
      'changed window attempted/completed accounting')
   call require(failed%calls(1)==actual_calls-2.and.failed%window_seconds(1)>=.004_rk, &
      'changed window actual calls/time retained')
   call require(all(failed%calls(2:)==0).and.failed%median_cpu==0.and.failed%low_cpu==0.and.failed%high_cpu==0, &
      'unattempted windows excluded and no failed aggregate statistics')
   actual_calls=0;mode=4
   call measure_batch(batch,values,failed,status)
   call require(status==numerical_failure.and.failed%limit_reason=='CALLBACK_FAILURE','callback failure rejected')
   call require(failed%window_count==1.and.failed%completed_windows==0.and.failed%calls(1)==1, &
      'failed callback counted as actual attempt')
   call require(failed%window_seconds(1)>=0.and.failed%wall_seconds(1)>=0,'failed window times retained')
   call require(failed%warm_recorded.and.failed%calibration_recorded,'failure retains warm/calibration evidence')
   actual_calls=0;mode=5
   call measure_batch(batch,values,failed,status)
   call require(status==numerical_failure.and.failed%limit_reason=='WARM_RESPONSE_FAILURE','warm failure rejected')
   call require(failed%warm_recorded.and..not.failed%calibration_recorded.and.failed%window_count==0, &
      'warm-only failure does not fabricate calibration/windows')
   call require(actual_calls==1.and.failed%completed_windows==0.and.sum(failed%calls)==0,'warm failure accounting')
   actual_calls=0;mode=6
   call measure_batch(batch,values,failed,status)
   call require(status==numerical_failure.and.failed%limit_reason=='CALIBRATION_RESPONSE_FAILURE', &
      'calibration failure rejected')
   call require(failed%warm_recorded.and.failed%calibration_recorded.and.failed%window_count==0, &
      'calibration failure records diagnostics without windows')
   call require(actual_calls==2.and.failed%calibration_cpu>=0.and.failed%calibration_wall>=0, &
      'calibration actual CPU/wall retained')
   print '(a)','PASS: actual nine-window query/total callbacks, CPU/call medians, changed responses and failures'
contains
   subroutine batch(v,st)
      real(rk),intent(out) :: v(:)
      integer,intent(out) :: st
      real(rk) :: sumwork,x
      integer :: n,p
      actual_calls=actual_calls+1;st=ok
      if((mode==4.and.actual_calls>=3).or.(mode==5.and.actual_calls==1).or.(mode==6.and.actual_calls==2)) then
         st=numerical_failure;v=0;return
      end if
      if(mode==2) then
         prep_calls=prep_calls+1;sumwork=0
         do n=1,2000
            x=real(n,rk)/2000;sumwork=sumwork+sqrt(1+x*x)
         end do
         work_sink=sumwork
      end if
      ! Recompute all entries on every call; no cached answer vector.
      do p=1,size(v)
         sumwork=0
         do n=1,800
            x=real(n+p,rk)/800;sumwork=sumwork+exp(-x*x)
         end do
         v(p)=sumwork/800
      end do
      if(mode==3.and.actual_calls>=3) v(1)=v(1)+1
   end subroutine
   subroutine check_measurement(m)
      type(measurement_t),intent(in) :: m
      call require(m%window_count==9.and.m%minimum_window==.004_rk,'nine four-ms windows')
      call require(m%requested_windows==9.and.m%completed_windows==9,'requested/attempted/completed windows')
      call require(m%warm_recorded.and.m%calibration_recorded,'warm/calibration actual calls recorded')
      call require(min(m%warm_cpu,m%warm_wall,m%calibration_cpu,m%calibration_wall)>=0, &
         'successful warm/calibration diagnostics')
      call require(all(m%window_seconds>=.004_rk).and.all(m%calls>0),'actual window duration and calls')
      call require(all(m%wall_seconds>=0),'diagnostic wall durations')
      call require(all(m%cpu_windows==m%window_seconds/real(m%calls,rk)),'CPU seconds/call readback')
      sorted=m%cpu_windows
      do i=2,9
         tmp=sorted(i);j=i-1
         do while(j>=1)
            if(sorted(j)<=tmp) exit
            sorted(j+1)=sorted(j);j=j-1
         end do
         sorted(j+1)=tmp
      end do
      call require(m%median_cpu==sorted(5).and.m%low_cpu==sorted(1).and.m%high_cpu==sorted(9),'median/min/max')
   end subroutine
   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(*),intent(in) :: label
      if(.not.condition) then
         print '(a)','FAIL: '//label
         error stop 1
      end if
   end subroutine
end program
