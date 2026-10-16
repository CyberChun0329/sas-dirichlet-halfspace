module paired_measure
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use replay_types, only: rk,ok,invalid,numerical_failure
   use replay_measure, only: measurement_t,batch_iface
   use unified_policy, only: policy_t,validate_policy
   implicit none
   private
   type,public :: pair_trace_t
      integer(int64) :: warm_start=0,warm_end=0,calibration_start=0,calibration_end=0
      integer(int64) :: window_start(9)=0,window_end(9)=0,rate=0,clockmax=0
   end type
   public :: measure_paired,pair_order
contains
   pure function pair_order(window,seed) result(order)
      integer,intent(in) :: window,seed
      integer :: order(2)
      order=[1,2]
      if(modulo(window+seed,2)==0) order=[2,1]
   end function

   ! One warm call and one calibration per side, followed by alternating actual
   ! CPU windows. Both sides remain in one task and share the same protocol.
   ! PREP, QUERY and TOTAL call this separately with genuinely different callbacks.
   subroutine measure_paired(left,right,left_values,right_values,measurements,statuses,p,formal,seed,trace)
      procedure(batch_iface) :: left,right
      real(rk),intent(out) :: left_values(:),right_values(:)
      type(measurement_t),intent(out) :: measurements(2)
      integer,intent(out) :: statuses(2)
      type(policy_t),intent(in) :: p
      logical,intent(in) :: formal
      integer,intent(in) :: seed
      type(pair_trace_t),optional,intent(out) :: trace(2)
      type(pair_trace_t) :: local_trace(2)
      real(rk) :: lb(size(left_values)),rb(size(right_values)),target
      real(rk) :: cpu,wall,side_wall(2),estimates(2)
      integer :: nw,side,phase,k,w,st,order(2)
      measurements=measurement_t();statuses=invalid;side_wall=0
      local_trace=pair_trace_t()
      if(present(trace)) trace=local_trace
      call validate_policy(p,st)
      if(st/=ok.or.min(size(left_values),size(right_values))<1) return
      if(formal.and..not.p%formal_allowed) then
         measurements%limit_reason='FORMAL_EXECUTION_NOT_RELEASED';return
      end if
      nw=p%coarse_windows;target=p%coarse_min_window
      if(formal) then
         nw=p%formal_windows;target=p%formal_min_window
      end if
      measurements%requested_windows=nw;measurements%minimum_window=target
      statuses=ok
      do side=1,2
         call system_clock(count_rate=local_trace(side)%rate,count_max=local_trace(side)%clockmax)
      end do
      ! Warm and calibration also rotate. They never contain a reference build.
      do phase=1,2
         order=pair_order(phase,seed)
         do k=1,2
            side=order(k)
            if(statuses(side)/=ok) cycle
            if(side==1.and.phase==1) then
               call timed_call(left,left_values,cpu,wall,st,local_trace(side)%warm_start,local_trace(side)%warm_end)
            else if(side==1) then
               call timed_call(left,left_values,cpu,wall,st, &
                  local_trace(side)%calibration_start,local_trace(side)%calibration_end)
            else if(phase==1) then
               call timed_call(right,right_values,cpu,wall,st,local_trace(side)%warm_start,local_trace(side)%warm_end)
            else
               call timed_call(right,right_values,cpu,wall,st, &
                  local_trace(side)%calibration_start,local_trace(side)%calibration_end)
            end if
            side_wall(side)=side_wall(side)+wall
            if(phase==1) then
               measurements(side)%warm_recorded=.true.
               measurements(side)%warm_cpu=cpu;measurements(side)%warm_wall=wall
               if(side==1) lb=left_values
               if(side==2) rb=right_values
            else
               measurements(side)%calibration_recorded=.true.
               measurements(side)%calibration_cpu=cpu;measurements(side)%calibration_wall=wall
            end if
            if(st/=ok.or.cpu<0.or.wall<0) then
               statuses(side)=numerical_failure;measurements(side)%limit_reason='WARM_OR_CALIBRATION_FAILURE'
            else if(max(cpu,wall)>p%callback_limit_seconds.or.side_wall(side)>p%candidate_timeout_seconds) then
               statuses(side)=numerical_failure;measurements(side)%limit_reason='CALLBACK_OR_CANDIDATE_TIME_LIMIT'
            else if(phase==2) then
               if(side==1) then
                  if(any(left_values/=lb)) statuses(side)=numerical_failure
               else
                  if(any(right_values/=rb)) statuses(side)=numerical_failure
               end if
               if(statuses(side)/=ok) measurements(side)%limit_reason='CALIBRATION_RESPONSE_CHANGED'
            end if
         end do
      end do
      estimates=max(measurements%calibration_cpu,1e-7_rk)
      do w=1,nw
         order=pair_order(w,seed)
         do k=1,2
            side=order(k)
            if(statuses(side)/=ok) cycle
            if(side==1) then
               call timed_window(left,lb,left_values,measurements(side),w,p,target,estimates(side),st, &
                  local_trace(side)%window_start(w),local_trace(side)%window_end(w),side_wall(side))
            else
               call timed_window(right,rb,right_values,measurements(side),w,p,target,estimates(side),st, &
                  local_trace(side)%window_start(w),local_trace(side)%window_end(w),side_wall(side))
            end if
            side_wall(side)=side_wall(side)+measurements(side)%wall_seconds(w)
            statuses(side)=st
            if(side_wall(side)>p%candidate_timeout_seconds) then
               statuses(side)=numerical_failure;measurements(side)%limit_reason='CANDIDATE_TIME_LIMIT'
            end if
         end do
      end do
      do side=1,2
         if(statuses(side)==ok.and.measurements(side)%completed_windows==nw) call finish(measurements(side))
      end do
      if(present(trace)) trace=local_trace
   end subroutine

   subroutine timed_call(fn,values,cpu,wall,status,clock_start,clock_end)
      procedure(batch_iface) :: fn
      real(rk),intent(out) :: values(:),cpu,wall
      integer,intent(out) :: status
      integer(int64),intent(out) :: clock_start,clock_end
      real(rk) :: c0,c1
      integer(int64) :: t0,t1,rate,clockmax
      call system_clock(t0,rate,clockmax);call cpu_time(c0)
      call fn(values,status)
      call cpu_time(c1);call system_clock(t1)
      clock_start=t0;clock_end=t1
      cpu=c1-c0;wall=-1
      if(rate>0) wall=clock_elapsed(t0,t1,rate,clockmax)
      if(.not.all(ieee_is_finite(values)).or..not.ieee_is_finite(cpu)) status=numerical_failure
   end subroutine

   subroutine timed_window(fn,baseline,values,m,w,p,target,estimate,status,clock_start,clock_end,prior_wall)
      procedure(batch_iface) :: fn
      real(rk),intent(in) :: baseline(:),target,prior_wall
      real(rk),intent(out) :: values(:)
      type(measurement_t),intent(inout) :: m
      integer,intent(in) :: w
      type(policy_t),intent(in) :: p
      real(rk),intent(inout) :: estimate
      integer,intent(out) :: status
      integer(int64),intent(out) :: clock_start,clock_end
      real(rk) :: c0,c1,elapsed,wall
      integer(int64) :: t0,t1,rate,clockmax,count,chunk,j
      integer :: st
      status=invalid;m%window_count=w;count=0
      call system_clock(t0,rate,clockmax)
      clock_start=t0;clock_end=t0
      if(rate<=0) then
         m%limit_reason='SYSTEM_CLOCK_UNAVAILABLE';return
      end if
      chunk=max(1_int64,int(min(real(p%call_limit,rk),target/max(estimate,1e-7_rk)),int64))
      call cpu_time(c0)
      do
         chunk=min(chunk,p%call_limit-count)
         if(chunk<=0) then
            m%limit_reason='BATCH_LIMIT';exit
         end if
         do j=1,chunk
            call fn(values,st);count=count+1
            ! Observe each returned callback before another is allowed to run.
            ! These identical guard clocks are part of both sides' window cost.
            call cpu_time(c1);call system_clock(t1)
            elapsed=c1-c0;wall=clock_elapsed(t0,t1,rate,clockmax)
            if(st/=ok) exit
            if(c0<0.or.elapsed<0.or..not.ieee_is_finite(elapsed)) exit
            if(max(elapsed,wall)>p%callback_limit_seconds.or. &
               prior_wall+wall>p%candidate_timeout_seconds) exit
         end do
         clock_end=t1
         m%calls(w)=count;m%window_seconds(w)=elapsed;m%wall_seconds(w)=wall
         if(count>0) m%cpu_windows(w)=elapsed/real(count,rk)
         if(st/=ok) then
            m%limit_reason='CALLBACK_FAILURE';exit
         end if
         if(c0<0.or.elapsed<0.or..not.ieee_is_finite(elapsed)) then
            m%limit_reason='CPU_TIME_UNAVAILABLE';exit
         end if
         if(max(elapsed,wall)>p%callback_limit_seconds) then
            m%limit_reason='CALLBACK_OR_WINDOW_TIME_LIMIT';exit
         end if
         if(prior_wall+wall>p%candidate_timeout_seconds) then
            m%limit_reason='CANDIDATE_TIME_LIMIT';exit
         end if
         if(elapsed>=target) then
            if(.not.all(ieee_is_finite(values)).or.any(values/=baseline)) then
               m%limit_reason='WINDOW_RESPONSE_CHANGED';exit
            end if
            m%completed_windows=w;estimate=m%cpu_windows(w);status=ok;return
         end if
         if(elapsed>0) then
            estimate=elapsed/real(count,rk)
            chunk=max(1_int64,int(min(real(p%call_limit-count,rk),(target-elapsed)/estimate+1),int64))
         else
            chunk=min(p%call_limit-count,2*chunk)
         end if
      end do
      status=numerical_failure
   end subroutine

   subroutine finish(m)
      type(measurement_t),intent(inout) :: m
      real(rk) :: a(9),b(9)
      integer :: n
      n=m%completed_windows;a(:n)=m%cpu_windows(:n)
      b(:n)=m%wall_seconds(:n)/real(m%calls(:n),rk)
      call sort(a(:n));call sort(b(:n))
      m%median_cpu=median(a(:n));m%median_wall=median(b(:n))
      m%low_cpu=a(1);m%high_cpu=a(n)
   end subroutine

   pure real(rk) function clock_elapsed(a,b,rate,clockmax) result(v)
      integer(int64),intent(in) :: a,b,rate,clockmax
      v=real(b-a,rk)/real(rate,rk)
      if(b<a) v=(real(clockmax-a,rk)+real(b,rk)+1)/real(rate,rk)
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
