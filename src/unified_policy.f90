module unified_policy
   use, intrinsic :: iso_fortran_env, only: real64,int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: policy_t,load_policy,candidate_nodes,validate_policy
   type policy_t
      integer :: coarse_nodes(12)=[2,4,6,8,12,16,24,32,0,0,0,0],coarse_count=8
      integer :: audit_nodes(12)=[2,4,6,8,12,16,24,32,36,40,44,48],audit_count=12
      integer :: max_refinement_rounds=3,formal_selection_passes=3,coarse_windows=1,formal_windows=9,max_workers=8
      real(real64) :: coarse_min_window=.001_real64,formal_min_window=.004_real64
      real(real64) :: callback_limit_seconds=30.0_real64,candidate_timeout_seconds=60.0_real64
      integer(int64) :: call_limit=10000000_int64
      integer :: reference_job_timeout_seconds=7200,job_timeout_seconds=1800,public_job_timeout_seconds=7200
      real(real64) :: minimum_free_gib=12.0_real64
      logical :: formal_allowed=.false.
   end type
contains
   subroutine validate_policy(p,status)
      type(policy_t),intent(in) :: p
      integer,intent(out) :: status
      integer :: i
      status=1
      if(p%coarse_count<2.or.p%coarse_count>12.or.p%audit_count<2.or.p%audit_count>12)return
      if(any(p%coarse_nodes(:p%coarse_count)<2).or.any(p%coarse_nodes(:p%coarse_count)>32))return
      if(any(p%audit_nodes(:p%audit_count)<2).or.any(p%audit_nodes(:p%audit_count)>48))return
      do i=2,p%coarse_count
         if(p%coarse_nodes(i)<=p%coarse_nodes(i-1))return
      end do
      do i=2,p%audit_count
         if(p%audit_nodes(i)<=p%audit_nodes(i-1))return
      end do
      if(p%coarse_nodes(1)/=2.or.p%coarse_nodes(p%coarse_count)/=32)return
      if(p%audit_nodes(1)/=2.or.p%audit_nodes(p%audit_count)/=48)return
      if(p%max_refinement_rounds<0.or.p%max_refinement_rounds>3)return
      if(p%formal_selection_passes<1.or.p%formal_selection_passes>3)return
      if(p%coarse_windows/=1.or.p%formal_windows/=9)return
      if(p%max_workers<1.or.p%max_workers>8.or.p%call_limit<1)return
      if(.not.all(ieee_is_finite([p%coarse_min_window,p%formal_min_window,p%callback_limit_seconds, &
         p%candidate_timeout_seconds,p%minimum_free_gib])))return
      if(p%coarse_min_window<.001_real64.or.p%formal_min_window<.004_real64)return
      if(p%callback_limit_seconds<=p%formal_min_window.or.p%callback_limit_seconds>30.0_real64.or. &
         p%candidate_timeout_seconds<p%callback_limit_seconds.or.p%candidate_timeout_seconds>60.0_real64)return
      if(p%call_limit>10000000_int64)return
      if(min(p%reference_job_timeout_seconds,p%job_timeout_seconds,p%public_job_timeout_seconds)<=0)return
      if(p%minimum_free_gib<12.0_real64)return
      status=0
   end subroutine

   subroutine load_policy(path,p,status)
      character(*),intent(in) :: path
      type(policy_t),intent(out) :: p
      integer,intent(out) :: status
      integer :: coarse_nodes(12),coarse_count,audit_nodes(12),audit_count,max_refinement_rounds,formal_selection_passes
      integer :: coarse_windows,formal_windows,max_workers,reference_job_timeout_seconds
      integer :: job_timeout_seconds,public_job_timeout_seconds,u,ios
      integer(int64) :: call_limit
      real(real64) :: coarse_min_window,formal_min_window,callback_limit_seconds,candidate_timeout_seconds,minimum_free_gib
      logical :: formal_allowed
      namelist /policy/ coarse_nodes,coarse_count,audit_nodes,audit_count,max_refinement_rounds,formal_selection_passes, &
         coarse_windows,formal_windows,coarse_min_window,formal_min_window,callback_limit_seconds, &
         call_limit,formal_allowed,candidate_timeout_seconds,reference_job_timeout_seconds, &
         job_timeout_seconds,public_job_timeout_seconds,max_workers,minimum_free_gib
      p=policy_t();status=1
      coarse_nodes=p%coarse_nodes;coarse_count=p%coarse_count;audit_nodes=p%audit_nodes;audit_count=p%audit_count
      max_refinement_rounds=p%max_refinement_rounds;coarse_windows=p%coarse_windows;formal_windows=p%formal_windows
      formal_selection_passes=p%formal_selection_passes
      coarse_min_window=p%coarse_min_window;formal_min_window=p%formal_min_window
      callback_limit_seconds=p%callback_limit_seconds;call_limit=p%call_limit;formal_allowed=p%formal_allowed
      candidate_timeout_seconds=p%candidate_timeout_seconds;reference_job_timeout_seconds=p%reference_job_timeout_seconds
      job_timeout_seconds=p%job_timeout_seconds;public_job_timeout_seconds=p%public_job_timeout_seconds
      max_workers=p%max_workers;minimum_free_gib=p%minimum_free_gib
      open(newunit=u,file=path,status='old',action='read',iostat=ios)
      if(ios/=0)return
      read(u,nml=policy,iostat=ios);close(u)
      if(ios/=0)return
      p%coarse_nodes=coarse_nodes;p%coarse_count=coarse_count;p%audit_nodes=audit_nodes;p%audit_count=audit_count
      p%max_refinement_rounds=max_refinement_rounds;p%coarse_windows=coarse_windows;p%formal_windows=formal_windows
      p%formal_selection_passes=formal_selection_passes
      p%coarse_min_window=coarse_min_window;p%formal_min_window=formal_min_window
      p%callback_limit_seconds=callback_limit_seconds;p%call_limit=call_limit;p%formal_allowed=formal_allowed
      p%candidate_timeout_seconds=candidate_timeout_seconds;p%reference_job_timeout_seconds=reference_job_timeout_seconds
      p%job_timeout_seconds=job_timeout_seconds;p%public_job_timeout_seconds=public_job_timeout_seconds
      p%max_workers=max_workers;p%minimum_free_gib=minimum_free_gib
      call validate_policy(p,status)
   end subroutine

   subroutine candidate_nodes(p,audit,allow_one,nodes,status)
      type(policy_t),intent(in) :: p
      logical,intent(in) :: audit,allow_one
      integer,allocatable,intent(out) :: nodes(:)
      integer,intent(out) :: status
      integer :: n,offset
      call validate_policy(p,status)
      if(status/=0)then
         allocate(nodes(0));return
      end if
      n=p%coarse_count
      if(audit)n=p%audit_count
      offset=0
      if(allow_one)offset=1
      allocate(nodes(n+offset))
      if(allow_one)nodes(1)=1
      if(audit)then
         nodes(offset+1:)=p%audit_nodes(:n)
      else
         nodes(offset+1:)=p%coarse_nodes(:n)
      end if
   end subroutine
end module
