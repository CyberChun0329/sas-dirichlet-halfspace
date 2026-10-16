! Explicit finite-list timing entry point. Native mathematics and timers are reused.
program targeted_driver
   use, intrinsic :: iso_fortran_env, only: int64,iostat_end
   use replay_types
   use replay_case_bank, only: case_spec_t,build_cases
   use replay_select, only: candidate_t,route_family,analyze_candidates
   use replay_measure, only: measurement_t
   use useful_screening, only: axis_cost_bounds
   use unified_policy, only: policy_t,load_policy
   use unified_engine, only: descriptor_t,prepared_t,make_descriptor,prepare_candidate,evaluate_candidate, &
      preparation_fingerprint,candidate_errors
   use unified_poly_reference, only: reference_t,read_reference
   use paired_measure, only: measure_paired,pair_order,pair_trace_t
   implicit none
   character(*),parameter :: task_header='job,family,case_index,repeats,reference,rows,values,keys,outdir,prep,budget,shared'
   character(*),parameter :: key_header='key,route,m,rank,nx,ny'
   character(*),parameter :: row_header='release,case,phase,axis,key,route,m,rank,nx,ny,EA,ER,query_mid,query_low,query_high,total_mid,total_low,total_high,numerical_valid,active,query_valid,total_valid,execution_state,candidate_valid,correctness_value_valid,correctness_numerical_status,correctness_finite,correctness_resource_state,correctness_resource_reason,correctness_status,correctness_state'
   character(*),parameter :: value_header='release,case,round,key,point,correctness_value,reference'
   character(*),parameter :: window_header='release,case,phase,round,pair,side,key,scope,window,order,calls,cpu_elapsed,wall_elapsed,cpu_per_batch,minimum_cpu_seconds,requested_windows,completed_windows,status,reason,clock_start,clock_end,clock_rate,clock_max'
   type correctness_t
      logical :: value_valid=.false.,finite=.false.,candidate_valid=.false.
      integer :: numerical_status=invalid,status=invalid,round=-1
      character(16) :: resource_state=''
      character(48) :: resource_reason=''
      character(32) :: state=''
   end type
   type(policy_t) :: policy
   type(case_spec_t),allocatable :: cases(:)
   type(case_spec_t) :: c
   type(reference_t) :: reference
   type(descriptor_t),allocatable :: bank(:)
   type(candidate_t),allocatable :: cpu(:),wall(:)
   type(correctness_t),allocatable :: correctness(:)
   type(prepared_t) :: prepared(2)
   type(measurement_t) :: measured(2)
   type(pair_trace_t) :: trace(2)
   real(rk),allocatable :: responses(:,:),lv(:),rv(:),caps(:)
   logical,allocatable :: row_seen(:,:),point_seen(:,:),eligible(:)
   character(4096),allocatable :: retained_rows(:,:)
   integer :: ci,repeats,ios,st,i,scope,pair_id,side,metric,setup(2),active_indices(2),statuses(2)
   integer :: urows,uwin,uprep,ucov,umeta,uprogress,ucorrect,first,last,ref_repeats
   integer :: preparation_attempts=0,scope_attempts(3)=0
   integer(int64) :: callback_calls(3)=0,anchor_before,anchor_after,anchor_rate,anchor_max
   character(1024) :: tasks_path,policy_path,reference_path,rows_path,values_path,keys_path,outdir
   character(128) :: release_id,job_id
   character(64) :: family,mode
   character(8) :: anchor_date
   character(10) :: anchor_time
   character(5) :: anchor_zone
   logical :: prep,budget,shared,audit,ref_audit,response_mismatch(2)
   call get_command_argument(1,mode);call get_command_argument(2,tasks_path)
   call get_command_argument(3,job_id);call get_command_argument(4,policy_path);call get_command_argument(5,release_id)
   if(command_argument_count()/=5.or.(mode/='--check'.and.mode/='--run')) &
      error stop 'targeted_driver --check|--run TASKS.csv JOB_ID POLICY.nml RELEASE_ID'
   if(min(len_trim(job_id),len_trim(release_id))==0) error stop 'empty job/release identity'
   call require_path(tasks_path);call require_path(policy_path)
   call load_policy(trim(policy_path),policy,st);call must(st,'policy')
   if(.not.policy%formal_allowed) error stop 'explicit formal timing policy required'
   call load_task()
   call build_cases(trim(family),cases,st);call must(st,'native case inputs')
   if(ci<1.or.ci>size(cases)) error stop 'task case index out of range'
   call repeated_case(ci,c)
   call load_reference()
   call load_keys()
   allocate(cpu(size(bank)),wall(size(bank)),correctness(size(bank)),responses(size(c%q,2),size(bank)), &
      row_seen(2,size(bank)),point_seen(size(c%q,2),size(bank)),eligible(size(bank)),retained_rows(2,size(bank)))
   responses=quiet_nan();row_seen=.false.;point_seen=.false.;retained_rows=''
   do i=1,size(bank)
      cpu(i)%key=bank(i)%key;cpu(i)%route=bank(i)%label
      cpu(i)%m=bank(i)%m;cpu(i)%rank=bank(i)%rank;cpu(i)%n=bank(i)%nx
      cpu(i)%nx=bank(i)%nx;cpu(i)%ny=bank(i)%ny
      cpu(i)%cpu=quiet_nan();cpu(i)%low=quiet_nan();cpu(i)%high=quiet_nan()
      cpu(i)%total=quiet_nan();cpu(i)%total_low=quiet_nan();cpu(i)%total_high=quiet_nan()
   end do
   call load_rows();call load_values()
   wall=cpu
   eligible=cpu%valid.and.reference%certified
   if(family=='deficit_original') then
      caps=[1e-4_rk,1e-6_rk,1e-8_rk]
   else
      caps=[1e-6_rk]
   end if
   if(mode=='--check') then
      print '(a)', 'CHECK_OK job='//trim(job_id)//' case='//trim(c%id)
      print '(a,i0,a,i0)', 'selected_keys=',size(bank),' eligible_keys=',count(eligible)
      print '(a)', 'coarse_rounds=0 fresh_correctness=0 generated_candidates=0 solver_calls=0 output_creation=0'
      stop
   end if
   call make_outputs()
   allocate(lv(size(c%q,2)),rv(size(c%q,2)))
   call measure_selected()
   call save_rows()
   call analyze_candidates(trim(c%id),trim(c%profile),c%rho,c%Fo,cpu,reference%floor,reference%afloor, &
      trim(outdir)//'/cpu',st,caps=caps);call must(st,'CPU finite-list selector')
   call analyze_candidates(trim(c%id),trim(c%profile),c%rho,c%Fo,wall,reference%floor,reference%afloor, &
      trim(outdir)//'/wall',st,caps=caps);call must(st,'WALL finite-list selector')
   write(ucov,'(*(g0,:,","))')trim(release_id),trim(job_id),trim(c%id),repeats,reference%certified, &
      size(bank),count(eligible),0,0,0,preparation_attempts,scope_attempts, &
      count(cpu%query_valid),count(cpu%total_valid),count(wall%query_valid),count(wall%total_valid),'FINITE_LISTED_POOL'
   write(umeta,*)'query_context_preparation_attempts=',preparation_attempts
   write(umeta,*)'prep_query_total_candidate_scope_attempts=',scope_attempts
   write(umeta,*)'prep_query_total_callback_calls=',callback_calls
   write(uprogress,'(a)')'COMPLETE FINITE_LISTED_POOL'
   close(urows);close(uwin);close(uprep);close(ucov);close(umeta);close(uprogress);close(ucorrect)
contains
   subroutine must(code,label)
      integer,intent(in) :: code
      character(*),intent(in) :: label
      if(code/=ok) then
         print *,trim(label),code
         error stop 'invalid targeted operation'
      end if
   end subroutine

   subroutine fields(line,f,nf)
      character(*),intent(in) :: line
      character(*),intent(out) :: f(:)
      integer,intent(out) :: nf
      integer :: a,b,l
      f='';nf=0;a=1;l=len_trim(line)
      if(index(line,'"')>0.or.index(line,achar(13))>0) error stop 'quoted/control CSV unsupported'
      do
         nf=nf+1;b=index(line(a:l),',')
         if(b==0) then
            if(nf<=size(f)) then
               if(l-a+1>len(f(nf))) error stop 'CSV field too long'
               f(nf)=line(a:l)
            end if
            exit
         end if
         b=a+b-1
         if(nf<=size(f)) then
            if(b-a>len(f(nf))) error stop 'CSV field too long'
            f(nf)=line(a:b-1)
         end if
         a=b+1
      end do
   end subroutine

   subroutine require_path(path)
      character(*),intent(in) :: path
      integer :: p,ch
      if(len_trim(path)==0) error stop 'empty relative path'
      if(path(1:1)=='/'.or.index('/'//trim(path)//'/','/../')>0.or.index(path,'//')>0) &
         error stop 'root-relative path required'
      do p=1,len_trim(path)
         ch=iachar(path(p:p))
         if((ch>=iachar('a').and.ch<=iachar('z')).or.(ch>=iachar('A').and.ch<=iachar('Z')).or. &
            (ch>=iachar('0').and.ch<=iachar('9')).or.index('_-./',path(p:p))>0) cycle
         error stop 'unsupported character in root-relative path'
      end do
   end subroutine

   logical function flag(value)
      character(*),intent(in) :: value
      if(value/='T'.and.value/='F') error stop 'flags must be literal T/F'
      flag=value=='T'
   end function

   subroutine load_task()
      character(4096) :: line
      character(1024) :: f(12)
      integer :: u,nf,matches
      open(newunit=u,file=trim(tasks_path),status='old',action='read',iostat=ios)
      if(ios/=0) error stop 'task file unavailable'
      read(u,'(a)',iostat=ios)line
      if(ios/=0.or.trim(line)/=task_header) error stop 'task header mismatch'
      matches=0
      do
         read(u,'(a)',iostat=ios)line
         if(ios==iostat_end) exit
         if(ios/=0) error stop 'task read failure'
         if(len_trim(line)==0) cycle
         call fields(line,f,nf)
         if(nf/=12) error stop 'task field count'
         if(f(1)/=job_id) cycle
         matches=matches+1
         if(matches/=1) error stop 'duplicate task job identity'
         if(len_trim(f(2))>len(family)) error stop 'task family too long'
         family=f(2)(:64)
         read(f(3),*,iostat=ios)ci
         if(ios/=0) error stop 'task case index'
         read(f(4),*,iostat=ios)repeats
         if(ios/=0.or..not.any(repeats==[1,10,100])) error stop 'task repeat count'
         reference_path=f(5);rows_path=f(6);values_path=f(7);keys_path=f(8);outdir=f(9)
         prep=flag(f(10));budget=flag(f(11));shared=flag(f(12))
      end do
      close(u)
      if(matches/=1) error stop 'task job not found'
      if(family/='paper'.and.family/='temperature_scale'.and.family/='deficit_original') &
         error stop 'supported families: paper, temperature_scale, deficit_original'
      call require_path(reference_path);call require_path(rows_path);call require_path(values_path)
      call require_path(keys_path);call require_path(outdir)
   end subroutine

   subroutine repeated_case(index,curr)
      integer,intent(in) :: index
      type(case_spec_t),intent(out) :: curr
      real(rk),allocatable :: batch(:,:)
      integer :: t,np
      curr=cases(index);np=size(curr%q,2)
      if(repeats==1) return
      allocate(batch(4,np*repeats))
      do t=1,repeats
         batch(:,(t-1)*np+1:t*np)=curr%q
      end do
      curr%q=batch
   end subroutine

   subroutine load_reference()
      character(128) :: magic,stored_release,stored_family
      type(case_spec_t) :: scan_case
      type(reference_t) :: scan_ref
      integer :: u,t
      open(newunit=u,file=trim(reference_path),status='old',action='read',iostat=ios)
      if(ios/=0) error stop 'existing native reference unavailable'
      read(u,'(a)',iostat=ios)magic
      if(ios/=0.or.magic/='UNIFIED_POLY_REFERENCE_V1') error stop 'native reference format'
      read(u,'(a)',iostat=ios)stored_release
      if(ios/=0.or.stored_release/=release_id) error stop 'native reference release'
      read(u,'(a)',iostat=ios)stored_family
      if(ios/=0.or.stored_family/=family) error stop 'native reference family'
      read(u,*,iostat=ios)first,last,ref_repeats,ref_audit
      if(ios/=0) error stop 'native reference batch header'
      if(first<1.or.last>size(cases).or.first>last.or.ref_repeats/=repeats.or.ci<first.or.ci>last) &
         error stop 'native reference batch does not contain target'
      audit=ref_audit
      if(audit.and..not.cases(ci)%audit) error stop 'target excluded by native audit batch'
      do t=first,ci
         if(audit.and..not.cases(t)%audit) cycle
         call repeated_case(t,scan_case)
         call read_reference(u,scan_case,scan_ref,st);call must(st,'native reference scan identity')
         if(t==ci) reference=scan_ref
      end do
      close(u)
      if(.not.allocated(reference%values)) error stop 'native reference target missing'
      if(.not.all(ieee_is_finite([reference%floor,reference%afloor])).or.reference%floor<=0.or.reference%afloor<0) &
         error stop 'invalid stored comparison floors'
   end subroutine

   subroutine load_keys()
      character(4096) :: line
      character(256) :: f(6)
      character(16) :: fam
      type(descriptor_t) :: d
      integer :: u,nf,route,m,r,nx,ny,ns,nd,nh
      allocate(bank(0))
      open(newunit=u,file=trim(keys_path),status='old',action='read',iostat=ios)
      if(ios/=0) error stop 'explicit keys unavailable'
      read(u,'(a)',iostat=ios)line
      if(ios/=0.or.trim(line)/=key_header) error stop 'key header mismatch'
      ns=0;nd=0;nh=0
      do
         read(u,'(a)',iostat=ios)line
         if(ios==iostat_end) exit
         if(ios/=0) error stop 'key read failure'
         if(len_trim(line)==0) cycle
         call fields(line,f,nf)
         if(nf/=6) error stop 'key field count'
         read(f(2),*,iostat=ios)route
         if(ios/=0) error stop 'key route'
         read(f(3),*,iostat=ios)m
         if(ios/=0) error stop 'key degree'
         read(f(4),*,iostat=ios)r
         if(ios/=0) error stop 'key rank'
         read(f(5),*,iostat=ios)nx
         if(ios/=0) error stop 'key nx'
         read(f(6),*,iostat=ios)ny
         if(ios/=0) error stop 'key ny'
         if(route<1.or.route>8.or.min(nx,ny)<1.or.max(nx,ny)>48.or.m>32.or.r>32) &
            error stop 'key descriptor outside native range'
         d=make_descriptor(route,m,r,nx,ny)
         if(d%key/=f(1).or.len_trim(d%label)==0) error stop 'key descriptor conflict'
         if(any(bank%key==d%key)) error stop 'duplicate explicit key'
         fam=route_family(d%label)
         select case(trim(fam))
         case('SAS')
            if(m<0.or.m>max(c%source%mx,c%source%my).or.r/=-1.or.nx/=ny) &
               error stop 'invalid or exact-zero-alias SAS descriptor'
            ns=ns+1
         case('DIRECT_2D');nd=nd+1
         case('SHARED_1D')
            if(.not.shared) error stop 'unreleased shared route'
            nh=nh+1
         case default;error stop 'unsupported route family'
         end select
         if(ns>merge(8,6,budget).or.nd>merge(16,12,budget).or.nh>6) error stop 'explicit key cap exceeded'
         bank=[bank,d]
      end do
      close(u)
      if(size(bank)==0) error stop 'empty explicit key list'
   end subroutine

   integer function key_index(key) result(idx)
      character(*),intent(in) :: key
      integer :: a
      idx=0
      do a=1,size(bank)
         if(bank(a)%key==key) then
            idx=a;return
         end if
      end do
   end function

   subroutine load_rows()
      character(4096) :: line
      character(256) :: f(31)
      type(correctness_t) :: rec
      integer :: u,nf,a,axis,degree,rank,nx,ny
      real(rk) :: ea,er
      logical :: valid,active,candidate_valid
      open(newunit=u,file=trim(rows_path),status='old',action='read',iostat=ios)
      if(ios/=0) error stop 'cached rows unavailable'
      read(u,'(a)',iostat=ios)line
      if(ios/=0.or.trim(line)/=row_header) error stop 'cached 31-field row header mismatch'
      do
         read(u,'(a)',iostat=ios)line
         if(ios==iostat_end) exit
         if(ios/=0) error stop 'cached row read failure'
         if(len_trim(line)==0) cycle
         call fields(line,f(:5),nf)
         if(nf<5.or.f(2)/=c%id) cycle
         a=key_index(f(5));if(a==0) cycle
         call fields(line,f,nf)
         if(nf/=31.or.f(1)/=release_id.or.f(3)/='coarse') error stop 'selected cached row identity'
         read(f(4),*,iostat=ios)axis
         if(ios/=0.or.axis<1.or.axis>2) error stop 'selected cached row axis'
         if(row_seen(axis,a)) error stop 'duplicate selected key/axis cache row'
         if(f(6)/=bank(a)%label) error stop 'selected route label conflict'
         read(f(7),*,iostat=ios)degree
         if(ios/=0) error stop 'cached row degree'
         read(f(8),*,iostat=ios)rank
         if(ios/=0) error stop 'cached row rank'
         read(f(9),*,iostat=ios)nx
         if(ios/=0) error stop 'cached row nx'
         read(f(10),*,iostat=ios)ny
         if(ios/=0) error stop 'cached row ny'
         if(any([degree,rank,nx,ny]/=[bank(a)%m,bank(a)%rank,bank(a)%nx,bank(a)%ny])) &
            error stop 'selected complete descriptor mismatch'
         read(f(11),*,iostat=ios)ea
         if(ios/=0) error stop 'cached EA'
         read(f(12),*,iostat=ios)er
         if(ios/=0) error stop 'cached ER'
         valid=flag(f(19));active=flag(f(20));candidate_valid=flag(f(24))
         if(valid.neqv.candidate_valid) error stop 'cached numerical/candidate validity conflict'
         rec=correctness_t();rec%candidate_valid=candidate_valid
         rec%value_valid=flag(f(25));rec%finite=flag(f(27))
         read(f(26),*,iostat=ios)rec%numerical_status
         if(ios/=0) error stop 'cached numerical status'
         read(f(30),*,iostat=ios)rec%status
         if(ios/=0) error stop 'cached correctness status'
         if(len_trim(f(28))>16.or.len_trim(f(29))>48.or.len_trim(f(31))>32) error stop 'cached status truncation'
         rec%resource_state=f(28)(:16);rec%resource_reason=f(29)(:48);rec%state=f(31)(:32)
         if(rec%value_valid.neqv.(rec%numerical_status==ok.and.rec%finite)) error stop 'cached value-valid conflict'
         if(valid.neqv.(rec%status==ok)) error stop 'cached overall validity conflict'
         if(valid.and.(.not.rec%value_valid.or.rec%resource_state/='ok'.or.rec%state/='OK')) &
            error stop 'cached valid row has HOLD status'
         if(any(row_seen(:,a))) then
            if(.not.same_scalar(ea,cpu(a)%EA).or..not.same_scalar(er,cpu(a)%ER).or. &
               (valid.neqv.cpu(a)%valid).or.(active.neqv.cpu(a)%active)) error stop 'cache axes disagree'
            if((rec%candidate_valid.neqv.correctness(a)%candidate_valid).or. &
               (rec%value_valid.neqv.correctness(a)%value_valid).or.(rec%finite.neqv.correctness(a)%finite).or. &
               rec%numerical_status/=correctness(a)%numerical_status.or.rec%status/=correctness(a)%status.or. &
               rec%resource_state/=correctness(a)%resource_state.or.rec%resource_reason/=correctness(a)%resource_reason.or. &
               rec%state/=correctness(a)%state) error stop 'cache correctness axes disagree'
         else
            cpu(a)%EA=ea;cpu(a)%ER=er;cpu(a)%valid=valid;cpu(a)%active=active;correctness(a)=rec
         end if
         retained_rows(axis,a)=line;row_seen(axis,a)=.true.
      end do
      close(u)
      if(.not.all(row_seen)) error stop 'selected cache key missing an axis'
   end subroutine

   subroutine load_values()
      character(4096) :: line
      character(256) :: f(7)
      integer :: u,nf,a,point,round
      real(rk) :: value,ref,ea,er
      logical :: active
      open(newunit=u,file=trim(values_path),status='old',action='read',iostat=ios)
      if(ios/=0) error stop 'cached values unavailable'
      read(u,'(a)',iostat=ios)line
      if(ios/=0.or.trim(line)/=value_header) error stop 'cached independent value header mismatch'
      do
         read(u,'(a)',iostat=ios)line
         if(ios==iostat_end) exit
         if(ios/=0) error stop 'cached value read failure'
         if(len_trim(line)==0) cycle
         call fields(line,f(:4),nf)
         if(nf<4.or.f(2)/=c%id) cycle
         a=key_index(f(4));if(a==0) cycle
         call fields(line,f,nf)
         if(nf/=7.or.f(1)/=release_id) error stop 'selected independent value identity'
         read(f(3),*,iostat=ios)round
         if(ios/=0.or.round<0) error stop 'selected independent value round'
         if(correctness(a)%round<0) correctness(a)%round=round
         if(round/=correctness(a)%round) error stop 'selected independent value round conflict'
         read(f(5),*,iostat=ios)point
         if(ios/=0.or.point<1.or.point>size(reference%values)) error stop 'selected independent value point'
         if(point_seen(point,a)) error stop 'duplicate selected independent value point'
         read(f(6),*,iostat=ios)value
         if(ios/=0) error stop 'selected independent numeric value'
         read(f(7),*,iostat=ios)ref
         if(ios/=0.or..not.same_scalar(ref,reference%values(point))) error stop 'selected reference bits mismatch'
         responses(point,a)=value;point_seen(point,a)=.true.
      end do
      close(u)
      if(.not.all(point_seen)) error stop 'selected independent numeric points incomplete'
      do a=1,size(bank)
         if(correctness(a)%finite.neqv.all(ieee_is_finite(responses(:,a)))) error stop 'cached finite flag conflict'
         if(.not.correctness(a)%value_valid) cycle
         call candidate_errors(responses(:,a),reference%values,ea,er,active)
         if(.not.same_scalar(ea,cpu(a)%EA).or..not.same_scalar(er,cpu(a)%ER).or.(active.neqv.cpu(a)%active)) &
            error stop 'cached error reduction conflict'
      end do
   end subroutine

   pure logical function same_scalar(a,b)
      real(rk),intent(in) :: a,b
      same_scalar=transfer(a,0_int64)==transfer(b,0_int64)
   end function

   pure logical function same_bits(a,b)
      real(rk),intent(in) :: a(:),b(:)
      integer :: t
      same_bits=.false.
      if(size(a)/=size(b)) return
      do t=1,size(a)
         if(.not.same_scalar(a(t),b(t))) return
      end do
      same_bits=.true.
   end function

   subroutine open_new(unit,name)
      integer,intent(out) :: unit
      character(*),intent(in) :: name
      open(newunit=unit,file=trim(outdir)//'/'//name,status='new',action='write',iostat=ios)
      if(ios/=0) error stop 'strict-new output file required'
   end subroutine

   subroutine make_outputs()
      logical :: exists
      integer :: exitstat,cmdstat,a,axis
      inquire(file=trim(outdir),exist=exists)
      if(exists) error stop 'strict-new output directory required'
      call execute_command_line('mkdir '//trim(outdir),exitstat=exitstat,cmdstat=cmdstat)
      if(cmdstat/=0.or.exitstat/=0) error stop 'cannot create strict-new output directory'
      call execute_command_line('mkdir '//trim(outdir)//'/cpu '//trim(outdir)//'/wall',exitstat=exitstat,cmdstat=cmdstat)
      if(cmdstat/=0.or.exitstat/=0) error stop 'cannot create selector output directories'
      call open_new(urows,'rows.csv');call open_new(uwin,'windows.csv');call open_new(uprep,'preparation.csv')
      call open_new(ucov,'coverage.csv');call open_new(umeta,'metadata.txt');call open_new(uprogress,'progress.log')
      call open_new(ucorrect,'cached_correctness_rows.csv')
      write(urows,'(a)')row_header;write(uwin,'(a)')window_header;write(ucorrect,'(a)')row_header
      write(uprep,'(a)')'release,case,phase,round,pair,key,cpu_median,wall_median,cpu_low,cpu_high,status,reason'
      write(ucov,'(a)')'release,job,case,repeat,reference_certified,selected_keys,eligible_keys,coarse_rounds,fresh_correctness,generated_candidates,query_context_preparations,prep_attempted,query_attempted,total_attempted,cpu_query_complete,cpu_total_complete,wall_query_complete,wall_total_complete,status'
      do axis=1,2
         do a=1,size(bank)
            write(ucorrect,'(a)')trim(retained_rows(axis,a))
         end do
      end do
      call system_clock(anchor_before,anchor_rate,anchor_max)
      call date_and_time(anchor_date,anchor_time,anchor_zone);call system_clock(anchor_after)
      write(umeta,'(a)')'release='//trim(release_id)
      write(umeta,'(a)')'job='//trim(job_id)
      write(umeta,'(a)')'mode=EXPLICIT_FINITE_LISTED_POOL'
      write(umeta,'(a)')'clock_anchor_date='//anchor_date//' time='//anchor_time//' zone='//anchor_zone
      write(umeta,*)'clock_anchor_before_after_rate_max=',anchor_before,anchor_after,anchor_rate,anchor_max
      write(umeta,'(a)')'reference='//trim(reference_path)
      write(umeta,'(a)')'cached_rows='//trim(rows_path)
      write(umeta,'(a)')'cached_values='//trim(values_path)
      write(umeta,'(a)')'explicit_keys='//trim(keys_path)
      write(umeta,'(a)')'coarse_rounds=0'
      write(umeta,'(a)')'fresh_correctness=0'
      write(umeta,'(a)')'generated_candidates=0'
      write(umeta,'(a)')'reference_generation=0'
      write(umeta,'(a)')'cost_protocol=native_paired_formal_nine_windows_explicit_one_to_one_v1'
      write(umeta,'(a)')'comparison_scope=finite_listed_pool; ALL_SPACE uses only measured non-SAS listed keys'
      write(umeta,'(a)')'reference_certified is the stored refinement/cross-check decision, not a true-error proof'
      write(umeta,'(a)')'cached_correctness_rows.csv retains historical costs for provenance only; rows.csv contains new costs'
      write(umeta,*)'prep_budget_shared=',prep,budget,shared
      write(umeta,*)'caps=',caps
      write(umeta,*)'formal_windows_minimum_seconds=',policy%formal_windows,policy%formal_min_window
      write(umeta,*)'callback_return_soft_limit_seconds=',policy%callback_limit_seconds
      write(umeta,*)'candidate_scope_soft_limit_seconds=',policy%candidate_timeout_seconds
      write(uprogress,'(a)')'START EXPLICIT_FINITE_LISTED_POOL'
      flush(umeta);flush(ucorrect);flush(uprogress)
   end subroutine

   subroutine measure_selected()
      integer,allocatable :: sas(:),space(:)
      integer :: a,p,s
      logical :: sas_mask(size(bank))
      do a=1,size(bank)
         sas_mask(a)=route_family(bank(a)%label)=='SAS'
      end do
      sas=pack([(a,a=1,size(bank))],eligible.and.sas_mask)
      space=pack([(a,a=1,size(bank))],eligible.and..not.sas_mask)
      pair_id=0
      do p=1,max(size(sas),size(space))
         active_indices=0
         if(p<=size(sas)) active_indices(1)=sas(p)
         if(p<=size(space)) active_indices(2)=space(p)
         setup=invalid
         do s=1,2
            a=active_indices(s);if(a==0) cycle
            preparation_attempts=preparation_attempts+1
            call prepare_candidate(c%source,bank(a),prepared(s),setup(s))
         end do
         pair_id=pair_id+1
         do scope=1,3
            if(scope==1.and..not.prep) cycle
            response_mismatch=.false.;lv=quiet_nan();rv=quiet_nan()
            scope_attempts(scope)=scope_attempts(scope)+count(active_indices>0)
            call measure_paired(left_callback,right_callback,lv,rv,measured,statuses,policy,.true.,pair_id+scope,trace)
            do side=1,2
               a=active_indices(side);if(a==0) cycle
               if(response_mismatch(side)) then
                  statuses(side)=numerical_failure;measured(side)%limit_reason='CORRECTNESS_RESPONSE_MISMATCH'
               end if
               call write_windows(side)
               if(scope==1) then
                  call write_preparation(a)
               else
                  metric=scope-1
                  call assign_cost(cpu(a),wall(a),measured(side),metric,statuses(side))
               end if
               write(uprogress,'(*(g0,:,","))')'SCOPE',pair_id,scope,trim(bank(a)%key),statuses(side), &
                  measured(side)%completed_windows,trim(measured(side)%limit_reason)
            end do
            flush(uwin);flush(uprep);flush(uprogress)
         end do
      end do
   end subroutine

   subroutine callback(which,v,s)
      integer,intent(in) :: which
      real(rk),intent(out) :: v(:)
      integer,intent(out) :: s
      type(prepared_t) :: fresh
      integer :: a
      ! An unmatched side has no descriptor, preparation, query or output row.
      v=0;s=invalid;a=active_indices(which)
      if(a==0) return
      callback_calls(scope)=callback_calls(scope)+1
      v=quiet_nan()
      select case(scope)
      case(1)
         call prepare_candidate(c%source,bank(a),fresh,s)
         if(s==ok) call preparation_fingerprint(fresh,v,s)
      case(2)
         s=setup(which)
         if(s==ok) call evaluate_candidate(prepared(which),c%q,c%observable,v,s)
      case(3)
         call prepare_candidate(c%source,bank(a),fresh,s)
         if(s==ok) call evaluate_candidate(fresh,c%q,c%observable,v,s)
      end select
      if(scope/=1.and.s==ok) then
         if(.not.same_bits(v,responses(:,a))) then
            s=numerical_failure;response_mismatch(which)=.true.
         end if
      end if
   end subroutine

   subroutine left_callback(v,s)
      real(rk),intent(out) :: v(:)
      integer,intent(out) :: s
      call callback(1,v,s)
   end subroutine

   subroutine right_callback(v,s)
      real(rk),intent(out) :: v(:)
      integer,intent(out) :: s
      call callback(2,v,s)
   end subroutine

   subroutine assign_cost(a,b,m,which,code)
      type(candidate_t),intent(inout) :: a,b
      type(measurement_t),intent(in) :: m
      integer,intent(in) :: which,code
      real(rk) :: mid,lo,hi
      integer :: axis,s
      do axis=1,2
         call axis_cost_bounds(m,axis,.false.,mid,lo,hi,s)
         if(code/=ok) s=invalid
         if(s/=ok) then
            mid=quiet_nan();lo=quiet_nan();hi=quiet_nan()
         end if
         if(axis==1) then
            call set_cost(a,which,mid,lo,hi,s)
         else
            call set_cost(b,which,mid,lo,hi,s)
         end if
      end do
   end subroutine

   subroutine set_cost(row,which,mid,lo,hi,code)
      type(candidate_t),intent(inout) :: row
      integer,intent(in) :: which,code
      real(rk),intent(in) :: mid,lo,hi
      if(which==1) then
         row%cpu=mid;row%low=lo;row%high=hi;row%query_valid=code==ok
      else
         row%total=mid;row%total_low=lo;row%total_high=hi;row%total_valid=code==ok
      end if
   end subroutine

   subroutine write_preparation(a)
      integer,intent(in) :: a
      real(rk) :: mid(2),lo(2),hi(2)
      integer :: s(2),axis
      do axis=1,2
         call axis_cost_bounds(measured(side),axis,.false.,mid(axis),lo(axis),hi(axis),s(axis))
         if(statuses(side)/=ok.or.s(axis)/=ok) then
            mid(axis)=quiet_nan();lo(axis)=quiet_nan();hi(axis)=quiet_nan()
         end if
      end do
      write(uprep,'(*(g0,:,","))')trim(release_id),trim(c%id),'formal',0,pair_id,trim(bank(a)%key), &
         mid(1),mid(2),lo(1),hi(1),statuses(side),trim(measured(side)%limit_reason)
   end subroutine

   subroutine write_windows(which)
      integer,intent(in) :: which
      integer :: a,w,order(2),position
      type(measurement_t) :: m
      a=active_indices(which);m=measured(which)
      if(m%warm_recorded) then
         order=pair_order(1,pair_id+scope);position=findloc(order,which,dim=1)
         call window_line(which,a,m,-1,position,1_int64,m%warm_cpu,m%warm_wall,0.0_rk)
      end if
      if(m%calibration_recorded) then
         order=pair_order(2,pair_id+scope);position=findloc(order,which,dim=1)
         call window_line(which,a,m,0,position,1_int64,m%calibration_cpu,m%calibration_wall,0.0_rk)
      end if
      do w=1,m%window_count
         order=pair_order(w,pair_id+scope);position=findloc(order,which,dim=1)
         call window_line(which,a,m,w,position,m%calls(w),m%window_seconds(w),m%wall_seconds(w),m%cpu_windows(w))
      end do
   end subroutine

   subroutine window_line(which,a,m,w,position,calls,cpu_elapsed,wall_elapsed,per_batch)
      integer,intent(in) :: which,a,w,position
      type(measurement_t),intent(in) :: m
      integer(int64),intent(in) :: calls
      real(rk),intent(in) :: cpu_elapsed,wall_elapsed,per_batch
      integer(int64) :: clock_start,clock_end
      character(5),parameter :: names(3)=[character(5)::'PREP','QUERY','TOTAL']
      clock_start=trace(which)%warm_start;clock_end=trace(which)%warm_end
      if(w==0) then
         clock_start=trace(which)%calibration_start;clock_end=trace(which)%calibration_end
      else if(w>0) then
         clock_start=trace(which)%window_start(w);clock_end=trace(which)%window_end(w)
      end if
      write(uwin,'(*(g0,:,","))')trim(release_id),trim(c%id),'formal',0,pair_id,which,trim(bank(a)%key), &
         trim(names(scope)),w,position,calls,cpu_elapsed,wall_elapsed,per_batch,m%minimum_window, &
         m%requested_windows,m%completed_windows,statuses(which),trim(m%limit_reason), &
         clock_start,clock_end,trace(which)%rate,trace(which)%clockmax
   end subroutine

   subroutine save_rows()
      type(candidate_t) :: r
      character(40) :: state
      integer :: a,axis
      do axis=1,2
         do a=1,size(bank)
            r=cpu(a);if(axis==2) r=wall(a)
            state='MEASURED'
            if(.not.r%query_valid.or..not.r%total_valid) state='COST_SCOPE_HOLD'
            if(.not.r%active) state='HELD_ABSOLUTE_ONLY'
            if(.not.r%valid) state=correctness(a)%state
            if(.not.reference%certified) state='HOLD_REFERENCE'
            write(urows,'(*(g0,:,","))')trim(release_id),trim(c%id),'formal',axis,trim(r%key),trim(r%route), &
               r%m,r%rank,r%nx,r%ny,r%EA,r%ER,r%cpu,r%low,r%high,r%total,r%total_low,r%total_high, &
               r%valid,r%active,r%query_valid,r%total_valid,trim(state),correctness(a)%candidate_valid, &
               correctness(a)%value_valid,correctness(a)%numerical_status,correctness(a)%finite, &
               trim(correctness(a)%resource_state),trim(correctness(a)%resource_reason),correctness(a)%status, &
               trim(correctness(a)%state)
         end do
      end do
   end subroutine
end program
