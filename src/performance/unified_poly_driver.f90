program unified_poly_driver
   use replay_types
   use replay_case_bank, only: case_spec_t,build_cases
   use unified_policy, only: policy_t,load_policy
   use unified_poly_reference
   implicit none
   type(policy_t) :: policy
   type(case_spec_t),allocatable :: cases(:)
   type(case_spec_t) :: c
   type(reference_t) :: reference
   integer :: first,last,repeats,ci,st,i,j,ios,uref,urefcsv
   character(1024) :: outdir,policy_path,arg,reference_path
   character(128) :: release_id
   character(64) :: family,mode
   logical :: audit
   call get_command_argument(1,family);call get_command_argument(2,outdir)
   call get_command_argument(3,policy_path);call get_command_argument(4,release_id)
   call get_command_argument(5,mode)
   if(len_trim(mode)==0) mode='reference'
   if(min(len_trim(family),len_trim(outdir),len_trim(policy_path),len_trim(release_id))==0) &
      error stop 'family outputdir policy_path release_id reference [first last] [audit] [repeatN]'
   if(mode/='reference') error stop 'only reference production is supported'
   call load_policy(trim(policy_path),policy,st);call must(st,'policy')
   audit=trim(family)=='paper_audit'
   if(audit)family='paper'
   call build_cases(trim(family),cases,st);call must(st,'case bank')
   first=1;last=size(cases);repeats=1
   call get_command_argument(6,arg)
   if(len_trim(arg)>0) then
      read(arg,*,iostat=ios)first
      if(ios/=0) error stop 'invalid first'
   end if
   call get_command_argument(7,arg)
   if(len_trim(arg)>0) then
      read(arg,*,iostat=ios)last
      if(ios/=0) error stop 'invalid last'
   end if
   do i=8,command_argument_count()
      call get_command_argument(i,arg)
      if(trim(arg)=='audit') then
         audit=.true.
      else if(index(arg,'repeat')==1) then
         read(arg(7:),*,iostat=ios)repeats
         if(ios/=0) error stop 'invalid repeat count'
      else
         error stop 'unknown optional argument'
      end if
   end do
   if(first<1.or.last>size(cases).or.first>last) error stop 'invalid case range'
   if(.not.any(repeats==[1,10,100])) error stop 'repeat count must be 1, 10 or 100'
   reference_path=trim(outdir)//'/references.dat'
   if(mode=='reference') then
      open(newunit=uref,file=trim(reference_path),status='new',action='write',iostat=ios)
      if(ios/=0) error stop 'fresh reference output required'
      write(uref,'(a)')'UNIFIED_POLY_REFERENCE_V1'
      write(uref,'(a)')trim(release_id);write(uref,'(a)')trim(family)
      write(uref,*)first,last,repeats,audit
      call open_new(urefcsv,'references.csv')
      write(urefcsv,'(a)')'release,case,repeat,point,x,y,z,t,observable,ref48,ref72,ref96,shared96,reference,order,EA,ER,crossEA,crossER,FR,FA,certified'
      do ci=first,last
         if(audit.and..not.cases(ci)%audit) cycle
         call current_case()
         call fresh_reference(c,reference,st);call must(st,'fresh spatial reference')
         call write_reference(uref,c,reference)
         do j=1,size(c%q,2)
            write(urefcsv,'(*(g0,:,","))')trim(release_id),trim(c%id),repeats,j,c%q(:,j),c%observable, &
               reference%v48(j),reference%v72(j),reference%v96(j),reference%shared(j),reference%values(j), &
               reference%order,reference%EA,reference%ER,reference%crossEA,reference%crossER, &
               reference%floor,reference%afloor,reference%certified
         end do
         flush(uref);flush(urefcsv)
      end do
      close(uref);close(urefcsv);stop
   end if
contains
   subroutine must(code,label)
      integer,intent(in) :: code
      character(*),intent(in) :: label
      if(code/=ok) then
         print *,trim(label),code;error stop 'invalid unified operation'
      end if
   end subroutine
   subroutine open_new(unit,name)
      integer,intent(out) :: unit
      character(*),intent(in) :: name
      open(newunit=unit,file=trim(outdir)//'/'//trim(name),status='new',action='write',iostat=ios)
      if(ios/=0) error stop 'output exists or directory unavailable'
   end subroutine
   subroutine current_case()
      real(rk),allocatable :: batch(:,:)
      integer :: n,t
      c=cases(ci);n=size(c%q,2)
      if(repeats>1) then
         allocate(batch(4,n*repeats))
         do t=1,repeats
            batch(:,(t-1)*n+1:t*n)=c%q
         end do
         c%q=batch
      end if
   end subroutine
end program
