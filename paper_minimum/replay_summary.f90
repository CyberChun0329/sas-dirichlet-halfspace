! Standalone, Fortran-only CSV analysis. Usage:
! replay_summary comparisons.csv budget_frontiers.csv NEW_report.md
program replay_summary
   use, intrinsic :: iso_fortran_env, only: real64,iostat_end,error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer,parameter :: rk=real64,maxfields=64
   type text_t
      character(:),allocatable :: value
   end type
   type tally_t
      character(256) :: label=''
      integer :: count=0
   end type
   type group_t
      character(64) :: family='',cost='',mode='',direction=''
      real(rk) :: cap=0
      integer :: rows=0,unmatched=0,unreached=0,floor_pairs=0,resolved_pairs=0,mixed_pairs=0
      integer :: sas_wins=0,space_wins=0,resolved_sas=0,resolved_space=0,ties=0,separated=0,veto=0
      integer :: space_dominates_rows=0,sas_dominates_rows=0
      integer :: positive_width=0,zero_width=0,sas10=0,sas100=0,space10=0,space100=0,finite_endpoints=0
      real(rk) :: ratio_min=huge(0.0_rk),ratio_max=-huge(0.0_rk)
      real(rk) :: sas_speed_min=huge(0.0_rk),sas_speed_max=-huge(0.0_rk)
      real(rk) :: space_speed_min=huge(0.0_rk),space_speed_max=-huge(0.0_rk)
      real(rk) :: start_min=huge(0.0_rk),end_max=-huge(0.0_rk)
      type(text_t),allocatable :: cases(:)
      type(tally_t),allocatable :: statuses(:),gains(:)
   end type
   type(group_t),allocatable :: comparisons(:),budgets(:)
   type(text_t),allocatable :: comparison_cases(:),budget_cases(:)
   character(4096) :: comparison_path,budget_path,report_path
   integer :: report,ios,i,comparison_rows,budget_rows,absolute_rows
   if(command_argument_count()/=3) call fail('usage: replay_summary comparisons.csv budget_frontiers.csv NEW_report.md')
   call get_command_argument(1,comparison_path);call get_command_argument(2,budget_path)
   call get_command_argument(3,report_path)
   allocate(comparisons(0),budgets(0),comparison_cases(0),budget_cases(0))
   comparison_rows=0;budget_rows=0;absolute_rows=0
   ! Validate and summarize both inputs before creating a report.
   call read_comparisons(trim(comparison_path))
   call read_budgets(trim(budget_path))
   open(newunit=report,file=trim(report_path),status='new',action='write',iostat=ios)
   if(ios/=0) call fail('report must be a new writable path: '//trim(report_path))
   write(report,'(a)') '# Fortran replay summary'
   write(report,'(a)') ''
   write(report,'(a)') 'Inputs: '//trim(comparison_path)//' and '//trim(budget_path)
   write(report,'(a,i0)') 'Comparison CSV rows: ',comparison_rows
   write(report,'(a,i0)') 'Distinct comparison case IDs: ',size(comparison_cases)
   write(report,'(a,i0)') 'Absolute-only evidence rows: ',absolute_rows
   write(report,'(a,i0)') 'Budget CSV rows: ',budget_rows
   write(report,'(a,i0)') 'Distinct budget case IDs: ',size(budget_cases)
   write(report,'(a)') ''
   write(report,'(a)') 'Counts describe exact saved rows and distinct case IDs within each group. '// &
      'Cases recur across families, modes, caps and budget directions; group counts must not be pooled as independent cases.'
   write(report,'(a)') 'No population confidence intervals are inferred. Speed ranges describe recorded cost ratios. '// &
      'A cost winner can remain at the reporting floor; it is not a resolved precision claim.'
   do i=1,size(comparisons)
      call write_comparison_group(report,comparisons(i))
   end do
   write(report,'(a)') ''
   write(report,'(a)') 'Budget gains count only intervals with end > start: 10-to-100 means >=10 and <100; '// &
      '100 means >=100. Finite terminal endpoints and zero-width rows are retained separately.'
   do i=1,size(budgets)
      call write_budget_group(report,budgets(i))
   end do
   close(report,iostat=ios)
   if(ios/=0) call fail('report close failed')
   print '(a,i0,a,i0)','SUMMARY_COMPLETE comparisons=',comparison_rows,' budget=',budget_rows
contains
   subroutine fail(message)
      character(*),intent(in) :: message
      write(error_unit,'(a)') trim(message)
      error stop 1
   end subroutine

   ! Preserve empty fields, including trailing empties. Quoted commas/doubled quotes are accepted.
   subroutine split_csv(line,fields,n)
      character(*),intent(in) :: line
      character(4096),intent(out) :: fields(maxfields)
      integer,intent(out) :: n
      integer :: p,q,k,length
      logical :: quoted,closed
      fields='';n=1;k=0;quoted=.false.;closed=.false.;length=len_trim(line);p=1
      do while(p<=length)
         if(line(p:p)=='"') then
            if(quoted) then
               if(p<length) then
                  if(line(p+1:p+1)=='"') then
                     k=k+1;fields(n)(k:k)='"';p=p+2;cycle
                  end if
               end if
               quoted=.false.;closed=.true.
            else
               if(k/=0.or.closed) call fail('malformed CSV quote')
               quoted=.true.
            end if
         else if(line(p:p)==','.and..not.quoted) then
            n=n+1
            if(n>maxfields) call fail('too many CSV fields')
            k=0;closed=.false.
         else
            if(closed.and.line(p:p)/=' ') call fail('characters after closing CSV quote')
            if(.not.closed) then
               k=k+1
               if(k>len(fields(n))) call fail('CSV field too long')
               fields(n)(k:k)=line(p:p)
            end if
         end if
         p=p+1
      end do
      if(quoted) call fail('unclosed CSV quote')
      do q=1,n
         fields(q)=adjustl(fields(q))
      end do
   end subroutine
   integer function column(header,n,name) result(index)
      character(*),intent(in) :: header(:),name
      integer,intent(in) :: n
      integer :: j
      index=0
      do j=1,n
         if(trim(header(j))==name) then
            if(index/=0) call fail('duplicate CSV header: '//name)
            index=j
         end if
      end do
      if(index==0) call fail('missing CSV header: '//name)
   end function
   real(rk) function number(text,required) result(v)
      character(*),intent(in) :: text
      logical,optional,intent(in) :: required
      integer :: st
      v=0
      if(len_trim(text)==0) then
         if(present(required)) then
            if(required) call fail('empty required numeric field')
         end if
         return
      end if
      read(text,*,iostat=st) v
      if(st/=0) call fail('invalid numeric CSV field: '//trim(text))
      if(.not.ieee_is_finite(v)) call fail('nonfinite CSV field: '//trim(text))
   end function
   logical function flag(text) result(v)
      character(*),intent(in) :: text
      select case(trim(text))
      case('T','TRUE','true');v=.true.
      case('F','FALSE','false');v=.false.
      case default;call fail('invalid logical CSV field: '//trim(text))
      end select
   end function
   function cost_label(text) result(label)
      character(*),intent(in) :: text
      character(64) :: label
      select case(trim(text))
      case('cpu','QUERY');label='QUERY'
      case('total','TOTAL');label='TOTAL'
      case('none');label='NONE'
      case default;call fail('unknown cost: '//trim(text))
      end select
   end function
   subroutine add_case(cases,id)
      type(text_t),allocatable,intent(inout) :: cases(:)
      character(*),intent(in) :: id
      type(text_t) :: entry
      integer :: j
      if(len_trim(id)==0) call fail('empty case ID')
      do j=1,size(cases)
         if(cases(j)%value==trim(id)) return
      end do
      entry%value=trim(id);cases=[cases,entry]
   end subroutine
   subroutine tally(list,label)
      type(tally_t),allocatable,intent(inout) :: list(:)
      character(*),intent(in) :: label
      type(tally_t) :: value
      integer :: j
      if(len_trim(label)==0) call fail('empty status/classification')
      do j=1,size(list)
         if(list(j)%label==label) then
            list(j)%count=list(j)%count+1;return
         end if
      end do
      value%label=label;value%count=1;list=[list,value]
   end subroutine
   integer function find_group(groups,family,cost,mode,cap,direction) result(g)
      type(group_t),allocatable,intent(inout) :: groups(:)
      character(*),intent(in) :: family,cost,mode,direction
      real(rk),intent(in) :: cap
      type(group_t) :: value
      integer :: j
      do j=1,size(groups)
         if(groups(j)%family==family.and.groups(j)%cost==cost.and.groups(j)%mode==mode &
            .and.groups(j)%cap==cap.and.groups(j)%direction==direction) then
            g=j;return
         end if
      end do
      value%family=family;value%cost=cost;value%mode=mode;value%cap=cap;value%direction=direction
      allocate(value%cases(0),value%statuses(0),value%gains(0));groups=[groups,value];g=size(groups)
   end function
   subroutine validate_numbers(header,fields,n)
      character(*),intent(in) :: header(:),fields(:)
      integer,intent(in) :: n
      real(rk) :: v
      integer :: j
      do j=1,n
         select case(trim(header(j)))
         case('rho','Fo','relative_floor','absolute_floor','cap','SAS_EA','SAS_ER','SPACE_EA','SPACE_ER', &
            'SPACE_over_SAS_cost','start','end','floor_capped_SPACE_over_SAS_ER')
            select case(trim(header(j)))
            case('SAS_EA','SAS_ER','SPACE_EA','SPACE_ER');v=number(fields(j))
            case default;v=number(fields(j),.true.)
            end select
            if(v<0) call fail('negative numeric CSV field: '//trim(header(j)))
         end select
      end do
   end subroutine
   subroutine next_row(u,line,fields,n,expected,done)
      integer,intent(in) :: u,expected
      character(*),intent(out) :: line
      character(4096),intent(out) :: fields(maxfields)
      integer,intent(out) :: n
      logical,intent(out) :: done
      integer :: st
      read(u,'(a)',iostat=st) line
      done=st==iostat_end
      if(done) return
      if(st/=0) call fail('CSV read failure')
      if(len_trim(line)==len(line)) call fail('CSV record exceeds buffer')
      if(len_trim(line)==0) call fail('empty CSV row')
      call split_csv(trim(line),fields,n)
      if(n/=expected) call fail('CSV column count mismatch')
   end subroutine
   subroutine open_csv(path,u,header,n)
      character(*),intent(in) :: path
      integer,intent(out) :: u,n
      character(4096),intent(out) :: header(maxfields)
      character(32768) :: line
      integer :: st
      open(newunit=u,file=path,status='old',action='read',iostat=st)
      if(st/=0) call fail('cannot open CSV: '//path)
      read(u,'(a)',iostat=st) line
      if(st/=0) call fail('missing CSV header: '//path)
      call split_csv(trim(line),header,n)
   end subroutine

   subroutine read_comparisons(path)
      character(*),intent(in) :: path
      character(4096) :: h(maxfields),f(maxfields)
      character(32768) :: line
      integer :: u,n,nrow,g,ic,ifa,ico,imo,ica,ist,iw,isep,ires,isk,ikk,ira,ikd,isd
      real(rk) :: cap,ratio
      logical :: done,separated,veto,resolved,pair,kd,sd
      call open_csv(path,u,h,n)
      ic=column(h,n,'case');ifa=column(h,n,'family');ico=column(h,n,'cost');imo=column(h,n,'mode')
      ica=column(h,n,'cap');ist=column(h,n,'status');iw=column(h,n,'winner');isep=column(h,n,'separated')
      ires=column(h,n,'resolution');isk=column(h,n,'SAS_key');ikk=column(h,n,'SPACE_key')
      ira=column(h,n,'SPACE_over_SAS_cost');ikd=column(h,n,'SPACE_dominates_SAS');isd=column(h,n,'SAS_dominates_SPACE')
      do
         call next_row(u,line,f,nrow,n,done)
         if(done) exit
         call validate_numbers(h,f,n);comparison_rows=comparison_rows+1
         call add_case(comparison_cases,trim(f(ic)))
         cap=number(f(ica),.true.);ratio=number(f(ira),.true.)
         separated=flag(f(isep));kd=flag(f(ikd));sd=flag(f(isd))
         veto=(kd.or.sd).and.f(iw)=='DOMINATED_OR_TIED'
         g=find_group(comparisons,trim(f(ifa)),trim(cost_label(f(ico))),trim(f(imo)),cap,'')
         comparisons(g)%rows=comparisons(g)%rows+1
         call add_case(comparisons(g)%cases,trim(f(ic)));call tally(comparisons(g)%statuses,trim(f(ist)))
         if(f(imo)=='ABSOLUTE') absolute_rows=absolute_rows+1
         if(f(ist)=='NO_ACTUAL_ERROR_MATCH') comparisons(g)%unmatched=comparisons(g)%unmatched+1
         if(f(ist)=='TARGET_NOT_BOTH_REACHED') comparisons(g)%unreached=comparisons(g)%unreached+1
         pair=len_trim(f(isk))>0.and.len_trim(f(ikk))>0
         resolved=pair.and.f(ires)=='BOTH_ABOVE_RELATIVE_FLOOR'
         if(pair) then
            if(ratio<=0) call fail('nonpositive selected-pair cost ratio')
            comparisons(g)%ratio_min=min(comparisons(g)%ratio_min,ratio)
            comparisons(g)%ratio_max=max(comparisons(g)%ratio_max,ratio)
            if(f(ires)=='BOTH_RELATIVE_FLOOR') comparisons(g)%floor_pairs=comparisons(g)%floor_pairs+1
            if(resolved) comparisons(g)%resolved_pairs=comparisons(g)%resolved_pairs+1
            if(f(ires)=='ONE_AT_RELATIVE_FLOOR') comparisons(g)%mixed_pairs=comparisons(g)%mixed_pairs+1
         end if
         if(separated) comparisons(g)%separated=comparisons(g)%separated+1
         if(separated.and.f(iw)/='SAS'.and.f(iw)/='SPACE') call fail('separated flag without a cost winner')
         if(veto) comparisons(g)%veto=comparisons(g)%veto+1
         if(kd) comparisons(g)%space_dominates_rows=comparisons(g)%space_dominates_rows+1
         if(sd) comparisons(g)%sas_dominates_rows=comparisons(g)%sas_dominates_rows+1
         select case(trim(f(iw)))
         case('SAS')
            if(.not.pair) call fail('winner without selected pair')
            comparisons(g)%sas_wins=comparisons(g)%sas_wins+1
            if(resolved) comparisons(g)%resolved_sas=comparisons(g)%resolved_sas+1
            comparisons(g)%sas_speed_min=min(comparisons(g)%sas_speed_min,ratio)
            comparisons(g)%sas_speed_max=max(comparisons(g)%sas_speed_max,ratio)
         case('SPACE')
            if(.not.pair) call fail('winner without selected pair')
            comparisons(g)%space_wins=comparisons(g)%space_wins+1
            if(resolved) comparisons(g)%resolved_space=comparisons(g)%resolved_space+1
            comparisons(g)%space_speed_min=min(comparisons(g)%space_speed_min,1/ratio)
            comparisons(g)%space_speed_max=max(comparisons(g)%space_speed_max,1/ratio)
         case('TIED','DOMINATED_OR_TIED')
            comparisons(g)%ties=comparisons(g)%ties+1
         case('NONE')
         case default
            call fail('unknown comparison winner: '//trim(f(iw)))
         end select
      end do
      close(u)
   end subroutine
   subroutine read_budgets(path)
      character(*),intent(in) :: path
      character(4096) :: h(maxfields),f(maxfields)
      character(32768) :: line
      integer :: u,n,nrow,g,ic,ifa,ico,idi,ista,iend,icl,iga,ife,izero,ires
      real(rk) :: start,finish
      logical :: done,zero,endpoint,resolved
      call open_csv(path,u,h,n)
      ic=column(h,n,'case');ifa=column(h,n,'family');ico=column(h,n,'cost');idi=column(h,n,'direction')
      ista=column(h,n,'start');iend=column(h,n,'end');icl=column(h,n,'classification');iga=column(h,n,'meaningful_gain')
      ife=column(h,n,'finite_bank_right_endpoint');izero=column(h,n,'zero_width');ires=column(h,n,'both_errors_resolved')
      do
         call next_row(u,line,f,nrow,n,done)
         if(done) exit
         call validate_numbers(h,f,n);budget_rows=budget_rows+1;call add_case(budget_cases,trim(f(ic)))
         start=number(f(ista),.true.);finish=number(f(iend),.true.);zero=flag(f(izero))
         endpoint=flag(f(ife));resolved=flag(f(ires))
         if(finish<start.or.(zero.neqv.(finish==start))) call fail('inconsistent budget interval/zero-width flag')
         select case(trim(f(idi)))
         case('nominal','SAS_conservative','SPACE_conservative')
         case default;call fail('unknown budget direction')
         end select
         select case(trim(f(iga)))
         case('NONE','SAS_10_to_100','SAS_100','SPACE_10_to_100','SPACE_100')
         case default;call fail('unknown meaningful gain (only >=10 tiers accepted)')
         end select
         g=find_group(budgets,trim(f(ifa)),trim(cost_label(f(ico))),'BUDGET',0.0_rk,trim(f(idi)))
         budgets(g)%rows=budgets(g)%rows+1;call add_case(budgets(g)%cases,trim(f(ic)))
         call tally(budgets(g)%statuses,trim(f(icl)));call tally(budgets(g)%gains,trim(f(iga)))
         budgets(g)%start_min=min(budgets(g)%start_min,start);budgets(g)%end_max=max(budgets(g)%end_max,finish)
         if(endpoint) budgets(g)%finite_endpoints=budgets(g)%finite_endpoints+1
         if(zero) then
            budgets(g)%zero_width=budgets(g)%zero_width+1
         else
            budgets(g)%positive_width=budgets(g)%positive_width+1
            select case(trim(f(iga)))
            case('SAS_10_to_100');budgets(g)%sas10=budgets(g)%sas10+1
            case('SAS_100');budgets(g)%sas100=budgets(g)%sas100+1
            case('SPACE_10_to_100');budgets(g)%space10=budgets(g)%space10+1
            case('SPACE_100');budgets(g)%space100=budgets(g)%space100+1
            end select
         end if
      end do
      close(u)
   end subroutine
   subroutine write_tallies(u,label,list)
      integer,intent(in) :: u
      character(*),intent(in) :: label
      type(tally_t),intent(in) :: list(:)
      integer :: j
      do j=1,size(list)
         write(u,'(a,a,a,i0)') label,trim(list(j)%label),': ',list(j)%count
      end do
   end subroutine
   subroutine write_range(u,label,lo,hi)
      integer,intent(in) :: u
      character(*),intent(in) :: label
      real(rk),intent(in) :: lo,hi
      if(lo>hi) then
         write(u,'(a)') label//': unavailable'
      else
         write(u,'(a,es15.7,a,es15.7)') label//': ',lo,' to ',hi
      end if
   end subroutine
   subroutine write_comparison_group(u,g)
      integer,intent(in) :: u
      type(group_t),intent(in) :: g
      write(u,'(a)') ''
      write(u,'(a,a,1x,a,1x,a,a,es15.7)') '## ',trim(g%family),trim(g%cost),trim(g%mode),' cap=',g%cap
      write(u,'(a,i0,a,i0)') 'rows=',g%rows,' distinct_case_ids=',size(g%cases)
      write(u,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'unmatched=',g%unmatched,' unreached=',g%unreached, &
         ' floor_pairs=',g%floor_pairs,' resolved_pairs=',g%resolved_pairs,' one_floor_pairs=',g%mixed_pairs
      write(u,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'SAS_cost_wins=',g%sas_wins,' SPACE_cost_wins=',g%space_wins, &
         ' resolved_SAS_wins=',g%resolved_sas,' resolved_SPACE_wins=',g%resolved_space,' ties_or_vetoed=',g%ties
      write(u,'(a,i0,a,i0)') 'separated_winners=',g%separated,' rows_with_dominance_veto=',g%veto
      write(u,'(a,i0,a,i0)') 'SPACE_dominates_SAS_rows=',g%space_dominates_rows, &
         ' SAS_dominates_SPACE_rows=',g%sas_dominates_rows
      call write_range(u,'SPACE/SAS cost ratio (selected pairs)',g%ratio_min,g%ratio_max)
      call write_range(u,'SAS winner speed ratio',g%sas_speed_min,g%sas_speed_max)
      call write_range(u,'SPACE winner speed ratio',g%space_speed_min,g%space_speed_max)
      call write_tallies(u,'status ',g%statuses)
   end subroutine
   subroutine write_budget_group(u,g)
      integer,intent(in) :: u
      type(group_t),intent(in) :: g
      write(u,'(a)') ''
      write(u,'(a)') '## '//trim(g%family)//' '//trim(g%cost)//' '//trim(g%direction)//' budget'
      write(u,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'rows=',g%rows,' distinct_case_ids=',size(g%cases), &
         ' positive_width_intervals=',g%positive_width,' zero_width_rows=',g%zero_width,' finite_endpoints=',g%finite_endpoints
      write(u,'(a,i0,a,i0,a,i0,a,i0)') 'positive_width_SAS_10_to_100=',g%sas10,' SAS_100=',g%sas100, &
         ' SPACE_10_to_100=',g%space10,' SPACE_100=',g%space100
      call write_range(u,'recorded finite bank budget span',g%start_min,g%end_max)
      call write_tallies(u,'classification ',g%statuses);call write_tallies(u,'all-width gain label ',g%gains)
   end subroutine
end program
