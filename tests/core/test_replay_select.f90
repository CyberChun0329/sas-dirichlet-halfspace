program test_replay_select
   use replay_types, only: rk,ok,invalid,quiet_nan
   use replay_select
   implicit none
   type(candidate_t) :: r(4),reordered(4)
   type(comparison_t) :: result,other
   type(budget_state_t) :: point
   type(budget_state_t),allocatable :: states(:)
   logical :: ss(4),kk(4),kept(4)
   integer :: metric,direction,sb,kb,n,i,status,unit,ios,query_only_matches,total_only_matches
   character(1024) :: outprefix
   character(4096) :: csvline
   ss=[.true.,.true.,.false.,.false.];kk=.not.ss
   call require(route_family('SPACE_RECT_ERFC')=='DIRECT_2D','rect family')
   call require(route_family('SPACE_SHARED_SOURCE')=='SHARED_1D','source shared family')
   call require(route_family('SPACE_SHARED_PROJECTED')=='SHARED_1D','projected shared family')
   call require(route_family('MAP_TRUTH')=='UNSUPPORTED','unsupported truth')
   ! Same selected keys, but at B=2 a newly eligible competitor vetoes the 10x gain.
   r(1)=row('s','SAS_LOW',1.0_rk,1e-4_rk,1e-5_rk)
   r(2)=row('d','SAS_LOW',3.0_rk,2e-2_rk,1e-4_rk)
   r(3)=row('k','SPACE_RECT_ERFC',1.0_rk,2e-3_rk,2e-5_rk)
   r(4)=row('q','SPACE_RECT_ERFC',2.0_rk,1e-2_rk,1e-6_rk)
   do metric=1,2
      do direction=1,3
         sb=1;kb=1
         if(direction==2) then
            sb=3;kb=2
         else if(direction==3) then
            sb=2;kb=3
         end if
         call budget_states(r,ss,kk,metric,sb,kb,1e-12_rk,1e-15_rk,states,n)
         call require(n==2,'event compression count')
         call require(states(1)%sas==1.and.states(1)%space==3,'first keys')
         call require(states(1)%gain=='SAS_10_to_100','10x gain tier')
         call require(states(1)%start==1.and.states(1)%finish==2,'EA veto boundary retained')
         call require(states(2)%sas==1.and.states(2)%space==3,'unchanged ER keys')
         call require(states(2)%gain=='NONE'.and.states(2)%finish==3,'finite ending')
         do i=1,3
            call budget_point(r,ss,kk,real(i,rk),metric,sb,kb,1e-12_rk,1e-15_rk,point)
            if(i==1) then
               call require(point%gain==states(1)%gain,'first event readback')
            else
               call require(point%gain==states(2)%gain,'remaining event readback')
            end if
         end do
      end do
   end do
   r(2)%valid=.false.;r(4)%valid=.false.
   call budget_states(r,ss,kk,1,1,1,1e-12_rk,1e-15_rk,states,n)
   call require(n==1.and.states(1)%start==states(1)%finish,'terminal zero width')
   call require(states(1)%gain=='NONE','zero width is not gain')
   r(2)%valid=.true.;r(4)%valid=.true.
   ! Conservative eligibility must use actual low/high bounds, and total is independent.
   r(1)%low=.5_rk;r(1)%high=1.5_rk;r(3)%low=.75_rk;r(3)%high=1.25_rk
   call budget_point(r,ss,kk,1.0_rk,1,3,2,1e-12_rk,1e-15_rk,point)
   call require(point%sas==0.and.point%space==3,'SAS conservative bounds')
   call budget_point(r,ss,kk,1.0_rk,1,2,3,1e-12_rk,1e-15_rk,point)
   call require(point%sas==1.and.point%space==0,'SPACE conservative bounds')
   r(1)%total=8;r(1)%total_low=7;r(1)%total_high=9
   call budget_point(r,ss,kk,1.0_rk,2,1,1,1e-12_rk,1e-15_rk,point)
   call require(point%sas==0.and.point%space==3,'direct total independent of query cost')
   ! Joint cap rejects the cheapest EA-ineligible row. Ordering does not pick winners.
   r(1)=row('s','SAS_LOW',1.0_rk,1e-5_rk,2e-4_rk)
   r(2)=row('t','SAS_LOW',2.0_rk,1e-5_rk,1e-5_rk)
   r(3)=row('k','SPACE_RECT_ERFC',3.0_rk,2e-5_rk,2e-5_rk)
   r(4)%valid=.false.
   call target_pair(r,ss,kk,1e-4_rk,1,1e-12_rk,result)
   call require(result%sas==2.and.result%space==3,'joint cap target')
   call matched_pair(r,ss,kk,1e-4_rk,1,1e-12_rk,1e-15_rk,result)
   call require(result%sas==2.and.result%winner=='SAS'.and.result%separated,'matched joint cap')
   reordered=r([2,1,4,3])
   call matched_pair(reordered,ss,kk,1e-4_rk,1,1e-12_rk,1e-15_rk,other)
   call require(reordered(other%sas)%key==r(result%sas)%key,'order independent SAS')
   call require(reordered(other%space)%key==r(result%space)%key,'order independent SPACE')
   ! A failed independent cost must not suppress the other valid cost.
   r(1)=row('s','SAS_LOW',1.0_rk,1e-5_rk,1e-5_rk);r(2)%valid=.false.
   r(3)=row('k','SPACE_RECT_ERFC',2.0_rk,2e-5_rk,2e-5_rk);r(4)%valid=.false.
   r(3)%total_valid=.false.;r(3)%total=quiet_nan();r(3)%total_low=quiet_nan();r(3)%total_high=quiet_nan()
   call matched_pair(r,ss,kk,1e-4_rk,1,1e-12_rk,1e-15_rk,result)
   call require(result%sas==1.and.result%space==3,'QUERY survives failed TOTAL')
   call matched_pair(r,ss,kk,1e-4_rk,2,1e-12_rk,1e-15_rk,result)
   call require(result%sas==0.and.result%space==0,'failed TOTAL excluded from match')
   call target_pair(r,ss,kk,1e-4_rk,1,1e-12_rk,result)
   call require(result%space==3,'QUERY target survives failed TOTAL')
   call target_pair(r,ss,kk,1e-4_rk,2,1e-12_rk,result)
   call require(result%space==0,'failed TOTAL excluded from target')
   call pareto_mask(r,kk,1,1e-12_rk,1e-15_rk,kept)
   call require(kept(3),'QUERY Pareto survives failed TOTAL')
   call budget_states(r,ss,kk,1,1,1,1e-12_rk,1e-15_rk,states,n)
   call require(states(n)%space==3,'QUERY cost event survives failed TOTAL')
   call budget_point(r,ss,kk,3.0_rk,2,1,1,1e-12_rk,1e-15_rk,point)
   call require(point%space==0,'failed TOTAL excluded from budget')
   r(3)=row('k','SPACE_RECT_ERFC',2.0_rk,2e-5_rk,2e-5_rk)
   r(3)%query_valid=.false.;r(3)%cpu=quiet_nan();r(3)%low=quiet_nan();r(3)%high=quiet_nan()
   call matched_pair(r,ss,kk,1e-4_rk,2,1e-12_rk,1e-15_rk,result)
   call require(result%sas==1.and.result%space==3,'TOTAL survives failed QUERY')
   call target_pair(r,ss,kk,1e-4_rk,1,1e-12_rk,result)
   call require(result%space==0,'failed QUERY excluded from target')
   call pareto_mask(r,kk,2,1e-12_rk,1e-15_rk,kept)
   call require(kept(3),'TOTAL Pareto survives failed QUERY')
   call budget_states(r,ss,kk,2,1,1,1e-12_rk,1e-15_rk,states,n)
   call require(states(n)%space==3,'TOTAL cost event survives failed QUERY')
   call budget_point(r,ss,kk,3.0_rk,1,1,1,1e-12_rk,1e-15_rk,point)
   call require(point%space==0,'failed QUERY excluded from budget')
   ! A SPACE point outside the matched error band still vetoes a claimed SAS win.
   r(1)=row('s','SAS_LOW',1.0_rk,1e-4_rk,1e-5_rk)
   r(2)%valid=.false.
   r(3)=row('k','SPACE_RECT_ERFC',2.0_rk,1e-4_rk,1e-5_rk)
   r(4)=row('q','SPACE_RECT_ERFC',.5_rk,1e-7_rk,1e-5_rk)
   call matched_pair(r,ss,kk,1e-3_rk,1,1e-12_rk,1e-15_rk,result)
   ! Pareto removes k, so no comparable pair remains; add EA tradeoff to retain it.
   r(3)%EA=1e-6_rk
   call matched_pair(r,ss,kk,1e-3_rk,1,1e-12_rk,1e-15_rk,result)
   call require(result%sas==1.and.result%space==3,'dominance witness pair')
   call require(result%space_dominates.and.result%winner=='DOMINATED_OR_TIED','opposite dominance veto')
   ! Floor-equivalent duplicates use key ties; both-floor keeps cost but no precision gain.
   r(1)=row('z','SAS_LOW',1.0_rk,1e-13_rk,1e-16_rk)
   r(2)=row('a','SAS_LOW',1.0_rk,5e-13_rk,5e-16_rk)
   r(3)=row('k','SPACE_RECT_ERFC',2.0_rk,3e-13_rk,3e-16_rk)
   r(4)%valid=.false.
   call pareto_mask(r,ss,1,1e-12_rk,1e-15_rk,kept)
   call require(.not.kept(1).and.kept(2),'floored Pareto key tie')
   call matched_pair(r,ss,kk,1e-10_rk,1,1e-12_rk,1e-15_rk,result)
   call require(result%status=='REPORTING_FLOOR_MATCH'.and.result%winner=='SAS','floor cost winner')
   call budget_point(r,ss,kk,2.0_rk,1,1,1,1e-12_rk,1e-15_rk,point)
   call require(point%classification=='BOTH_RELATIVE_FLOOR'.and.point%gain=='NONE','floor suppresses precision')
   r(1)%active=.false.;r(2)%active=.false.;r(3)%active=.false.
   call target_pair(r,ss,kk,1e-4_rk,1,1e-12_rk,result)
   call require(result%sas==0.and.result%space==0,'weak reference no ranking')
   ! Optional caller-owned temporary directory exercises CSV new/append behavior.
   call get_command_argument(1,outprefix)
   if(len_trim(outprefix)>0) then
      call analyze_candidates('fixture','U',1.0_rk,.1_rk,r,1e-12_rk,1e-15_rk,trim(outprefix),status)
      call require(status==ok,'new CSV files')
      call analyze_candidates('fixture','U',1.0_rk,.1_rk,r,1e-12_rk,1e-15_rk,trim(outprefix),status)
      call require(status==invalid,'existing CSV protected')
      call analyze_candidates('append','U',1.0_rk,.1_rk,r,1e-12_rk,1e-15_rk,trim(outprefix),status,append=.true.)
      call require(status==ok,'append CSV files')
      r(1)=row('s','SAS_LOW',1.0_rk,1e-5_rk,1e-5_rk);r(2)%valid=.false.
      r(3)=row('k','SPACE_RECT_ERFC',2.0_rk,2e-5_rk,2e-5_rk);r(4)%valid=.false.
      r(3)%total_valid=.false.;r(3)%total=quiet_nan();r(3)%total_low=quiet_nan();r(3)%total_high=quiet_nan()
      call analyze_candidates('query_only','U',1.0_rk,.1_rk,r,1e-12_rk,1e-15_rk,trim(outprefix),status, &
         caps=[1e-4_rk],append=.true.)
      call require(status==ok,'QUERY-only analysis')
      r(3)=row('k','SPACE_RECT_ERFC',2.0_rk,2e-5_rk,2e-5_rk)
      r(3)%query_valid=.false.;r(3)%cpu=quiet_nan();r(3)%low=quiet_nan();r(3)%high=quiet_nan()
      call analyze_candidates('total_only','U',1.0_rk,.1_rk,r,1e-12_rk,1e-15_rk,trim(outprefix),status, &
         caps=[1e-4_rk],append=.true.)
      call require(status==ok,'TOTAL-only analysis')
      query_only_matches=0;total_only_matches=0
      open(newunit=unit,file=trim(outprefix)//'/comparisons.csv',status='old',action='read')
      do
         read(unit,'(a)',iostat=ios) csvline
         if(ios/=0) exit
         if(index(csvline,'query_only,')==1.and.index(csvline,'DIRECT_2D,cpu,ACTUAL_ERROR,')>0 &
            .and.index(csvline,',ACTUAL_ERROR_MATCH,')>0) query_only_matches=query_only_matches+1
         if(index(csvline,'total_only,')==1.and.index(csvline,'DIRECT_2D,total,ACTUAL_ERROR,')>0 &
            .and.index(csvline,',ACTUAL_ERROR_MATCH,')>0) total_only_matches=total_only_matches+1
         if(index(csvline,'query_only,')==1.and.index(csvline,'DIRECT_2D,total,ACTUAL_ERROR,')>0) &
            call require(index(csvline,',NO_ACTUAL_ERROR_MATCH,')>0,'failed TOTAL absent from CSV match')
         if(index(csvline,'total_only,')==1.and.index(csvline,'DIRECT_2D,cpu,ACTUAL_ERROR,')>0) &
            call require(index(csvline,',NO_ACTUAL_ERROR_MATCH,')>0,'failed QUERY absent from CSV match')
      end do
      close(unit)
      call require(query_only_matches==1.and.total_only_matches==1,'analysis metric masks preserve valid counterpart')
   end if
   print '(a)','PASS: joint caps, order, Pareto floors, dominance veto, all cost events and 10x EA boundary'
contains
   function row(key,route,c,er,ea) result(v)
      character(*),intent(in) :: key,route
      real(rk),intent(in) :: c,er,ea
      type(candidate_t) :: v
      v%key=key;v%route=route;v%cpu=c;v%low=c;v%high=c
      ! Fixture same budget events across both costs and all directions.
      v%total=c;v%total_low=c;v%total_high=c
      v%EA=ea;v%ER=er;v%valid=.true.;v%active=.true.;v%query_valid=.true.;v%total_valid=.true.
   end function
   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(*),intent(in) :: label
      if(.not.condition) then
         print '(a)','FAIL: '//label
         error stop 1
      end if
   end subroutine
end program
