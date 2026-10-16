module replay_select
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use replay_types, only: rk,ok,invalid
   implicit none
   private
   public :: candidate_t,comparison_t,budget_state_t,analyze_candidates
   public :: pareto_mask,matched_pair,target_pair,budget_point,budget_states,route_family
   type candidate_t
      character(128) :: key='',route=''
      integer :: n=0,nx=0,ny=0,m=0,rank=0
      real(rk) :: EA=0,ER=0,cpu=0,low=0,high=0,total=0,total_low=0,total_high=0
      logical :: valid=.false.,active=.false.,query_valid=.false.,total_valid=.false.
   end type
   type comparison_t
      integer :: sas=0,space=0
      character(40) :: status='NO_ACTUAL_ERROR_MATCH',winner='NONE',resolution='UNAVAILABLE'
      logical :: separated=.false.,space_dominates=.false.,sas_dominates=.false.
      real(rk) :: ratio=0
   end type
   type budget_state_t
      integer :: sas=0,space=0
      real(rk) :: start=0,finish=0,ratio=0
      character(32) :: classification='UNAVAILABLE',gain='NONE'
      logical :: resolved=.false.
   end type
contains
   pure function route_family(route) result(family)
      character(*),intent(in) :: route
      character(16) :: family
      family='UNSUPPORTED'
      if(index(route,'SAS')==1) then
         family='SAS';return
      end if
      select case(trim(route))
      case('SPACE_POLY','SPACE_RECT_ERFC','SPACE_POLY_FACTOR')
         family='DIRECT_2D'
      case('SPACE_SHARED_SOURCE','SPACE_SHARED_PROJECTED')
         family='SHARED_1D'
      end select
   end function

   pure logical function usable(r,metric)
      type(candidate_t),intent(in) :: r
      integer,intent(in) :: metric
      usable=r%valid.and.r%active
      if(.not.usable) return
      usable=all(ieee_is_finite([r%EA,r%ER]))
      if(.not.usable) return
      usable=min(r%EA,r%ER)>=0
      if(.not.usable) return
      select case(metric)
      case(1)
         usable=r%query_valid.and.all(ieee_is_finite([r%cpu,r%low,r%high]))
         if(usable) usable=min(r%cpu,r%low)>0.and.r%low<=r%cpu.and.r%cpu<=r%high
      case(2)
         usable=r%total_valid.and.all(ieee_is_finite([r%total,r%total_low,r%total_high]))
         if(usable) usable=min(r%total,r%total_low)>0.and.r%total_low<=r%total.and.r%total<=r%total_high
      case default
         usable=.false.
      end select
   end function

   pure real(rk) function cost(r,metric,bound) result(v)
      type(candidate_t),intent(in) :: r
      integer,intent(in) :: metric,bound
      v=huge(1.0_rk)
      if(metric==1) then
         select case(bound)
         case(1);v=r%cpu
         case(2);v=r%low
         case(3);v=r%high
         end select
      else
         select case(bound)
         case(1);v=r%total
         case(2);v=r%total_low
         case(3);v=r%total_high
         end select
      end if
   end function

   pure subroutine pareto_mask(rows,eligible,metric,floor,afloor,kept)
      type(candidate_t),intent(in) :: rows(:)
      logical,intent(in) :: eligible(:)
      integer,intent(in) :: metric
      real(rk),intent(in) :: floor,afloor
      logical,intent(out) :: kept(size(rows))
      integer :: i,j
      real(rk) :: a,e,c,qa,qe,qc
      kept=.false.
      do i=1,size(rows)
         if(.not.eligible(i).or..not.usable(rows(i),metric)) cycle
         a=max(rows(i)%EA,afloor);e=max(rows(i)%ER,floor);c=cost(rows(i),metric,1)
         kept(i)=.true.
         do j=1,size(rows)
            if(i==j.or..not.eligible(j).or..not.usable(rows(j),metric)) cycle
            qa=max(rows(j)%EA,afloor);qe=max(rows(j)%ER,floor);qc=cost(rows(j),metric,1)
            if(qc<=c.and.qa<=a.and.qe<=e) then
               if(qc<c.or.qa<a.or.qe<e.or.rows(j)%key<rows(i)%key) then
                  kept(i)=.false.;exit
               end if
            end if
         end do
      end do
   end subroutine

   pure logical function dominates(a,b,metric,floor,afloor)
      type(candidate_t),intent(in) :: a,b
      integer,intent(in) :: metric
      real(rk),intent(in) :: floor,afloor
      dominates=cost(a,metric,1)<=cost(b,metric,1).and.a%EA<=b%EA+afloor.and.a%ER<=b%ER+floor
   end function

   pure subroutine matched_pair(rows,ss,kk,cap,metric,floor,afloor,result)
      type(candidate_t),intent(in) :: rows(:)
      logical,intent(in) :: ss(:),kk(:)
      real(rk),intent(in) :: cap,floor,afloor
      integer,intent(in) :: metric
      type(comparison_t),intent(out) :: result
      logical :: ps(size(rows)),pk(size(rows)),match,atfloor
      real(rk) :: sumcost,bestcost
      integer :: i,j,s,k
      result=comparison_t();bestcost=huge(1.0_rk)
      call pareto_mask(rows,ss,metric,floor,afloor,ps)
      call pareto_mask(rows,kk,metric,floor,afloor,pk)
      do i=1,size(rows)
         if(.not.ps(i)) cycle
         do j=1,size(rows)
            if(.not.pk(j)) cycle
            if(max(rows(i)%EA,rows(i)%ER,rows(j)%EA,rows(j)%ER)>cap) cycle
            atfloor=max(rows(i)%ER,rows(j)%ER)<=floor
            match=atfloor
            if(min(rows(i)%ER,rows(j)%ER)>floor) &
               match=max(rows(i)%ER,rows(j)%ER)/min(rows(i)%ER,rows(j)%ER)<=10
            if(.not.match) cycle
            sumcost=cost(rows(i),metric,1)+cost(rows(j),metric,1)
            if(sumcost>bestcost) cycle
            if(sumcost==bestcost.and.result%sas>0) then
               if(rows(i)%key>rows(result%sas)%key) cycle
               if(rows(i)%key==rows(result%sas)%key.and.rows(j)%key>=rows(result%space)%key) cycle
            end if
            result%sas=i;result%space=j;bestcost=sumcost
         end do
      end do
      s=result%sas;k=result%space
      if(s==0) return
      result%status='ACTUAL_ERROR_MATCH'
      if(max(rows(s)%ER,rows(k)%ER)<=floor) result%status='REPORTING_FLOOR_MATCH'
      do i=1,size(rows)
         if(pk(i)) result%space_dominates=result%space_dominates.or.dominates(rows(i),rows(s),metric,floor,afloor)
         if(ps(i)) result%sas_dominates=result%sas_dominates.or.dominates(rows(i),rows(k),metric,floor,afloor)
      end do
      call cost_winner(rows,metric,result)
      if((result%winner=='SAS'.and.result%space_dominates).or. &
         (result%winner=='SPACE'.and.result%sas_dominates).or.result%winner=='TIED') then
         result%winner='DOMINATED_OR_TIED';result%separated=.false.
      end if
      call relative_resolution(rows,floor,result)
   end subroutine

   pure subroutine target_pair(rows,ss,kk,cap,metric,floor,result)
      type(candidate_t),intent(in) :: rows(:)
      logical,intent(in) :: ss(:),kk(:)
      real(rk),intent(in) :: cap,floor
      integer,intent(in) :: metric
      type(comparison_t),intent(out) :: result
      integer :: i,side,pick
      result=comparison_t();result%status='TARGET_NOT_BOTH_REACHED'
      do side=1,2
         pick=0
         do i=1,size(rows)
            if(side==1.and..not.ss(i)) cycle
            if(side==2.and..not.kk(i)) cycle
            if(.not.usable(rows(i),metric)) cycle
            if(max(rows(i)%EA,rows(i)%ER)>cap) cycle
            if(pick>0) then
               if(cost(rows(i),metric,1)>cost(rows(pick),metric,1)) cycle
               if(cost(rows(i),metric,1)==cost(rows(pick),metric,1).and.rows(i)%key>=rows(pick)%key) cycle
            end if
            pick=i
         end do
         if(side==1) result%sas=pick
         if(side==2) result%space=pick
      end do
      if(result%sas==0.or.result%space==0) return
      result%status='TARGET_BOTH_REACHED'
      if(cap<=floor) result%status='TARGET_FLOOR_BOTH_REACHED'
      call cost_winner(rows,metric,result);call relative_resolution(rows,floor,result)
   end subroutine

   pure subroutine cost_winner(rows,metric,r)
      type(candidate_t),intent(in) :: rows(:)
      integer,intent(in) :: metric
      type(comparison_t),intent(inout) :: r
      real(rk) :: s,k
      s=cost(rows(r%sas),metric,1);k=cost(rows(r%space),metric,1);r%ratio=k/s
      r%winner='TIED'
      if(s<k) then
         r%winner='SAS';r%separated=cost(rows(r%sas),metric,3)<cost(rows(r%space),metric,2)
      else if(k<s) then
         r%winner='SPACE';r%separated=cost(rows(r%space),metric,3)<cost(rows(r%sas),metric,2)
      end if
   end subroutine
   pure subroutine relative_resolution(rows,floor,r)
      type(candidate_t),intent(in) :: rows(:)
      real(rk),intent(in) :: floor
      type(comparison_t),intent(inout) :: r
      r%resolution='ONE_AT_RELATIVE_FLOOR'
      if(max(rows(r%sas)%ER,rows(r%space)%ER)<=floor) r%resolution='BOTH_RELATIVE_FLOOR'
      if(min(rows(r%sas)%ER,rows(r%space)%ER)>floor) r%resolution='BOTH_ABOVE_RELATIVE_FLOOR'
   end subroutine

   pure integer function best(rows,eligible,B,metric,bound,floor) result(pick)
      type(candidate_t),intent(in) :: rows(:)
      logical,intent(in) :: eligible(:)
      real(rk),intent(in) :: B,floor
      integer,intent(in) :: metric,bound
      real(rk) :: error
      integer :: i
      error=huge(1.0_rk);pick=0
      do i=1,size(rows)
         if(.not.eligible(i).or..not.usable(rows(i),metric)) cycle
         if(cost(rows(i),metric,bound)<=B) error=min(error,rows(i)%ER)
      end do
      do i=1,size(rows)
         if(.not.eligible(i).or..not.usable(rows(i),metric)) cycle
         if(cost(rows(i),metric,bound)>B.or.rows(i)%ER>max(floor,error+floor)) cycle
         if(pick>0) then
            if(cost(rows(i),metric,bound)>cost(rows(pick),metric,bound)) cycle
            if(cost(rows(i),metric,bound)==cost(rows(pick),metric,bound).and.rows(i)%key>=rows(pick)%key) cycle
         end if
         pick=i
      end do
   end function

   pure subroutine budget_point(rows,ss,kk,B,metric,sbound,kbound,floor,afloor,state)
      type(candidate_t),intent(in) :: rows(:)
      logical,intent(in) :: ss(:),kk(:)
      real(rk),intent(in) :: B,floor,afloor
      integer,intent(in) :: metric,sbound,kbound
      type(budget_state_t),intent(out) :: state
      real(rk) :: sa,ka,ratio
      logical :: sj,kj
      integer :: i,s,k
      state=budget_state_t();state%start=B;state%finish=B
      s=best(rows,ss,B,metric,sbound,floor);k=best(rows,kk,B,metric,kbound,floor)
      state%sas=s;state%space=k
      if(s==0.or.k==0) return
      sa=huge(1.0_rk);ka=sa
      do i=1,size(rows)
         if(.not.usable(rows(i),metric)) cycle
         if(ss(i).and.cost(rows(i),metric,sbound)<=B) sa=min(sa,rows(i)%EA)
         if(kk(i).and.cost(rows(i),metric,kbound)<=B) ka=min(ka,rows(i)%EA)
      end do
      sj=rows(s)%EA<=ka+afloor;kj=rows(k)%EA<=sa+afloor
      ratio=max(rows(k)%ER,floor)/max(rows(s)%ER,floor);state%ratio=ratio
      state%resolved=min(rows(s)%ER,rows(k)%ER)>floor
      state%classification='NEITHER_100'
      if(ratio>=100.and.sj) then
         state%classification='SAS_100';state%gain='SAS_100'
      else if(ratio<=.01_rk.and.kj) then
         state%classification='SPACE_100';state%gain='SPACE_100'
      else if(ratio>=10.and.sj) then
         state%gain='SAS_10_to_100'
      else if(ratio<=.1_rk.and.kj) then
         state%gain='SPACE_10_to_100'
      end if
      if(max(rows(s)%ER,rows(k)%ER)<=floor) then
         state%classification='BOTH_RELATIVE_FLOOR';state%gain='NONE'
      end if
   end subroutine

   pure subroutine budget_states(rows,ss,kk,metric,sbound,kbound,floor,afloor,states,nstates)
      type(candidate_t),intent(in) :: rows(:)
      logical,intent(in) :: ss(:),kk(:)
      integer,intent(in) :: metric,sbound,kbound
      real(rk),intent(in) :: floor,afloor
      type(budget_state_t),allocatable,intent(out) :: states(:)
      integer,intent(out) :: nstates
      real(rk) :: events(3*size(rows)),v,previous_event
      type(budget_state_t) :: state
      integer :: i,j,b,n
      n=0;nstates=0
      do i=1,size(rows)
         if(.not.(ss(i).or.kk(i)).or..not.usable(rows(i),metric)) cycle
         do b=1,3
            n=n+1;events(n)=cost(rows(i),metric,b)
         end do
      end do
      allocate(states(n))
      do i=2,n
         v=events(i);j=i-1
         do while(j>=1)
            if(events(j)<=v) exit
            events(j+1)=events(j);j=j-1
         end do
         events(j+1)=v
      end do
      previous_event=-huge(1.0_rk)
      do i=1,n
         if(i>1) then
            if(events(i)==previous_event) cycle
         end if
         previous_event=events(i)
         ! Evaluate all events, including new competitors which only change an EA veto.
         call budget_point(rows,ss,kk,events(i),metric,sbound,kbound,floor,afloor,state)
         if(nstates>0) then
            if(state%sas==states(nstates)%sas.and.state%space==states(nstates)%space &
               .and.state%classification==states(nstates)%classification.and.state%gain==states(nstates)%gain) cycle
            states(nstates)%finish=events(i)
         end if
         nstates=nstates+1;states(nstates)=state
      end do
      if(nstates>0) then
         states(nstates)%finish=events(n)
         ! Terminal event has zero width; retain identity, suppress any gain claim.
         if(states(nstates)%finish==states(nstates)%start) states(nstates)%gain='NONE'
      end if
   end subroutine

   ! The caller supplies an existing cell directory. New files are never overwritten;
   ! append=.true. adds a new case to existing files, without repeating headers.
   subroutine analyze_candidates(caseid,sourceprofile,rho,Fo,rows,floor,afloor,outprefix,status,caps,append)
      character(*),intent(in) :: caseid,sourceprofile,outprefix
      real(rk),intent(in) :: rho,Fo,floor,afloor
      type(candidate_t),intent(in) :: rows(:)
      integer,intent(out) :: status
      real(rk),optional,intent(in) :: caps(:)
      logical,optional,intent(in) :: append
      real(rk),allocatable :: cc(:)
      type(comparison_t) :: result
      type(budget_state_t),allocatable :: states(:)
      character(16),parameter :: families(3)=[character(16)::'DIRECT_2D','SHARED_1D','ALL_SPACE']
      character(20),parameter :: directions(3)=[character(20)::'nominal','SAS_conservative','SPACE_conservative']
      character(5),parameter :: metrics(2)=[character(5)::'cpu','total']
      character(16) :: fam
      character(12) :: mode
      character(8) :: openstatus
      logical :: ss(size(rows)),kk(size(rows)),add
      integer :: u,v,ios,f,c,m,d,i,ns,sb,kb,md
      status=invalid;add=.false.
      if(present(append)) add=append
      if(.not.all(ieee_is_finite([rho,Fo,floor,afloor]))) return
      if(floor<=0.or.afloor<0) return
      cc=[1e-4_rk,1e-6_rk,1e-8_rk,1e-10_rk]
      if(present(caps)) cc=caps
      if(any(.not.ieee_is_finite(cc)).or.any(cc<=0)) return
      openstatus='new'
      if(add) openstatus='old'
      open(newunit=u,file=trim(outprefix)//'/comparisons.csv',status=trim(openstatus), &
         position=merge('append','rewind',add),action='write',iostat=ios)
      if(ios/=0) return
      open(newunit=v,file=trim(outprefix)//'/budget_frontiers.csv',status=trim(openstatus), &
         position=merge('append','rewind',add),action='write',iostat=ios)
      if(ios/=0) then
         close(u);return
      end if
      if(.not.add) then
         write(u,'(a)',iostat=ios) 'case,profile,rho,Fo,relative_floor,absolute_floor,family,cost,mode,cap,'// &
            'status,winner,separated,resolution,SAS_key,SPACE_key,SAS_EA,SAS_ER,SPACE_EA,SPACE_ER,'// &
            'SPACE_over_SAS_cost,SPACE_dominates_SAS,SAS_dominates_SPACE'
         if(ios/=0) goto 900
         write(v,'(a)',iostat=ios) 'case,profile,rho,Fo,relative_floor,absolute_floor,family,cost,direction,'// &
            'start,end,classification,meaningful_gain,SAS_key,SPACE_key,SAS_EA,SAS_ER,SPACE_EA,SPACE_ER,'// &
            'floor_capped_SPACE_over_SAS_ER,both_errors_resolved,finite_bank_right_endpoint,zero_width'
         if(ios/=0) goto 900
      end if
      do f=1,3
         do m=1,2
            do i=1,size(rows)
               fam=route_family(rows(i)%route)
               ss(i)=fam=='SAS'.and.usable(rows(i),m)
               kk(i)=(fam==families(f).or.(f==3.and.(fam=='DIRECT_2D'.or.fam=='SHARED_1D'))).and.usable(rows(i),m)
            end do
            do c=1,size(cc)
               do md=1,2
                  mode='ACTUAL_ERROR'
                  if(md==1) call matched_pair(rows,ss,kk,cc(c),m,floor,afloor,result)
                  if(md==2) then
                     mode='TARGET';call target_pair(rows,ss,kk,cc(c),m,floor,result)
                  end if
                  if(.not.any(ss).and..not.any(kk)) result%status='HELD_ABSOLUTE_ONLY_OR_UNAVAILABLE'
                  call write_comparison(u,caseid,sourceprofile,rho,Fo,floor,afloor,families(f),metrics(m), &
                     mode,cc(c),rows,result,ios)
                  if(ios/=0) goto 900
               end do
            end do
            do d=1,3
               sb=1;kb=1
               if(d==2) then
                  sb=3;kb=2
               else if(d==3) then
                  sb=2;kb=3
               end if
               call budget_states(rows,ss,kk,m,sb,kb,floor,afloor,states,ns)
               do i=1,ns
                  call write_budget(v,caseid,sourceprofile,rho,Fo,floor,afloor,families(f),metrics(m), &
                     directions(d),rows,states(i),i==ns,ios)
                  if(ios/=0) goto 900
               end do
            end do
         end do
      end do
      ! Explicit absolute-only evidence for weak-reference rows; no relative rank is fabricated.
      do i=1,size(rows)
         if(.not.rows(i)%valid.or.rows(i)%active) cycle
         result=comparison_t();result%status='HELD_ABSOLUTE_ONLY';result%resolution='ER_UNAVAILABLE'
         fam=route_family(rows(i)%route)
         if(fam=='SAS') then
            result%sas=i
         else
            result%space=i
         end if
         call write_comparison(u,caseid,sourceprofile,rho,Fo,floor,afloor,fam,'none','ABSOLUTE', &
            0.0_rk,rows,result,ios)
         if(ios/=0) goto 900
      end do
      status=ok
900   close(u,iostat=ios)
      if(ios/=0) status=invalid
      close(v,iostat=ios)
      if(ios/=0) status=invalid
   end subroutine

   subroutine selected_fields(rows,s,k,sk,kk,ea,er)
      type(candidate_t),intent(in) :: rows(:)
      integer,intent(in) :: s,k
      character(128),intent(out) :: sk,kk
      character(26),intent(out) :: ea(2),er(2)
      integer :: i,j
      sk='';kk='';ea='';er=''
      if(s>0) sk=rows(s)%key
      if(k>0) kk=rows(k)%key
      do j=1,2
         i=s
         if(j==2) i=k
         if(i==0) cycle
         write(ea(j),'(es25.17e3)') rows(i)%EA
         if(rows(i)%active) write(er(j),'(es25.17e3)') rows(i)%ER
      end do
   end subroutine
   subroutine write_comparison(u,id,profile,rho,Fo,floor,afloor,family,metric,mode,cap,rows,r,ios)
      integer,intent(in) :: u
      character(*),intent(in) :: id,profile,family,metric,mode
      real(rk),intent(in) :: rho,Fo,floor,afloor,cap
      type(candidate_t),intent(in) :: rows(:)
      type(comparison_t),intent(in) :: r
      integer,intent(out) :: ios
      character(128) :: sk,kk
      character(26) :: ea(2),er(2)
      call selected_fields(rows,r%sas,r%space,sk,kk,ea,er)
      write(u,'(*(g0,:,","))',iostat=ios) trim(id),trim(profile),rho,Fo,floor,afloor,trim(family), &
         trim(metric),trim(mode),cap,trim(r%status),trim(r%winner),r%separated,trim(r%resolution), &
         trim(sk),trim(kk),trim(adjustl(ea(1))),trim(adjustl(er(1))),trim(adjustl(ea(2))),trim(adjustl(er(2))), &
         r%ratio,r%space_dominates,r%sas_dominates
   end subroutine
   subroutine write_budget(u,id,profile,rho,Fo,floor,afloor,family,metric,direction,rows,r,last,ios)
      integer,intent(in) :: u
      character(*),intent(in) :: id,profile,family,metric,direction
      real(rk),intent(in) :: rho,Fo,floor,afloor
      type(candidate_t),intent(in) :: rows(:)
      type(budget_state_t),intent(in) :: r
      logical,intent(in) :: last
      integer,intent(out) :: ios
      character(128) :: sk,kk
      character(26) :: ea(2),er(2)
      call selected_fields(rows,r%sas,r%space,sk,kk,ea,er)
      write(u,'(*(g0,:,","))',iostat=ios) trim(id),trim(profile),rho,Fo,floor,afloor,trim(family), &
         trim(metric),trim(direction),r%start,r%finish,trim(r%classification),trim(r%gain),trim(sk),trim(kk), &
         trim(adjustl(ea(1))),trim(adjustl(er(1))),trim(adjustl(ea(2))),trim(adjustl(er(2))), &
         r%ratio,r%resolved,last,r%start==r%finish
   end subroutine
end module
