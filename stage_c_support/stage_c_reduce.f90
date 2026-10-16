program stage_c_reduce
  use iso_fortran_env, only: real64
  use ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  integer,parameter :: r=real64,nm=34
  real(r),parameter :: case_abs=1e-15_r,case_rel=1e-12_r,root_rel=1e-9_r,ref_rel=1e-11_r
  real(r),parameter :: periodic_limit=1e-9_r,coord_limit=8e-15_r
  real(r),parameter :: duties(3)=[.1_r,.4_r,.7_r]
  character(1024) :: indir,outdir
  integer :: checks,holds
  call get_command_argument(1,indir);call get_command_argument(2,outdir)
  if(min(len_trim(indir),len_trim(outdir))==0)error stop 'stage_c_reduce FRESH_RESULTS_DIR NEW_OUTPUT_DIR'
  holds=0
  call preflight()
  open(newunit=checks,file=trim(outdir)//'/refinement_checks.csv',status='new',action='write')
  write(checks,'(a)')'source,metric,observed,limit,state'
  call statuses()
  call case1()
  call gaussian_refinements()
  call roots()
  call periodic()
  close(checks)
  call finish()
contains
  subroutine preflight()
    integer :: i
    character(32) :: tag
    call need('case1_profiles.csv');call need('case3_table.csv')
    call need('case2_profiles.csv');call need('case2_convergence.csv')
    call need('periodic_refinement_union.csv')
    call need('validation_summary.txt');call need('periodic_summary.txt')
    do i=1,8
      write(tag,'(i0)')i
      call need('case3_grid_'//trim(tag)//'.csv')
      call need('periodic_grid_'//trim(tag)//'.csv')
    end do
  end subroutine
  subroutine need(name)
    character(*),intent(in) :: name
    logical :: exists
    inquire(file=trim(indir)//'/'//name,exist=exists)
    if(.not.exists)then
      print '(a)','Missing required fresh input: '//name
      error stop 'run both native aggregates first; no input substitution'
    end if
  end subroutine
  subroutine csv(name,u)
    character(*),intent(in) :: name
    integer,intent(out) :: u
    character(4096) :: header
    open(newunit=u,file=trim(indir)//'/'//name,status='old',action='read')
    read(u,'(a)')header
  end subroutine
  subroutine check(source,key,value,limit)
    character(*),intent(in) :: source,key
    real(r),intent(in) :: value,limit
    character(4) :: state
    state='PASS'
    if(.not.ieee_is_finite(value).or.value>limit)then
      holds=holds+1;state='HOLD'
    end if
    write(checks,'(a,",",a,2(",",es25.16),",",a)')source,key,value,limit,state
  end subroutine
  subroutine observe(source,key,value)
    character(*),intent(in) :: source,key
    real(r),intent(in) :: value
    write(checks,'(a,",",a,",",es25.16,",,OBSERVE")')source,key,value
  end subroutine
  subroutine statuses()
    integer :: i
    character(16) :: tag
    call status_one('core_summary.txt');call status_one('case1_summary.txt');call status_one('gaussian_summary.txt')
    call status_one('precision_summary.txt');call status_one('table_s2_summary.txt');call status_one('validation_summary.txt')
    call status_one('periodic_unit.txt');call status_one('periodic_summary.txt')
    do i=1,8
      write(tag,'(i0)')i
      call status_one('roots_'//trim(tag)//'_summary.txt');call status_one('periodic_grid_'//trim(tag)//'_summary.txt')
    end do
    do i=1,4
      write(tag,'(i0)')i
      call status_one('periodic_invariants_'//trim(tag)//'.txt')
    end do
  end subroutine
    subroutine status_one(name)
      character(*),intent(in) :: name
      character(4096) :: line
      integer :: u,ios
      logical :: passed,failed
      call need(name);open(newunit=u,file=trim(indir)//'/'//name,status='old',action='read')
      passed=.false.;failed=.false.
      do
        read(u,'(a)',iostat=ios)line
        if(ios/=0)exit
        if(index(line,'status=PASS')>0.or.index(line,'status: PASS')>0)passed=.true.
        if(index(line,'status=FAIL')>0.or.index(line,'status=HOLD')>0)failed=.true.
      end do
      close(u)
      call check(name,'missing_PASS_or_explicit_failure',merge(0.0_r,1.0_r,passed.and..not.failed),0.0_r)
    end subroutine
  subroutine case1()
    real(r) :: a(11),delta(4),reference(4),ab(4,4),rel(4,4),mixed(4,4),ratio
    integer :: u,out,ios,n,g,k,groups(2),count(4),row(4,4),j
    character(24),parameter :: labels(4)=[character(24)::'GL200_GL400','GL400_partition120','SAS_partition120','partition80_120']
    character(16) :: tag
    call csv('case1_profiles.csv',u)
    ab=0;rel=0;mixed=0;count=0;row=0;n=0
    do
      read(u,*,iostat=ios)a
      if(ios<0)exit
      if(ios/=0)error stop 'invalid Case1 CSV'
      n=n+1;if(n>150)error stop 'extra Case1 row'
      g=(n-1)/50+1;j=mod(n-1,50)
      if(maxval(abs(a(:3)-[2.0_r*(g-1),0.0_r,.1_r+1.4_r*j/49]))>coord_limit)error stop 'Case1 input identity'
      if(a(3)<.5_r-coord_limit)cycle
      groups=[g,4]
      delta=abs([a(8)-a(9),a(9)-a(6),a(4)-a(6),a(5)-a(6)])
      reference=[a(9),a(6),a(6),a(6)]
      if(.not.all(ieee_is_finite(delta)))then
        call check('case1_profiles.csv','nonfinite_value',1.0_r,0.0_r);cycle
      end if
      do j=1,2
        g=groups(j);count(g)=count(g)+1
        do k=1,4
          ratio=delta(k)/(case_abs+case_rel*abs(reference(k)))
          ab(k,g)=max(ab(k,g),delta(k));rel(k,g)=max(rel(k,g),delta(k)/max(abs(reference(k)),tiny(1.0_r)))
          if(ratio>mixed(k,g))then
            mixed(k,g)=ratio;row(k,g)=n
          end if
        end do
      end do
    end do
    close(u);if(n/=150)error stop 'incomplete Case1 bank'
    open(newunit=out,file=trim(outdir)//'/case1_z_ge_half.csv',status='new',action='write')
    write(out,'(a)')'line_index_1to3_or4_all,metric,points,max_absolute,max_relative,max_over_mixed_limit,max_mixed_source_row,atol,rtol'
    do g=1,4
      write(tag,'(i0)')g
      call check('case1_profiles.csv:line='//trim(tag),'selected_count',abs(real(count(g)-merge(36,108,g<4),r)),0.0_r)
      do k=1,4
        write(out,'(i0,",",a,",",i0,3(",",es25.16),",",i0,2(",",es25.16))')g,trim(labels(k)),count(g),ab(k,g),rel(k,g),mixed(k,g),row(k,g),case_abs,case_rel
        call check('case1_profiles.csv:line='//trim(tag),trim(labels(k))//'_mixed_ratio',mixed(k,g),1.0_r)
      end do
    end do
    close(out)
  end subroutine
  subroutine gaussian_refinements()
    integer :: u,ios,m,n
    real(r) :: a(9),refine,eval,tolerance,coeff
    refine=0;eval=0;tolerance=0;coeff=0;n=0
    call csv('case2_profiles.csv',u)
    do
      read(u,*,iostat=ios)m,a
      if(ios<0)exit
      if(ios/=0)error stop 'invalid Gaussian profile CSV'
      n=n+1
      if(.not.all(ieee_is_finite(a)))call check('case2_profiles.csv','nonfinite_profile',1.0_r,0.0_r)
      refine=max(refine,abs(a(5)*a(4))/(case_abs+case_rel*abs(a(4))))
      eval=max(eval,abs(a(6)-a(7))/(case_abs+1e-11_r*abs(a(7))))
    end do
    close(u);if(n/=600)error stop 'Gaussian profile coverage must be six M x 100'
    call check('case2_profiles.csv','Gaussian80_120_mixed_ratio',refine,1.0_r)
    call check('case2_profiles.csv','same_polynomial_evaluation_mixed_ratio',eval,1.0_r)
    call csv('case2_convergence.csv',u);n=0
    do
      read(u,*,iostat=ios)m,a
      if(ios<0)exit
      if(ios/=0)error stop 'invalid Gaussian convergence CSV'
      n=n+1;if(m/=2*n)error stop 'Gaussian convergence order identity'
      if(.not.all(ieee_is_finite(a)))call check('case2_convergence.csv','nonfinite_convergence',1.0_r,0.0_r)
      tolerance=max(tolerance,abs(a(8)));coeff=max(coeff,abs(a(9)))
      call check('case2_convergence.csv','same_polynomial_transient_mixed_ratio',abs(a(3)-a(4))/(case_abs+1e-11_r*abs(a(4))),1.0_r)
    end do
    close(u);if(n/=6)error stop 'six Gaussian convergence orders required'
    call check('case2_convergence.csv','SAS_tolerance_change_abs',tolerance,1e-11_r)
    call check('case2_convergence.csv','stored_cosine512_1024_change_abs',coeff,1e-12_r)
  end subroutine
  subroutine roots()
    real(r) :: z,v(8),table(2,4,3),table_sas(2,4,3),a(6),zdepth(4),reduction,other,maxerr,maxref
    integer :: u,out,i,k,j,n,part,node,ios,p,seen(401,2),ts(2,4,3)
    character(16) :: tag
    character :: profile
    zdepth=[.1_r,.5_r,1.0_r,2.0_r];table=0;table_sas=0;ts=0;seen=0;maxerr=0;maxref=0;n=0
    open(newunit=out,file=trim(outdir)//'/case3_z4_ratios.csv',status='new',action='write')
    write(out,'(a)')'profile,source_chunk,node,z,SAS_t90_over_z2,reference64_t90_over_z2,reference96_t90_over_z2'
    do part=1,8
      write(tag,'(i0)')part;call csv('case3_grid_'//trim(tag)//'.csv',u)
      do
        read(u,*,iostat=ios)node,z,profile,v
        if(ios<0)exit
        if(ios/=0.or.node<0.or.node>400)error stop 'invalid Case3 grid row'
        p=1;if(profile=='Q')p=2
        if(profile/='Q'.and.profile/='U')error stop 'unknown root profile'
        if(abs(z-(.1_r+.00975_r*node))>coord_limit.or.mod(node,8)/=part-1)error stop 'Case3 node identity'
        if(seen(node+1,p)/=0)error stop 'duplicate Case3 root'
        seen(node+1,p)=1;n=n+1
        call check('case3_grid_'//trim(tag)//'.csv','finite_positive_roots',merge(0.0_r,1.0_r,all(ieee_is_finite(v(:3))).and.all(v(:3)>0)),0.0_r)
        maxerr=max(maxerr,abs(v(1)-v(3))/(case_abs+root_rel*abs(v(3))))
        maxref=max(maxref,abs(v(2)-v(3))/(case_abs+ref_rel*abs(v(3))))
        if(node==400)write(out,'(a,2(",",i0),4(",",es25.16))')profile,part,node,z,v(1)/z**2,v(2)/z**2,v(3)/z**2
      end do
      close(u)
    end do
    close(out);if(n/=802.or.any(seen/=1))error stop 'incomplete 401x2 root bank'
    call check('case3_grid_1..8.csv','SAS_reference96_mixed_ratio',maxerr,1.0_r)
    call check('case3_grid_1..8.csv','reference64_96_mixed_ratio',maxref,1.0_r)
    call csv('case3_table.csv',u);n=0
    do
      read(u,*,iostat=ios)i,profile,a
      if(ios<0)exit
      if(ios/=0.or.i<1.or.i>3)error stop 'invalid table root row'
      p=1;if(profile=='Q')p=2
      if(profile/='Q'.and.profile/='U')error stop 'unknown table profile'
      k=minloc(abs(zdepth-a(1)),dim=1)
      if(abs(zdepth(k)-a(1))>coord_limit.or.ts(p,k,i)/=0)error stop 'table root identity'
      ts(p,k,i)=1;table_sas(p,k,i)=a(2);table(p,k,i)=a(4);n=n+1
      call check('case3_table.csv','SAS_reference96_mixed_ratio',abs(a(2)-a(4))/(case_abs+root_rel*abs(a(4))),1.0_r)
      call check('case3_table.csv','reference64_96_mixed_ratio',abs(a(3)-a(4))/(case_abs+ref_rel*abs(a(4))),1.0_r)
    end do
    close(u);if(n/=24.or.any(ts/=1).or.any(table<=0).or.any(table_sas<=0))error stop 'incomplete/invalid 24-root table'
    open(newunit=out,file=trim(outdir)//'/case3_six_reductions.csv',status='new',action='write')
    write(out,'(a)')'geometry,z,U_SAS,Q_SAS,U_reference96,Q_reference96,reduction_SAS_percent,reduction_reference96_percent'
    do i=1,3
      do j=1,2
        k=1;if(j==2)k=4
        reduction=100*(1-table(2,k,i)/table(1,k,i));other=100*(1-table_sas(2,k,i)/table_sas(1,k,i))
        write(out,'(i0,7(",",es25.16))')i,zdepth(k),table_sas(:,k,i),table(:,k,i),other,reduction
        call check('case3_table.csv','reduction_percent_consistency',abs(other-reduction),200*root_rel+root_rel*abs(reduction))
      end do
    end do
    close(out)
  end subroutine
  subroutine periodic()
    real(r) :: rows(nm,96),v(nm),L,d,z
    integer :: seen(96),union_seen(96),u,out,i,k,p,depth,j,ios,n,index
    character :: profile
    character(16) :: tag
    seen=0;union_seen=0;rows=ieee_value(0.0_r,ieee_quiet_nan)
    do i=1,8
      write(tag,'(i0)')i;call csv('periodic_grid_'//trim(tag)//'.csv',u)
      do k=1,12
        read(u,*,iostat=ios)j,L,profile,d,z,v
        if(ios/=0.or.j/=i-1)error stop 'incomplete periodic grid'
        p=1;if(k>6)p=2
        depth=mod(k-1,2)+1;index=12*(i-1)+k
        if(profile/=merge('U','Q',p==1).or.abs(L-(.2_r+4.0_r*j/35))>coord_limit.or. &
          abs(d-duties((mod(k-1,6))/2+1))>coord_limit)error stop 'periodic parameter identity'
        if(abs(z-merge(.1_r,.5_r,depth==1))>coord_limit)error stop 'periodic depth identity'
        if(v(33)/=merge(1.0_r,0.0_r,j==7.or.(j==0.and.profile=='Q'.and. &
          abs(d-.1_r)<coord_limit)))error stop 'fixed periodic refinement identity'
        rows(:,index)=v;seen(index)=1
      end do
      read(u,*,iostat=ios)
      if(ios==0)error stop 'extra periodic grid row'
      close(u)
    end do
    if(any(seen/=1))error stop 'missing periodic depth'
    if(count(rows(33,:)==1)/=14)error stop 'fixed 14 refinement controls missing'
    call csv('periodic_refinement_union.csv',u);n=0
    do
      read(u,*,iostat=ios)j,L,profile,d,z,v
      if(ios<0)exit
      if(ios/=0.or.j<0.or.j>7)error stop 'invalid refinement union'
      p=1;if(profile=='Q')p=2
      if(profile/='U'.and.profile/='Q')error stop 'invalid union profile'
      k=nint((d-.1_r)/.3_r)+1;depth=1;if(abs(z-.5_r)<coord_limit)depth=2
      if(k<1.or.k>3.or.abs(z-merge(.1_r,.5_r,depth==1))>coord_limit)error stop 'union parameter identity'
      if(abs(d-duties(k))>coord_limit.or.abs(L-(.2_r+4.0_r*j/35))>coord_limit)error stop 'union geometry/duty identity'
      index=12*j+6*(p-1)+2*(k-1)+depth
      if(union_seen(index)/=0.or.v(33)/=1.or.rows(33,index)/=1)error stop 'duplicate/nonselected union row'
      call check('periodic_refinement_union.csv','base_swing_identity',abs(v(5)-rows(5,index)),1e-13_r)
      union_seen(index)=1;n=n+1;rows(25:33,index)=v(25:33)
    end do
    close(u)
    if(n/=14.or.any(union_seen/=nint(rows(33,:))))error stop 'refinement union must contain the fixed 14 depths'
    call check('periodic_refinement_union.csv','selected_union_membership',real(count(rows(33,:)==1)-n,r),0.0_r)
    open(newunit=out,file=trim(outdir)//'/periodic_s5.csv',status='new',action='write')
    write(out,'(a)')'grid_source,source_row,j,L,profile,d,z,unrounded_swing'
    do k=1,12
      index=84+k;p=1;if(k>6)p=2
      profile=merge('U','Q',p==1);d=duties((mod(k-1,6))/2+1);z=merge(.1_r,.5_r,mod(k,2)==1)
      write(out,'(a,",",i0,",",i0,",",es25.16,",",a,3(",",es25.16))')'periodic_grid_8.csv',k,7,1.0_r,profile,d,z,rows(5,index)
    end do
    close(out)
    call periodic_checks(rows)
    call compare_summary('refinement_union_depths=',real(n,r),0.0_r)
    call compare_summary('failed_depth_acceptance=',0.0_r,0.0_r)
    call summary_truth('all_chunk_acceptance=')
  end subroutine
  subroutine periodic_checks(v)
    real(r),intent(in) :: v(nm,96)
    integer :: i
    character(16) :: tag
    character(64) :: source
    do i=1,96
      write(tag,'(i0)')i;source='periodic_grid/union:depth='//trim(tag)
      call check(trim(source),'reference_value_abs',v(15,i),periodic_limit)
      call check(trim(source),'reference_swing_abs',v(16,i),periodic_limit)
      call check(trim(source),'independent_stationary_swing_abs',v(20,i),periodic_limit)
      call check(trim(source),'phase128_512_swing_abs',v(24,i),periodic_limit)
      call check(trim(source),'two_root_brackets',abs(v(32,i)-2),0.0_r)
      call check(trim(source),'swing_above_own_steady',max(0.0_r,v(5,i)-v(7,i)),periodic_limit)
      if(v(33,i)/=1)cycle
      call check(trim(source),'doubled_quadrature_swing_abs',v(25,i),periodic_limit)
      call check(trim(source),'history400_800_over_bound',v(26,i)/(v(27,i)+periodic_limit),1.0_r)
      call check(trim(source),'spatial48_SAS_swing_abs',v(28,i),periodic_limit)
      call check(trim(source),'spatial32_48_swing_abs',v(29,i),periodic_limit)
    end do
    call observe('periodic_grid/union','max_stationary_phase_difference_no_new_limit',maxval(v(21,:)))
    call observe('periodic_grid/union','max_derivative_residual_no_new_limit',maxval(v(23,:)))
    call observe('periodic_grid/union','max_analytic_tail_temperature',maxval(v(8,:)))
    call observe('periodic_grid/union','max_analytic_tail_swing',maxval(v(9,:)))
  end subroutine
  subroutine compare_summary(key,expected,tolerance)
    character(*),intent(in) :: key
    real(r),intent(in) :: expected,tolerance
    real(r) :: value
    integer :: u,ios,p
    character(4096) :: line
    logical :: found
    found=.false.;open(newunit=u,file=trim(indir)//'/periodic_summary.txt',status='old',action='read')
    do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      p=index(line,key);if(p==0)cycle
      read(line(p+len(key):),*,iostat=ios)value
      if(ios/=0)error stop 'malformed periodic summary scalar'
      found=.true.;call check('periodic_summary.txt',key,abs(value-expected),tolerance)
    end do
    close(u)
    if(.not.found)error stop 'missing periodic summary scalar'
  end subroutine
  subroutine summary_truth(key)
    character(*),intent(in) :: key
    integer :: u,ios,p
    logical :: value,found
    character(4096) :: line
    found=.false.;open(newunit=u,file=trim(indir)//'/periodic_summary.txt',status='old',action='read')
    do
      read(u,'(a)',iostat=ios)line
      if(ios/=0)exit
      p=index(line,key);if(p==0)cycle
      read(line(p+len(key):),*,iostat=ios)value
      if(ios/=0)error stop 'malformed periodic summary logical'
      found=.true.;call check('periodic_summary.txt',key,merge(0.0_r,1.0_r,value),0.0_r)
    end do
    close(u)
    if(.not.found)error stop 'missing periodic summary logical'
  end subroutine
  subroutine finish()
    integer :: u
    open(newunit=u,file=trim(outdir)//'/stage_c_support_summary.txt',status='new',action='write')
    write(u,'(a)')'Pure Fortran reduction of caller-selected fresh results; no physical solver or performance callback invoked.'
    write(u,'(a)')'Fixed mixed tolerances: Case1 1e-15+1e-12|ref|; root agreement 1e-15+1e-9|ref|; root refinement 1e-15+1e-11|ref|.'
    write(u,'(a)')'Periodic absolute tolerance 1e-9; history refinement uses the existing bound+1e-9. OBSERVE rows have no new numerical acceptance limit.'
    write(u,'(a)')'Gaussian reference80/120 mixed tolerance 1e-15+1e-12|ref|; same-polynomial 1e-15+1e-11|ref|; SAS tolerance change 1e-11 absolute; stored cosine change 1e-12 absolute.'
    write(u,'(a,i0)')'HOLD_checks=',holds
    if(holds==0)then
      write(u,'(a)')'status=PASS'
    else
      write(u,'(a)')'status=HOLD'
    end if
    close(u)
    if(holds/=0)error stop 'stage-C support HOLD; retain outputs and original fresh evidence'
  end subroutine
end program
