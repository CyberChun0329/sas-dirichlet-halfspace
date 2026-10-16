program periodic_driver
  use periodic_validation
  use validation_sas, only: rk,qk,sas_center_step
  use validation_reference, only: ref_periodic_context,ref_prepare_periodic,ref_periodic_value, &
    ref_periodic_derivative_value
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  character(32) :: mode,arg
  integer :: index,ios
  call get_command_argument(1,mode)
  select case(trim(mode))
  case('unit')
    call unit_checks()
  case('grid','invariants')
    call get_command_argument(2,arg);read(arg,*,iostat=ios) index
    if(ios/=0) error stop 'usage: periodic_driver grid <1..8> | invariants <1..4>'
    if(trim(mode)=='grid') then
      if(index<1.or.index>8) error stop 'grid chunk must be 1..8'
      call grid_chunk(index)
    else
      if(index<1.or.index>4) error stop 'invariant case must be 1..4'
      call invariant_case(index)
    end if
  case('aggregate')
    call aggregate()
  case default
    error stop 'usage: periodic_driver unit | grid <1..8> | invariants <1..4> | aggregate'
  end select
contains
  subroutine new_output(path,u)
    character(*), intent(in) :: path
    integer, intent(out) :: u
    integer :: status
    open(newunit=u,file=path,status='new',action='write',iostat=status)
    if(status/=0) then
      write(*,'(A)') 'Cannot create new output (existing files are preserved): '//trim(path)
      error stop 2
    end if
  end subroutine

  subroutine unit_checks()
    type(periodic_context) :: ctx
    type(ref_periodic_context) :: reference
    type(extrema_result) :: naive,fixed
    type(tail_result) :: bounds
    real(rk) :: z,S,analytic,v,left,right,max_steady,max_step,max_endpoint,max_switch,max_derivative
    real(rk) :: max_zero,full_min,full_max,age,step,phase
    real(qk) :: value,remainder,refine,next,remnext,refnext,zeta_s,a,direct,edge,recurrence,max_rec,max_ref
    character :: profile
    integer :: u,k,p,depth,status,failures
    logical :: okay
    call new_output('results/periodic_unit.txt',u)
    failures=0;max_steady=0;max_step=0;max_endpoint=0;max_switch=0;max_derivative=0
    max_rec=0;max_ref=0
    write(u,'(A)') 'Fortran periodic unit checks.'
    write(u,'(A)') 'S-T is the instantaneous direct integral, not the cumulative deficit D.'
    do p=1,2
      profile='U';if(p==2) profile='Q'
      do depth=1,2
        z=0.1_rk;if(depth==2) z=0.5_rk
        call periodic_prepare(ctx,z,1.0_rk,profile)
        call ref_prepare_periodic(reference,z,1.0_rk,profile,64)
        S=sum(reference%weight)
        max_steady=max(max_steady,abs(S-ctx%steady))
        if(profile=='U') then
          analytic=2/acos(-1.0_rk)*atan(1/(z*sqrt(2+z*z)))
          max_steady=max(max_steady,abs(analytic-ctx%steady))
        end if
        if(periodic_deficit(ctx,0.0_rk)/=ctx%steady.or. &
          periodic_deficit(ctx,-1.0_rk)/=ctx%steady) failures=failures+1
        do k=1,4
          select case(k)
          case(1);age=0.02_rk
          case(2);age=0.5_rk
          case(3);age=2.0_rk
          case(4);age=20.0_rk
          end select
          step=sas_center_step(age,z,1.0_rk,1.0_rk,profile,2e-12_rk,status)
          if(status/=0) failures=failures+1
          max_step=max(max_step,abs(step+periodic_deficit(ctx,age)-ctx%steady))
        end do
        bounds=tail_bounds(z,1.0_rk,profile,0.4_rk)
        max_endpoint=max(max_endpoint,abs(periodic_value(ctx,0.0_rk,0.4_rk)- &
          periodic_value(ctx,period,0.4_rk)))
        if(abs(periodic_value(ctx,0.0_rk,0.4_rk)-periodic_value(ctx,period,0.4_rk))> &
          bounds%swing+1e-12_rk) failures=failures+1
        left=periodic_value(ctx,0.2_rk-1e-10_rk,0.4_rk)
        right=periodic_value(ctx,0.2_rk+1e-10_rk,0.4_rk)
        max_switch=max(max_switch,abs(left-right))
        v=periodic_derivative(ctx,0.037_rk,0.4_rk)
        max_derivative=max(max_derivative,abs(v-ref_periodic_derivative_value(reference,0.037_rk,0.4_rk,400)))
        write(u,'(A,1X,A,1X,F5.2,3(1X,ES24.16E3))') 'source',profile,z,ctx%steady,S,bounds%temperature
      end do
    end do
    do k=1,3
      zeta_s=real(k,qk)+1.5_qk;a=399.6_qk
      call hurwitz_tail(zeta_s,a,value,remainder,refine)
      call hurwitz_tail(zeta_s,a+1,next,remnext,refnext)
      recurrence=abs(value-next-a**(-zeta_s));max_rec=max(max_rec,recurrence);max_ref=max(max_ref,refine)
      direct=0
      do p=0,19999
        direct=direct+(a+p)**(-zeta_s)
      end do
      edge=a+20000
      okay=value>=direct+edge**(1-zeta_s)/(zeta_s-1)-1e-30_qk.and. &
        value<=direct+edge**(1-zeta_s)/(zeta_s-1)+edge**(-zeta_s)+1e-30_qk
      if(.not.okay.or.recurrence>1e-30_qk.or.refine>1e-30_qk) failures=failures+1
      write(u,'(A,F4.1,4(1X,ES25.16E3),1X,L1)') 'zeta ',real(zeta_s,rk),value,remainder,refine,recurrence,okay
    end do
    call periodic_prepare(ctx,0.5_rk,1.0_rk,'U')
    max_zero=0;full_min=huge(1.0_rk);full_max=-huge(1.0_rk)
    do k=0,16
      phase=period*k/16
      max_zero=max(max_zero,abs(periodic_value(ctx,phase,0.0_rk)))
      v=periodic_value(ctx,phase,1.0_rk);full_min=min(full_min,v);full_max=max(full_max,v)
    end do
    bounds=tail_bounds(0.5_rk,1.0_rk,'U',1.0_rk)
    if(max_zero>=1e-12_rk.or.full_max-full_min>bounds%swing+1e-9_rk) failures=failures+1
    call periodic_prepare(ctx,0.1_rk,1.0_rk,'Q')
    call periodic_extrema(ctx,0.4_rk,400,128,naive,.false.)
    call periodic_extrema(ctx,0.4_rk,400,128,fixed)
    if(fixed%swing<=naive%swing+1e-8_rk.or.fixed%tmin<=0.or.fixed%tmin>=period/128) failures=failures+1
    if(max_steady>1e-10_rk.or.max_step>1e-9_rk.or.max_switch>1e-8_rk.or.max_derivative>1e-9_rk) &
      failures=failures+1
    write(u,'(A,ES24.16E3)') 'max_steady_spatial64_or_uniform_analytic_difference=',max_steady
    write(u,'(A,ES24.16E3)') 'max_step_plus_instantaneous_deficit_error=',max_step
    write(u,'(A,ES24.16E3)') 'max_endpoint_difference=',max_endpoint
    write(u,'(A,ES24.16E3)') 'max_switch_neighbour_difference_2e-10_phase=',max_switch
    write(u,'(A,ES24.16E3)') 'max_independent_derivative_difference=',max_derivative
    write(u,'(A,ES24.16E3)') 'zero_duty_max=',max_zero
    write(u,'(A,ES24.16E3)') 'full_duty_finite_history_swing=',full_max-full_min
    write(u,'(A,ES24.16E3)') 'full_duty_tail_swing_bound=',bounds%swing
    write(u,'(A,ES24.16E3)') 'endpoint_search_swing_correction=',fixed%swing-naive%swing
    write(u,'(A,2(1X,ES24.16E3))') 'naive_min_phase_value=',naive%tmin,naive%minimum
    write(u,'(A,2(1X,ES24.16E3))') 'searched_min_phase_value=',fixed%tmin,fixed%minimum
    write(u,'(A,I0)') 'failures=',failures
    if(failures==0) then
      write(u,'(A)') 'status=PASS'
    else
      write(u,'(A)') 'status=FAIL'
    end if
    close(u)
    if(failures/=0) error stop 'periodic unit checks failed; inspect results/periodic_unit.txt'
    write(*,'(A)') 'periodic unit PASS: results/periodic_unit.txt'
  end subroutine

  subroutine summarize(u,rows,n)
    integer, intent(in) :: u,n
    real(rk), intent(in) :: rows(nmetrics,n)
    integer :: i,selected,root_count,failures
    real(rk) :: max_order,max_history,max_space48,max_spacechange,max_history_ratio
    logical :: good
    failures=0;selected=0;root_count=0
    max_order=0;max_history=0;max_space48=0;max_spacechange=0;max_history_ratio=0
    do i=1,n
      good=all(ieee_is_finite(rows(1:m_phase_change,i)))
      if(.not.good) then
        failures=failures+1;cycle
      end if
      root_count=root_count+nint(rows(m_root_brackets,i))
      good=rows(m_root_brackets,i)==2.and.rows(m_ref_value_diff,i)<1e-9_rk.and. &
        rows(m_ref_swing_diff,i)<1e-9_rk.and.rows(m_root_swing_diff,i)<1e-9_rk.and. &
        rows(m_phase_change,i)<1e-9_rk.and.rows(m_min,i)>=-1e-10_rk.and. &
        rows(m_max,i)<=rows(m_steady,i)+1e-10_rk
      if(.not.good) failures=failures+1
      if(rows(m_selected,i)==1) then
        selected=selected+1
        if(.not.all(ieee_is_finite(rows(m_order_change:m_history_swing,i)))) then
          failures=failures+1;cycle
        end if
        max_order=max(max_order,rows(m_order_change,i));max_history=max(max_history,rows(m_history_change,i))
        max_space48=max(max_space48,rows(m_spatial48_diff,i))
        max_spacechange=max(max_spacechange,rows(m_spatial48_change,i))
        max_history_ratio=max(max_history_ratio,rows(m_history_change,i)/rows(m_history_bound,i))
        if(rows(m_order_change,i)>=1e-9_rk.or.rows(m_spatial48_diff,i)>=1e-9_rk.or. &
          rows(m_history_change,i)>=rows(m_history_bound,i)+1e-9_rk) failures=failures+1
      end if
    end do
    write(u,'(A,I0)') 'depth_count=',n
    write(u,'(A,I0)') 'independent_derivative_root_count=',root_count
    write(u,'(A,I0)') 'selected_refinement_count=',selected
    write(u,'(A,ES24.16E3)') 'max_reference_temperature_difference=',maxval(rows(m_ref_value_diff,:))
    write(u,'(A,ES24.16E3)') 'max_reference_swing_difference=',maxval(rows(m_ref_swing_diff,:))
    write(u,'(A,ES24.16E3)') 'max_stationary_swing_difference=',maxval(rows(m_root_swing_diff,:))
    write(u,'(A,ES24.16E3)') 'max_stationary_phase_difference=',maxval(rows(m_root_phase_diff,:))
    write(u,'(A,ES24.16E3)') 'max_independent_derivative_difference=',maxval(rows(m_derivative_diff,:))
    write(u,'(A,ES24.16E3)') 'max_derivative_root_residual=',maxval(rows(m_derivative_residual,:))
    write(u,'(A,ES24.16E3)') 'max_phase128_to512_swing_change=',maxval(rows(m_phase_change,:))
    write(u,'(A,ES24.16E3)') 'max_order64_16_to128_32_swing_change=',max_order
    write(u,'(A,ES24.16E3)') 'max_history400_to800_swing_change=',max_history
    write(u,'(A,ES24.16E3)') 'max_history_change_divided_by_analytic_bound=',max_history_ratio
    write(u,'(A,ES24.16E3)') 'max_spatial48_SAS_swing_difference=',max_space48
    write(u,'(A,ES24.16E3)') 'max_spatial32_to48_swing_change=',max_spacechange
    write(u,'(A,ES24.16E3)') 'max_tail_temperature=',maxval(rows(m_tailt,:))
    write(u,'(A,ES24.16E3)') 'max_tail_swing=',maxval(rows(m_tails,:))
    write(u,'(A,ES24.16E3)') 'max_zeta_EM_remainder=',maxval(rows(m_zeta_rem,:))
    write(u,'(A,ES24.16E3)') 'max_zeta_16_to32_refinement=',maxval(rows(m_zeta_ref,:))
    write(u,'(A,I0)') 'failed_depth_acceptance=',failures
    if(failures==0) then
      write(u,'(A)') 'status=PASS'
    else
      write(u,'(A)') 'status=FAIL'
    end if
    write(u,'(A)') 'Extrema are numerical estimates, without interval certification.'
  end subroutine

  subroutine grid_chunk(chunk)
    integer, intent(in) :: chunk
    real(rk) :: L,d,z,v(nmetrics),rows(nmetrics,12)
    real(rk), parameter :: duties(3)=[0.1_rk,0.4_rk,0.7_rk]
    character :: profile
    character(128) :: path
    integer :: u,su,p,k,depth,count,j,metric
    j=chunk-1;L=0.2_rk+4.0_rk*j/35
    write(path,'(A,I0,A)') 'results/periodic_grid_',chunk,'.csv';call new_output(trim(path),u)
    write(path,'(A,I0,A)') 'results/periodic_grid_',chunk,'_summary.txt';call new_output(trim(path),su)
    write(u,'(A)') periodic_header;count=0
    do p=1,2
      profile='U';if(p==2) profile='Q'
      do k=1,3
        d=duties(k)
        do depth=1,2
          z=0.1_rk;if(depth==2) z=0.5_rk
          call periodic_depth(j,L,profile,d,z,v);count=count+1;rows(:,count)=v
          write(u,'(I0,A,ES25.16E3,A,A,A,ES25.16E3,A,ES25.16E3,*(A,ES25.16E3))') &
            j,',',L,',',profile,',',d,',',z,(',',v(metric),metric=1,nmetrics)
          flush(u)
          write(*,'(A,I0,1X,A,2(1X,F4.1),1X,ES12.4E3)') 'periodic depth complete j=',j,profile,d,z,v(m_swing)
        end do
      end do
    end do
    call summarize(su,rows,count);close(u);close(su)
  end subroutine

  subroutine invariant_case(which)
    integer, intent(in) :: which
    type(periodic_context) :: ctx,fullctx
    type(tail_result) :: bound,other,fullbound
    real(rk) :: z,L,d,mean2048,mean4096,mean_error,mean_change,complement,phase,v,c
    real(rk) :: waves(4096),max_zero,fullmin,fullmax
    character :: profile
    character(128) :: path
    integer :: u,k,failures
    select case(which)
    case(1);z=0.1_rk;L=1;profile='U';d=0.1_rk
    case(2);z=0.5_rk;L=1;profile='Q';d=0.4_rk
    case(3);z=0.1_rk;L=0.2_rk;profile='Q';d=0.7_rk
    case(4);z=0.5_rk;L=1;profile='U';d=0.7_rk
    end select
    write(path,'(A,I0,A)') 'results/periodic_invariants_',which,'.txt';call new_output(trim(path),u)
    call periodic_prepare(ctx,z,L,profile);bound=tail_bounds(z,L,profile,d);other=tail_bounds(z,L,profile,1-d)
    do k=0,4095
      waves(k+1)=periodic_value(ctx,period*k/4096,d)
    end do
    mean4096=real(sum(real(waves,qk))/4096,rk)
    mean2048=real(sum(real(waves(1:4096:2),qk))/2048,rk)
    mean_error=abs(mean2048-d*ctx%steady);mean_change=abs(mean4096-mean2048);complement=0
    do k=0,2047
      phase=period*k/2048
      v=periodic_value(ctx,modulo(phase+d*period,period),d)
      c=periodic_value(ctx,phase,1-d);complement=max(complement,abs(c+v-ctx%steady))
    end do
    failures=0
    if(mean_error>bound%temperature+2*mean_change+1e-9_rk) failures=failures+1
    if(complement>bound%temperature+other%temperature+1e-9_rk) failures=failures+1
    write(u,'(A,ES24.16E3)') 'z=',z
    write(u,'(A,ES24.16E3)') 'L=',L
    write(u,'(A,A)') 'profile=',profile
    write(u,'(A,ES24.16E3)') 'd=',d
    write(u,'(A,I0)') 'history=',400
    write(u,'(A,ES24.16E3)') 'mean2048=',mean2048
    write(u,'(A,ES24.16E3)') 'mean4096=',mean4096
    write(u,'(A,ES24.16E3)') 'steady=',ctx%steady
    write(u,'(A,ES24.16E3)') 'mean_error=',mean_error
    write(u,'(A,ES24.16E3)') 'mean_resolution_change=',mean_change
    write(u,'(A,ES24.16E3)') 'mean_temperature_tail_bound=',bound%temperature
    write(u,'(A,ES24.16E3)') 'complement_error=',complement
    write(u,'(A,ES24.16E3)') 'complement_temperature_tail_bound=',bound%temperature+other%temperature
    if(which==1) then
      call periodic_prepare(fullctx,0.5_rk,1.0_rk,'U');max_zero=0;fullmin=huge(1.0_rk);fullmax=-huge(1.0_rk)
      do k=0,256
        phase=period*k/256;max_zero=max(max_zero,abs(periodic_value(fullctx,phase,0.0_rk)))
        v=periodic_value(fullctx,phase,1.0_rk);fullmin=min(fullmin,v);fullmax=max(fullmax,v)
      end do
      fullbound=tail_bounds(0.5_rk,1.0_rk,'U',1.0_rk)
      if(max_zero>=1e-12_rk.or.fullmax-fullmin>fullbound%swing+1e-9_rk) failures=failures+1
      write(u,'(A,ES24.16E3)') 'zero_duty_max=',max_zero
      write(u,'(A,ES24.16E3)') 'full_duty_finite_history_swing=',fullmax-fullmin
      write(u,'(A,ES24.16E3)') 'full_duty_swing_tail_bound=',fullbound%swing
    end if
    write(u,'(A,I0)') 'failures=',failures
    if(failures==0) then
      write(u,'(A)') 'status=PASS'
    else
      write(u,'(A)') 'status=FAIL'
    end if
    close(u)
    if(failures/=0) error stop 'periodic invariant acceptance failed'
  end subroutine

  subroutine aggregate()
    real(rk) :: allrows(nmetrics,96),v(nmetrics),L,d,z,geometry(3,48)
    real(rk) :: read_L,read_d,read_z
    character :: profile,read_profile
    character(512) :: path
    character(4096) :: line
    integer :: chunk,u,out,su,j,read_j,k,row,design_count,ios,depth,metric
    logical :: all_chunks_pass
    all_chunks_pass=.true.;design_count=0
    read_j=-1;read_L=0;read_d=0;read_z=0;read_profile='?'
    do chunk=1,8
      write(path,'(A,I0,A)') 'results/periodic_grid_',chunk,'_summary.txt'
      open(newunit=u,file=trim(path),status='old',action='read')
      do
        read(u,'(A)',iostat=ios) line
        if(ios/=0) exit
        if(trim(line)=='status=FAIL') all_chunks_pass=.false.
      end do
      close(u)
      write(path,'(A,I0,A)') 'results/periodic_grid_',chunk,'.csv'
      open(newunit=u,file=trim(path),status='old',action='read');read(u,'(A)') line
      if(trim(line)/=periodic_header) error stop 'periodic aggregate: CSV header mismatch'
      do k=1,12
        read(u,'(A)',iostat=ios) line
        if(ios/=0) error stop 'periodic aggregate: incomplete chunk'
        read(line,*,iostat=ios) j,L,profile,d,z,v
        if(ios/=0.or.j/=chunk-1) error stop 'periodic aggregate: invalid chunk row'
        if(v(m_selected)/=merge(1.0_rk,0.0_rk, &
          j==7.or.(j==0.and.profile=='Q'.and.abs(d-0.1_rk)<1e-14_rk))) &
          error stop 'periodic aggregate: fixed refinement membership mismatch'
        row=(chunk-1)*12+k;allrows(:,row)=v
        if(mod(k,2)==1) then
          read_j=j;read_L=L;read_profile=profile;read_d=d;read_z=z
        else
          if(j/=read_j.or.L/=read_L.or.profile/=read_profile.or.d/=read_d.or. &
            read_z/=0.1_rk.or.z/=0.5_rk) error stop 'periodic aggregate: invalid depth pair'
          design_count=design_count+1
          geometry(:,design_count)=[real(j,rk),L,d]
        end if
      end do
      read(u,'(A)',iostat=ios) line
      if(ios==0) error stop 'periodic aggregate: extra chunk rows'
      close(u)
    end do
    if(design_count/=48) error stop 'periodic aggregate: wrong design count'
    if(count(allrows(m_selected,:)==1)/=14)error stop 'fixed refinement membership mismatch'
    call new_output('results/periodic_refinement_union.csv',out)
    write(out,'(A)')periodic_header
    do k=1,48
      j=nint(geometry(1,k));L=geometry(2,k);d=geometry(3,k)
      profile='U';if(mod(k-1,6)>=3)profile='Q'
      do depth=1,2
        row=2*(k-1)+depth
        if(allrows(m_selected,row)/=1)cycle
        z=.1_rk;if(depth==2)z=.5_rk
        write(out,'(I0,A,ES25.16E3,A,A,A,ES25.16E3,A,ES25.16E3,*(A,ES25.16E3))') &
          j,',',L,',',profile,',',d,',',z,(',',allrows(metric,row),metric=1,nmetrics)
      end do
    end do
    close(out)
    call new_output('results/periodic_summary.txt',su)
    call summarize(su,allrows,96)
    write(su,'(A,I0)') 'design_count=',design_count
    write(su,'(A,I0)') 'original_control_refinement_depths=',14
    write(su,'(A,I0)') 'refinement_union_depths=',count(allrows(m_selected,:)==1)
    write(su,'(A,ES24.16E3)') 'endpoint_search_correction=',allrows(m_endpoint_correction,93)
    write(su,'(A,L1)') 'all_chunk_acceptance=',all_chunks_pass
    close(su)
    if(.not.all_chunks_pass) error stop 'periodic aggregate acceptance failed'
  end subroutine
end program
