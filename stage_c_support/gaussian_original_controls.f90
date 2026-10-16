program gaussian_original_controls
   use replay_types
   use source_projection, only: project_source
   use sas_general, only: sas_general_t,prepare_sas_general,evaluate_sas_general
   use spatial_routes, only: spatial_t,prepare_spatial,evaluate_spatial,evaluate_shared
   implicit none
   type(source_t) :: original,projected,scaled,scaled_poly
   type(sas_general_t) :: base_sas,scaled_sas
   type(spatial_t) :: base_space,scaled_space,projected_space
   type(gauss_rule_t) :: rule
   integer,parameter :: degrees(6)=[2,4,6,8,10,12]
   real(rk),parameter :: tolerance=2.e-10_rk
   real(rk) :: q(4),qq(4),sv,pv,gv,vv,gg,scale,difference,max_difference
   integer :: unit,summary,k,t,setup_status,holds,rows
   logical :: identity,accepted
   character(1024) :: outdir
   character(64) :: stage
   ! Exact original eight control definitions, not candidate or timing runs.
   ! Uses the imported numerical modules; both output files must be new.
   if(command_argument_count()/=1) error stop 'one existing new output directory required'
   call get_command_argument(1,outdir)
   open(newunit=unit,file=trim(outdir)//'/controls.csv',status='new',action='write')
   open(newunit=summary,file=trim(outdir)//'/summary.txt',status='new',action='write')
   write(unit,'(a)') 'control,M,n,scale,SAS,independent_polynomial,Gaussian,SAS_original,Gaussian_original,difference,coefficient_identity,accepted,status,stage'
   holds=0;rows=0;max_difference=0
   call init_rule(rule,64,t)
   if(t==ok) call init_gaussian(original,[2._rk,2._rk],1._rk,1._rk,1._rk,t)
   setup_status=t
   q=[.5_rk,.5_rk,1._rk,1._rk]
   do k=1,6
      t=setup_status;sv=quiet_nan();pv=quiet_nan();stage='SHARED_PREPARATION'
      if(t==ok) then
         stage='PROJECT_SOURCE';call project_source(original,degrees(k),projected,t)
      end if
      if(t==ok) then
         stage='SAS_PREPARATION';call prepare_sas_general(projected,base_sas,t)
      end if
      if(t==ok) then
         stage='POLYNOMIAL_SPACE_PREPARATION';call prepare_spatial(projected,projected_space,t)
      end if
      if(t==ok) then
         stage='SAS_EVALUATION';sv=evaluate_sas_general(base_sas,q,rule,OBS_T,t)
      end if
      if(t==ok) then
         stage='POLYNOMIAL_EVALUATION';pv=evaluate_spatial(projected_space,q,rule,rule,OBS_T,t)
      end if
      difference=quiet_nan();accepted=.false.
      if(all(ieee_is_finite([sv,pv]))) then
         difference=abs(sv-pv);accepted=t==ok.and.difference<=tolerance
      end if
      call count_result(accepted,difference)
      write(unit,'(*(g0,:,","))') 'CaseII',degrees(k),64,1._rk,sv,pv,0,0,0, &
         difference,t==ok,accepted,t,trim(stage)
   end do
   call init_gaussian(original,[2._rk,2._rk],1._rk,2._rk,1._rk,t)
   if(t==ok) call project_source(original,8,projected,t)
   if(t==ok) call prepare_sas_general(projected,base_sas,t)
   if(t==ok) call prepare_spatial(original,base_space,t)
   if(t==ok) call init_rule(rule,32,t)
   q=[.9_rk,1.8_rk,.1_rk,.1_rk];sv=quiet_nan();gv=quiet_nan()
   if(t==ok) sv=evaluate_sas_general(base_sas,q,rule,OBS_T,t)
   if(t==ok) gv=evaluate_shared(base_space,q,rule,OBS_T,t)
   setup_status=t
   do k=1,2
      scale=merge(.1_rk,10._rk,k==1)
      t=setup_status;vv=quiet_nan();gg=quiet_nan();stage='SIMILARITY_BASE_PREPARATION'
      if(t==ok) then
         stage='SCALED_SOURCE_PREPARATION';call init_gaussian(scaled,[2._rk,2._rk],scale,2*scale,1._rk,t)
      end if
      if(t==ok) then
         stage='PROJECT_SOURCE';call project_source(scaled,8,scaled_poly,t)
      end if
      if(t==ok) then
         stage='SAS_PREPARATION';call prepare_sas_general(scaled_poly,scaled_sas,t)
      end if
      if(t==ok) then
         stage='GAUSSIAN_SPACE_PREPARATION';call prepare_spatial(scaled,scaled_space,t)
      end if
      qq=q;qq(:3)=q(:3)*scale;qq(4)=q(4)*scale*scale
      if(t==ok) then
         stage='SAS_EVALUATION';vv=evaluate_sas_general(scaled_sas,qq,rule,OBS_T,t)
      end if
      if(t==ok) then
         stage='GAUSSIAN_EVALUATION';gg=evaluate_shared(scaled_space,qq,rule,OBS_T,t)
      end if
      identity=.false.;accepted=.false.;difference=quiet_nan()
      if(t==ok) identity=all(projected%c==scaled_poly%c)
      if(all(ieee_is_finite([vv,sv,gg,gv]))) then
         difference=max(abs(vv-sv),abs(gg-gv))
         accepted=t==ok.and.difference<=tolerance.and.identity
      end if
      call count_result(accepted,difference)
      write(unit,'(*(g0,:,","))') 'M8_similarity',8,32,scale,vv,0,gg,sv,gv, &
         difference,identity,accepted,t,trim(stage)
   end do
   close(unit)
   write(summary,'(a,i0)') 'original_control_rows=',rows
   write(summary,'(a,i0)') 'HOLD_controls=',holds
   write(summary,'(a,es26.17e3)') 'fixed_absolute_tolerance=',tolerance
   write(summary,'(a,es26.17e3)') 'maximum_finite_difference=',max_difference
   write(summary,'(a)') 'scope=original six n64 polynomial controls and two n32 dual-route scale/identity controls'
   write(summary,'(a)') 'performance_timing=not_executed'
   close(summary)
   if(holds/=0.or.rows/=8) stop 1
   print '(a)', 'PASS original eight Gaussian controls'
contains
   subroutine count_result(passed,delta)
      logical,intent(in) :: passed
      real(rk),intent(in) :: delta
      rows=rows+1
      if(.not.passed) holds=holds+1
      if(ieee_is_finite(delta)) max_difference=max(max_difference,delta)
   end subroutine
end program
