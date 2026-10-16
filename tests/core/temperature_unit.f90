program unit_driver
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use temperature_common
   use prepared_sas
   use prepared_space
   implicit none
   type(patch_source_t) :: source, other, empty
   type(prepared_sas_t) :: ps, ps_other
   type(prepared_space_t) :: pd, pd_other
   type(gauss_rule_t) :: rule, r32, r48, r12, r24, bad
   type(query_work_t) :: ws, wd, wd_deep
   real(rk) :: c(0:2,0:2), d(0:2,0:2), query(4), moved(4)
   real(rk) :: vs, vd, v1, v2, v3, parity(2)
   integer :: status, ss, sd, n, k, checks
   integer, parameter :: orders(6) = [1,2,8,32,48,72]
   checks = 0
   do k = 1, size(orders)
      n = orders(k)
      call init_rule(rule,n,status)
      call require(status == ok,'GL construction')
      call close(sum(rule%w),2.0_rk,2.0e-13_rk,'GL mass')
      call close(sum(rule%w*rule%x),0.0_rk,2.0e-13_rk,'GL odd moment')
      if (n >= 2) call close(sum(rule%w*rule%x**2),2.0_rk/3.0_rk,2.0e-13_rk,'GL quadratic')
      if (n >= 3) call close(sum(rule%w*rule%x**4),2.0_rk/5.0_rk,2.0e-13_rk,'GL quartic')
   end do
   call init_rule(bad,0,status)
   call require(status == invalid,'reject zero order')
   call init_rule(r32,32,status)
   call require(status == ok,'rule32')
   call init_rule(r48,48,status)
   call require(status == ok,'rule48')
   call init_rule(r12,12,status)
   call init_rule(r24,24,status)

   c = 0.0_rk
   c(0,0) = 1.0_rk
   call setup(c,1.0_rk,2.0_rk,1.0_rk,source,ps,pd)
   query = [0.2_rk,-0.3_rk,0.1_rk,1.0_rk]
   vs = evaluate_sas(ps,query,r32,ss,ws)
   vd = evaluate_space(pd,query,r48,r48,sd,wd)
   call require(ss == ok .and. sd == ok,'uniform evaluation')
   call require(vs > 0.0_rk .and. vs < 1.0_rk,'uniform bound')
   call close(vs,vd,3.0e-12_rk,'uniform routes')
   call require(ws%v_nodes > 0 .and. ws%x_nodes == 0 .and. ws%y_nodes == 0,'analytic axes')
   call require(wd%kernel_evals == wd%x_nodes*wd%y_nodes,'actual tensor work')

   c = 0.0_rk
   c(0,0) = 1.0_rk/1.55_rk
   c(0,1) = 0.2_rk/1.55_rk
   c(1,0) = 0.15_rk/1.55_rk
   c(1,1) = 0.2_rk/1.55_rk
   call setup(c,1.0_rk,0.75_rk,1.0_rk,source,ps,pd)
   query = [0.2_rk,-0.3_rk,0.1_rk,0.75_rk]
   vs = evaluate_sas(ps,query,r48,ss)
   vd = evaluate_space(pd,query,r48,r48,sd)
   call close(vs,vd,3.0e-12_rk,'mixed odd source')
   call setup(c,1.0_rk,0.75_rk,2.0_rk,other,ps_other,pd_other)
   moved = query; moved(4) = query(4)/2.0_rk
   v1 = evaluate_sas(ps_other,moved,r48,ss)
   v2 = evaluate_space(pd_other,moved,r48,r48,sd)
   call close(v1,vs,1.0e-12_rk,'diffusivity SAS')
   call close(v2,vd,1.0e-12_rk,'diffusivity space')

   call init_source(other,c,1.0_rk,0.75_rk,1.0_rk,status,0.3_rk,-0.2_rk)
   call require(status == ok,'translated source')
   call prepare_sas(other,ps_other,status)
   call prepare_space(other,pd_other,status)
   moved = query; moved(1:2) = query(1:2)+[0.3_rk,-0.2_rk]
   v1 = evaluate_sas(ps_other,moved,r48,ss)
   v2 = evaluate_space(pd_other,moved,r48,r48,sd)
   call close(v1,vs,1.0e-12_rk,'translation SAS')
   call close(v2,vd,1.0e-12_rk,'translation space')

   call setup(c,3.0_rk,2.25_rk,1.0_rk,other,ps_other,pd_other)
   moved(1:3) = 3.0_rk*query(1:3); moved(4) = 9.0_rk*query(4)
   v1 = evaluate_sas(ps_other,moved,r48,ss)
   v2 = evaluate_space(pd_other,moved,r48,r48,sd)
   call close(v1,vs,2.0e-12_rk,'geometry scaling SAS')
   call close(v2,vd,2.0e-12_rk,'geometry scaling space')

   d = 0.0_rk
   d(0,2) = -0.3_rk; d(2,1) = 0.13_rk; d(2,2) = -0.2_rk
   v1 = vs
   call setup(d,1.0_rk,0.75_rk,1.0_rk,other,ps_other,pd_other)
   v2 = evaluate_sas(ps_other,query,r48,ss)
   call setup(c+d,1.0_rk,0.75_rk,1.0_rk,other,ps_other,pd_other)
   v3 = evaluate_sas(ps_other,query,r48,ss)
   call close(v3,v1+v2,1.0e-12_rk,'signed linearity SAS')
   v1 = vd
   call setup(d,1.0_rk,0.75_rk,1.0_rk,other,ps_other,pd_other)
   v2 = evaluate_space(pd_other,query,r48,r48,sd)
   call setup(c+d,1.0_rk,0.75_rk,1.0_rk,other,ps_other,pd_other)
   v3 = evaluate_space(pd_other,query,r48,r48,sd)
   call close(v3,v1+v2,1.0e-12_rk,'signed linearity space')

   c = 0.0_rk; c(0,1) = 1.0_rk
   call setup(c,1.0_rk,0.75_rk,1.0_rk,source,ps,pd)
   query = [0.0_rk,0.0_rk,0.1_rk,1.0_rk]
   call close(evaluate_sas(ps,query,r48,ss),0.0_rk,1.0e-13_rk,'odd central zero SAS')
   call close(evaluate_space(pd,query,r48,r48,sd),0.0_rk,1.0e-13_rk,'odd central zero space')
   query(1) = 0.6_rk
   parity(1) = evaluate_sas(ps,query,r48,ss)
   parity(2) = evaluate_space(pd,query,r48,r48,sd)
   query(1) = -0.6_rk
   v1 = evaluate_sas(ps,query,r48,ss)
   v2 = evaluate_space(pd,query,r48,r48,sd)
   call require(v1 < 0.0_rk .and. v2 < 0.0_rk,'negative response not clipped')
   call close(v1,-parity(1),1.0e-12_rk,'odd reflection SAS')
   call close(v2,-parity(2),1.0e-12_rk,'odd reflection space')

   c = 0.0_rk
   call setup(c,1.0_rk,0.75_rk,1.0_rk,source,ps,pd)
   call close(evaluate_sas(ps,query,r48,ss),0.0_rk,0.0_rk,'zero source SAS')
   call close(evaluate_space(pd,query,r48,r48,sd),0.0_rk,0.0_rk,'zero source space')
   c(0,0) = 1.0_rk
   call setup(c,1.0_rk,0.75_rk,1.0_rk,source,ps,pd)
   query = [0.2_rk,-0.3_rk,0.1_rk,1.0_rk]
   vd = evaluate_space(pd,query,r12,r24,sd,wd)
   call require(sd == ok .and. mod(wd%x_nodes,12) == 0 .and. mod(wd%y_nodes,24) == 0,'independent axis orders')
   moved = query; moved(3) = 0.8_rk
   vd = evaluate_space(pd,moved,r12,r24,sd,wd_deep)
   call require(wd%x_nodes /= wd_deep%x_nodes .or. wd%y_nodes /= wd_deep%y_nodes,'query-dependent partition')

   moved = query; moved(3) = 0.0_rk
   vs = evaluate_sas(ps,moved,r48,ss)
   vd = evaluate_space(pd,moved,r48,r48,sd)
   call require(ss == invalid .and. sd == invalid,'reject surface query')
   call require(.not.ieee_is_finite(vs) .and. .not.ieee_is_finite(vd),'failed query NaN')
   moved = query; moved(4) = 0.0_rk
   vs = evaluate_sas(ps,moved,r48,ss)
   vd = evaluate_space(pd,moved,r48,r48,sd)
   call require(ss == invalid .and. sd == invalid,'reject zero time')
   moved = query; moved(1) = quiet_nan()
   vs = evaluate_sas(ps,moved,r48,ss)
   vd = evaluate_space(pd,moved,r48,r48,sd)
   call require(ss == invalid .and. sd == invalid,'reject NaN coordinate')
   vs = evaluate_sas(ps,query,bad,ss)
   vd = evaluate_space(pd,query,bad,r48,sd)
   call require(ss == invalid .and. sd == invalid,'reject unready rules')
   call init_source(empty,c,-1.0_rk,1.0_rk,1.0_rk,status)
   call require(status == invalid .and. .not.empty%ready,'reject geometry')
   c(0,0) = quiet_nan()
   call init_source(empty,c,1.0_rk,1.0_rk,1.0_rk,status)
   call require(status == invalid .and. .not.empty%ready,'reject NaN coefficient')
   c = 0.0_rk; c(0,0) = 1.0_rk
   call init_source(empty,c,1.0_rk,1.0_rk,1.0_rk,status)
   empty%alpha = quiet_nan()
   call prepare_sas(empty,ps_other,ss)
   call prepare_space(empty,pd_other,sd)
   call require(ss == invalid .and. sd == invalid,'revalidate externally corrupted source')

   ! Public coefficients may change after init_source; derived degrees must not
   ! be copied back from the caller after either preparation has normalized them.
   c = 0.0_rk
   call init_source(source,c,1.0_rk,0.75_rk,1.0_rk,status)
   call require(status == ok,'mutable coefficient source initialization')
   source%c(0,1) = 1.0_rk
   call prepare_sas(source,ps,ss)
   call prepare_space(source,pd,sd)
   call require(ss == ok .and. sd == ok,'recompute degrees after coefficient change')
   call setup(source%c,1.0_rk,0.75_rk,1.0_rk,other,ps_other,pd_other)
   query = [0.9_rk,0.675_rk,0.01875_rk,0.5625_rk]
   vs = evaluate_sas(ps,query,r48,ss)
   vd = evaluate_space(pd,query,r48,r48,sd)
   call require(ss == ok .and. sd == ok,'mutated coefficient evaluation')
   call require(vs > 0.0_rk .and. vd > 0.0_rk,'retain new odd term')
   v1 = evaluate_sas(ps_other,query,r48,ss)
   v2 = evaluate_space(pd_other,query,r48,r48,sd)
   call require(ss == ok .and. sd == ok,'fresh coefficient evaluation')
   call close(vs,v1,1.0e-12_rk,'mutated versus fresh SAS')
   call close(vd,v2,1.0e-12_rk,'mutated versus fresh space')
   call close(vs,vd,3.0e-12_rk,'mutated coefficient routes')
   write(*,'(a,i0)') 'UNIT_PASS,',checks
contains
   subroutine setup(c,a,b,alpha,source,sas,space)
      real(rk), intent(in) :: c(0:2,0:2),a,b,alpha
      type(patch_source_t), intent(out) :: source
      type(prepared_sas_t), intent(out) :: sas
      type(prepared_space_t), intent(out) :: space
      integer :: st
      call init_source(source,c,a,b,alpha,st)
      call require(st == ok,'source preparation')
      call prepare_sas(source,sas,st)
      call require(st == ok,'SAS preparation')
      call prepare_space(source,space,st)
      call require(st == ok,'space preparation')
   end subroutine
   subroutine require(condition,label)
      logical, intent(in) :: condition
      character(*), intent(in) :: label
      checks = checks+1
      if (.not.condition) then
         write(*,'(a)') 'UNIT_FAIL: '//label
         error stop 1
      end if
   end subroutine
   subroutine close(a,b,tol,label)
      real(rk), intent(in) :: a,b,tol
      character(*), intent(in) :: label
      call require(ieee_is_finite(a) .and. ieee_is_finite(b),label//' finite')
      call require(abs(a-b) <= tol*(1.0_rk+abs(b)),label)
   end subroutine
end program
