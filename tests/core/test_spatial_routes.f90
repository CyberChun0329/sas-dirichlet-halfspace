program test_spatial_routes
   use replay_types
   use spatial_routes
   use gaussian_axes, only: monomial_axis
   use temperature_common, only: patch_source_t,init_source
   use prepared_space, only: prepared_space_t,prepare_space,evaluate_space
   implicit none
   type(source_t) :: s,s2
   type(spatial_t) :: p,p2
   type(patch_source_t) :: oldsource
   type(prepared_space_t) :: oldprepared
   type(gauss_rule_t) :: r24,r32,r48
   real(rk) :: c(0:2,0:2),hc(0:32,0:32),rect(0:4,0:6),query(4),qq(4,2)
   real(rk) :: v,v2,v3,d,dp,dm,steady,expected,h
   real(rk) :: moments6(0:6),moments2(0:2)
   real(rk),allocatable :: w(:,:,:)
   integer :: st,failures,i,j,obs
   failures=0
   call init_rule(r24,24,st);call assert_status('rule24',st)
   call init_rule(r32,32,st);call assert_status('rule32',st)
   call init_rule(r48,48,st);call assert_status('rule48',st)
   call profile_source('N2',1.0_rk,2.0_rk,s,st);call assert_status('N2',st)
   call prepare_spatial(s,p,st);call assert_status('prepare N2',st)
   c=s%c
   call init_source(oldsource,c,1.0_rk,2.0_rk,1.0_rk,st)
   call prepare_space(oldsource,oldprepared,st);call assert_status('old prepare',st)
   query=[.3_rk,.4_rk,.2_rk,.75_rk]
   v=evaluate_spatial(p,query,r24,r32,OBS_T,st);call assert_status('2D N2',st)
   v2=evaluate_space(oldprepared,query,r24,r32,st);call assert_status('direct N2',st)
   call close('direct low-order actual2D',v,v2,2e-14_rk)
   v2=evaluate_shared(p,query,r48,OBS_T,st);call assert_status('shared N2',st)
   call close('shared/2D N2',v,v2,2e-10_rk)
   call check(v>0 .and. v<1.1_rk,'positive bounded temperature')
   d=evaluate_spatial(p,query,r32,r32,OBS_D,st);call assert_status('deficit',st)
   call check(d>=0,'nonnegative accumulated deficit')
   v2=evaluate_shared(p,query,r48,OBS_D,st);call assert_status('shared deficit',st)
   call close('shared/2D physical deficit',d,v2,2e-10_rk)
   h=1e-4_rk
   query(4)=.75_rk+h;dp=evaluate_spatial(p,query,r32,r32,OBS_D,st)
   query(4)=.75_rk-h;dm=evaluate_spatial(p,query,r32,r32,OBS_D,st)
   query(4)=.75_rk
   steady=evaluate_steady(p,query,r32,r32,st);call assert_status('steady N2',st)
   call close('dD/dtp=steady-T',(dp-dm)/(2*h),steady-v,3e-9_rk)
   do obs=OBS_T,OBS_D
      v=evaluate_spatial(p,query,r32,r32,obs,st)
      v2=evaluate_spatial_complementary(p,query,r32,r32,obs,st)
      call assert_status('complementary',st)
      call close('fresh steady minus correction',v,v2,2e-14_rk)
   end do
   query(4)=0
   v=evaluate_spatial(p,query,r24,r24,OBS_T,st)
   call check(st==invalid .and. .not.ieee_is_finite(v),'t0 explicitly unsupported')
   query(4)=.75_rk
   s2=s;s2%alpha=2
   call prepare_spatial(s2,p2,st);call assert_status('alpha2',st)
   query(4)=.375_rk
   v=evaluate_spatial(p2,query,r32,r32,OBS_T,st)
   d=evaluate_spatial(p2,query,r32,r32,OBS_D,st)
   query(4)=.75_rk
   v2=evaluate_spatial(p,query,r32,r32,OBS_T,st)
   v3=evaluate_spatial(p,query,r32,r32,OBS_D,st)
   call close('alpha*time temperature',v,v2,2e-14_rk)
   call close('physical time deficit',d,v3/2,2e-14_rk)
   s2=s;s2%cx=3;s2%cy=-4
   call prepare_spatial(s2,p2,st);call assert_status('translated source',st)
   query(1)=query(1)+3;query(2)=query(2)-4
   v=evaluate_spatial(p2,query,r32,r32,OBS_T,st)
   call close('center translation',v,v2,2e-14_rk)
   query=[0.0_rk,0.0_rk,.2_rk,1e8_rk]
   call profile_source('U',1.0_rk,2.0_rk,s,st);call prepare_spatial(s,p,st)
   steady=evaluate_steady(p,query,r32,r32,st)
   expected=2/pi*atan(2/(query(3)*sqrt(1+4+query(3)**2)))
   call close('independent rectangle solid angle',steady,expected,2e-13_rk)
   v=evaluate_spatial(p,query,r32,r32,OBS_T,st)
   call close('long-time steady limit',v,steady,2e-12_rk)
   d=evaluate_spatial(p,query,r32,r32,OBS_D,st)
   v2=evaluate_shared(p,query,r48,OBS_D,st)
   call close('stable small-s deficit',d,v2,3e-10_rk)
   query=[.4_rk,.8_rk,.3_rk,.2_rk]
   call init_gaussian(s,[2.0_rk,4.0_rk],1.0_rk,2.0_rk,1.0_rk,st)
   call prepare_spatial(s,p,st)
   do obs=OBS_T,OBS_D
      v=evaluate_spatial(p,query,r32,r32,obs,st);call assert_status('true Gaussian 2D',st)
      v2=evaluate_shared(p,query,r48,obs,st);call assert_status('true Gaussian shared',st)
      call close('true Gaussian shared/2D',v,v2,2e-10_rk)
   end do
   call monomial_axis(3.0_rk,.4_rk,1.0_rk,6,moments6,st)
   call assert_status('high monomial axis',st)
   call monomial_axis(3.0_rk,.4_rk,1.0_rk,2,moments2,st)
   do i=0,2
      call close('high axis preserves low physical moments',moments6(i),moments2(i),2e-14_rk)
   end do
   ! Complete non-square coefficient matrix, including all high-order edges.
   do i=0,6
      do j=0,4
         rect(j,i)=(-1.0_rk)**(i+j)/(1+i+j)**3
      end do
   end do
   call init_poly(s,rect,1.0_rk,2.0_rk,1.0_rk,st);call prepare_spatial(s,p,st)
   call check(spatial_precision(p)==128,'high-degree intermediates visible')
   v=evaluate_spatial(p,query,r32,r48,OBS_T,st)
   v2=evaluate_shared(p,query,r48,OBS_T,st);call assert_status('rectangular shared',st)
   call close('complete rectangular tensor',v,v2,2e-10_rk)
   qq(:,1)=query;qq(:,2)=[-.7_rk,1.5_rk,.4_rk,1.0_rk]
   call space_operator(p,qq,r32,r48,OBS_T,w,st);call assert_status('polynomial operator',st)
   do i=1,2
      v=evaluate_spatial(p,qq(:,i),r32,r48,OBS_T,st)
      call close('polynomial operator applied',sum(w(:,:,i)*s%c),v,2e-14_rk)
   end do
   hc=0;hc(0,0)=1;hc(32,32)=.1_rk
   call init_poly(s,hc,1.0_rk,2.0_rk,1.0_rk,st);call prepare_spatial(s,p,st)
   v=evaluate_spatial(p,query,r48,r48,OBS_T,st)
   v2=evaluate_shared(p,query,r48,OBS_T,st);call assert_status('M32 shared',st)
   call close('M32 full tensor precision',v,v2,2e-8_rk)
   if(failures/=0) then
      write(*,'(a,i0)') 'FAIL spatial checks: ',failures
      stop 1
   end if
   write(*,'(a)') 'PASS spatial routes functional checks'
contains
   subroutine check(condition,label)
      logical,intent(in) :: condition
      character(*),intent(in) :: label
      if(.not.condition) then
         failures=failures+1;write(*,'(a)') 'FAIL '//label
      end if
   end subroutine
   subroutine assert_status(label,status)
      character(*),intent(in) :: label
      integer,intent(in) :: status
      call check(status==ok,label//' status')
   end subroutine
   subroutine close(label,a,b,tol)
      character(*),intent(in) :: label
      real(rk),intent(in) :: a,b,tol
      logical :: pass
      pass=ieee_is_finite(a).and.ieee_is_finite(b)
      if(pass) pass=abs(a-b)<=tol*max(1.0_rk,abs(a),abs(b))
      call check(pass,label)
      if(.not.pass) write(*,'(2es25.16,a,es12.3)') a,b,' tolerance ',tol
   end subroutine
end program
