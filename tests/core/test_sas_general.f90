program test_sas_general
   use replay_types
   use gaussian_axes
   use sas_general
   use temperature_common, only: patch_source_t,init_source
   use prepared_sas, only: prepared_sas_t,prepare_sas,evaluate_sas
   use spatial_routes,only: spatial_t,prepare_spatial,evaluate_spatial,evaluate_steady
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(gauss_rule_t) :: reference,rule,bad_rule
   type(source_t) :: source,other,scaled,unsupported
   type(sas_general_t) :: prepared,second,third,unprepared
   type(patch_source_t) :: direct_source
   type(prepared_sas_t) :: direct
   real(rk) :: values(0:32),ref(0:32),negative(0:32),monomials(0:32),bound(0:32)
   real(qk) :: high(0:32)
   real(rk) :: u,x,a,b,term,xx,w,v1,v2,v3,expected,profile
   real(rk) :: query(4),qscaled(4),C(0:4,0:6),Cother(0:4,0:6),Wop(0:4,0:6),c32(0:32,0:32)
   real(rk),parameter :: us(7)=[.001_rk,.1_rk,.5_rk,1.0_rk,3.0_rk,12.0_rk,30.0_rk]
   real(rk),parameter :: xs(7)=[0.0_rk,.3_rk,1.0_rk,1.7_rk,-2.4_rk,20.0_rk,2.e5_rk]
   character(2),parameter :: names(6)=[character(2)::'U','Q','N1','N2','B1','B2']
   integer :: status,failures,checks,j,k,iu,ix,n,p,q,observable
   failures=0;checks=0
   call init_rule(reference,512,status)
   call require(status==ok,'external axis GL rule')
   call init_rule(rule,32,status)
   call require(status==ok,'time GL rule')

   ! Exact zero-u values include the nonzero even Chebyshev moments.
   a=1.7_rk
   call cheb_axis(0.0_rk,.4_rk,a,32,values,status,high,bound)
   call require(status==ok,'zero-u Chebyshev status')
   do j=0,32
      expected=0
      if(mod(j,2)==0) expected=2*a/real(1-j*j,rk)
      call close(values(j),expected,4.e-16_rk,'zero-u Chebyshev moment')
   end do
   call monomial_axis(0.0_rk,.4_rk,a,32,monomials,status,high)
   call require(status==ok,'zero-u physical monomial status')
   do j=0,32
      expected=0
      if(mod(j,2)==0) expected=2*a**(j+1)/real(j+1,rk)
      call close(monomials(j),expected,3.e-15_rk,'zero-u physical monomial')
   end do
   call require(axis_precision(2)==64 .and. axis_precision(32)==128,'visible axis precision')

   ! Isolated reference quadrature is exclusively a test oracle. Production
   ! axes never call it. Each tensor Chebyshev degree 0..32 is checked.
   a=1.2_rk
   do iu=1,size(us)
      u=us(iu)
      do ix=1,size(xs)
         x=xs(ix)
         call axis_reference(u,x,a,ref)
         call cheb_axis(u,x,a,32,values,status,high,bound)
         call require(status==ok,'high axis status')
         do j=0,32
            call close(values(j),ref(j),2.e-13_rk,'Chebyshev axis versus external GL')
            call require(bound(j)>=0 .and. ieee_is_finite(bound(j)),'series truncation bound')
         end do
         call cheb_axis(u,-x,a,32,negative,status)
         call require(status==ok,'reflected axis status')
         do j=0,32
            call close(negative(j),(-1.0_rk)**j*values(j),4.e-15_rk,'axis reflection parity')
         end do
      end do
   end do
   do n=0,32
      call cheb_axis(.9_rk,.4_rk,a,n,values(:n),status)
      call require(status==ok,'individual degree status')
      call axis_reference(.9_rk,.4_rk,a,ref)
      call close(values(n),ref(n),2.e-13_rk,'individual degrees 0..32')
   end do
   call cheb_axis(2.3_rk,0.0_rk,a,32,values,status)
   do j=1,32,2
      call close(values(j),0.0_rk,2.e-30_rk,'centered odd Chebyshev parity')
   end do

   a=2.3_rk;u=.8_rk;x=3.2_rk
   call monomial_axis(u,x,a,32,monomials,status,high)
   call require(status==ok,'physical monomial outside status')
   ref=0
   do k=1,reference%n
      xx=a*reference%x(k);w=a*reference%w(k)*exp(-(u*(xx-x))**2)
      term=1
      do j=0,32
         ref(j)=ref(j)+w*term;term=term*xx
      end do
   end do
   do j=0,32
      call close(monomials(j),ref(j),3.e-13_rk,'physical monomial versus external GL')
   end do

   a=.7_rk;u=2.4_rk;x=1.3_rk;b=1.8_rk
   profile=gaussian_profile_axis(u,x,a,b,status)
   call require(status==ok,'true Gaussian axis status')
   expected=0
   do k=1,reference%n
      xx=a*reference%x(k)
      expected=expected+a*reference%w(k)*exp(-(u*(xx-x))**2-(b*xx/a)**2)
   end do
   call close(profile,expected,3.e-14_rk,'true finite Gaussian axis versus GL')
   a=1;u=1.e-6_rk;x=1.e6_rk;b=1.e-6_rk
   profile=gaussian_profile_axis(u,x,a,b,status)
   expected=0
   do k=1,reference%n
      xx=a*reference%x(k)
      expected=expected+a*reference%w(k)*exp(-(u*(xx-x))**2-(b*xx/a)**2)
   end do
   call require(status==ok,'narrow far Gaussian profile status')
   call close(profile,expected,3.e-14_rk,'narrow far Gaussian profile avoids erf cancellation')

   ! Compare the low-degree temperature call with its direct evaluator.
   do n=1,size(names)
      call profile_source(names(n),1.0_rk,2.0_rk,source,status)
      source%alpha=2.7_rk;source%cx=.4_rk;source%cy=-.8_rk
      call prepare_sas_general(source,prepared,status)
      call require(status==ok,'low source preparation')
      call init_source(direct_source,source%c,source%a,source%b,source%alpha,status,source%cx,source%cy)
      call prepare_sas(direct_source,direct,status)
      do ix=1,3
         query=[real(ix-2,rk)*2.3_rk,.7_rk,.3_rk,.2_rk]
         v1=evaluate_sas_general(prepared,query,rule,OBS_T,status)
         call require(status==ok,'low temperature status')
         v2=evaluate_sas(direct,query,rule,status)
         call require(status==ok,'direct temperature status')
         call require(v1==v2,'low binary64 baseline bit agreement')
         v3=evaluate_sas_general(prepared,query,rule,OBS_D,status)
         call require(status==ok .and. v3>0,'positive profile deficit')
      end do
   end do

   do p=0,6
      do q=0,4
         C(q,p)=(-1.0_rk)**(p+q)/real(1+p+q,rk)**3
      end do
   end do
   call init_poly(source,C,1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(source,prepared,status)
   call init_poly(other,-.7_rk*C,1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(other,second,status)
   call init_poly(scaled,1.3_rk*C,1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(scaled,third,status)
   query=[.4_rk,.8_rk,.3_rk,.2_rk]
   do observable=OBS_T,OBS_D
      v1=evaluate_sas_general(prepared,query,rule,observable,status)
      call require(status==ok,'signed rectangular source status')
      v2=evaluate_sas_general(second,query,rule,observable,status)
      v3=evaluate_sas_general(third,query,rule,observable,status)
      call close(v2,-.7_rk*v1,2.e-14_rk,'signed scaling')
      call close(v3,v1-v2/2.333333333333333333_rk,2.e-14_rk,'signed linearity')
      call sas_features(prepared,query,rule,observable,Wop,status)
      call require(status==ok,'feature status')
      call close(sum(C*Wop),v1,3.e-14_rk,'full rectangular tensor feature contraction')
   end do
   v1=evaluate_sas_general(prepared,query,rule,OBS_T,status)
   v2=evaluate_sas_complementary(prepared,query,rule,status)
   call require(status==ok,'complementary route status')
   call close(v2,v1,2.e-12_rk,'direct versus complementary intervals')

   Cother=0;Cother(0,0)=-.3_rk;Cother(4,5)=.2_rk;Cother(3,6)=-.5_rk
   call init_poly(other,Cother,1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(other,second,status)
   call init_poly(other,C+Cother,1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(other,third,status)
   do observable=OBS_T,OBS_D
      v1=evaluate_sas_general(prepared,query,rule,observable,status)
      v2=evaluate_sas_general(second,query,rule,observable,status)
      v3=evaluate_sas_general(third,query,rule,observable,status)
      call close(v3,v1+v2,3.e-14_rk,'independent signed tensor superposition')
   end do

   scaled=source;scaled%alpha=3.7_rk
   call prepare_sas_general(scaled,second,status)
   qscaled=query;qscaled(4)=query(4)*scaled%alpha
   v1=evaluate_sas_general(second,query,rule,OBS_T,status)
   v2=evaluate_sas_general(prepared,qscaled,rule,OBS_T,status)
   call close(v1,v2,3.e-15_rk,'diffusivity time scaling T')
   v1=evaluate_sas_general(second,query,rule,OBS_D,status)
   v2=evaluate_sas_general(prepared,qscaled,rule,OBS_D,status)/scaled%alpha
   call close(v1,v2,3.e-15_rk,'diffusivity time scaling D')

   c32=0;c32(0,0)=1;c32(32,32)=.1_rk
   call init_poly(source,c32,1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(source,prepared,status)
   query=[2.3_rk,-2.4_rk,.12_rk,.7_rk]
   v1=evaluate_sas_general(prepared,query,rule,OBS_T,status)
   call require(status==ok .and. v1>0,'degree32 positive outside temperature')
   v2=evaluate_sas_general(prepared,query,rule,OBS_D,status)
   call require(status==ok .and. v2>0,'degree32 positive outside deficit')

   call init_gaussian(source,[1.4_rk,2.1_rk],1.0_rk,2.0_rk,1.0_rk,status)
   call prepare_sas_general(source,prepared,status)
   v1=evaluate_sas_general(prepared,query,rule,OBS_T,status)
   call require(status==ok .and. v1>0,'true Gaussian temperature control')
   v2=evaluate_sas_general(prepared,query,rule,OBS_D,status)
   call require(status==ok .and. v2>0,'true Gaussian deficit control')

   call cheb_axis(-1.0_rk,0.0_rk,1.0_rk,32,values,status)
   call require(status==invalid .and. .not.ieee_is_finite(values(0)),'negative-u rejected')
   call monomial_axis(1.0_rk,0.0_rk,0.0_rk,32,values,status)
   call require(status==invalid .and. .not.ieee_is_finite(values(0)),'zero width rejected')
   call monomial_axis(0.0_rk,0.0_rk,1.e20_rk,32,values,status)
   call require(status==numerical_failure .and. .not.ieee_is_finite(values(0)),'unrepresentable physical moment fails')
   call cheb_axis(0.0_rk,0.0_rk,1.e308_rk,0,values(:0),status)
   call require(status==numerical_failure .and. .not.ieee_is_finite(values(0)),'unrepresentable zero-u area fails')
   profile=gaussian_profile_axis(0.0_rk,0.0_rk,1.e308_rk,1.e-6_rk,status)
   call require(status==numerical_failure .and. .not.ieee_is_finite(profile),'unrepresentable Gaussian area fails')
   v1=evaluate_sas_general(unprepared,query,rule,OBS_T,status)
   call require(status==invalid .and. .not.ieee_is_finite(v1),'unprepared rejected')
   v1=evaluate_sas_general(prepared,query,bad_rule,OBS_T,status)
   call require(status==invalid .and. .not.ieee_is_finite(v1),'unready rule rejected')
   v1=evaluate_sas_general(prepared,query,rule,99,status)
   call require(status==invalid .and. .not.ieee_is_finite(v1),'invalid observable rejected')
   query(3)=0
   v1=evaluate_sas_general(prepared,query,rule,OBS_T,status)
   call require(status==invalid .and. .not.ieee_is_finite(v1),'zero depth rejected')
   unsupported=source;unsupported%kind=-1
   call prepare_sas_general(unsupported,prepared,status)
   call require(status==invalid,'unsupported source rejected')

   call low_deficit_checks()

   write(*,'(a,i0,a,i0)') 'sas_general unit checks=',checks,' failures=',failures
   if(failures/=0) error stop 1
contains
   subroutine low_deficit_checks()
      type(source_t) :: s,translated,different,summed
      type(sas_general_t) :: ps,pt,pd,psum
      type(spatial_t) :: space
      type(gauss_rule_t) :: r96
      real(rk) :: q0(4),qt(4),qp(4),csmall(0:2,0:2),wlow(0:2,0:2)
      real(rk) :: ds,dr,dold,base,changed,added,h,derivative,steady,Tstep,dages(4)
      real(rk) :: max_absolute,max_relative,max_operator_difference,derivative_difference
      real(rk),parameter :: ages(2)=[1.0_rk,100.0_rk]
      character(2),parameter :: profiles(3)=[character(2)::'U','Q','N2']
      integer :: ip,ia,point,k,st
      call init_rule(r96,96,st);call require(st==ok,'low deficit independent rule96')
      max_absolute=0;max_relative=0;max_operator_difference=0;derivative_difference=0
      do ip=1,3
         call profile_source(profiles(ip),1.0_rk,2.0_rk,s,st)
         call prepare_sas_general(s,ps,st);call require(st==ok,'low deficit source preparation')
         call prepare_spatial(s,space,st);call require(st==ok,'low deficit independent preparation')
         do ia=1,2
            do point=1,4
               q0=[0.0_rk,0.0_rk,.1_rk,ages(ia)]
               if(point==2) q0(3)=1.5_rk
               if(point==3) q0(:3)=[.9_rk,1.8_rk,.1_rk]
               if(point==4) q0(:3)=[1.2_rk,0.0_rk,.5_rk]
               ds=evaluate_sas_general(ps,q0,rule,OBS_D,st)
               call require(st==ok,'binary64 low deficit response')
               dr=evaluate_spatial(space,q0,r96,r96,OBS_D,st)
               call require(st==ok,'independent 2D low deficit response')
               call low_close(ds,dr,'low D versus independent SPACE2D96')
               max_absolute=max(max_absolute,abs(ds-dr))
               if(abs(dr)>=1.e-6_rk) max_relative=max(max_relative,abs(ds-dr)/abs(dr))
               call sas_features(ps,q0,rule,OBS_D,wlow,st)
               call require(st==ok.and.all(ieee_is_finite(wlow)),'low D feature operator initialized')
               dold=sum(s%c*wlow)
               call low_close(ds,dold,'binary64 low D versus unchanged qk feature diagnostic')
               max_operator_difference=max(max_operator_difference,abs(ds-dold))
            end do
         end do
         q0=[0.0_rk,0.0_rk,.5_rk,1.0_rk]
         base=evaluate_sas_general(ps,q0,rule,OBS_D,st)
         translated=s;translated%cx=.4_rk;translated%cy=-.8_rk
         call prepare_sas_general(translated,pt,st)
         qt=q0;qt(:2)=q0(:2)+[translated%cx,translated%cy]
         changed=evaluate_sas_general(pt,qt,rule,OBS_D,st)
         call require(st==ok,'translated low D status')
         call low_close(changed,base,'low D translation invariance')
         translated=s;translated%alpha=3.7_rk
         call prepare_sas_general(translated,pt,st)
         changed=evaluate_sas_general(pt,q0,rule,OBS_D,st)
         qt=q0;qt(4)=q0(4)*translated%alpha
         ds=evaluate_sas_general(ps,qt,rule,OBS_D,st)/translated%alpha
         call low_close(changed,ds,'low D physical alpha scaling')

         ! Fourth-order age differentiation limits subtraction/truncation
         ! noise while testing dD/dtp=steady-T(tp), in physical time.
         h=q0(4)*5.e-4_rk
         do k=1,4
            qp=q0
            select case(k)
            case(1);qp(4)=q0(4)-2*h
            case(2);qp(4)=q0(4)-h
            case(3);qp(4)=q0(4)+h
            case(4);qp(4)=q0(4)+2*h
            end select
            dages(k)=evaluate_sas_general(ps,qp,rule,OBS_D,st)
            call require(st==ok,'low D derivative sample')
         end do
         derivative=(dages(1)-8*dages(2)+8*dages(3)-dages(4))/(12*h)
         steady=evaluate_steady(space,q0,r96,r96,st)
         Tstep=evaluate_sas_general(ps,q0,rule,OBS_T,st)
         call low_close(derivative,steady-Tstep,'physical low D age derivative')
         derivative_difference=max(derivative_difference,abs(derivative-(steady-Tstep)))
      end do

      ! Nine genuinely signed low-degree coefficients, including mixed terms.
      csmall(0,:)=[-.2_rk,.4_rk,.15_rk]
      csmall(1,:)=[.3_rk,-.5_rk,.25_rk]
      csmall(2,:)=[-.1_rk,.2_rk,-.3_rk]
      call init_poly(different,csmall,1.0_rk,2.0_rk,1.0_rk,st)
      call prepare_sas_general(different,pd,st);call prepare_spatial(different,space,st)
      q0=[.9_rk,1.8_rk,.1_rk,1.0_rk]
      changed=evaluate_sas_general(pd,q0,rule,OBS_D,st)
      dr=evaluate_spatial(space,q0,r96,r96,OBS_D,st)
      call low_close(changed,dr,'signed full low tensor versus SPACE2D96')
      call init_poly(summed,s%c+csmall,1.0_rk,2.0_rk,1.0_rk,st)
      call prepare_sas_general(summed,psum,st)
      base=evaluate_sas_general(ps,q0,rule,OBS_D,st)
      added=evaluate_sas_general(psum,q0,rule,OBS_D,st)
      call low_close(added,base+changed,'signed independent low tensor linearity')
      call sas_features(pd,q0,rule,OBS_D,wlow,st)
      call require(st==ok.and.all(ieee_is_finite(wlow)),'signed qk features initialized')
      call low_close(changed,sum(csmall*wlow),'signed binary64 D and qk features')
      write(*,'(a,4es24.15)') 'lowD maxabs/maxrel/qk-feature-diff/derivative-diff=', &
         max_absolute,max_relative,max_operator_difference,derivative_difference
   end subroutine

   subroutine low_close(actual,desired,label)
      real(rk),intent(in) :: actual,desired
      character(*),intent(in) :: label
      logical :: finite
      finite=ieee_is_finite(actual).and.ieee_is_finite(desired)
      call require(finite,'finite '//label)
      if(.not.finite) return
      call require(abs(actual-desired)<=1.e-12_rk+2.e-11_rk*abs(desired),label)
      if(abs(actual-desired)>1.e-12_rk+2.e-11_rk*abs(desired)) &
         write(*,'(3es24.15)') actual,desired,actual-desired
   end subroutine

   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(*),intent(in) :: label
      checks=checks+1
      if(condition) return
      failures=failures+1
      if(failures<=20) write(*,'(a)') 'FAIL '//label
   end subroutine

   subroutine close(actual,desired,tolerance,label)
      real(rk),intent(in) :: actual,desired,tolerance
      character(*),intent(in) :: label
      call require(ieee_is_finite(actual) .and. ieee_is_finite(desired) .and. &
         abs(actual-desired)<=tolerance*max(1.0_rk,abs(desired)),label)
      if(failures>0 .and. failures<=20 .and. abs(actual-desired)>tolerance*max(1.0_rk,abs(desired))) &
         write(*,'(3es24.15)') actual,desired,actual-desired
   end subroutine

   subroutine axis_reference(u,x,a,result)
      real(rk),intent(in) :: u,x,a
      real(rk),intent(out) :: result(0:32)
      real(rk) :: tx(0:32),weight,s
      integer :: i,j
      result=0
      do i=1,reference%n
         s=reference%x(i);tx(0)=1;tx(1)=s
         do j=2,32
            tx(j)=2*s*tx(j-1)-tx(j-2)
         end do
         weight=a*reference%w(i)*exp(-(u*(a*s-x))**2)
         result=result+weight*tx
      end do
   end subroutine
end program
