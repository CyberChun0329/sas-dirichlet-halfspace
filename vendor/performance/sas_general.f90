module sas_general
   use replay_types, only: rk,qk,source_t,SRC_POLY,SRC_GAUSSIAN,OBS_T,OBS_D, &
      ok,invalid,numerical_failure,quiet_nan,gauss_rule_t,valid_rule,mapped_rule, &
      init_poly,check_query,valid_geometry,pi,sqrtpi,add_compensated
   use temperature_common, only: patch_source_t,init_source
   use prepared_sas, only: prepared_sas_t,prepare_sas,evaluate_sas
   use gaussian_axes, only: cheb_axis,gaussian_profile_axis
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   type,public :: sas_general_t
      private
      type(source_t) :: source
      type(prepared_sas_t) :: low
      logical :: ready=.false.,low_ready=.false.
   end type
   public :: prepare_sas_general,evaluate_sas_general,sas_features,evaluate_sas_complementary
contains
   subroutine prepare_sas_general(source,prepared,status)
      type(source_t),intent(in) :: source
      type(sas_general_t),intent(out) :: prepared
      integer,intent(out) :: status
      type(patch_source_t) :: low_source
      real(rk) :: c(0:2,0:2)
      prepared=sas_general_t();status=invalid
      if(.not.source%ready .or. .not.valid_geometry(source)) return
      select case(source%kind)
      case(SRC_POLY)
         if(.not.allocated(source%c)) return
         if(any(lbound(source%c)/=0)) return
         call init_poly(prepared%source,source%c,source%a,source%b,source%alpha, &
            status,source%cx,source%cy)
         if(status/=ok) return
         if(max(prepared%source%mx,prepared%source%my)<=2) then
            c=0
            c(:prepared%source%my,:prepared%source%mx)= &
               prepared%source%c(:prepared%source%my,:prepared%source%mx)
            call init_source(low_source,c,source%a,source%b,source%alpha,status,source%cx,source%cy)
            if(status/=ok) return
            call prepare_sas(low_source,prepared%low,status)
            if(status/=ok) return
            prepared%low_ready=.true.
         end if
      case(SRC_GAUSSIAN)
         if(.not.all(ieee_is_finite(source%beta)) .or. any(source%beta<=0)) return
         prepared%source=source;prepared%source%mx=0;prepared%source%my=0
      case default
         ! A sampled map must first be projected to its declared polynomial.
         return
      end select
      prepared%ready=.true.;status=ok
   end subroutine

   function evaluate_sas_general(prepared,query,rule,observable,status) result(value)
      type(sas_general_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rule
      integer,intent(in) :: observable
      integer,intent(out) :: status
      real(rk) :: value
      value=quiet_nan();status=invalid
      if(.not.prepared%ready) return
      if(observable==OBS_T .and. prepared%low_ready) then
         ! Use the low-degree binary64 temperature evaluator.
         value=evaluate_sas(prepared%low,query,rule,status)
         return
      end if
      call integral(prepared,query,rule,observable,.false.,value,status)
   end function

   subroutine sas_features(prepared,query,rule,observable,W,status)
      type(sas_general_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rule
      integer,intent(in) :: observable
      real(rk),intent(out) :: W(0:,0:)
      integer,intent(out) :: status
      real(rk) :: value
      W=quiet_nan();status=invalid
      if(.not.prepared%ready) return
      if(prepared%source%kind/=SRC_POLY) return
      if(size(W,1)<prepared%source%my+1 .or. size(W,2)<prepared%source%mx+1) return
      call integral(prepared,query,rule,observable,.false.,value,status,W)
   end subroutine

   function evaluate_sas_complementary(prepared,query,rule,status) result(value)
      ! Independent interval route: steady v-integral minus [0,v0]. The
      ! steady integral uses the same explicit cutoff 12, not a cached final
      ! response or an unreported elementary/extended-precision steady value.
      type(sas_general_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rule
      integer,intent(out) :: status
      real(rk) :: value
      call integral(prepared,query,rule,OBS_T,.true.,value,status)
   end function

   subroutine integral(prepared,query,rule,observable,complementary,value,status,W)
      type(sas_general_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rule
      integer,intent(in) :: observable
      logical,intent(in) :: complementary
      real(rk),intent(out) :: value
      integer,intent(out) :: status
      real(rk),optional,intent(out) :: W(0:,0:)
      real(rk) :: x,y,z,dt,v0,h,candidates(9),lower,upper,sign_region,u
      real(rk) :: ax(0:32),ay(0:32)
      real(rk),allocatable :: v(:),weights(:)
      real(qk) :: axq(0:32),ayq(0:32),term,contraction,total,correction,next
      real(qk) :: feature(0:32,0:32),prefactor,weightq,spanq,prefactorq
      integer :: region,regions,i,p,q,mx,my
      logical :: high_degree
      value=quiet_nan();status=invalid
      if(present(W)) W=quiet_nan()
      if(.not.prepared%ready .or. .not.valid_rule(rule)) return
      if(observable/=OBS_T .and. observable/=OBS_D) return
      call check_query(prepared%source,query,x,y,z,dt,status)
      if(status/=ok) return
      if(observable==OBS_D.and.prepared%low_ready.and..not.present(W)) then
         call low_deficit(prepared%source,query,rule,x,y,z,dt,value,status)
         return
      end if
      mx=prepared%source%mx;my=prepared%source%my
      high_degree=max(mx,my)>2
      if(z>=24*sqrt(dt)) then
         v0=12
         if(observable==OBS_T) then
            value=0
            if(present(W)) W=0
            return
         end if
      else
         v0=z/(2*sqrt(dt))
      end if
      spanq=max(real(prepared%source%a,qk)+abs(real(x,qk)), &
         real(prepared%source%b,qk)+abs(real(y,qk)))
      ! A geometry knot beyond 12 cannot split [0,12]. This bound only
      ! represents the integration partition; it never clips a response.
      h=real(min(12.0_qk,real(z,qk)/spanq),rk)
      candidates=[h/4,h,4*h,16*h,.25_rk,1.0_rk,3.0_rk,6.0_rk,v0]
      prefactorq=2/(real(pi,qk)*real(sqrtpi,qk))
      total=0;correction=0;feature=0
      regions=1
      if(observable==OBS_D .or. complementary) regions=2
      do region=1,regions
         sign_region=1
         if(complementary) then
            lower=0;upper=12
            if(region==2) then
               upper=v0;sign_region=-1
            end if
         else if(observable==OBS_D) then
            if(region==1) then
               lower=0;upper=v0
            else
               lower=v0;upper=12
            end if
         else
            lower=v0;upper=12
         end if
         call mapped_rule(lower,upper,candidates,rule,v,weights,status)
         if(status/=ok) return
         do i=1,size(v)
            u=v(i)/z
            if(.not.ieee_is_finite(u)) then
               status=numerical_failure;return
            end if
            if(prepared%source%kind==SRC_GAUSSIAN) then
               ax(0)=gaussian_profile_axis(u,x,prepared%source%a,prepared%source%beta(1),status)
               if(status/=ok) return
               ay(0)=gaussian_profile_axis(u,y,prepared%source%b,prepared%source%beta(2),status)
               if(status/=ok) return
               axq(0)=real(ax(0),qk);ayq(0)=real(ay(0),qk)
            else if(high_degree) then
               ! Optional outputs preserve the analytic binary128 values until
               ! the complete tensor contraction has been accumulated.
               call cheb_axis(u,x,prepared%source%a,mx,ax(:mx),status,axq(:mx))
               if(status/=ok) return
               call cheb_axis(u,y,prepared%source%b,my,ay(:my),status,ayq(:my))
               if(status/=ok) return
            else
               call cheb_axis(u,x,prepared%source%a,mx,ax(:mx),status)
               if(status/=ok) return
               call cheb_axis(u,y,prepared%source%b,my,ay(:my),status)
               if(status/=ok) return
               axq(:mx)=real(ax(:mx),qk);ayq(:my)=real(ay(:my),qk)
            end if
            weightq=real(weights(i)*sign_region,qk)*exp(-real(v(i),qk)**2)*prefactorq
            if(observable==OBS_D .and. region==2) then
               prefactor=weightq/(4*real(prepared%source%alpha,qk))
            else
               axq(:mx)=real(u,qk)*axq(:mx);ayq(:my)=real(u,qk)*ayq(:my)
               prefactor=weightq
               if(observable==OBS_D) prefactor=prefactor*real(query(4),qk)
            end if
            contraction=0
            do p=0,mx
               do q=0,my
                  term=prefactor*ayq(q)*axq(p)
                  if(prepared%source%kind==SRC_POLY) then
                     feature(q,p)=feature(q,p)+term
                     contraction=contraction+real(prepared%source%c(q,p),qk)*term
                  else
                     contraction=term
                  end if
               end do
            end do
            if(.not.ieee_is_finite(contraction)) then
               status=numerical_failure;return
            end if
            next=total+contraction
            if(abs(total)>=abs(contraction)) then
               correction=correction+((total-next)+contraction)
            else
               correction=correction+((contraction-next)+total)
            end if
            total=next
         end do
      end do
      if(abs(total+correction)>real(huge(1.0_rk),qk)) then
         status=numerical_failure;return
      end if
      value=real(total+correction,rk)
      if(present(W)) then
         if(maxval(abs(feature(:my,:mx)))>real(huge(1.0_rk),qk)) then
            value=quiet_nan();status=numerical_failure;return
         end if
         W=0;W(:my,:mx)=real(feature(:my,:mx),rk)
      end if
      if(.not.ieee_is_finite(value)) then
         value=quiet_nan();status=numerical_failure;return
      end if
      if(present(W)) then
         if(.not.all(ieee_is_finite(W))) then
            W=quiet_nan();value=quiet_nan();status=numerical_failure;return
         end if
      end if
      status=ok
   end subroutine

   subroutine low_deficit(source,query,rule,x,y,z,dt,value,status)
      ! Binary64 direct deficit for the complete low-degree polynomial. The
      ! source-dependent tensor and v integral both use compensated binary64
      ! sums. Feature/operator requests continue through the existing qk path.
      type(source_t),intent(in) :: source
      real(rk),intent(in) :: query(4),x,y,z,dt
      type(gauss_rule_t),intent(in) :: rule
      real(rk),intent(out) :: value
      integer,intent(out) :: status
      real(rk) :: v0,span,h,lo,hi,u,scale,weight,term
      real(rk) :: ax(0:2),ay(0:2),contraction,tensor_correction,total,correction
      real(rk) :: cuts(9)
      real(rk),allocatable :: v(:),weights(:)
      integer :: region,i,p,q
      value=quiet_nan();status=invalid
      if(source%kind/=SRC_POLY.or.max(source%mx,source%my)>2) return
      if(z>=24*sqrt(dt)) then
         v0=12
      else
         v0=z/(2*sqrt(dt))
      end if
      span=max(source%a+abs(x),source%b+abs(y));h=12
      if(z/12<span) h=z/span
      cuts=[h/4,h,4*h,16*h,.25_rk,1.0_rk,3.0_rk,6.0_rk,v0]
      total=0;correction=0
      do region=1,2
         if(region==1) then
            lo=0;hi=v0
         else
            lo=v0;hi=12
         end if
         call mapped_rule(lo,hi,cuts,rule,v,weights,status)
         if(status/=ok) return
         do i=1,size(v)
            u=v(i)/z
            if(.not.ieee_is_finite(u)) then
               status=numerical_failure;return
            end if
            call cheb_axis(u,x,source%a,source%mx,ax(:source%mx),status)
            if(status/=ok) return
            call cheb_axis(u,y,source%b,source%my,ay(:source%my),status)
            if(status/=ok) return
            weight=weights(i)*(2/(pi*sqrtpi))*exp(-v(i)*v(i))
            if(region==1) then
               ! Here min(tp,z^2/(4*alpha*v^2))=physical tp. Balanced
               ! factors avoid forming u^2 or multiplying a large tp first.
               scale=u*sqrt(query(4))
               ax(:source%mx)=scale*ax(:source%mx)
               ay(:source%my)=scale*ay(:source%my)
            else
               ! Above v0, u^2*min(...)=1/(4*alpha). Divide separately
               ! rather than squaring z/v or forming a possibly large 4*alpha.
               weight=(weight/4)/source%alpha
            end if
            contraction=0;tensor_correction=0
            do p=0,source%mx
               do q=0,source%my
                  term=source%c(q,p)*ax(p)*ay(q)
                  if(.not.ieee_is_finite(term)) then
                     status=numerical_failure;return
                  end if
                  call add_compensated(contraction,tensor_correction,term)
               end do
            end do
            term=weight*(contraction+tensor_correction)
            if(.not.ieee_is_finite(term)) then
               status=numerical_failure;return
            end if
            call add_compensated(total,correction,term)
         end do
      end do
      value=total+correction
      if(.not.ieee_is_finite(value)) then
         value=quiet_nan();status=numerical_failure;return
      end if
      status=ok
   end subroutine
end module
