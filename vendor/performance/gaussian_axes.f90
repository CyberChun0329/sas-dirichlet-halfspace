module gaussian_axes
   ! Analytic finite-interval Gaussian axes. No spatial quadrature is used.
   ! u is the inverse Gaussian length: exp[-u**2*(xi-position)**2].
   use replay_types, only: rk,qk,ok,invalid,numerical_failure,quiet_nan
   use temperature_common, only: sqrtpi,stable_expm1
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   integer,parameter,public :: axis_series_max_terms=160
   real(qk),parameter,public :: axis_series_relative_tolerance=epsilon(1.0_qk)/16
   real(qk),parameter,public :: axis_recurrence_relative_budget=1.e-26_qk
   public :: cheb_axis,monomial_axis,gaussian_profile_axis,axis_precision
   public :: try_recurrence_axis
contains
   pure integer function axis_precision(M) result(bits)
      integer,intent(in) :: M
      bits=64
      if(M>2) bits=128
   end function

   subroutine cheb_axis(u,position,width,M,values,status,high_values,truncation_bound)
      real(rk),intent(in) :: u,position,width
      integer,intent(in) :: M
      real(rk),intent(out) :: values(0:M)
      integer,intent(out) :: status
      real(qk),optional,intent(out) :: high_values(0:M)
      real(rk),optional,intent(out) :: truncation_bound(0:M)
      real(qk) :: high(0:M),bound
      real(rk) :: low(0:2)
      values=quiet_nan();status=invalid
      if(present(high_values)) high_values=real(quiet_nan(),qk)
      if(present(truncation_bound)) truncation_bound=quiet_nan()
      if(.not.valid_axis(u,position,width,M)) return
      if(M<=2) then
         call low_axis(u,position,width,M,low,status)
         if(status/=ok) return
         values=low(0:M)
         ! The optional quad result is explicitly requested extended arithmetic;
         ! the mandatory low-degree output uses binary64 arithmetic.
         if(present(high_values) .or. present(truncation_bound)) then
            call analytic_axis(u,position,width,M,.true.,high,bound,status)
            if(status/=ok) return
         end if
      else
         call analytic_axis(u,position,width,M,.true.,high,bound,status)
         if(status/=ok) return
         if(maxval(abs(high))>real(huge(1.0_rk),qk)) then
            status=numerical_failure;return
         end if
         values=real(high,rk)
      end if
      if(.not.all(ieee_is_finite(values))) then
         values=quiet_nan();status=numerical_failure;return
      end if
      if(present(high_values)) high_values=high
      if(present(truncation_bound)) truncation_bound=real(bound,rk)
      status=ok
   end subroutine

   subroutine monomial_axis(u,position,width,M,values,status,high_values,truncation_bound)
      real(rk),intent(in) :: u,position,width
      integer,intent(in) :: M
      real(rk),intent(out) :: values(0:M)
      integer,intent(out) :: status
      real(qk),optional,intent(out) :: high_values(0:M)
      real(rk),optional,intent(out) :: truncation_bound(0:M)
      real(qk) :: high(0:M),bound,scale
      real(rk) :: low(0:2)
      integer :: j
      values=quiet_nan();status=invalid
      if(present(high_values)) high_values=real(quiet_nan(),qk)
      if(present(truncation_bound)) truncation_bound=quiet_nan()
      if(.not.valid_axis(u,position,width,M)) return
      if(M<=2) then
         call low_axis(u,position,width,M,low,status)
         if(status/=ok) return
         values(0)=low(0)
         if(M>=1) then
            scale=real(width,qk)*real(low(1),qk)
            if(abs(scale)>real(huge(1.0_rk),qk)) then
               values=quiet_nan();status=numerical_failure;return
            end if
            values(1)=width*low(1)
         end if
         if(M>=2) then
            scale=real(width,qk)**2*((real(low(2),qk)+real(low(0),qk))/2)
            if(abs(scale)>real(huge(1.0_rk),qk)) then
               values=quiet_nan();status=numerical_failure;return
            end if
            values(2)=width*(width*((low(2)+low(0))/2))
         end if
      end if
      if(M>2 .or. present(high_values) .or. present(truncation_bound)) then
         call analytic_axis(u,position,width,M,.false.,high,bound,status)
         if(status/=ok) return
         if(maxval(abs(high))>real(huge(1.0_rk),qk)) then
            values=quiet_nan();status=numerical_failure;return
         end if
         if(M>2) values=real(high,rk)
         if(present(high_values)) high_values=high
         if(present(truncation_bound)) then
            scale=1
            do j=0,M
               truncation_bound(j)=real(bound*scale,rk)
               scale=scale*real(width,qk)
            end do
         end if
      end if
      if(.not.all(ieee_is_finite(values))) then
         values=quiet_nan();status=numerical_failure;return
      end if
      status=ok
   end subroutine

   function gaussian_profile_axis(u,position,width,beta,status) result(value)
      real(rk),intent(in) :: u,position,width,beta
      integer,intent(out) :: status
      real(rk) :: value
      real(qk) :: Uq,Xq,Aq,Bq,lambda,center,lo,hi,diff,prefactor,resultq
      real(qk) :: beta_u,delta_u,linear,quadratic,previous,current,following,total,tail
      integer :: n,contractive
      value=quiet_nan();status=invalid
      if(.not.valid_axis(u,position,width,0)) return
      if(.not.ieee_is_finite(beta) .or. beta<0) return
      Uq=real(u,qk);Xq=real(position,qk);Aq=real(width,qk);Bq=real(beta,qk)
      beta_u=Uq*Aq;delta_u=Uq*Xq
      linear=2*beta_u*delta_u;quadratic=-(beta_u*beta_u+Bq*Bq)
      if(abs(linear)+abs(quadratic)<=2) then
         ! A narrow shifted interval must not subtract equal rounded erf
         ! values. Integrate the convergent Gaussian-profile series exactly.
         previous=1;current=linear;total=2
         contractive=max(4,ceiling(2*(abs(linear)+2*abs(quadratic))))
         do n=2,axis_series_max_terms
            following=(linear*current+2*quadratic*previous)/real(n,qk)
            if(mod(n,2)==0) total=total+2*following/real(n+1,qk)
            previous=current;current=following
            tail=4*max(abs(previous),abs(current))
            if(n>=contractive .and. tail<=axis_series_relative_tolerance) exit
         end do
         if(n>axis_series_max_terms) then
            status=numerical_failure;return
         end if
         resultq=Aq*exp(-delta_u*delta_u)*total
      else
         lambda=sqrt(Uq*Uq+(Bq/Aq)**2)
         center=(Uq/lambda)**2*Xq
         lo=lambda*(-Aq-center);hi=lambda*(Aq-center)
         if(lo>=0) then
            diff=erfc(lo)-erfc(hi)
         else if(hi<=0) then
            diff=erfc(-hi)-erfc(-lo)
         else
            diff=erf(hi)-erf(lo)
         end if
         ! Complete the square without subtracting two nearly equal squares.
         prefactor=exp(-((Uq/lambda)*(Bq/Aq)*Xq)**2)
         resultq=prefactor*sqrt(acos(-1.0_qk))*diff/(2*lambda)
      end if
      if(abs(resultq)>real(huge(1.0_rk),qk)) then
         status=numerical_failure;return
      end if
      value=real(resultq,rk)
      if(.not.ieee_is_finite(value)) then
         value=quiet_nan();status=numerical_failure;return
      end if
      status=ok
   end function

   pure logical function valid_axis(u,position,width,M)
      real(rk),intent(in) :: u,position,width
      integer,intent(in) :: M
      valid_axis=all(ieee_is_finite([u,position,width]))
      if(valid_axis) valid_axis=u>=0 .and. width>0 .and. M>=0 .and. M<=32
   end function

   subroutine analytic_axis(u,position,width,M,chebyshev,values,bound,status)
      real(rk),intent(in) :: u,position,width
      integer,intent(in) :: M
      logical,intent(in) :: chebyshev
      real(qk),intent(out) :: values(0:M),bound
      integer,intent(out) :: status
      real(qk) :: a,beta,delta,variation,anchor,scale,tail_log_scale,local(0:M),error,rounding(0:M)
      logical :: accepted
      integer :: j
      values=0;bound=0;status=numerical_failure
      a=real(width,qk);beta=a*real(u,qk);delta=real(position,qk)*real(u,qk)
      tail_log_scale=log(2*a)
      if(.not.chebyshev) tail_log_scale=tail_log_scale+max(0.0_qk,real(M,qk)*log(a))
      if(tail_log_scale>log(huge(1.0_qk))-2) return
      if(u==0) then
         do j=0,M
            if(mod(j,2)/=0) cycle
            if(chebyshev) then
               values(j)=2*a/real(1-j*j,qk)
            else
               values(j)=2*a**(j+1)/real(j+1,qk)
            end if
         end do
         status=ok;return
      end if
      call try_recurrence_axis(u,position,width,M,chebyshev,values,rounding,accepted)
      if(accepted) then
         ! Finite-endpoint recurrences have no series truncation. Their
         ! internal roundoff estimate is not a rigorous error certificate.
         bound=0;status=ok;return
      end if
      variation=beta*beta+2*abs(beta*delta)
      if(variation<=2) then
         call panel(0.0_qk,1.0_qk,-delta,beta,M,chebyshev,local,error,status)
         if(status/=ok) return
         values=a*local;bound=a*error
      else if(abs(delta)<=beta) then
         anchor=real(position,qk)/a
         call segment(0.0_qk,beta+delta,anchor,-1,local,error,status)
         if(status/=ok) return
         values=local;bound=error
         call segment(0.0_qk,beta-delta,anchor,1,local,error,status)
         if(status/=ok) return
         values=values+local;bound=bound+error
      else if(delta>beta) then
         call segment(delta-beta,delta+beta,1.0_qk,-1,values,bound,status)
         if(status/=ok) return
      else
         call segment(-delta-beta,-delta+beta,-1.0_qk,1,values,bound,status)
         if(status/=ok) return
      end if
      if(.not.chebyshev) then
         scale=1
         do j=0,M
            values(j)=values(j)*scale
            scale=scale*a
         end do
      end if
      if(.not.all(ieee_is_finite(values))) return
      status=ok
   contains
      subroutine segment(rmin,rmax,origin,direction,total,error_total,istat)
         real(qk),intent(in) :: rmin,rmax,origin
         integer,intent(in) :: direction
         real(qk),intent(out) :: total(0:M),error_total
         integer,intent(out) :: istat
         real(qk) :: r,dr,mid,h,rcenter,part(0:M),part_error,correction(0:M),next(0:M),omitted
         integer :: count
         total=0;error_total=0;correction=0;r=rmin;istat=numerical_failure
         count=0
         do while(r<rmax)
            ! The entire remaining Gaussian integral has this analytic bound
            ! because every normalized basis function has magnitude <=1 on
            ! the patch. Termination is at binary128 rounding resolution,
            ! with the omitted bound returned explicitly, not a spatial
            ! accuracy cutoff or a repair of a computed signed value.
            omitted=(a/beta)*(sqrt(acos(-1.0_qk))/2)*erfc(r)
            if(maxval(abs(total))>0) then
               if(omitted<=epsilon(1.0_qk)*maxval(abs(total))/32) then
                  error_total=error_total+omitted
                  exit
               end if
            end if
            ! This only omits a remainder whose full physical absolute bound
            ! underflows binary128. There is no fixed spatial/tail cutoff.
            if(exp(tail_log_scale-r*r)==0) exit
            dr=min(rmax-r,2/(1+r))
            if(dr<=0 .or. r+dr==r) return
            rcenter=r+dr/2
            mid=origin+real(direction,qk)*((r-rmin)+dr/2)/beta
            h=dr/(2*beta)
            ! Physical s increases with the local polynomial coordinate. On
            ! the left/outside-right branch its positive Gaussian distance
            ! decreases, so the Gaussian linear coefficient reverses sign.
            call panel(mid,h,rcenter,real(direction,qk)*dr/2,M,chebyshev,part,part_error,istat)
            if(istat/=ok) return
            part=a*part
            next=total+part
            where(abs(total)>=abs(part))
               correction=correction+((total-next)+part)
            elsewhere
               correction=correction+((part-next)+total)
            end where
            total=next;error_total=error_total+a*part_error
            r=r+dr;count=count+1
            if(count>100000) then
               istat=numerical_failure;return
            end if
         end do
         total=total+correction;istat=ok
      end subroutine
   end subroutine

   subroutine try_recurrence_axis(u,position,width,M,chebyshev,values,rounding,accepted)
      ! Exact finite-interval Gaussian moments, with a conservative qk
      ! roundoff-policy gate. Both consuming routes use the same primitive.
      ! Rejection changes the arithmetic route, never the mathematical input
      ! or a computed signed moment; the original panels remain the fallback.
      real(rk),intent(in) :: u,position,width
      integer,intent(in) :: M
      logical,intent(in) :: chebyshev
      real(qk),intent(out) :: values(0:M),rounding(0:M)
      logical,intent(out) :: accepted
      real(qk) :: a,Uq,beta,delta,q,lo,hi,em,ep,d,s,rf1,rf2,fac,t,gamma,eq,et,eb,e0
      real(qk) :: c(0:M),e(0:M),h(0:M),f(0:M),budget,term1,term2,term3,boundary
      real(qk) :: force,Psum,Qsum,r,rpower,scale,nextscale,logscale
      integer :: k
      values=real(quiet_nan(),qk);rounding=values;accepted=.false.
      if(.not.valid_axis(u,position,width,M) .or. u==0) return
      a=real(width,qk);Uq=real(u,qk);beta=a*Uq;delta=real(position,qk)*Uq
      q=real(position,qk)/a
      ! These are branch-selection limits, not interval/source cutoffs.
      if(beta<.5_qk .or. abs(q)>2) return
      logscale=log(2*a)
      if(.not.chebyshev) logscale=logscale+max(0.0_qk,real(M,qk)*log(a))
      if(logscale>log(huge(1.0_qk))-2) return
      lo=-Uq*(a+real(position,qk));hi=Uq*(a-real(position,qk))
      em=exp(-lo*lo);ep=exp(-hi*hi)
      if(lo>=0) then
         rf1=erfc(lo);rf2=erfc(hi)
      else if(hi<=0) then
         rf1=erfc(-hi);rf2=erfc(-lo)
      else
         rf1=erf(hi);rf2=-erf(-lo)
      end if
      fac=sqrt(acos(-1.0_qk))/(2*beta)
      c=0;e=0;h=0;f=0;c(0)=fac*(rf1-rf2)
      if(.not.ieee_is_finite(c(0)) .or. c(0)<=tiny(1.0_qk)/axis_recurrence_relative_budget) return
      budget=axis_recurrence_relative_budget*c(0)
      ! Assumes small-ulp qk elementary-function errors; argument sensitivity
      ! and both erf/erfc operands are included. This is an estimate, not an
      ! intrinsic-library accuracy proof or a signed-contraction certificate.
      gamma=64*epsilon(1.0_qk)/(1-64*epsilon(1.0_qk))
      eb=gamma*(em+ep+2*(abs(lo)*em+abs(hi)*ep)*(beta+abs(delta)))
      e0=gamma*(abs(c(0))+abs(fac)*(abs(rf1)+abs(rf2))+(beta+abs(delta))*(em+ep)/beta)
      e(0)=e0
      if(.not.ieee_is_finite(e0) .or. e0>budget) return
      t=1/(beta*beta);eq=gamma*abs(q);et=gamma*t;s=em+ep
      if(delta>=0) then
         d=ep*axis_expm1(-4*beta*delta)
      else
         d=-em*axis_expm1(4*beta*delta)
      end if
      if(M>=1) then
         term1=q*c(0);term2=.5_qk*t*d;c(1)=term1+term2
         e(1)=abs(q)*e0+.5_qk*t*eb+abs(c(0))*eq+.5_qk*abs(d)*et+ &
            gamma*(abs(term1)+abs(term2))
         if(.not.ieee_is_finite(e(1)) .or. e(1)>budget) return
         if(abs(c(1))>c(0)+e(1)) return
      end if
      if(chebyshev) then
         h(0)=c(0);f(0)=e0
         if(M>=1) then
            h(1)=2*c(1);f(1)=2*e(1)+gamma*2*abs(c(1))
         end if
         r=1
         if(abs(q)>1) r=(abs(q)+sqrt((abs(q)-1)*(abs(q)+1)))*(1+gamma)
         Psum=0;Qsum=0;rpower=1
         do k=1,M-1
            boundary=s
            if(mod(k,2)==0) boundary=-d
            term1=2*q*c(k);term2=-c(k-1);term3=t*(real(k,qk)*h(k-1)-boundary)
            c(k+1)=term1+term2+term3
            force=real(k,qk)*t*f(k-1)+t*eb+2*(abs(c(k))+e(k))*eq+ &
               (abs(real(k,qk)*h(k-1)-boundary)+real(k,qk)*f(k-1)+eb)*et+ &
               gamma*(abs(term1)+abs(term2)+t*(real(k,qk)*abs(h(k-1))+abs(boundary)))
            Psum=r*Psum+force;Qsum=r*Qsum+Psum
            e(k+1)=real(k+1,qk)*r*rpower*e(1)+real(k,qk)*rpower*e0+Qsum
            if(.not.ieee_is_finite(e(k+1)) .or. e(k+1)>budget) return
            if(.not.ieee_is_finite(c(k+1)) .or. abs(c(k+1))>c(0)+e(k+1)) return
            h(k+1)=2*c(k+1)+h(k-1)
            f(k+1)=2*e(k+1)+f(k-1)+gamma*(2*abs(c(k+1))+abs(h(k-1)))
            rpower=rpower*r
         end do
      else
         do k=2,M
            boundary=-s
            if(mod(k-1,2)==0) boundary=d
            term1=q*c(k-1);term2=real(k-1,qk)*t*c(k-2)/2;term3=t*boundary/2
            c(k)=term1+term2+term3
            e(k)=abs(q)*e(k-1)+real(k-1,qk)*t*e(k-2)/2+t*eb/2+ &
               (abs(c(k-1))+e(k-1))*eq+ &
               (real(k-1,qk)*(abs(c(k-2))+e(k-2))+abs(boundary)+eb)*et/2+ &
               gamma*(abs(term1)+abs(term2)+abs(term3))
            if(.not.ieee_is_finite(e(k)) .or. e(k)>budget) return
            if(.not.ieee_is_finite(c(k)) .or. abs(c(k))>c(0)+e(k)) return
         end do
      end if
      scale=a
      do k=0,M
         values(k)=scale*c(k);rounding(k)=scale*e(k)+gamma*abs(values(k))
         if(k<M .and. .not.chebyshev) then
            nextscale=scale*a
            if(.not.ieee_is_finite(nextscale)) return
            scale=nextscale
         end if
      end do
      accepted=all(ieee_is_finite(values)).and.all(ieee_is_finite(rounding))
   end subroutine

   pure real(qk) function axis_expm1(x) result(value)
      real(qk),intent(in) :: x
      real(qk) :: term,next
      integer :: k
      if(abs(x)>.25_qk) then
         value=exp(x)-1;return
      end if
      value=x;term=x
      do k=2,80
         term=term*x/real(k,qk);next=value+term
         value=next
         if(abs(term)<=epsilon(1.0_qk)*abs(value)/16) exit
      end do
   end function

   subroutine panel(mid,h,gauss_mid,gauss_half,M,chebyshev,values,bound,status)
      real(qk),intent(in) :: mid,h,gauss_mid,gauss_half
      integer,intent(in) :: M
      logical,intent(in) :: chebyshev
      real(qk),intent(out) :: values(0:M),bound
      integer,intent(out) :: status
      real(qk) :: coefficients(0:M,0:M),moments(0:M),A,B,previous,current,following
      real(qk) :: term,tail,scale
      integer :: n,j,p,contractive
      coefficients=0;coefficients(0,0)=1
      if(M>0) then
         coefficients(1,0)=mid;coefficients(1,1)=h
      end if
      do p=2,M
         if(chebyshev) then
            coefficients(p,:)=2*mid*coefficients(p-1,:)-coefficients(p-2,:)
            coefficients(p,1:)=coefficients(p,1:)+2*h*coefficients(p-1,:M-1)
         else
            coefficients(p,:)=mid*coefficients(p-1,:)
            coefficients(p,1:)=coefficients(p,1:)+h*coefficients(p-1,:M-1)
         end if
      end do
      A=-2*gauss_mid*gauss_half;B=-gauss_half*gauss_half
      contractive=max(4,ceiling(2*(abs(A)+2*abs(B))))
      previous=1;current=A;moments=0
      do j=0,M
         if(mod(j,2)==0) moments(j)=2/real(j+1,qk)
         if(mod(j,2)==1) moments(j)=2*A/real(j+2,qk)
      end do
      status=numerical_failure;tail=huge(1.0_qk)
      do n=2,axis_series_max_terms
         following=(A*current+2*B*previous)/real(n,qk)
         do j=0,M
            if(mod(j+n,2)==0) then
               term=2*following/real(j+n+1,qk)
               moments(j)=moments(j)+term
            end if
         end do
         previous=current;current=following
         ! Beyond contractive, each two-coefficient envelope contracts by
         ! at least 1/2. Thus 4*envelope bounds the remaining absolute series.
         if(n>=contractive) then
            tail=4*max(abs(previous),abs(current))
            if(tail<=axis_series_relative_tolerance) exit
         end if
      end do
      if(n>axis_series_max_terms) return
      scale=h*exp(-gauss_mid*gauss_mid)
      values=matmul(coefficients,moments)*scale
      ! |T_m(s)|<=1 and |s**m|<=1 on the real patch. This is a series
      ! truncation bound, NOT a certificate for binary128 rounding or final
      ! cancellation of a signed source contraction.
      bound=2*abs(scale)*tail
      if(.not.all(ieee_is_finite(values))) return
      status=ok
   end subroutine

   subroutine low_axis(u,position,width,M,values,status)
      ! Same dedicated degree<=2 binary64 formulas as prepared_sas. The
      ! source/query temperature baseline itself calls that frozen module.
      real(rk),intent(in) :: u,position,width
      integer,intent(in) :: M
      real(rk),intent(out) :: values(0:2)
      integer,intent(out) :: status
      real(rk) :: beta,delta,A,B,variation,lo,hi,difference,k0,k1,k2,e0,e1
      real(rk) :: change,factor,shift,ib,previous,current,following,term
      real(rk) :: total(0:2),correction(0:2),moments(0:2,0:24)
      integer :: n,p
      logical :: small,huge_change
      values=0;status=numerical_failure
      if(real(width,qk)*real(u,qk)>real(huge(1.0_rk),qk)) return
      if(abs(real(position,qk))*real(u,qk)>real(huge(1.0_rk),qk)) return
      if(u==0 .and. width>huge(1.0_rk)/2) return
      beta=width*u;delta=position*u
      if(.not.all(ieee_is_finite([beta,delta]))) return
      small=.false.
      if(beta<=.5_rk) then
         A=2*beta*delta;B=-beta*beta;variation=abs(A)+abs(B)
         small=variation<=.25_rk
      end if
      if(M==0) small=small .and. (abs(delta)>1 .or. abs(position)/width>128 .or. u==0)
      if(small) then
         moments=0
         do n=0,24
            if(mod(n,2)==0) then
               moments(0,n)=2/real(n+1,rk)
               moments(2,n)=4/real(n+3,rk)-2/real(n+1,rk)
            else
               moments(1,n)=2/real(n+2,rk)
            end if
         end do
         previous=1;current=A;total=moments(:,0)+moments(:,1)*A;correction=0
         do n=2,24
            following=(A*current+2*B*previous)/real(n,rk)
            do p=0,M
               term=moments(p,n)*following
               call low_add(total(p),correction(p),term)
            end do
            previous=current;current=following
            if(n>=4 .and. max(abs(previous),abs(current))<1.e-20_rk) exit
         end do
         values(:M)=width*low_gaussian(delta)*(total(:M)+correction(:M))
      else
         if(beta==0 .or. beta>huge(1.0_rk)-abs(delta)) return
         lo=-beta-delta;hi=beta-delta
         if(max(abs(lo),abs(hi))<=1) then
            difference=erf(hi)-erf(lo)
         else if(lo>=0) then
            difference=erfc(lo)-erfc(hi)
         else if(hi<=0) then
            difference=erfc(-hi)-erfc(-lo)
         else
            difference=erf(hi)-erf(lo)
         end if
         k0=sqrtpi*difference/2;factor=width/beta;values(0)=factor*k0
         if(M>0) then
            e0=low_gaussian(lo);e1=low_gaussian(hi)
            if(k0==0 .and. e0==0 .and. e1==0) then
               status=ok;return
            end if
            huge_change=.false.
            if(beta>.25_rk) huge_change=abs(delta)>(huge(1.0_rk)/4)/beta
            if(huge_change) then
               if(delta<=0) then
                  k1=e0/2
               else
                  k1=-e1/2
               end if
            else
               change=-4*(beta*delta)
               if(change>=0) then
                  k1=-e0*stable_expm1(-change)/2
               else
                  k1=e1*stable_expm1(change)/2
               end if
            end if
            ib=1/beta;shift=position/width
            values(1)=factor*(shift*k0+k1*ib)
            if(M==2) then
               if(abs(shift)>sqrt(huge(1.0_rk))/4) return
               k2=(k0-hi*e1+lo*e0)/2
               values(2)=factor*((2*shift*shift-1)*k0+4*shift*k1*ib+2*k2*ib*ib)
            end if
         end if
      end if
      if(.not.all(ieee_is_finite(values(:M)))) return
      status=ok
   end subroutine

   pure real(rk) function low_gaussian(x) result(value)
      real(rk),intent(in) :: x
      if(abs(x)>sqrt(huge(1.0_rk))) then
         value=0
      else
         value=exp(-x*x)
      end if
   end function

   pure subroutine low_add(total,correction,term)
      real(rk),intent(inout) :: total,correction
      real(rk),intent(in) :: term
      real(rk) :: next
      next=total+term
      if(abs(total)>=abs(term)) then
         correction=correction+((total-next)+term)
      else
         correction=correction+((term-next)+total)
      end if
      total=next
   end subroutine
end module
