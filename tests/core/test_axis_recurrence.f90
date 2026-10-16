program test_axis_recurrence
   use replay_types
   use gaussian_axes
   implicit none
   real(rk),parameter :: widths(3)=[.2_rk,1.0_rk,2.3_rk]
   real(rk),parameter :: betas(6)=[.51_rk,1.0_rk,3.0_rk,10.0_rk,30.0_rk,100.0_rk]
   real(rk),parameter :: shifts(7)=[0.0_rk,.3_rk,.9_rk,-.9999999999999999_rk,1.0_rk,1.2_rk,-1.2_rk]
   real(qk) :: nodes(128),weights(128),fast(0:32),estimate(0:32),oracle(0:32),mass
   real(qk) :: high(0:32),scale,err,max_fast_error,max_output_error
   real(rk) :: values(0:32),bound(0:32),a,u,x
   integer :: ia,ib,iq,basis,m,j,st,checks,failed,accepted_count(2),fallback_count(2)
   logical :: accepted,cheb
   checks=0;failed=0;accepted_count=0;fallback_count=0;max_fast_error=0;max_output_error=0
   call quad_rule(nodes,weights)
   do ia=1,size(widths)
      a=widths(ia)
      do ib=1,size(betas)
         u=betas(ib)/a
         do iq=1,size(shifts)
            x=shifts(iq)*a
            do basis=1,2
               cheb=basis==1
               call quadrature_oracle(u,x,a,32,cheb,oracle,mass)
               call try_recurrence_axis(u,x,a,32,cheb,fast,estimate,accepted)
               if(accepted) then
                  accepted_count(basis)=accepted_count(basis)+1
                  scale=real(a,qk)
                  do j=0,32
                     err=abs(fast(j)-oracle(j))/(scale*mass)
                     max_fast_error=max(max_fast_error,err)
                     call require(err<=2.e-25_qk,'accepted recurrence mass-scaled quad oracle')
                     call require(estimate(j)>=0 .and. ieee_is_finite(estimate(j)),'finite roundoff estimate')
                     if(.not.cheb) scale=scale*real(a,qk)
                  end do
               else
                  fallback_count(basis)=fallback_count(basis)+1
               end if
               if(cheb) then
                  call cheb_axis(u,x,a,32,values,st,high,bound)
               else
                  call monomial_axis(u,x,a,32,values,st,high,bound)
               end if
               call require(st==ok,'production axis status including fallback')
               scale=real(a,qk)
               do j=0,32
                  ! Positive Gaussian mass is the meaningful scale even when
                  ! an individual signed moment is zero or a tiny tail.
                  err=abs(real(values(j),qk)-oracle(j))/(scale*mass)
                  max_output_error=max(max_output_error,err)
                  call require(err<=2.e-12_qk,'binary64 mass-scaled output versus quad oracle')
                  call require(bound(j)>=0 .and. ieee_is_finite(bound(j)),'truncation-only bound')
                  if(.not.cheb) scale=scale*real(a,qk)
               end do
            end do
         end do
      end do
   end do
   call try_recurrence_axis(100.0_rk,.9_rk,1.0_rk,32,.true.,fast,estimate,accepted)
   call require(accepted,'full degree32 interior fast branch')
   call try_recurrence_axis(.1_rk,.3_rk,1.0_rk,32,.true.,fast,estimate,accepted)
   call require(.not.accepted,'small beta retains panels')
   call try_recurrence_axis(1.0_rk,20.0_rk,1.0_rk,32,.true.,fast,estimate,accepted)
   call require(.not.accepted,'far observer retains panels')
   do m=0,32
      do basis=1,2
         cheb=basis==1
         call quadrature_oracle(100.0_rk,0.0_rk,1.0_rk,m,cheb,oracle,mass)
         call try_recurrence_axis(100.0_rk,0.0_rk,1.0_rk,m,cheb,fast(:m),estimate(:m),accepted)
         call require(accepted,'all degree prefixes supported')
         do j=0,m
            call require(abs(fast(j)-oracle(j))<=2.e-25_qk*mass,'prefix moments')
            if(mod(j,2)==1) call require(fast(j)==0,'exact centered odd parity')
         end do
      end do
   end do
   call require(all(accepted_count>0).and.all(fallback_count>0),'both bases cross gate')
   print '(a,i0,a,i0)','AXIS_RECURRENCE checks=',checks,' failures=',failed
   print '(a,2i6,a,2i6)','accepted Cheb/monomial=',accepted_count,' fallback=',fallback_count
   print '(a,es16.7,a,es16.7)','max_fast_mass_scaled_error=',max_fast_error, &
      ' max_binary64_mass_scaled_error=',max_output_error
   if(failed/=0) error stop 1
contains
   subroutine require(condition,message)
      logical,intent(in) :: condition
      character(*),intent(in) :: message
      checks=checks+1
      if(condition) return
      failed=failed+1
      if(failed<=8) print '(a)','FAIL '//message
   end subroutine

   subroutine quad_rule(xx,ww)
      ! Independent test quadrature with qk nodes/weights, not binary64
      ! quadrature promoted to qk and not a production spatial fallback.
      real(qk),intent(out) :: xx(:),ww(:)
      real(qk) :: z,p0,p1,p2,dp,step,pie
      integer :: n,i,k,iteration
      n=size(xx);pie=acos(-1.0_qk)
      do i=1,n/2
         z=cos(pie*(real(i,qk)-.25_qk)/(real(n,qk)+.5_qk))
         do iteration=1,80
            p0=1;p1=z
            do k=2,n
               p2=(real(2*k-1,qk)*z*p1-real(k-1,qk)*p0)/real(k,qk)
               p0=p1;p1=p2
            end do
            dp=real(n,qk)*(z*p1-p0)/(z*z-1);step=p1/dp;z=z-step
            if(abs(step)<=8*epsilon(1.0_qk)) exit
         end do
         if(iteration>80) error stop 'quad oracle root did not converge'
         p0=1;p1=z
         do k=2,n
            p2=(real(2*k-1,qk)*z*p1-real(k-1,qk)*p0)/real(k,qk)
            p0=p1;p1=p2
         end do
         dp=real(n,qk)*(z*p1-p0)/(z*z-1)
         xx(i)=-z;xx(n+1-i)=z;ww(i)=2/((1-z*z)*dp*dp);ww(n+1-i)=ww(i)
      end do
   end subroutine

   subroutine quadrature_oracle(ur,xr,ar,degree,chebyshev,result,positive_mass)
      real(rk),intent(in) :: ur,xr,ar
      integer,intent(in) :: degree
      logical,intent(in) :: chebyshev
      real(qk),intent(out) :: result(0:32),positive_mass
      real(qk) :: beta,delta,lo,hi,nearest,cut,half,mid,r,s,weight,c0,c1,c2,scale
      integer :: i,k
      beta=real(ur,qk)*real(ar,qk);delta=real(ur,qk)*real(xr,qk)
      lo=real(ur,qk)*(-real(ar,qk)-real(xr,qk));hi=real(ur,qk)*(real(ar,qk)-real(xr,qk))
      nearest=0
      if(lo>0) nearest=lo
      if(hi<0) nearest=-hi
      ! Oracle only: omitted positive tail has exponential ratio<=exp(-100),
      ! well below the tested qk budget, including representable tiny tails.
      cut=sqrt(nearest*nearest+100);lo=max(lo,-cut);hi=min(hi,cut)
      half=(hi-lo)/2;mid=(hi+lo)/2;result=0
      do i=1,size(nodes)
         r=mid+half*nodes(i);s=(r+delta)/beta
         weight=half*weights(i)*exp(-r*r)/beta
         c0=1;result(0)=result(0)+weight
         if(degree>=1) then
            c1=s;result(1)=result(1)+weight*c1
         end if
         do k=2,degree
            if(chebyshev) then
               c2=2*s*c1-c0
            else
               c2=s*c1
            end if
            result(k)=result(k)+weight*c2;c0=c1;c1=c2
         end do
      end do
      positive_mass=result(0);scale=real(ar,qk)
      do k=0,degree
         result(k)=scale*result(k)
         if(.not.chebyshev) scale=scale*real(ar,qk)
      end do
   end subroutine
end program
