module validation_sas
   ! Accuracy verification only; performance measurements use separate modules.
   use iso_fortran_env, only: real64,real128
   use three_dimensional_analysis, only: integrate_gk15,lower_incomplete_gamma_indexed
   use gaussian_axes, only: monomial_axis
   implicit none
   private
   integer,parameter,public :: rk=real64,qk=real128
   real(rk),parameter :: pi=acos(-1.0_rk),sp=sqrt(pi)
   public :: sas_temperature,sas_center_step,center_integrand,sas_uniform_deficit,sas_gamma_product
contains
   function sas_temperature(x,y,z,t,bounds,cx,cy,eps,status) result(value)
      real(rk),intent(in) :: x,y,z,t,bounds(4),cx(0:),cy(0:),eps
      integer,intent(out) :: status
      real(rk) :: value,error,lower,ax,ay,xc,yc
      real(qk) :: px(0:size(cx)-1),py(0:size(cy)-1)
      integer :: i,j,axis_status
      ax=(bounds(2)-bounds(1))/2;ay=(bounds(4)-bounds(3))/2
      xc=(bounds(1)+bounds(2))/2;yc=(bounds(3)+bounds(4))/2
      px=0;py=0
      do i=0,size(cx)-1
         do j=0,i
            px(j)=px(j)+real(cx(i),qk)*binom(i,j)*real(xc,qk)**(i-j)
         end do
      end do
      do i=0,size(cy)-1
         do j=0,i
            py(j)=py(j)+real(cy(i),qk)*binom(i,j)*real(yc,qk)**(i-j)
         end do
      end do
      value=0;status=0;axis_status=0
      if(t<=0) return
      lower=z/(2*sqrt(t))
      if(lower>=28) return
      call integrate_gk15(integrand,lower,28.0_rk,eps,eps,value,error,status)
      if(axis_status/=0) status=axis_status
   contains
      function integrand(v) result(f)
         real(rk),intent(in) :: v
         real(rk) :: f,axv(0:size(cx)-1),ayv(0:size(cy)-1)
         real(qk) :: aq(0:size(cx)-1),bq(0:size(cy)-1),aa,bb
         integer :: s1,s2
         call monomial_axis(v/z,x-xc,ax,size(cx)-1,axv,s1,aq)
         call monomial_axis(v/z,y-yc,ay,size(cy)-1,ayv,s2,bq)
         if(s1/=0.or.s2/=0) then
            axis_status=1;f=0;return
         end if
         aa=sum(px*aq);bb=sum(py*bq)
         f=real(2*real(v,qk)**2*exp(-real(v,qk)**2)*aa*bb/(acos(-1.0_qk)**1.5_qk*real(z,qk)**2),rk)
      end function
   end function

   pure function binom(n,k) result(v)
      integer,intent(in) :: n,k
      real(qk) :: v
      integer :: i
      v=1
      do i=1,k
         v=v*real(n-i+1,qk)/real(i,qk)
      end do
   end function

   pure function qaxis(r,a) result(v)
      real(rk),intent(in) :: r,a
      real(rk) :: v,term,s
      integer :: k
      if(r<0.8_rk) then
         term=1;s=1.0_rk/3
         do k=1,22
            term=term*(-r*r)/k
            s=s+term/((2*k+1.0_rk)*(2*k+3.0_rk))
         end do
         v=4*a*s
      else
         v=a*(sp*erf(r)/r*(1-0.5_rk/(r*r))+exp(-r*r)/(r*r))
      end if
   end function

   pure function center_integrand(v,z,a,b,profile) result(f)
      real(rk),intent(in) :: v,z,a,b
      character,intent(in) :: profile
      real(rk) :: f,u
      u=v/z
      if(profile=='U') then
         f=2/sp*erf(a*u)*erf(b*u)*exp(-v*v)
      else
         f=2/(pi*sp*z*z)*v*v*exp(-v*v)*qaxis(a*u,a)*qaxis(b*u,b)
      end if
   end function

   function sas_center_step(t,z,a,b,profile,eps,status) result(value)
      real(rk),intent(in) :: t,z,a,b,eps
      character,intent(in) :: profile
      integer,intent(out) :: status
      real(rk) :: value,error,lower
      value=0;status=0
      if(t<=0) return
      lower=z/(2*sqrt(t))
      if(lower>=10) return
      call integrate_gk15(f,lower,10.0_rk,eps,eps,value,error,status)
   contains
      function f(v) result(vv)
         real(rk),intent(in) :: v
         real(rk) :: vv
         vv=center_integrand(v,z,a,b,profile)
      end function
   end function

   function sas_uniform_deficit(x,y,z,tp,alpha,bounds,eps,status) result(value)
      real(rk),intent(in) :: x,y,z,tp,alpha,bounds(4),eps
      integer,intent(out) :: status
      real(rk) :: value,v0,part1,part2,err,L
      integer :: s1,s2
      L=min(bounds(2)-bounds(1),bounds(4)-bounds(3))/2
      v0=z/sqrt(4*alpha*tp)
      call integrate_gk15(g,0.0_rk,min(v0,28.0_rk),eps,eps,part1,err,s1)
      part2=0;s2=0
      if(v0<28) call integrate_gk15(h,v0,28.0_rk,eps,eps,part2,err,s2)
      value=alpha*tp/(2*sp*L*L)*part1+z*z/(8*sp*L*L)*part2
      status=max(s1,s2)
   contains
      function g(v) result(f)
         real(rk),intent(in) :: v
         real(rk) :: f
         f=erfdiff((bounds(1)-x)*v/z,(bounds(2)-x)*v/z)* &
            erfdiff((bounds(3)-y)*v/z,(bounds(4)-y)*v/z)*exp(-v*v)
      end function
      function h(v) result(f)
         real(rk),intent(in) :: v
         real(rk) :: f
         f=g(v)/(v*v)
      end function
   end function

   pure function erfdiff(lo,hi) result(v)
      real(rk),intent(in) :: lo,hi
      real(rk) :: v
      if(lo>0) then
         v=erfc(lo)-erfc(hi)
      else if(hi<0) then
         v=erfc(-hi)-erfc(-lo)
      else
         v=erf(hi)-erf(lo)
      end if
   end function

   function sas_gamma_product(A,B,C,m,n,eps,status) result(value)
      real(rk),intent(in) :: A,B,C,eps
      integer,intent(in) :: m,n
      integer,intent(out) :: status
      real(rk) :: value,error
      integer :: gamma_status
      gamma_status=0
      call integrate_gk15(f,0.0_rk,28.0_rk,eps,eps,value,error,status)
      if(gamma_status/=0) status=gamma_status
   contains
      function f(v) result(fv)
         real(rk),intent(in) :: v
         real(rk) :: fv,u,g1,g2
         integer :: s1,s2
         u=v/C
         g1=lower_incomplete_gamma_indexed(m+1,(A*u)**2,s1)
         g2=lower_incomplete_gamma_indexed(n+1,(B*u)**2,s2)
         gamma_status=max(gamma_status,s1,s2)
         fv=g1*g2*exp(-v*v)/u**(m+n)/C
      end function
   end function
end module
