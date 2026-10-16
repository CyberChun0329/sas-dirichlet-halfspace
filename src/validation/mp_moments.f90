module mp_moments
   ! Manuscript line/J recurrence; no spatial quadrature or spatial fallback.
   use mp_arithmetic
   implicit none
   private
   public :: qk,mp_steady,mp_gaussian_coeff,mp_moment_table,mp_gaussian_coeff_mp
contains
   subroutine line_moments(m,lo,hi,rho,a)
      integer,intent(in) :: m
      type(mp_real),intent(in) :: lo,hi,rho
      type(mp_real),intent(out) :: a(0:m)
      type(mp_real) :: rlo,rhi
      integer :: j
      rlo=sqrt(lo*lo+rho*rho);rhi=sqrt(hi*hi+rho*rho)
      a(0)=asinh(hi/rho)-asinh(lo/rho)
      if(m>=1) a(1)=rhi-rlo
      do j=2,m
         a(j)=(hi**(j-1)*rhi-lo**(j-1)*rlo-(j-1)*rho*rho*a(j-2))/j
      end do
   end subroutine

   subroutine mp_moment_table(x,y,z,bounds,mx,my,table)
      ! Inputs/outputs belong to the current mp_arithmetic arena. This
      ! routine does not reset it, so two precisions can coexist for subtraction.
      type(mp_real),intent(in) :: x,y,z,bounds(4)
      integer,intent(in) :: mx,my
      type(mp_real),intent(out) :: table(0:mx,0:my)
      type(mp_real) :: u(2),v(2),g(0:mx,2),h(0:my,2),a,b,pi2
      integer :: i,j,m,n,si,sj
      if(mx<0.or.my<0) error stop 'Negative moment degree'
      if(mp_compare(z,mp_integer(0))<=0) error stop 'Nonpositive moment z'
      if(mp_compare(bounds(2),bounds(1))<=0.or.mp_compare(bounds(4),bounds(3))<=0) &
         error stop 'Invalid moment bounds'
      u=[bounds(1)-x,bounds(2)-x];v=[bounds(3)-y,bounds(4)-y]
      pi2=2*mp_pi()
      do i=1,2
         call line_moments(mx,u(1),u(2),sqrt(v(i)**2+z*z),g(:,i))
         call line_moments(my,v(1),v(2),sqrt(u(i)**2+z*z),h(:,i))
      end do
      table=mp_integer(0)
      do i=1,2
         si=2*i-3
         do j=1,2
            sj=2*j-3
            table(0,0)=table(0,0)+si*sj*atan(u(i)*v(j)/(z*sqrt(u(i)**2+v(j)**2+z*z)))/pi2
            if(mx>=1) table(1,0)=table(1,0)-z/pi2*si*sj*asinh(v(j)/sqrt(u(i)**2+z*z))
         end do
      end do
      do m=2,mx
         a=-v(1)*g(m-2,1)+v(2)*g(m-2,2)
         table(m,0)=z/pi2*a-z*z*table(m-2,0)
      end do
      do m=0,mx
         if(my>=1) table(m,1)=-z/pi2*(-g(m,1)+g(m,2))
         do n=2,my
            a=-u(1)**(m+1)*h(n-2,1)+u(2)**(m+1)*h(n-2,2)
            b=-v(1)**(n-1)*g(m,1)+v(2)**(n-1)*g(m,2)
            table(m,n)=((n-1)*z*a-m*z*b)/(pi2*(m+n-1))-(n-1)*z*z*table(m,n-2)/(m+n-1)
         end do
      end do
   end subroutine

   function binomial(n,k) result(v)
      integer,intent(in) :: n,k
      type(mp_real) :: v
      integer :: i
      v=mp_integer(1)
      do i=1,min(k,n-k)
         v=v*(n-i+1)/i
      end do
   end function

   real(qk) function mp_steady(x,y,z,bounds,cx,cy,bits,shifted) result(value)
      ! qk arguments retain their binary128 input precision through a 37 digit
      ! decimal conversion. Physical-power coefficients are the default.
      ! shifted=.true. contracts cx(k)*cy(l)*J(k,l) without binomial shifting.
      ! This wrapper resets the arena on entry AND exit; do not retain handles
      ! across it or call it concurrently. Only the final value is rounded qk.
      real(qk),intent(in) :: x,y,z,bounds(4),cx(0:),cy(0:)
      integer,intent(in) :: bits
      logical,optional,intent(in) :: shifted
      type(mp_real) :: xx,yy,zz,bd(4),a(0:ubound(cx,1)),b(0:ubound(cy,1))
      type(mp_real) :: j(0:ubound(cx,1),0:ubound(cy,1)),total
      integer :: i,k,l,m,n,mx,my
      logical :: direct
      call mp_begin(bits)
      mx=ubound(cx,1);my=ubound(cy,1)
      xx=mp_from_qk(x);yy=mp_from_qk(y);zz=mp_from_qk(z)
      do i=1,4
         bd(i)=mp_from_qk(bounds(i))
      end do
      call mp_moment_table(xx,yy,zz,bd,mx,my,j)
      direct=.false.;if(present(shifted)) direct=shifted
      a=mp_integer(0);b=mp_integer(0)
      do k=0,mx
         if(direct) then
            a(k)=mp_from_qk(cx(k))
         else
            do m=k,mx
               a(k)=a(k)+mp_from_qk(cx(m))*binomial(m,k)*xx**(m-k)
            end do
         end if
      end do
      do l=0,my
         if(direct) then
            b(l)=mp_from_qk(cy(l))
         else
            do n=l,my
               b(l)=b(l)+mp_from_qk(cy(n))*binomial(n,l)*yy**(n-l)
            end do
         end if
      end do
      total=mp_integer(0)
      do k=0,mx
         do l=0,my
            total=total+a(k)*b(l)*j(k,l)
         end do
      end do
      if(.not.mp_finite(total)) error stop 'Nonfinite steady recurrence'
      value=mp_to_qk(total)
      call mp_reset()
   end function

   function bessel_i_at_two(k) result(s)
      integer,intent(in) :: k
      type(mp_real) :: s,term,tolerance
      integer :: j
      term=mp_integer(1)
      do j=1,k
         term=term/j
      end do
      s=term;tolerance=mp_integer(2)**(-mp_bits()-16)
      do j=1,10000
         term=term/j/(j+k)
         s=s+term
         if(mp_compare(abs(term),abs(s)*tolerance)<0) return
      end do
      error stop 'Bessel series did not converge'
   end function

   subroutine mp_gaussian_coeff_mp(m,cx)
      ! exp(-4 X^2) = exp(-2)[I0(2)+2 sum (-1)^k Ik(2) T2k(X)].
      ! Bessel Ik(2) is the convergent sum 1/(j! (j+k)!).
      ! No fitting or binary64 constants enter coefficient construction.
      integer,intent(in) :: m
      type(mp_real),intent(out) :: cx(0:m)
      type(mp_real) :: t(0:m,0:m),co,scale
      integer :: j,k
      if(m<0) error stop 'Negative Gaussian order'
      t=mp_integer(0);t(0,0)=mp_integer(1)
      if(m>=1) t(1,1)=mp_integer(1)
      do k=2,m
         do j=1,k
            t(j,k)=2*t(j-1,k-1)
         end do
         do j=0,k-2
            t(j,k)=t(j,k)-t(j,k-2)
         end do
      end do
      scale=exp(mp_integer(-2));cx=mp_integer(0)
      do k=0,m/2
         co=scale*bessel_i_at_two(k)
         if(k>0) co=2*(-1)**k*co
         do j=0,2*k
            cx(j)=cx(j)+co*t(j,2*k)
         end do
      end do
   end subroutine

   subroutine mp_gaussian_coeff(m,bits,cx)
      integer,intent(in) :: m,bits
      real(qk),intent(out) :: cx(0:m)
      type(mp_real) :: high(0:m)
      integer :: j
      call mp_begin(bits)
      call mp_gaussian_coeff_mp(m,high)
      do j=0,m
         cx(j)=mp_to_qk(high(j))
      end do
      call mp_reset()
   end subroutine
end module
