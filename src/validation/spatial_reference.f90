! Independent spatial Poisson/heat-kernel reference. No SAS axis moments are used.
module validation_reference
  use iso_fortran_env, only: real64, real128
  implicit none
  private
  integer, parameter, public :: rk=real64
  real(rk), parameter :: pi=acos(-1.0_rk)
  type, public :: ref_periodic_context
    real(rk), allocatable :: r2(:), weight(:)
  end type
  public :: ref_rule, ref_cheb, ref_temperature, ref_moment, ref_deficit, ref_t90
  public :: ref_periodic, ref_periodic_derivative, ref_pgamma, ref_qgamma
  public :: ref_prepare_periodic, ref_periodic_value, ref_periodic_derivative_value
  public :: ref_operator
contains
  subroutine ref_rule(n,x,w)
    integer, intent(in) :: n
    real(rk), intent(out) :: x(n),w(n)
    real(real128) :: z,old,p0,p1,p2,dp,qpi
    integer :: i,j,it
    if(n<1) error stop 'ref_rule: positive order required'
    qpi=acos(-1.0_real128)
    do i=1,(n+1)/2
      z=cos(qpi*(real(i,real128)-0.25_real128)/(real(n,real128)+0.5_real128))
      do it=1,100
        p0=1.0_real128; p1=z
        do j=2,n
          p2=((2*j-1)*z*p1-(j-1)*p0)/j
          p0=p1; p1=p2
        end do
        dp=n*(z*p1-p0)/(z*z-1.0_real128)
        old=z; z=old-p1/dp
        if(abs(z-old)<=8*epsilon(z)) exit
      end do
      if(it>100) error stop 'ref_rule: Newton failed'
      ! Re-evaluate the derivative at the accepted node.
      p0=1.0_real128; p1=z
      do j=2,n
        p2=((2*j-1)*z*p1-(j-1)*p0)/j
        p0=p1; p1=p2
      end do
      dp=n*(z*p1-p0)/(z*z-1.0_real128)
      x(i)=real(-z,rk); x(n+1-i)=real(z,rk)
      w(i)=real(2.0_real128/((1-z*z)*dp*dp),rk); w(n+1-i)=w(i)
    end do
  end subroutine

  subroutine ref_cheb(m,n,c)
    integer, intent(in) :: m,n
    real(rk), intent(out) :: c(0:m)
    real(real128) :: theta,s,qpi
    integer :: j,k
    if(m<0.or.n<1) error stop 'ref_cheb: invalid order'
    qpi=acos(-1.0_real128)
    do k=0,m
      if(mod(k,2)==1) then
        c(k)=0; cycle
      end if
      s=0
      do j=1,n
        theta=(real(j,real128)-0.5_real128)*qpi/n
        s=s+exp(-4*cos(theta)**2)*cos(k*theta)
      end do
      if(k>0) s=2*s
      c(k)=real(s/n,rk)
    end do
  end subroutine

  elemental function ref_pgamma(q) result(p)
    real(rk), intent(in) :: q
    real(rk) :: p,term,s
    integer :: k
    if(q<=0) then
      p=0
    else if(q<2) then
      ! gamma(a,q)/Gamma(a) = exp(-q) q**a/Gamma(a+1) times
      ! [1+q/(a+1)+q**2/((a+1)(a+2))+...], a=3/2.
      term=1; s=1
      do k=1,200
        term=term*q/(1.5_rk+k); s=s+term
        if(term<=epsilon(s)*s) exit
      end do
      p=exp(-q)*q*sqrt(q)/(0.75_rk*sqrt(pi))*s
    else
      p=1-ref_qgamma(q)
    end if
  end function

  elemental function ref_qgamma(q) result(v)
    real(rk), intent(in) :: q
    real(rk) :: v
    if(q<=0) then
      v=1
    else
      v=erfc(sqrt(q))+2*sqrt(q/pi)*exp(-q)
    end if
  end function

  subroutine axis(lo,hi,n,cut,split,x,w)
    real(rk), intent(in) :: lo,hi,cut
    integer, intent(in) :: n
    logical, intent(in) :: split
    real(rk), allocatable, intent(out) :: x(:),w(:)
    real(rk) :: q(n),v(n),edges(3),a,b
    integer :: count,j,k
    if(hi<=lo) error stop 'reference: invalid bounds'
    count=1; edges(1)=lo; edges(2)=hi
    if(split.and.cut>lo.and.cut<hi) then
      count=2; edges(2)=cut; edges(3)=hi
    end if
    allocate(x(n*count),w(n*count)); call ref_rule(n,q,v)
    do j=1,count
      a=edges(j); b=edges(j+1); k=(j-1)*n
      x(k+1:k+n)=(a+b)/2+(b-a)*q/2; w(k+1:k+n)=(b-a)*v/2
    end do
  end subroutine

  subroutine geometry(x,y,z,bounds,n,split,xx,yy,r2,weight)
    real(rk), intent(in) :: x,y,z,bounds(4)
    integer, intent(in) :: n
    logical, intent(in) :: split
    real(rk), allocatable, intent(out) :: xx(:),yy(:),r2(:),weight(:)
    real(rk), allocatable :: q(:),r(:),w(:),v(:)
    integer :: i,j,k,nt
    if(z<=0) error stop 'reference: positive interior depth required'
    call axis(bounds(1),bounds(2),n,x,split,q,w)
    call axis(bounds(3),bounds(4),n,y,split,r,v)
    nt=size(q)*size(r); allocate(xx(nt),yy(nt),r2(nt),weight(nt)); k=0
    do j=1,size(r)
      do i=1,size(q)
        k=k+1; xx(k)=q(i); yy(k)=r(j)
        r2(k)=(q(i)-x)**2+(r(j)-y)**2+z*z
        weight(k)=z/(2*pi)*w(i)*v(j)/(r2(k)*sqrt(r2(k)))
      end do
    end do
  end subroutine

  pure function polynomial(x,c) result(v)
    real(rk), intent(in) :: x(:)
    real(rk), intent(in) :: c(0:)
    real(rk) :: v(size(x))
    integer :: k
    v=c(ubound(c,1))
    do k=ubound(c,1)-1,0,-1
      v=v*x+c(k)
    end do
  end function

  pure function chebyshev_value(x,c) result(v)
    real(rk), intent(in) :: x(:),c(0:)
    real(rk) :: v(size(x)),b0(size(x)),b1(size(x)),b2(size(x))
    integer :: k
    b1=0;b2=0
    do k=ubound(c,1),1,-1
      b0=2*x*b1-b2+c(k);b2=b1;b1=b0
    end do
    v=x*b1-b2+c(0)
  end function

  pure function compensated_sum(x) result(v)
    real(rk),intent(in) :: x(:)
    real(rk) :: v,c,t
    integer :: i
    v=0;c=0
    do i=1,size(x)
      t=v+x(i)
      if(abs(v)>=abs(x(i)))then
        c=c+((v-t)+x(i))
      else
        c=c+((x(i)-t)+v)
      end if
      v=t
    end do
    v=v+c
  end function

  function ref_temperature(x,y,z,bounds,cx,cy,t,n,split,gaussian,chebyshev,normalized) result(v)
    real(rk), intent(in) :: x,y,z,bounds(4),cx(0:),cy(0:),t
    integer, intent(in) :: n
    logical, intent(in) :: split
    logical, intent(in), optional :: gaussian,chebyshev,normalized
    real(rk) :: v
    real(rk), allocatable :: xx(:),yy(:),r2(:),w(:)
    call geometry(x,y,z,bounds,n,split,xx,yy,r2,w)
    if(present(normalized)) then
      if(normalized) then
        w=w*polynomial((2*xx-bounds(1)-bounds(2))/(bounds(2)-bounds(1)),cx)* &
          polynomial((2*yy-bounds(3)-bounds(4))/(bounds(4)-bounds(3)),cy)
      else
        w=w*polynomial(xx,cx)*polynomial(yy,cy)
      end if
    else if(present(chebyshev)) then
      if(chebyshev) then
        w=w*chebyshev_value(xx,cx)*chebyshev_value(yy,cy)
      else
        w=w*polynomial(xx,cx)*polynomial(yy,cy)
      end if
    else if(present(gaussian)) then
      if(gaussian) then
        w=w*exp(-4*(xx*xx+yy*yy))
      else
        w=w*polynomial(xx,cx)*polynomial(yy,cy)
      end if
    else
      w=w*polynomial(xx,cx)*polynomial(yy,cy)
    end if
    if(t==0) then
      v=0
    else if(t<0) then
      v=compensated_sum(w)
    else
      v=compensated_sum(w*ref_qgamma(r2/(4*t)))
    end if
  end function

  function ref_moment(x,y,z,bounds,m,n,order,split,shifted) result(v)
    ! The normalized Poisson moment includes z/(2*pi), matching manuscript J.
    real(rk), intent(in) :: x,y,z,bounds(4)
    integer, intent(in) :: m,n,order
    logical, intent(in) :: split,shifted
    real(rk) :: v
    real(rk), allocatable :: xx(:),yy(:),r2(:),w(:)
    if(min(m,n)<0) error stop 'ref_moment: negative exponent'
    call geometry(x,y,z,bounds,order,split,xx,yy,r2,w)
    if(shifted) then
      xx=xx-x; yy=yy-y
    end if
    v=compensated_sum(w*xx**m*yy**n)
  end function

  subroutine ref_operator(x,y,z,bounds,t,order,epsilon,perturbation,mass)
    real(rk), intent(in) :: x,y,z,bounds(4),t,epsilon
    integer, intent(in) :: order
    real(rk), intent(out) :: perturbation,mass
    real(rk), allocatable :: xx(:),yy(:),r2(:),w(:)
    call geometry(x,y,z,bounds,order,.true.,xx,yy,r2,w)
    if(t==0) then
      perturbation=0; mass=0; return
    else if(t>0) then
      w=w*ref_qgamma(r2/(4*t))
    end if
    mass=compensated_sum(w); perturbation=epsilon*compensated_sum(w*sin(3*xx-2*yy))
  end subroutine

  function ref_deficit(x,y,z,tp,alpha,bounds,n,split) result(v)
    real(rk), intent(in) :: x,y,z,tp,alpha,bounds(4)
    integer, intent(in) :: n
    logical, intent(in) :: split
    real(rk) :: v,L
    real(rk), allocatable :: xx(:),yy(:),r2(:),w(:),b(:),q(:)
    if(tp<0.or.alpha<=0) error stop 'ref_deficit: invalid time/diffusivity'
    if(tp==0) then
      v=0; return
    end if
    call geometry(x,y,z,bounds,n,split,xx,yy,r2,w)
    L=min(bounds(2)-bounds(1),bounds(4)-bounds(3))/2
    b=r2/(4*alpha); q=b/tp
    v=alpha/L**2*compensated_sum(w*(tp*ref_pgamma(q)+2*b*erfc(sqrt(q))))
  end function

  function ref_t90(z,a,b,profile,order) result(t)
    real(rk), intent(in) :: z,a,b
    character(*), intent(in) :: profile
    integer, intent(in) :: order
    real(rk) :: t,lo,hi,target,v
    real(rk), allocatable :: xx(:),yy(:),r2(:),w(:)
    integer :: k
    call geometry(0.0_rk,0.0_rk,z,[-a,a,-b,b],order,.true.,xx,yy,r2,w)
    if(profile=='Q') then
      w=w*(1-xx**2/a**2)*(1-yy**2/b**2)
    else if(profile/='U') then
      error stop 'ref_t90: profile must be U or Q'
    end if
    target=0.9_rk*compensated_sum(w); lo=0; hi=z*z+a*a+b*b
    do k=1,100
      if(compensated_sum(w*ref_qgamma(r2/(4*hi)))>=target) exit
      hi=2*hi
    end do
    if(k>100) error stop 'ref_t90: no bracket'
    do k=1,100
      t=(lo+hi)/2; v=compensated_sum(w*ref_qgamma(r2/(4*t)))
      if(v<target) then
        lo=t
      else
        hi=t
      end if
      if(hi-lo<=8*epsilon(t)*max(hi,tiny(t))) exit
    end do
    t=(lo+hi)/2
  end function

  subroutine ref_prepare_periodic(ctx,z,L,profile,order)
    type(ref_periodic_context), intent(out) :: ctx
    real(rk), intent(in) :: z,L
    character(*), intent(in) :: profile
    integer, intent(in) :: order
    real(rk), allocatable :: xx(:),yy(:)
    if(L<=0) error stop 'ref_periodic: positive L required'
    call geometry(0.0_rk,0.0_rk,z,[0.0_rk,L,0.0_rk,L],order,.false.,xx,yy,ctx%r2,ctx%weight)
    ctx%weight=4*ctx%weight
    if(profile=='Q') then
      ctx%weight=ctx%weight*(1-xx**2/L**2)*(1-yy**2/L**2)
    else if(profile/='U') then
      error stop 'ref_periodic: profile must be U or Q'
    end if
  end subroutine

  function ref_periodic_value(ctx,phase,duty,n_history) result(v)
    type(ref_periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: phase,duty
    integer, intent(in) :: n_history
    real(rk) :: v,a,b,qa,qb,s,comp,term,y,tmp
    integer :: i,k
    if(duty<0.or.duty>1.or.n_history<1) error stop 'ref_periodic: invalid history/duty'
    v=0
    if(duty==0) return
    do i=1,size(ctx%r2)
      s=0; comp=0
      do k=0,n_history-1
        a=phase+0.5_rk*k; b=a-0.5_rk*duty
        if(a<=0) cycle
        qa=ctx%r2(i)/(4*a)
        if(b<=0) then
          term=ref_qgamma(qa)
        else
          qb=ctx%r2(i)/(4*b)
          if(qa<2) then
            ! Late pulses: avoid subtracting two Q values close to unity.
            term=ref_pgamma(qb)-ref_pgamma(qa)
          else
            term=ref_qgamma(qa)-ref_qgamma(qb)
          end if
        end if
        y=term-comp; tmp=s+y; comp=(tmp-s)-y; s=tmp
      end do
      v=v+ctx%weight(i)*s
    end do
  end function

  function ref_periodic_derivative_value(ctx,phase,duty,n_history) result(v)
    type(ref_periodic_context), intent(in) :: ctx
    real(rk), intent(in) :: phase,duty
    integer, intent(in) :: n_history
    real(rk) :: v,a,b,s,qa,qb,ha,hb
    integer :: i,k
    if(duty<0.or.duty>1.or.n_history<1) error stop 'ref_periodic: invalid history/duty'
    v=0
    if(duty==0) return
    do i=1,size(ctx%r2)
      s=0
      do k=0,n_history-1
        a=phase+0.5_rk*k; b=a-0.5_rk*duty
        ha=0; hb=0
        if(a>0) then
          qa=ctx%r2(i)/(4*a)
          ha=2/sqrt(pi)*qa*sqrt(qa)*exp(-qa)/a
        end if
        if(b>0) then
          qb=ctx%r2(i)/(4*b)
          hb=2/sqrt(pi)*qb*sqrt(qb)*exp(-qb)/b
        end if
        s=s+(ha-hb)
      end do
      v=v+ctx%weight(i)*s
    end do
  end function

  function ref_periodic(phase,z,L,profile,duty,n_history,order) result(v)
    real(rk), intent(in) :: phase,z,L,duty
    character(*), intent(in) :: profile
    integer, intent(in) :: n_history,order
    real(rk) :: v
    type(ref_periodic_context) :: ctx
    call ref_prepare_periodic(ctx,z,L,profile,order)
    v=ref_periodic_value(ctx,phase,duty,n_history)
  end function

  function ref_periodic_derivative(phase,z,L,profile,duty,n_history,order) result(v)
    real(rk), intent(in) :: phase,z,L,duty
    character(*), intent(in) :: profile
    integer, intent(in) :: n_history,order
    real(rk) :: v
    type(ref_periodic_context) :: ctx
    call ref_prepare_periodic(ctx,z,L,profile,order)
    v=ref_periodic_derivative_value(ctx,phase,duty,n_history)
  end function
end module
