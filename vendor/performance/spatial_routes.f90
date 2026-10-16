module spatial_routes
   use replay_types
   use gaussian_axes, only: monomial_axis,gaussian_profile_axis
   implicit none
   private
   type,public :: spatial_t
      private
      type(source_t) :: source
      real(rk),allocatable :: power(:,:)
      real(qk),allocatable :: physical_power(:,:)
      logical :: ready=.false.,low_order=.false.
   end type
   public :: prepare_spatial,evaluate_spatial,evaluate_shared,evaluate_steady
   public :: evaluate_spatial_complementary,space_operator,spatial_precision
   public :: spatial_nodes,radial_kernel
contains
   subroutine prepare_spatial(source,prepared,status)
      type(source_t),intent(in) :: source
      type(spatial_t),intent(out) :: prepared
      integer,intent(out) :: status
      real(qk),allocatable :: tx(:,:),ty(:,:),c(:,:)
      integer :: p,q,i,j,mx,my
      prepared=spatial_t();status=invalid
      if(.not.source%ready .or. .not.valid_geometry(source)) return
      select case(source%kind)
      case(SRC_POLY)
         if(.not.allocated(source%c)) return
         call init_poly(prepared%source,source%c,source%a,source%b,source%alpha,status,source%cx,source%cy)
      case(SRC_GAUSSIAN)
         call init_gaussian(prepared%source,source%beta,source%a,source%b,source%alpha,status)
         prepared%source%cx=source%cx;prepared%source%cy=source%cy
      case default
         return
      end select
      if(status/=ok) return
      if(source%kind==SRC_POLY) then
         mx=prepared%source%mx;my=prepared%source%my
         prepared%low_order=max(mx,my)<=2
         allocate(tx(0:mx,0:mx),ty(0:my,0:my),c(0:my,0:mx))
         call power_basis(mx,tx);call power_basis(my,ty);c=0
         ! Full tensor conversion, retaining every supplied nonzero coefficient.
         do p=0,mx
            do q=0,my
               do i=0,mx
                  do j=0,my
                     c(j,i)=c(j,i)+real(prepared%source%c(q,p),qk)*tx(i,p)*ty(j,q)
                  end do
               end do
            end do
         end do
         allocate(prepared%power(0:my,0:mx),prepared%physical_power(0:my,0:mx))
         prepared%power=real(c,rk)
         do p=0,mx
            do q=0,my
               prepared%physical_power(q,p)=c(q,p)/(real(source%a,qk)**p*real(source%b,qk)**q)
            end do
         end do
         if(.not.all(ieee_is_finite(prepared%physical_power))) then
            status=numerical_failure;return
         end if
      end if
      prepared%ready=.true.;status=ok
   end subroutine

   subroutine power_basis(m,t)
      integer,intent(in) :: m
      real(qk),intent(out) :: t(0:m,0:m)
      integer :: j
      t=0;t(0,0)=1
      if(m>=1) t(1,1)=1
      do j=2,m
         t(1:j,j)=2*t(0:j-1,j-1)
         t(:,j)=t(:,j)-t(:,j-2)
      end do
   end subroutine

   integer function spatial_precision(prepared) result(bits)
      type(spatial_t),intent(in) :: prepared
      ! Precision of analytic/shared intermediates; actual 2D quadrature is rk.
      bits=64
      if(prepared%source%kind==SRC_POLY .and. .not.prepared%low_order) bits=128
      if(prepared%source%kind==SRC_GAUSSIAN) bits=128
   end function

   function evaluate_spatial(prepared,query,rx,ry,observable,status) result(value)
      type(spatial_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rx,ry
      integer,intent(in) :: observable
      integer,intent(out) :: status
      real(rk) :: value
      value=integrate_spatial(prepared,query,rx,ry,observable,status)
   end function

   function evaluate_steady(prepared,query,rx,ry,status) result(value)
      type(spatial_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rx,ry
      integer,intent(out) :: status
      real(rk) :: value
      value=integrate_spatial(prepared,query,rx,ry,0,status)
   end function

   function evaluate_spatial_complementary(prepared,query,rx,ry,observable,status) result(value)
      type(spatial_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rx,ry
      integer,intent(in) :: observable
      integer,intent(out) :: status
      real(rk) :: value,steady,correction
      value=quiet_nan();status=invalid
      if(observable/=OBS_T .and. observable/=OBS_D) return
      ! Fresh steady and complementary integrals; no stored query answers.
      steady=integrate_spatial(prepared,query,rx,ry,0,status)
      if(status/=ok) return
      correction=integrate_spatial(prepared,query,rx,ry,observable+2,status)
      if(status/=ok) return
      if(observable==OBS_D) steady=steady*query(4)
      value=steady-correction
      if(.not.ieee_is_finite(value)) status=numerical_failure
   end function

   function integrate_spatial(prepared,query,rx,ry,mode,status) result(value)
      type(spatial_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rx,ry
      integer,intent(in) :: mode
      integer,intent(out) :: status
      real(rk) :: value,x,y,z,dt,f,k,total,correction,term
      real(rk),allocatable :: xn(:),wx(:),yn(:),wy(:)
      integer :: ix,iy
      value=quiet_nan();status=invalid
      if(mode<0 .or. mode>4 .or. .not.prepared%ready) return
      call check_query(prepared%source,query,x,y,z,dt,status)
      if(status/=ok) return
      call spatial_nodes(prepared%source,x,y,z,dt,rx,ry,xn,wx,yn,wy,status)
      if(status/=ok) return
      total=0;correction=0
      do iy=1,size(yn)
         do ix=1,size(xn)
            f=prepared_value(prepared,xn(ix),yn(iy))
            k=radial_kernel((xn(ix)-x)**2+(yn(iy)-y)**2+z*z,z,dt,prepared%source%alpha,mode)
            term=wx(ix)*wy(iy)*f*k
            if(.not.ieee_is_finite(term)) then
               status=numerical_failure;return
            end if
            call add_compensated(total,correction,term)
         end do
      end do
      value=total+correction;status=ok
      if(.not.ieee_is_finite(value)) status=numerical_failure
   end function

   function prepared_value(p,x,y) result(f)
      type(spatial_t),intent(in) :: p
      real(rk),intent(in) :: x,y
      real(rk) :: f,row,Xr,Yr
      integer :: i,j
      if(p%source%kind==SRC_POLY .and. p%low_order) then
         Xr=x/p%source%a;Yr=y/p%source%b;f=0
         do j=p%source%my,0,-1
            row=0
            do i=p%source%mx,0,-1
               row=row*Xr+p%power(j,i)
            end do
            f=f*Yr+row
         end do
      else
         ! High degrees use the complete Chebyshev series directly in 2D.
         f=source_value(p%source,x+p%source%cx,y+p%source%cy)
      end if
   end function

   subroutine spatial_nodes(s,x,y,z,dt,rx,ry,xn,wx,yn,wy,status)
      type(source_t),intent(in) :: s
      real(rk),intent(in) :: x,y,z,dt
      type(gauss_rule_t),intent(in) :: rx,ry
      real(rk),allocatable,intent(out) :: xn(:),wx(:),yn(:),wy(:)
      integer,intent(out) :: status
      real(rk),allocatable :: cx(:),cy(:)
      real(rk) :: bx(9),by(9)
      bx=[x-z/4,x+z/4,x-z,x+z,x-4*z,x+4*z,x-2*sqrt(dt),x+2*sqrt(dt),x]
      by=[y-z/4,y+z/4,y-z,y+z,y-4*z,y+4*z,y-2*sqrt(dt),y+2*sqrt(dt),y]
      cx=bx;cy=by
      call mapped_rule(-s%a,s%a,cx,rx,xn,wx,status)
      if(status/=ok) return
      call mapped_rule(-s%b,s%b,cy,ry,yn,wy,status)
   end subroutine

   pure function lower_gamma_three_halves(s) result(p)
      real(rk),intent(in) :: s
      real(rk) :: p,term,total
      integer :: j
      if(s<1) then
         ! P(3/2,s^2)=2/sqrt(pi)*sum (-1)^j s^(2j+3)/(j!*(j+3/2)).
         term=2*s**3/(1.5_rk*sqrtpi);total=term
         do j=1,80
            term=term*(-s*s)/real(j,rk)*(real(j,rk)+.5_rk)/(real(j,rk)+1.5_rk)
            total=total+term
            if(abs(term)<=epsilon(total)*abs(total)) exit
         end do
         p=total
      else if(s>27) then
         p=1
      else
         p=erf(s)-2*s*exp(-s*s)/sqrtpi
      end if
   end function

   pure function radial_kernel(r2,z,dt,alpha,mode) result(k)
      real(rk),intent(in) :: r2,z,dt,alpha
      integer,intent(in) :: mode
      real(rk) :: k,r,s,poisson,p,q,e,c
      k=quiet_nan()
      r=sqrt(r2);poisson=z/(2*pi*r2*r)
      if(mode==0) then
         k=poisson;return
      end if
      s=r/(2*sqrt(dt));p=lower_gamma_three_halves(s)
      if(s>27) then
         e=0;q=0
      else
         e=erfc(s);q=e+2*s*exp(-s*s)/sqrtpi
      end if
      c=z*e/(4*pi*alpha*r)
      select case(mode)
      case(OBS_T)
         k=poisson*q
      case(OBS_D)
         k=(dt/alpha)*poisson*p+c
      case(3)
         k=poisson*p
      case(4)
         k=(dt/alpha)*poisson*q-c
      end select
   end function

   function evaluate_shared(prepared,query,rule,observable,status) result(value)
      type(spatial_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rule
      integer,intent(in) :: observable
      integer,intent(out) :: status
      real(rk) :: value,x,y,z,dt,v0,lo,h,span,u,term,total,correction,factor,axg,ayg
      real(rk),allocatable :: v(:),w(:),ax(:),ay(:)
      real(qk),allocatable :: axq(:),ayq(:)
      real(qk) :: cq,totalq,termq
      integer :: i,p,q,mx,my
      value=quiet_nan();status=invalid
      if(.not.prepared%ready .or. .not.valid_rule(rule)) return
      if(observable/=OBS_T .and. observable/=OBS_D) return
      call check_query(prepared%source,query,x,y,z,dt,status)
      if(status/=ok) return
      v0=z/(2*sqrt(dt));lo=0
      if(observable==OBS_T) lo=v0
      if(lo>=12) then
         value=0;status=ok;return
      end if
      span=max(prepared%source%a+abs(x),prepared%source%b+abs(y));h=z/span
      call mapped_rule(lo,12.0_rk,[h/4,h,4*h,16*h,.25_rk,1.0_rk,3.0_rk,6.0_rk,v0],rule,v,w,status)
      if(status/=ok) return
      mx=prepared%source%mx;my=prepared%source%my
      allocate(ax(0:mx),ay(0:my),axq(0:mx),ayq(0:my))
      total=0;correction=0;totalq=0
      do i=1,size(v)
         u=v(i)/z
         factor=w(i)*2/(pi*sqrtpi)*exp(-v(i)*v(i))
         if(observable==OBS_D) factor=factor*min(query(4),z*z/(4*prepared%source%alpha*v(i)*v(i)))
         if(prepared%source%kind==SRC_GAUSSIAN) then
            axg=gaussian_profile_axis(u,x,prepared%source%a,prepared%source%beta(1),status)
            if(status/=ok) return
            ayg=gaussian_profile_axis(u,y,prepared%source%b,prepared%source%beta(2),status)
            if(status/=ok) return
            term=factor*(u*axg)*(u*ayg)
         else if(prepared%low_order) then
            call monomial_axis(u,x,prepared%source%a,mx,ax,status)
            if(status/=ok) return
            call monomial_axis(u,y,prepared%source%b,my,ay,status)
            if(status/=ok) return
            term=0
            do p=0,mx
               do q=0,my
                  term=term+real(prepared%physical_power(q,p),rk)*(u*ax(p))*(u*ay(q))
               end do
            end do
            term=factor*term
         else
            call monomial_axis(u,x,prepared%source%a,mx,ax,status,axq)
            if(status/=ok) return
            call monomial_axis(u,y,prepared%source%b,my,ay,status,ayq)
            if(status/=ok) return
            cq=0
            do p=0,mx
               do q=0,my
                  cq=cq+prepared%physical_power(q,p)*(real(u,qk)*axq(p))*(real(u,qk)*ayq(q))
               end do
            end do
            termq=real(factor,qk)*cq;totalq=totalq+termq
            if(.not.ieee_is_finite(termq)) then
               status=numerical_failure;return
            end if
            cycle
         end if
         if(.not.ieee_is_finite(term)) then
            status=numerical_failure;return
         end if
         call add_compensated(total,correction,term)
      end do
      value=total+correction+real(totalq,rk);status=ok
      if(.not.ieee_is_finite(value)) status=numerical_failure
      ! Omitted v>12 tail: <=max|f|*erfc(12) for T; multiply by physical tp for D.
   end function

   subroutine space_operator(prepared,queries,rx,ry,observable,W,status)
      type(spatial_t),intent(in) :: prepared
      real(rk),intent(in) :: queries(:,:)
      type(gauss_rule_t),intent(in) :: rx,ry
      integer,intent(in) :: observable
      real(rk),allocatable,intent(out) :: W(:,:,:)
      integer,intent(out) :: status
      real(rk),allocatable :: xn(:),wx(:),yn(:),wy(:),comp(:,:,:)
      real(rk) :: x,y,z,dt,k,tx(0:32),ty(0:32)
      integer :: nq,j,ix,iy,p,q
      status=invalid
      if(.not.prepared%ready .or. size(queries,1)/=4) return
      if(observable/=OBS_T .and. observable/=OBS_D) return
      nq=size(queries,2)
      select case(prepared%source%kind)
      case(SRC_POLY)
         allocate(W(0:ubound(prepared%source%c,1),0:ubound(prepared%source%c,2),nq))
      case default
         return
      end select
      allocate(comp(0:ubound(W,1),0:ubound(W,2),nq));W=0;comp=0
      do j=1,nq
         call check_query(prepared%source,queries(:,j),x,y,z,dt,status)
         if(status/=ok) return
         call spatial_nodes(prepared%source,x,y,z,dt,rx,ry,xn,wx,yn,wy,status)
         if(status/=ok) return
         do iy=1,size(yn)
            do ix=1,size(xn)
               k=wx(ix)*wy(iy)*radial_kernel((xn(ix)-x)**2+(yn(iy)-y)**2+z*z, &
                                            z,dt,prepared%source%alpha,observable)
               call cheb_values(xn(ix)/prepared%source%a,ubound(W,2),tx)
               call cheb_values(yn(iy)/prepared%source%b,ubound(W,1),ty)
               do p=0,ubound(W,2)
                  do q=0,ubound(W,1)
                     call add_compensated(W(q,p,j),comp(q,p,j),k*tx(p)*ty(q))
                  end do
               end do
            end do
         end do
      end do
      W=W+comp;status=ok
      if(.not.all(ieee_is_finite(W))) status=numerical_failure
   end subroutine

   subroutine cheb_values(x,m,t)
      real(rk),intent(in) :: x
      integer,intent(in) :: m
      real(rk),intent(out) :: t(0:32)
      integer :: i
      t=0;t(0)=1
      if(m>0) t(1)=x
      do i=2,m
         t(i)=2*x*t(i-1)-t(i-2)
      end do
   end subroutine
end module
