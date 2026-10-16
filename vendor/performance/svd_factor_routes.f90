module svd_factor_routes
   use replay_types
   use spatial_routes,only: spatial_nodes,radial_kernel
   implicit none
   private
   public :: poly_factor_t,prepare_poly_factor,evaluate_poly_factor
   type poly_factor_t
      type(source_t) :: source
      ! Binary64 power coefficients, degree by complete thin-SVD component.
      real(rk),allocatable :: px(:,:),py(:,:)
      integer :: rank=0
      logical :: ready=.false.
   end type
   interface
      subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
         import rk
         character,intent(in) :: jobu,jobvt
         integer,intent(in) :: m,n,lda,ldu,ldvt,lwork
         real(rk),intent(inout) :: a(lda,*),work(*)
         real(rk),intent(out) :: s(*),u(ldu,*),vt(ldvt,*)
         integer,intent(out) :: info
      end subroutine
   end interface
contains
   subroutine thin_factors(matrix,yfactor,xfactor,status)
      real(rk),intent(in) :: matrix(:,:)
      real(rk),allocatable,intent(out) :: yfactor(:,:),xfactor(:,:)
      integer,intent(out) :: status
      real(rk),allocatable :: copy(:,:),u(:,:),vt(:,:),s(:),work(:)
      real(rk) :: workspace(1)
      integer :: m,n,k,j,info,lwork
      status=invalid;m=size(matrix,1);n=size(matrix,2);k=min(m,n)
      if(k<1 .or. .not.all(ieee_is_finite(matrix))) return
      allocate(copy(m,n),u(m,k),vt(k,n),s(k));copy=matrix
      call dgesvd('S','S',m,n,copy,m,s,u,m,vt,k,workspace,-1,info)
      if(info/=0 .or. .not.ieee_is_finite(workspace(1))) then
         status=numerical_failure;return
      end if
      lwork=max(1,ceiling(workspace(1)));allocate(work(lwork));copy=matrix
      call dgesvd('S','S',m,n,copy,m,s,u,m,vt,k,work,lwork,info)
      if(info/=0) then
         status=numerical_failure;return
      end if
      allocate(yfactor(0:m-1,k),xfactor(0:n-1,k))
      do j=1,k
         yfactor(:,j)=u(:,j)*s(j)
         xfactor(:,j)=vt(j,:)
      end do
      status=ok
      if(.not.all(ieee_is_finite(yfactor)) .or. .not.all(ieee_is_finite(xfactor))) status=numerical_failure
   end subroutine

   subroutine factor_powers(factors,power)
      real(rk),intent(in) :: factors(0:,:)
      real(rk),allocatable,intent(out) :: power(:,:)
      real(rk) :: c0(0:size(factors,1)-1),c1(0:size(factors,1)-1),temp(0:size(factors,1)-1)
      integer :: m,j,k
      m=size(factors,1)-1;allocate(power(0:m,size(factors,2)));power=0
      do k=1,size(factors,2)
         if(m<=1) then
            power(:,k)=factors(:,k)
         else
            ! Chebyshev-to-power Clenshaw recurrence, entirely binary64.
            c0=0;c1=0;c0(0)=factors(m-1,k);c1(0)=factors(m,k)
            do j=m-2,0,-1
               temp=c0;c0=-c1;c0(0)=c0(0)+factors(j,k)
               c1(1:m)=2*c1(0:m-1);c1(0)=0;c1=c1+temp
            end do
            c0(1:m)=c0(1:m)+c1(0:m-1);power(:,k)=c0
         end if
      end do
   end subroutine

   subroutine prepare_poly_factor(source,prepared,status)
      type(source_t),intent(in) :: source
      type(poly_factor_t),intent(out) :: prepared
      integer,intent(out) :: status
      real(rk),allocatable :: yf(:,:),xf(:,:)
      prepared=poly_factor_t();status=invalid
      if(.not.source%ready .or. source%kind/=SRC_POLY .or. .not.valid_geometry(source)) return
      if(.not.allocated(source%c)) return
      ! Use the supplied rectangular matrix, including trailing zero components.
      call thin_factors(source%c,yf,xf,status)
      if(status/=ok) return
      prepared%source=source;prepared%rank=size(xf,2)
      call factor_powers(xf,prepared%px);call factor_powers(yf,prepared%py)
      if(.not.all(ieee_is_finite(prepared%px)) .or. .not.all(ieee_is_finite(prepared%py))) then
         status=numerical_failure;return
      end if
      prepared%ready=.true.;status=ok
   end subroutine

   subroutine axis_factors(nodes,width,power,values)
      real(rk),intent(in) :: nodes(:),width,power(0:,:)
      real(rk),allocatable,intent(out) :: values(:,:)
      real(rk) :: t,v
      integer :: i,j,k
      allocate(values(size(nodes),size(power,2)))
      do k=1,size(power,2)
         do i=1,size(nodes)
            t=nodes(i)/width;v=0
            do j=size(power,1)-1,0,-1
               v=v*t+power(j,k)
            end do
            values(i,k)=v
         end do
      end do
   end subroutine

   function evaluate_poly_factor(prepared,query,rx,ry,observable,status) result(value)
      type(poly_factor_t),intent(in) :: prepared
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: rx,ry
      integer,intent(in) :: observable
      integer,intent(out) :: status
      real(rk) :: value,x,y,z,dt,total,correction,term,f,kernel
      real(rk),allocatable :: xn(:),wx(:),yn(:),wy(:),fx(:,:),fy(:,:)
      integer :: ix,iy
      value=quiet_nan();status=invalid
      if(.not.prepared%ready .or. (observable/=OBS_T .and. observable/=OBS_D)) return
      call check_query(prepared%source,query,x,y,z,dt,status)
      if(status/=ok) return
      call spatial_nodes(prepared%source,x,y,z,dt,rx,ry,xn,wx,yn,wy,status)
      if(status/=ok) return
      call axis_factors(xn,prepared%source%a,prepared%px,fx)
      call axis_factors(yn,prepared%source%b,prepared%py,fy)
      total=0;correction=0
      do iy=1,size(yn)
         do ix=1,size(xn)
            f=sum(fy(iy,:)*fx(ix,:))
            kernel=radial_kernel((xn(ix)-x)**2+(yn(iy)-y)**2+z*z,z,dt,prepared%source%alpha,observable)
            term=wx(ix)*wy(iy)*f*kernel
            if(.not.ieee_is_finite(term)) then
               status=numerical_failure;return
            end if
            call add_compensated(total,correction,term)
         end do
      end do
      value=total+correction;status=ok
      if(.not.ieee_is_finite(value)) status=numerical_failure
   end function

end module
