module replay_types
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use temperature_common, only: rk,ok,invalid,numerical_failure,pi,sqrtpi,quiet_nan, &
      gauss_rule_t,init_rule,mapped_rule,valid_rule
   implicit none
   integer,parameter :: qk=selected_real_kind(30,300)
   integer,parameter :: SRC_POLY=1,SRC_GAUSSIAN=2,OBS_T=1,OBS_D=2
   type source_t
      integer :: kind=0,mx=0,my=0
      real(rk) :: a=0,b=0,alpha=1,cx=0,cy=0,beta(2)=[1.0_rk,1.0_rk]
      real(rk),allocatable :: c(:,:)
      logical :: ready=.false.
   end type
contains
   subroutine init_poly(s,c,a,b,alpha,status,cx,cy)
      type(source_t),intent(out) :: s
      real(rk),intent(in) :: c(0:,0:),a,b,alpha
      integer,intent(out) :: status
      real(rk),optional,intent(in) :: cx,cy
      integer :: p,q
      s=source_t();status=invalid
      if(size(c,1)<1 .or. size(c,2)<1 .or. maxval(shape(c))>33) return
      if(.not.all(ieee_is_finite(c))) return
      s%kind=SRC_POLY;s%a=a;s%b=b;s%alpha=alpha
      if(present(cx)) s%cx=cx
      if(present(cy)) s%cy=cy
      allocate(s%c(0:size(c,1)-1,0:size(c,2)-1));s%c=c
      do p=0,ubound(s%c,2)
         do q=0,ubound(s%c,1)
            if(c(q,p)/=0) then
               s%mx=max(s%mx,p);s%my=max(s%my,q)
            end if
         end do
      end do
      if(.not.valid_geometry(s)) return
      s%ready=.true.;status=ok
   end subroutine

   subroutine init_gaussian(s,beta,a,b,alpha,status)
      type(source_t),intent(out) :: s
      real(rk),intent(in) :: beta(2),a,b,alpha
      integer,intent(out) :: status
      s=source_t();status=invalid
      if(.not.all(ieee_is_finite(beta)) .or. any(beta<=0)) return
      s%kind=SRC_GAUSSIAN;s%beta=beta;s%a=a;s%b=b;s%alpha=alpha
      if(.not.valid_geometry(s)) return
      s%ready=.true.;status=ok
   end subroutine

   pure logical function valid_geometry(s)
      type(source_t),intent(in) :: s
      valid_geometry=all(ieee_is_finite([s%a,s%b,s%alpha,s%cx,s%cy]))
      if(valid_geometry) valid_geometry=min(s%a,s%b,s%alpha)>0
   end function

   subroutine check_query(s,query,x,y,z,dt,status)
      type(source_t),intent(in) :: s
      real(rk),intent(in) :: query(4)
      real(rk),intent(out) :: x,y,z,dt
      integer,intent(out) :: status
      status=invalid;x=quiet_nan();y=x;z=x;dt=x
      if(.not.s%ready .or. .not.valid_geometry(s)) return
      if(.not.all(ieee_is_finite(query)) .or. min(query(3),query(4))<=0) return
      x=query(1)-s%cx;y=query(2)-s%cy;z=query(3)
      if(s%alpha>1) then
         if(query(4)>huge(dt)/s%alpha) return
      end if
      dt=query(4)*s%alpha
      if(.not.all(ieee_is_finite([x,y,z,dt])) .or. dt<tiny(dt)) return
      status=ok
   end subroutine

   pure function source_value(s,x,y) result(v)
      type(source_t),intent(in) :: s
      real(rk),intent(in) :: x,y
      real(rk) :: v,Xr,Yr,tx(0:32),ty(0:32)
      integer :: j,p,q
      v=0;Xr=x-s%cx;Yr=y-s%cy
      if(abs(Xr)>s%a .or. abs(Yr)>s%b) return
      select case(s%kind)
      case(SRC_POLY)
         Xr=Xr/s%a;Yr=Yr/s%b;tx=0;ty=0;tx(0)=1;ty(0)=1
         if(s%mx>0) tx(1)=Xr
         if(s%my>0) ty(1)=Yr
         do j=2,s%mx
            tx(j)=2*Xr*tx(j-1)-tx(j-2)
         end do
         do j=2,s%my
            ty(j)=2*Yr*ty(j-1)-ty(j-2)
         end do
         do p=0,s%mx
            do q=0,s%my
               v=v+s%c(q,p)*tx(p)*ty(q)
            end do
         end do
      case(SRC_GAUSSIAN)
         v=exp(-(s%beta(1)*Xr/s%a)**2-(s%beta(2)*Yr/s%b)**2)
      end select
   end function

   subroutine profile_source(name,a,b,s,status)
      character(*),intent(in) :: name
      real(rk),intent(in) :: a,b
      type(source_t),intent(out) :: s
      integer,intent(out) :: status
      real(rk) :: c(0:2,0:2),Bmat(0:2,0:2),g(0:2,0:2)
      integer :: p,q,i,j
      c=0
      select case(trim(name))
      case('U')
         c(0,0)=1
      case('Q')
         c(0,0)=.25_rk;c(0,2)=-.25_rk;c(2,0)=-.25_rk;c(2,2)=.25_rk
      case('N1')
         c(0,0)=1/1.55_rk;c(0,1)=.2_rk/1.55_rk;c(1,0)=.15_rk/1.55_rk;c(1,1)=.2_rk/1.55_rk
      case('N2')
         c(0,0)=.775_rk;c(0,2)=-.125_rk;c(2,0)=-.075_rk;c(1,1)=.25_rk;c(2,2)=.025_rk
      case('B1','B2')
         if(trim(name)=='B1') then
            Bmat(0,:)=[1.0_rk,.45_rk,.25_rk];Bmat(1,:)=[.6_rk,.35_rk,.8_rk];Bmat(2,:)=[.2_rk,.7_rk,.9_rk]
         else
            Bmat(0,:)=[.3_rk,.9_rk,1.0_rk];Bmat(1,:)=[.8_rk,.25_rk,.45_rk];Bmat(2,:)=[.7_rk,.6_rk,.2_rk]
         end if
         g(0,:)=[.375_rk,-.5_rk,.125_rk];g(1,:)=[.25_rk,0.0_rk,-.25_rk];g(2,:)=[.375_rk,.5_rk,.125_rk]
         do p=0,2
            do q=0,2
               do i=0,2
                  do j=0,2
                     c(q,p)=c(q,p)+Bmat(i,j)*g(i,q)*g(j,p)
                  end do
               end do
            end do
         end do
      case default
         s=source_t();status=invalid;return
      end select
      call init_poly(s,c,a,b,1.0_rk,status)
   end subroutine

   subroutine add_compensated(total,correction,term)
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
