module temperature_common
   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: iso_c_binding, only: c_double
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   integer, parameter, public :: rk = real64
   integer, parameter, public :: ok = 0, invalid = 1, numerical_failure = 2
   real(rk), parameter, public :: pi = acos(-1.0_rk), sqrtpi = sqrt(pi)

   ! C(y degree,x degree) is an ordinary, unprimed tensor Chebyshev series.
   ! Coefficients represent the prescribed temperature, not heat flux.
   type, public :: patch_source_t
      real(rk) :: c(0:2,0:2) = 0.0_rk
      real(rk) :: a = 0.0_rk, b = 0.0_rk, alpha = 1.0_rk
      real(rk) :: cx = 0.0_rk, cy = 0.0_rk
      integer :: mx = 0, my = 0
      logical :: ready = .false.
   end type

   ! Rules are caller-owned; source preparation does not hide node costs.
   type, public :: gauss_rule_t
      integer :: n = 0
      real(rk), allocatable :: x(:), w(:)
      logical :: ready = .false.
   end type

   type, public :: query_work_t
      integer :: v_nodes = 0, x_nodes = 0, y_nodes = 0, kernel_evals = 0
   end type

   public :: init_source, init_rule, valid_rule, validate_query, mapped_rule
   public :: quiet_nan, stable_expm1

   interface
      pure function c_expm1(x) bind(C, name="expm1") result(value)
         import c_double
         real(c_double), value :: x
         real(c_double) :: value
      end function
   end interface
contains
   pure function quiet_nan() result(value)
      real(rk) :: value
      value = ieee_value(0.0_rk, ieee_quiet_nan)
   end function

   pure function stable_expm1(x) result(value)
      real(rk), intent(in) :: x
      real(rk) :: value
      value = real(c_expm1(real(x,c_double)),rk)
   end function

   subroutine init_source(source,c,a,b,alpha,status,cx,cy)
      type(patch_source_t), intent(out) :: source
      real(rk), intent(in) :: c(0:2,0:2), a, b, alpha
      integer, intent(out) :: status
      real(rk), intent(in), optional :: cx, cy
      real(rk) :: center_x, center_y, limit
      integer :: ix, iy
      source = patch_source_t()
      status = invalid
      center_x = 0.0_rk
      center_y = 0.0_rk
      if (present(cx)) center_x = cx
      if (present(cy)) center_y = cy
      if (.not.all(ieee_is_finite(c))) return
      if (.not.all(ieee_is_finite([a,b,alpha,center_x,center_y]))) return
      if (a <= 0.0_rk .or. b <= 0.0_rk .or. alpha <= 0.0_rk) return
      ! Keep squared distances representable; this is not a physical scale gate.
      limit = sqrt(huge(1.0_rk))/16.0_rk
      if (max(a,b,abs(center_x),abs(center_y)) > limit) return
      if (min(a,b) < sqrt(tiny(1.0_rk))) return
      source%c = c
      source%a = a
      source%b = b
      source%alpha = alpha
      source%cx = center_x
      source%cy = center_y
      do ix = 0, 2
         do iy = 0, 2
            if (c(iy,ix) /= 0.0_rk) then
               source%mx = max(source%mx,ix)
               source%my = max(source%my,iy)
            end if
         end do
      end do
      source%ready = .true.
      status = ok
   end subroutine

   subroutine validate_query(source,query,x,y,z,dt,status)
      type(patch_source_t), intent(in) :: source
      real(rk), intent(in) :: query(4)
      real(rk), intent(out) :: x, y, z, dt
      integer, intent(out) :: status
      real(rk) :: limit
      status = invalid
      x = quiet_nan(); y = x; z = x; dt = x
      if (.not.source%ready) return
      if (.not.all(ieee_is_finite(query))) return
      if (query(3) <= 0.0_rk .or. query(4) <= 0.0_rk) return
      limit = sqrt(huge(1.0_rk))/16.0_rk
      if (maxval(abs(query(1:3))) > limit) return
      if (query(3) < sqrt(tiny(1.0_rk))) return
      if (source%alpha > 1.0_rk) then
         if (query(4) > huge(1.0_rk)/source%alpha) return
      end if
      x = query(1)-source%cx
      y = query(2)-source%cy
      z = query(3)
      dt = source%alpha*query(4)
      if (.not.ieee_is_finite(dt)) return
      if (dt < tiny(1.0_rk)) return
      status = ok
   end subroutine

   subroutine init_rule(rule,n,status)
      type(gauss_rule_t), intent(out) :: rule
      integer, intent(in) :: n
      integer, intent(out) :: status
      real(rk) :: z, old, p0, p1, p2, derivative
      integer :: i, iteration, half
      status = invalid
      if (n < 1 .or. n > 512) return
      allocate(rule%x(n),rule%w(n))
      rule%n = n
      half = (n+1)/2
      do i = 1, half
         z = cos(pi*(real(i,rk)-0.25_rk)/(real(n,rk)+0.5_rk))
         do iteration = 1, 100
            call legendre(n,z,p1,p0)
            derivative = real(n,rk)*(z*p1-p0)/(z*z-1.0_rk)
            old = z
            z = old-p1/derivative
            if (abs(z-old) <= 2.0_rk*epsilon(z)) exit
         end do
         if (iteration > 100) then
            status = numerical_failure
            return
         end if
         ! Evaluate the derivative at the accepted root, not the previous iterate.
         call legendre(n,z,p1,p0)
         derivative = real(n,rk)*(z*p1-p0)/(z*z-1.0_rk)
         p2 = 2.0_rk/((1.0_rk-z*z)*derivative*derivative)
         rule%x(i) = -z
         rule%x(n+1-i) = z
         rule%w(i) = p2
         rule%w(n+1-i) = p2
      end do
      rule%ready = .true.
      if (.not.valid_rule(rule)) then
         rule%ready = .false.
         status = numerical_failure
         return
      end if
      status = ok
   contains
      pure subroutine legendre(order,z,pn,previous)
         integer, intent(in) :: order
         real(rk), intent(in) :: z
         real(rk), intent(out) :: pn, previous
         real(rk) :: following
         integer :: j
         previous = 1.0_rk
         pn = z
         do j = 2, order
            following = ((2.0_rk*real(j,rk)-1.0_rk)*z*pn &
                        -(real(j,rk)-1.0_rk)*previous)/real(j,rk)
            previous = pn
            pn = following
         end do
      end subroutine
   end subroutine

   pure logical function valid_rule(rule)
      type(gauss_rule_t), intent(in) :: rule
      valid_rule = .false.
      if (.not.rule%ready .or. rule%n < 1) return
      if (.not.allocated(rule%x)) return
      if (.not.allocated(rule%w)) return
      if (size(rule%x) /= rule%n .or. size(rule%w) /= rule%n) return
      if (.not.all(ieee_is_finite(rule%x))) return
      if (.not.all(ieee_is_finite(rule%w))) return
      if (any(rule%w <= 0.0_rk) .or. any(abs(rule%x) >= 1.0_rk)) return
      valid_rule = .true.
   end function

   subroutine mapped_rule(lo,hi,candidates,rule,x,w,status)
      real(rk), intent(in) :: lo, hi, candidates(:)
      type(gauss_rule_t), intent(in) :: rule
      real(rk), allocatable, intent(out) :: x(:), w(:)
      integer, intent(out) :: status
      real(rk) :: cuts(size(candidates)+2), key, mid, half
      integer :: count, i, j, unique, offset
      status = invalid
      if (.not.valid_rule(rule)) return
      if (.not.all(ieee_is_finite([lo,hi]))) return
      if (.not.all(ieee_is_finite(candidates))) return
      if (hi < lo) return
      if (hi == lo) then
         allocate(x(0),w(0))
         status = ok
         return
      end if
      count = 1
      cuts(count) = lo
      do i = 1, size(candidates)
         if (candidates(i) > lo .and. candidates(i) < hi) then
            count = count+1
            cuts(count) = candidates(i)
         end if
      end do
      count = count+1
      cuts(count) = hi
      do i = 2, count
         key = cuts(i)
         j = i-1
         do while (j >= 1)
            if (cuts(j) <= key) exit
            cuts(j+1) = cuts(j)
            j = j-1
         end do
         cuts(j+1) = key
      end do
      unique = 1
      do i = 2, count
         if (cuts(i) /= cuts(unique)) then
            unique = unique+1
            cuts(unique) = cuts(i)
         end if
      end do
      allocate(x((unique-1)*rule%n),w((unique-1)*rule%n))
      offset = 0
      do i = 1, unique-1
         half = (cuts(i+1)-cuts(i))/2.0_rk
         mid = cuts(i)+half
         x(offset+1:offset+rule%n) = mid+half*rule%x
         w(offset+1:offset+rule%n) = half*rule%w
         offset = offset+rule%n
      end do
      status = ok
   end subroutine
end module
