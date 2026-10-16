module prepared_space
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use temperature_common, only: rk, pi, sqrtpi, ok, invalid, numerical_failure, &
      patch_source_t, gauss_rule_t, query_work_t, validate_query, mapped_rule, &
      valid_rule, quiet_nan, init_source
   implicit none
   private

   type, public :: prepared_space_t
      private
      type(patch_source_t) :: patch_source
      real(rk) :: coeff(0:2,0:2) = 0.0_rk
      logical :: ready = .false.
   end type

   public :: prepare_space, evaluate_space

contains

   subroutine prepare_space(source, prepared, status)
      type(patch_source_t), intent(in) :: source
      type(prepared_space_t), intent(out) :: prepared
      integer, intent(out) :: status
      real(rk) :: x_power(0:2,0:2)

      status = invalid
      if (.not.source%ready) return
      call init_source(prepared%patch_source,source%c,source%a,source%b,source%alpha, &
                       status,source%cx,source%cy)
      if (status /= ok) return

      ! Exact full tensor conversion T0=1,T1=X,T2=2X^2-1.
      ! Avoid a generic SVD or multiplication by structural zeros at degree<=2.
      x_power(:,0) = source%c(:,0)-source%c(:,2)
      x_power(:,1) = source%c(:,1)
      x_power(:,2) = 2.0_rk*source%c(:,2)
      prepared%coeff(0,:) = x_power(0,:)-x_power(2,:)
      prepared%coeff(1,:) = x_power(1,:)
      prepared%coeff(2,:) = 2.0_rk*x_power(2,:)
      if (.not.all(ieee_is_finite(prepared%coeff))) then
         prepared%coeff = 0.0_rk
         status = numerical_failure
         return
      end if
      prepared%ready = .true.
      status = ok
   end subroutine

   function evaluate_space(prepared, query, rule_x, rule_y, status, work) result(value)
      type(prepared_space_t), intent(in) :: prepared
      real(rk), intent(in) :: query(4)
      type(gauss_rule_t), intent(in) :: rule_x, rule_y
      integer, intent(out) :: status
      type(query_work_t), intent(out), optional :: work
      real(rk) :: value
      real(rk) :: x, y, z, dt, sqrt_dt
      real(rk) :: cuts_x(9), cuts_y(9)
      real(rk), allocatable :: x_nodes(:), wx(:), y_nodes(:), wy(:), x_poly(:,:)
      real(rk) :: f, dx, dy, r2, radius, s, tail, kernel
      real(rk) :: term, sum_value, correction, next_sum
      type(query_work_t) :: counts
      integer :: axis_status, ix, iy, degree

      value = quiet_nan()
      status = invalid
      counts = query_work_t()
      if (present(work)) work = counts
      if (.not.prepared%ready) return
      if (.not.valid_rule(rule_x) .or. .not.valid_rule(rule_y)) return
      call validate_query(prepared%patch_source, query, x, y, z, dt, status)
      if (status /= ok) return

      sqrt_dt = sqrt(dt)
      cuts_x = [x-z/4.0_rk, x+z/4.0_rk, x-z, x+z, x-4.0_rk*z, &
                x+4.0_rk*z, x-2.0_rk*sqrt_dt, x+2.0_rk*sqrt_dt, x]
      cuts_y = [y-z/4.0_rk, y+z/4.0_rk, y-z, y+z, y-4.0_rk*z, &
                y+4.0_rk*z, y-2.0_rk*sqrt_dt, y+2.0_rk*sqrt_dt, y]
      call mapped_rule(-prepared%patch_source%a, prepared%patch_source%a, &
                       cuts_x, rule_x, x_nodes, wx, axis_status)
      if (axis_status /= ok) then
         status = axis_status
         return
      end if
      call mapped_rule(-prepared%patch_source%b, prepared%patch_source%b, &
                       cuts_y, rule_y, y_nodes, wy, axis_status)
      if (axis_status /= ok) then
         status = axis_status
         return
      end if

      counts%x_nodes = size(x_nodes)
      counts%y_nodes = size(y_nodes)
      counts%kernel_evals = counts%x_nodes*counts%y_nodes
      allocate(x_poly(0:prepared%patch_source%my,counts%x_nodes))
      do ix = 1, counts%x_nodes
         dx = x_nodes(ix)/prepared%patch_source%a
         do degree = 0, prepared%patch_source%my
            select case (prepared%patch_source%mx)
            case (0)
               x_poly(degree,ix) = prepared%coeff(degree,0)
            case (1)
               x_poly(degree,ix) = prepared%coeff(degree,1)*dx+prepared%coeff(degree,0)
            case (2)
               x_poly(degree,ix) = (prepared%coeff(degree,2)*dx + &
                                    prepared%coeff(degree,1))*dx+prepared%coeff(degree,0)
            case default
               status = invalid
               return
            end select
         end do
      end do

      sum_value = 0.0_rk
      correction = 0.0_rk
      do iy = 1, counts%y_nodes
         dy = y_nodes(iy)/prepared%patch_source%b
         do ix = 1, counts%x_nodes
            select case (prepared%patch_source%my)
            case (0)
               f = x_poly(0,ix)
            case (1)
               f = x_poly(1,ix)*dy+x_poly(0,ix)
            case (2)
               f = (x_poly(2,ix)*dy+x_poly(1,ix))*dy+x_poly(0,ix)
            case default
               status = invalid
               return
            end select

            r2 = (x_nodes(ix)-x)**2 + (y_nodes(iy)-y)**2 + z*z
            if (.not.ieee_is_finite(r2) .or. r2 <= 0.0_rk .or. &
                .not.ieee_is_finite(f)) then
               status = numerical_failure
               if (present(work)) work = counts
               return
            end if
            radius = sqrt(r2)
            s = radius/(2.0_rk*sqrt_dt)
            if (.not.ieee_is_finite(s)) then
               tail = 0.0_rk
            else if (s > 27.0_rk) then
               tail = 0.0_rk
            else
               tail = erfc(s) + 2.0_rk*s/sqrtpi*exp(-s*s)
            end if
            kernel = (z/r2)/radius/(2.0_rk*pi)*tail
            term = wx(ix)*wy(iy)*f*kernel
            if (.not.ieee_is_finite(kernel) .or. .not.ieee_is_finite(term)) then
               status = numerical_failure
               if (present(work)) work = counts
               return
            end if

            ! Neumaier compensated summation for signed source polynomials.
            next_sum = sum_value + term
            if (abs(sum_value) >= abs(term)) then
               correction = correction + ((sum_value-next_sum)+term)
            else
               correction = correction + ((term-next_sum)+sum_value)
            end if
            sum_value = next_sum
         end do
      end do

      value = sum_value + correction
      if (.not.ieee_is_finite(value)) then
         value = quiet_nan()
         status = numerical_failure
      else
         status = ok
      end if
      if (present(work)) work = counts
   end function

end module
