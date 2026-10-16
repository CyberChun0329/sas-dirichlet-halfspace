program precision_driver
   use mp_arithmetic
   use mp_moments
   implicit none
   integer,parameter :: degrees(2,4)=reshape([0,0,3,2,8,8,12,12],[2,4])
   character(8),parameter :: geometry(7,3)=reshape([character(8) :: &
      '0.17','-0.21','0.43','-0.8','1.2','-1.7','0.9', &
      '0','0','100','-1','1','-1','1', &
      '1','-1','0.08','-1','1','-1','1'],[7,3])
   type(mp_real) :: x,y,z,bd(4),lo(0:12,0:12),hi(0:12,0:12),change,maximum
   type(mp_real) :: a,b,c,d,epsilon,relative,pi_error,elementary_error
   character(1024) :: output_dir
   character(:),allocatable :: max_text,cancel_text,pi_text,elementary_text
   real(qk) :: coefficients(0:12),physical,direct,xx,yy,bdq(4),closed_form
   integer :: g,i,k,m,n,u,summary,imax,gmax,cycles,ios
   logical :: exists
   output_dir='results'
   if(command_argument_count()>0) call get_command_argument(1,output_dir)
   ! Refuse existing evidence before any numerical work. STATUS='new' below
   ! also protects against a file appearing after this preflight check.
   inquire(file=trim(output_dir)//'/precision.csv',exist=exists)
   if(exists) error stop 'precision.csv already exists; fresh output directory required'
   inquire(file=trim(output_dir)//'/precision_summary.txt',exist=exists)
   if(exists) error stop 'precision_summary.txt already exists; fresh output directory required'
   call mp_begin(320)
   if(mp_abi_bytes()/=32) error stop 'Unexpected MPFR LP64 structure size'
   if(mp_version()/='4.2.2') error stop 'Unexpected loaded MPFR version'
   epsilon=mp('1e-70');a=(mp_integer(1)+epsilon)-1
   if(mp_precision_of(a)/=320) error stop 'MPFR precision ABI check failed'
   relative=abs(a-epsilon)/epsilon
   cancel_text=mp_string(relative,40)
   if(mp_compare(relative,mp('1e-24'))>=0) error stop 'MPFR cancellation unit check failed'
   a=mp_pi()
   b=mp('3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798')
   pi_error=abs(a-b);pi_text=mp_string(pi_error,40)
   if(mp_compare(pi_error,mp('1e-94'))>=0) error stop 'MPFR pi unit check failed'
   a=abs(4*atan(mp_integer(1))-mp_pi())
   b=abs(sqrt(mp_integer(2))**2-2)
   c=abs(asinh(mp('0.75'))-log(mp_integer(2)))
   d=abs(exp(log(mp_integer(2)))-2)
   elementary_error=a+b+c+d;elementary_text=mp_string(elementary_error,40)
   if(mp_compare(elementary_error,mp('1e-94'))>=0) error stop 'MPFR elementary unit check failed'
   ! Assignment aliases an immutable value, including across arena growth.
   a=mp('1.23456789012345678901234567890123456789');b=a;c=mp_integer(0)
   do i=1,4096
      c=c+1
   end do
   if(mp_compare(a,b)/=0.or.mp_compare(c,mp_integer(4096))/=0) error stop 'Arena growth/alias check failed'
   do cycles=1,128
      call mp_begin(224)
      a=mp('0.125');b=a*a;c=b+1
      if(mp_precision_of(c)/=224) error stop 'MPFR precision lifecycle check failed'
      if(mp_compare(c,mp('1.015625'))/=0) error stop 'Arena lifecycle value check failed'
      call mp_reset()
      if(mp_live_count()/=0) error stop 'Arena lifecycle reset failed'
   end do
   ! This also verifies the qk wrapper and coefficient API, not high precision
   ! subtraction. The precision acceptance below uses MP objects throughout.
   coefficients=0;coefficients(3)=1
   physical=mp_steady(0.0_qk,0.0_qk,0.5_qk,[-1.0_qk,1.0_qk,-1.0_qk,1.0_qk], &
      coefficients(:3),[1.0_qk],224)
   direct=mp_steady(0.0_qk,0.0_qk,0.5_qk,[-1.0_qk,1.0_qk,-1.0_qk,1.0_qk], &
      coefficients(:3),[1.0_qk],320,.true.)
   if(abs(physical)>1e-60_qk.or.abs(direct)>1e-60_qk) error stop 'Odd symmetric moment check failed'
   physical=mp_steady(0.0_qk,0.0_qk,0.5_qk,[-1.0_qk,1.0_qk,-1.0_qk,1.0_qk], &
      [1.0_qk],[1.0_qk],224)
   closed_form=2/acos(-1.0_qk)*atan(4.0_qk/3)
   if(abs(physical-closed_form)>1e-32_qk) error stop 'Uniform solid angle check failed'
   xx=0.17_qk;yy=-0.21_qk;bdq=[-0.8_qk,1.2_qk,-1.7_qk,0.9_qk]
   physical=mp_steady(xx,yy,0.43_qk,bdq,[0.0_qk,0.0_qk,0.0_qk,1.0_qk],[0.0_qk,0.0_qk,1.0_qk],320)
   direct=mp_steady(xx,yy,0.43_qk,bdq,[xx**3,3*xx**2,3*xx,1.0_qk],[yy**2,2*yy,1.0_qk],320,.true.)
   if(abs(physical-direct)>1e-32_qk) error stop 'Physical/direct binomial shift check failed'
   call mp_gaussian_coeff(12,320,coefficients)
   if(maxval(abs(coefficients(1:12:2)))>0.0_qk) error stop 'Gaussian odd coefficient check failed'
   if(coefficients(0)<=0.999_qk.or.coefficients(0)>=1.0_qk) &
      error stop 'Gaussian coefficient constant bound check failed'
   open(newunit=u,file=trim(output_dir)//'/precision.csv',status='new',action='write',iostat=ios)
   if(ios/=0) error stop 'Cannot create new precision.csv; existing evidence preserved'
   write(u,'(A)') 'geometry,m,n,low_bits,high_bits,low_value,high_value,scaled_change'
   max_text='0';imax=1;gmax=1
   do g=1,3
      call mp_begin(224)
      call load_geometry(g,x,y,z,bd)
      call mp_moment_table(x,y,z,bd,12,12,lo)
      call mp_set_precision(320)
      call load_geometry(g,x,y,z,bd)
      call mp_moment_table(x,y,z,bd,12,12,hi)
      maximum=mp(max_text)
      do k=1,4
         m=degrees(1,k);n=degrees(2,k)
         change=abs(lo(m,n)-hi(m,n))/(1+abs(hi(m,n)))
         if(mp_compare(change,maximum)>0) then
            maximum=change;gmax=g;imax=k
         end if
         write(u,'(I0,",",I0,",",I0,",224,320,",A,",",A,",",A)') &
            g,m,n,mp_string(lo(m,n),75),mp_string(hi(m,n),105),mp_string(change,45)
      end do
      max_text=mp_string(maximum,105)
   end do
   close(u)
   if(mp_compare(maximum,mp('1e-20'))>=0) error stop 'Fresh grid precision refinement failed'
   open(newunit=summary,file=trim(output_dir)//'/precision_summary.txt',status='new',action='write',iostat=ios)
   if(ios/=0) error stop 'Cannot create new precision_summary.txt; existing evidence preserved'
   write(summary,'(A)') 'status: PASS'
   write(summary,'(A)') 'arithmetic: loaded MPFR '//mp_version()//' via Fortran ISO_C_BINDING; RNDN; LP64 structure 32 bytes'
   write(summary,'(A)') 'precision_ABI: mpfr_get_prec confirms actual 224/320 bit values PASS'
   write(summary,'(A)') 'precision: 224 bits (67.43 decimal digits) and 320 bits (96.33 decimal digits)'
   write(summary,'(A)') 'fresh_grid_count: 12; 3 decimal-string geometries x J(0,0), J(3,2), J(8,8), J(12,12)'
   write(summary,'(A)') 'difference: subtract unrounded MPFR values at 320 bits, divide by 1+abs(high)'
   write(summary,'(A)') 'maximum_scaled_change: '//max_text
   write(summary,'(A,I0,A,I0,A,I0)') 'maximum_location: geometry=',gmax,', m=',degrees(1,imax),', n=',degrees(2,imax)
   write(summary,'(A)') 'cancellation_1_plus_1e_minus70_minus_1_relative_error: '//cancel_text
   write(summary,'(A)') 'pi_literal_absolute_error: '//pi_text
   write(summary,'(A)') 'elementary_identity_error_sum: '//elementary_text
   write(summary,'(A)') 'arena: immutable copy across dynamic growth; 128 clear/reinitialize cycles PASS'
   write(summary,'(A)') 'API_smoke: odd symmetry; uniform solid angle; nonzero-coordinate T(3,2) binomial shift; Gaussian M12 parity/constant bound PASS'
   write(summary,'(A)') 'scope: MPFR unit checks and the 12-point precision-refinement grid'
   close(summary)
   print '(A)', 'precision PASS; maximum scaled change = '//max_text
   call mp_reset()
contains
   subroutine load_geometry(g,x,y,z,bd)
      integer,intent(in) :: g
      type(mp_real),intent(out) :: x,y,z,bd(4)
      integer :: j
      x=mp(trim(geometry(1,g)));y=mp(trim(geometry(2,g)));z=mp(trim(geometry(3,g)))
      do j=1,4
         bd(j)=mp(trim(geometry(j+3,g)))
      end do
   end subroutine
end program
