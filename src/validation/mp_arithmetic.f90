module mp_arithmetic
   ! Minimal MPFR 4.2.x binding for the local LP64 ABI (macOS/Homebrew).
   ! mp_real is an immutable arena index. Intrinsic assignment copies the
   ! index, never MPFR ownership. reset is legal only between independent
   ! problems, after every handle from the old problem has been discarded.
   ! This module has one arena and is deliberately not thread safe/reentrant.
   use, intrinsic :: iso_c_binding
   implicit none
   private
   integer, parameter, public :: qk=selected_real_kind(30,300)
   integer(c_int), parameter :: rndn=0_c_int
   type, bind(C) :: mpfr_c
      integer(c_long) :: precision
      integer(c_int) :: sign
      integer(c_long) :: exponent
      type(c_ptr) :: limbs
   end type
   type, public :: mp_real
      private
      integer :: slot=0
   end type
   type(mpfr_c), allocatable, target, save :: arena(:)
   integer, save :: used=0, working_bits=320
   public :: mp_begin,mp_reset,mp_set_precision,mp_live_count,mp_bits,mp_abi_bytes,mp_precision_of,mp_version
   public :: mp,mp_integer,mp_from_qk,mp_pi,mp_string,mp_to_qk,mp_finite,mp_compare
   public :: operator(+),operator(-),operator(*),operator(/),operator(**)
   public :: abs,sqrt,asinh,atan,exp,log
   interface mp
      module procedure mp_decimal
   end interface
   interface operator(+)
      module procedure add_mp,add_mi,add_im
   end interface
   interface operator(-)
      module procedure sub_mp,neg_mp,sub_mi,sub_im
   end interface
   interface operator(*)
      module procedure mul_mp,mul_mi,mul_im
   end interface
   interface operator(/)
      module procedure div_mp,div_mi,div_im
   end interface
   interface operator(**)
      module procedure pow_mi
   end interface
   interface abs
      module procedure abs_mp
   end interface
   interface sqrt
      module procedure sqrt_mp
   end interface
   interface asinh
      module procedure asinh_mp
   end interface
   interface atan
      module procedure atan_mp
   end interface
   interface exp
      module procedure exp_mp
   end interface
   interface log
      module procedure log_mp
   end interface
   interface
      subroutine c_init(a,bits) bind(C,name='mpfr_init2')
         import mpfr_c,c_long
         type(mpfr_c) :: a
         integer(c_long),value :: bits
      end subroutine
      subroutine c_clear(a) bind(C,name='mpfr_clear')
         import mpfr_c
         type(mpfr_c) :: a
      end subroutine
      function c_set_str(a,s,base,rnd) bind(C,name='mpfr_set_str') result(r)
         import mpfr_c,c_char,c_int
         type(mpfr_c) :: a
         character(c_char) :: s(*)
         integer(c_int),value :: base,rnd
         integer(c_int) :: r
      end function
      function c_set_si(a,n,rnd) bind(C,name='mpfr_set_si') result(r)
         import mpfr_c,c_long,c_int
         type(mpfr_c) :: a
         integer(c_long),value :: n
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_get_str(s,e,base,n,a,rnd) bind(C,name='mpfr_get_str') result(p)
         import mpfr_c,c_ptr,c_long,c_int,c_size_t
         type(c_ptr),value :: s
         integer(c_long) :: e
         integer(c_int),value :: base,rnd
         integer(c_size_t),value :: n
         type(mpfr_c),intent(in) :: a
         type(c_ptr) :: p
      end function
      subroutine c_free_str(s) bind(C,name='mpfr_free_str')
         import c_ptr
         type(c_ptr),value :: s
      end subroutine
      function c_add(c,a,b,rnd) bind(C,name='mpfr_add') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a,b
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_sub(c,a,b,rnd) bind(C,name='mpfr_sub') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a,b
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_mul(c,a,b,rnd) bind(C,name='mpfr_mul') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a,b
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_div(c,a,b,rnd) bind(C,name='mpfr_div') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a,b
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_pow(c,a,n,rnd) bind(C,name='mpfr_pow_si') result(r)
         import mpfr_c,c_int,c_long
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_long),value :: n
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_neg(c,a,rnd) bind(C,name='mpfr_neg') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_abs(c,a,rnd) bind(C,name='mpfr_abs') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_sqrt(c,a,rnd) bind(C,name='mpfr_sqrt') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_asinh(c,a,rnd) bind(C,name='mpfr_asinh') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_atan(c,a,rnd) bind(C,name='mpfr_atan') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_exp(c,a,rnd) bind(C,name='mpfr_exp') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_log(c,a,rnd) bind(C,name='mpfr_log') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         type(mpfr_c),intent(in) :: a
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_pi(c,rnd) bind(C,name='mpfr_const_pi') result(r)
         import mpfr_c,c_int
         type(mpfr_c) :: c
         integer(c_int),value :: rnd
         integer(c_int) :: r
      end function
      function c_cmp(a,b) bind(C,name='mpfr_cmp') result(r)
         import mpfr_c,c_int
         type(mpfr_c),intent(in) :: a,b
         integer(c_int) :: r
      end function
      function c_number(a) bind(C,name='mpfr_number_p') result(r)
         import mpfr_c,c_int
         type(mpfr_c),intent(in) :: a
         integer(c_int) :: r
      end function
      function c_get_prec(a) bind(C,name='mpfr_get_prec') result(r)
         import mpfr_c,c_long
         type(mpfr_c),intent(in) :: a
         integer(c_long) :: r
      end function
      function c_get_version() bind(C,name='mpfr_get_version') result(p)
         import c_ptr
         type(c_ptr) :: p
      end function
   end interface
contains
   subroutine mp_begin(bits)
      integer,intent(in) :: bits
      call mp_reset()
      call mp_set_precision(bits)
   end subroutine
   subroutine mp_set_precision(bits)
      integer,intent(in) :: bits
      if(bits<64) error stop 'MPFR precision must be at least 64 bits'
      working_bits=bits
   end subroutine
   subroutine mp_reset()
      integer :: i
      do i=1,used
         call c_clear(arena(i))
      end do
      used=0
   end subroutine
   integer function mp_live_count()
      mp_live_count=used
   end function
   integer function mp_bits()
      mp_bits=working_bits
   end function
   integer function mp_abi_bytes()
      type(mpfr_c) :: a
      mp_abi_bytes=int(c_sizeof(a))
   end function
   integer function mp_precision_of(a)
      type(mp_real),intent(in) :: a
      mp_precision_of=int(c_get_prec(arena(a%slot)))
   end function
   function mp_version() result(s)
      character(:),allocatable :: s
      character(c_char),pointer :: chars(:)
      integer :: i
      call c_f_pointer(c_get_version(),chars,[32])
      s=''
      do i=1,32
         if(chars(i)==c_null_char) return
         s=s//chars(i)
      end do
      error stop 'MPFR version string exceeded buffer'
   end function
   function fresh() result(a)
      type(mp_real) :: a
      type(mpfr_c),allocatable,target :: larger(:)
      if(.not.allocated(arena)) allocate(arena(256))
      if(used==size(arena)) then
         allocate(larger(2*size(arena)))
         ! Transfer the C structures, not MPFR values. No destructor is
         ! associated with the Fortran storage; clear occurs only in reset.
         larger(:used)=arena(:used)
         call move_alloc(larger,arena)
      end if
      used=used+1;a%slot=used
      call c_init(arena(used),int(working_bits,c_long))
   end function
   function mp_decimal(s) result(a)
      character(*),intent(in) :: s
      type(mp_real) :: a
      integer(c_int) :: r
      a=fresh()
      r=c_set_str(arena(a%slot),trim(s)//c_null_char,10_c_int,rndn)
      if(r/=0) error stop 'Invalid MPFR decimal input'
   end function
   function mp_integer(n) result(a)
      integer,intent(in) :: n
      type(mp_real) :: a
      integer(c_int) :: r
      a=fresh();r=c_set_si(arena(a%slot),int(n,c_long),rndn)
   end function
   function mp_from_qk(x) result(a)
      real(qk),intent(in) :: x
      type(mp_real) :: a
      character(64) :: s
      write(s,'(ES48.36E4)') x
      a=mp(trim(adjustl(s)))
   end function
   function mp_string(a,digits) result(s)
      type(mp_real),intent(in) :: a
      integer,optional,intent(in) :: digits
      character(:),allocatable :: s,t
      character(c_char),pointer :: chars(:)
      character(32) :: es
      type(c_ptr) :: p
      integer(c_long) :: e
      integer :: n,i,j
      n=105;if(present(digits)) n=digits
      if(.not.mp_finite(a)) error stop 'Nonfinite MPFR conversion'
      p=c_get_str(c_null_ptr,e,10_c_int,int(n,c_size_t),arena(a%slot),rndn)
      if(.not.c_associated(p)) error stop 'MPFR string allocation failed'
      call c_f_pointer(p,chars,[n+3])
      t=''
      do i=1,n+2
         if(chars(i)==c_null_char) exit
         t=t//chars(i)
      end do
      call c_free_str(p)
      j=1;if(t(1:1)=='-') j=2
      write(es,'(I0)') e-1
      s=t(:j)//'.'//t(j+1:)//'E'//trim(es)
   end function
   real(qk) function mp_to_qk(a) result(x)
      type(mp_real),intent(in) :: a
      character(:),allocatable :: s
      s=mp_string(a,40)
      read(s,*) x
   end function
   logical function mp_finite(a)
      type(mp_real),intent(in) :: a
      mp_finite=c_number(arena(a%slot))/=0
   end function
   integer function mp_compare(a,b)
      type(mp_real),intent(in) :: a,b
      mp_compare=int(c_cmp(arena(a%slot),arena(b%slot)))
   end function
   function mp_pi() result(c)
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_pi(arena(c%slot),rndn)
   end function
   function add_mp(a,b) result(c)
      type(mp_real),intent(in) :: a,b
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_add(arena(c%slot),arena(a%slot),arena(b%slot),rndn)
   end function
   function sub_mp(a,b) result(c)
      type(mp_real),intent(in) :: a,b
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_sub(arena(c%slot),arena(a%slot),arena(b%slot),rndn)
   end function
   function mul_mp(a,b) result(c)
      type(mp_real),intent(in) :: a,b
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_mul(arena(c%slot),arena(a%slot),arena(b%slot),rndn)
   end function
   function div_mp(a,b) result(c)
      type(mp_real),intent(in) :: a,b
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_div(arena(c%slot),arena(a%slot),arena(b%slot),rndn)
   end function
   function pow_mi(a,n) result(c)
      type(mp_real),intent(in) :: a
      integer,intent(in) :: n
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_pow(arena(c%slot),arena(a%slot),int(n,c_long),rndn)
   end function
   function neg_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_neg(arena(c%slot),arena(a%slot),rndn)
   end function
   function abs_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_abs(arena(c%slot),arena(a%slot),rndn)
   end function
   function sqrt_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_sqrt(arena(c%slot),arena(a%slot),rndn)
   end function
   function asinh_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_asinh(arena(c%slot),arena(a%slot),rndn)
   end function
   function atan_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_atan(arena(c%slot),arena(a%slot),rndn)
   end function
   function exp_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_exp(arena(c%slot),arena(a%slot),rndn)
   end function
   function log_mp(a) result(c)
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      integer(c_int) :: r
      c=fresh();r=c_log(arena(c%slot),arena(a%slot),rndn)
   end function
   function add_mi(a,n) result(c)
      type(mp_real),intent(in) :: a
      integer,intent(in) :: n
      type(mp_real) :: c
      c=a+mp_integer(n)
   end function
   function add_im(n,a) result(c)
      integer,intent(in) :: n
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      c=mp_integer(n)+a
   end function
   function sub_mi(a,n) result(c)
      type(mp_real),intent(in) :: a
      integer,intent(in) :: n
      type(mp_real) :: c
      c=a-mp_integer(n)
   end function
   function sub_im(n,a) result(c)
      integer,intent(in) :: n
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      c=mp_integer(n)-a
   end function
   function mul_mi(a,n) result(c)
      type(mp_real),intent(in) :: a
      integer,intent(in) :: n
      type(mp_real) :: c
      c=a*mp_integer(n)
   end function
   function mul_im(n,a) result(c)
      integer,intent(in) :: n
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      c=mp_integer(n)*a
   end function
   function div_mi(a,n) result(c)
      type(mp_real),intent(in) :: a
      integer,intent(in) :: n
      type(mp_real) :: c
      c=a/mp_integer(n)
   end function
   function div_im(n,a) result(c)
      integer,intent(in) :: n
      type(mp_real),intent(in) :: a
      type(mp_real) :: c
      c=mp_integer(n)/a
   end function
end module
