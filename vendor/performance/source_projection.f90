module source_projection
   use replay_types
   implicit none
   private
   public :: project_source, truncate_poly_rank, truncate_matrix_rank
   ! Standard LP64 LAPACK DGESVD, linked with macOS -framework Accelerate.
   ! This is the mathematical SVD primitive, not a performance-study variant.
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
   subroutine project_source(source,M,polynomial,status)
      type(source_t),intent(in) :: source
      integer,intent(in) :: M
      type(source_t),intent(out) :: polynomial
      integer,intent(out) :: status
      real(rk),allocatable :: c(:,:)
      real(rk) :: ax(0:M),ay(0:M)
      integer :: i,j
      polynomial=source_t();status=invalid
      if(.not.source%ready .or. M<0 .or. M>32) return
      allocate(c(0:M,0:M));c=0
      select case(source%kind)
      case(SRC_POLY)
         do i=0,min(M,ubound(source%c,2))
            do j=0,min(M,ubound(source%c,1))
               c(j,i)=source%c(j,i)
            end do
         end do
      case(SRC_GAUSSIAN)
         call gaussian_coefficients(source%beta(1),M,ax,status)
         if(status/=ok) return
         call gaussian_coefficients(source%beta(2),M,ay,status)
         if(status/=ok) return
         do i=0,M
            c(:,i)=ay*ax(i)
         end do
      case default
         return
      end select
      call init_poly(polynomial,c,source%a,source%b,source%alpha,status,source%cx,source%cy)
   end subroutine

   subroutine gaussian_coefficients(beta,M,c,status)
      real(rk),intent(in) :: beta
      integer,intent(in) :: M
      real(rk),intent(out) :: c(0:M)
      integer,intent(out) :: status
      real(qk) :: t,logterm,term,total
      integer :: n,j
      status=numerical_failure;c=0;t=real(beta,qk)**2/2
      ! Sum positive terms for exp(-t) I_n(t), without subtractive cancellation.
      ! Log terms avoid overflow of the unscaled Bessel function.
      if(t>10000) return
      if(t==0) then
         c(0)=1;status=ok;return
      end if
      do n=0,M/2
         logterm=real(n,qk)*log(t/2)-log_gamma(real(n+1,qk));total=0
         do j=0,20000
            term=exp(logterm-t);total=total+term
            if(j>int(t) .and. term<=epsilon(total)*total) exit
            logterm=logterm+2*log(t/2)-log(real(j+1,qk))-log(real(j+n+1,qk))
         end do
         if(j>20000) return
         c(2*n)=real(total,rk)
         if(n>0) c(2*n)=2*(-1.0_rk)**n*c(2*n)
      end do
      status=ok
   end subroutine

   subroutine truncate_matrix_rank(A,rank,B,status)
      real(rk),intent(in) :: A(0:,0:)
      integer,intent(in) :: rank
      real(rk),allocatable,intent(out) :: B(:,:)
      integer,intent(out) :: status
      real(rk),allocatable :: copy(:,:),u(:,:),vt(:,:),s(:),work(:)
      real(rk) :: query(1)
      integer :: m,n,k,r,info,lwork,j
      status=invalid;m=size(A,1);n=size(A,2);k=min(m,n)
      if(min(m,n)<1 .or. .not.all(ieee_is_finite(A))) return
      if(rank>k) return
      allocate(B(0:m-1,0:n-1));B=0
      if(rank<0 .or. rank==k) then
         B=A;status=ok;return
      end if
      if(rank==0) then
         status=ok;return
      end if
      r=rank;allocate(copy(m,n),u(m,k),vt(k,n),s(k));copy=A
      call dgesvd('S','S',m,n,copy,m,s,u,m,vt,k,query,-1,info)
      if(info/=0 .or. .not.ieee_is_finite(query(1))) then
         status=numerical_failure;return
      end if
      lwork=max(1,ceiling(query(1)));allocate(work(lwork));copy=A
      call dgesvd('S','S',m,n,copy,m,s,u,m,vt,k,work,lwork,info)
      if(info/=0) then
         status=numerical_failure;return
      end if
      do j=1,r
         B=B+spread(u(:,j)*s(j),2,n)*spread(vt(j,:),1,m)
      end do
      status=ok
      if(.not.all(ieee_is_finite(B))) status=numerical_failure
   end subroutine

   subroutine truncate_poly_rank(poly,rank,lowpoly,status)
      type(source_t),intent(in) :: poly
      integer,intent(in) :: rank
      type(source_t),intent(out) :: lowpoly
      integer,intent(out) :: status
      real(rk),allocatable :: c(:,:)
      lowpoly=source_t();status=invalid
      if(.not.poly%ready .or. poly%kind/=SRC_POLY) return
      call truncate_matrix_rank(poly%c,rank,c,status)
      if(status/=ok) return
      call init_poly(lowpoly,c,poly%a,poly%b,poly%alpha,status,poly%cx,poly%cy)
   end subroutine
end module
