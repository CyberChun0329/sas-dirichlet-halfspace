program test_source_projection
   use replay_types
   use source_projection
   implicit none
   type(source_t) :: gs,p,low
   type(gauss_rule_t) :: space
   real(rk) :: theta,weight,value,expected,cdiag(0:2,0:2)
   integer :: status,m,n,i
   call init_rule(space,96,status);call require(status==ok,'GL96 preparation')

   ! Independent angular quadrature checks the exact Gaussian coefficient formula.
   call init_gaussian(gs,[2.1_rk,1.3_rk],1.0_rk,.7_rk,1.0_rk,status)
   call project_source(gs,12,p,status);call require(status==ok,'Gaussian Bessel projection')
   do m=0,12
      do n=0,12
         value=0
         do i=1,space%n
            theta=pi/2*(space%x(i)+1);weight=space%w(i)/2
            value=value+weight*exp(-(2.1_rk*cos(theta))**2)*cos(m*theta)
         end do
         if(m>0) value=2*value
         expected=0
         do i=1,space%n
            theta=pi/2*(space%x(i)+1);weight=space%w(i)/2
            expected=expected+weight*exp(-(1.3_rk*cos(theta))**2)*cos(n*theta)
         end do
         if(n>0) expected=2*expected
         call near(p%c(n,m),value*expected,3e-14_rk,'Bessel/angular projection')
      end do
   end do

   cdiag=0;cdiag(0,0)=4;cdiag(1,1)=2;cdiag(2,2)=-1
   call init_poly(p,cdiag,1.0_rk,.7_rk,1.0_rk,status)
   call truncate_poly_rank(p,1,low,status);call require(status==ok,'signed diagonal SVD')
   cdiag=0;cdiag(0,0)=4
   call require(maxval(abs(low%c-cdiag))<3e-15_rk,'SVD retains largest singular component')
   call truncate_poly_rank(p,0,low,status)
   call require(status==ok.and.all(low%c==0),'polynomial rank zero')
   print '(a)','PASS Gaussian projection and polynomial rank truncation'
contains
   subroutine require(condition,label)
      logical,intent(in) :: condition
      character(*),intent(in) :: label
      if(.not.condition) then
         print '(a)','FAIL '//label
         error stop 1
      end if
   end subroutine
   subroutine near(got,wanted,tolerance,label)
      real(rk),intent(in) :: got,wanted,tolerance
      character(*),intent(in) :: label
      if(.not.ieee_is_finite(got) .or. abs(got-wanted)>tolerance) then
         print '(a,3es24.15)','FAIL '//label//': ',got,wanted,abs(got-wanted)
         error stop 1
      end if
   end subroutine
end program
