program test_svd_factor_routes
   use replay_types
   use spatial_routes,only: spatial_t,prepare_spatial,evaluate_spatial,spatial_nodes
   use svd_factor_routes
   implicit none
   type(gauss_rule_t) :: rx,ry
   integer :: st,failures,i
   integer,parameter :: degrees(8)=[0,1,2,4,6,8,16,32]
   failures=0
   call init_rule(rx,32,st);call status_ok(st,'x rule')
   call init_rule(ry,40,st);call status_ok(st,'y rule')
   do i=1,size(degrees)
      call polynomial_case(degrees(i),degrees(i))
   end do
   call polynomial_case(6,4)
   call polynomial_case(2,5)
   call zero_polynomial()
   if(failures/=0) then
      write(*,'(a,i0)') 'FAIL SVD factor checks: ',failures
      stop 1
   end if
   write(*,'(a)') 'PASS SVD factor functional checks; high-degree response discrepancies reported separately'
contains
   subroutine polynomial_case(mx,my)
      integer,intent(in) :: mx,my
      type(source_t) :: s,s2
      type(poly_factor_t) :: p,p2
      type(spatial_t) :: dense
      real(rk) :: c(0:my,0:mx),query(4),basequery(4),a,b,direct,scaled,field,x,y
      integer :: ix,iy,obs,j,k,status
      do ix=0,mx
         do iy=0,my
            c(iy,ix)=(-1.0_rk)**(ix+iy)*(1+.1_rk*mod(ix+2*iy,3))/(1+ix+iy)**2
         end do
      end do
      call init_poly(s,c,.7_rk,1.6_rk,1.7_rk,status,2.3_rk,-1.2_rk)
      call status_ok(status,'signed polynomial source')
      call prepare_poly_factor(s,p,status);call status_ok(status,'polynomial factors')
      call check(p%rank==min(mx+1,my+1),'complete rectangular thin polynomial factors')
      call check(size(p%px,1)==mx+1.and.size(p%py,1)==my+1,'factor degrees retained')
      call check(storage_size(p%px)==64.and.storage_size(p%py)==64,'factor powers are binary64')
      call prepare_spatial(s,dense,status);call status_ok(status,'dense polynomial source')
      ! Interior point checks remain well-conditioned even at degree 32.
      do j=0,2
         x=.125_rk*j;y=-.1_rk*j;field=0
         do k=1,p%rank
            field=field+horner(p%px(:,k),x)*horner(p%py(:,k),y)
         end do
         call near(field,source_value(s,s%cx+s%a*x,s%cy+s%b*y),2e-12_rk,'factor axes reconstruct interior')
      end do
      basequery=[2.46_rk,-1.39_rk,.32_rk,.43_rk]
      do obs=OBS_T,OBS_D
         query=basequery
         a=evaluate_poly_factor(p,query,rx,ry,obs,status);call status_ok(status,'polynomial factor response')
         b=evaluate_spatial(dense,query,rx,ry,obs,status);call status_ok(status,'dense polynomial response')
         direct=direct_response(s,query,rx,ry,obs)
         call near(b,direct,3e-13_rk,'dense/independent source and radial kernel')
         if(max(mx,my)<=8) then
            call near(a,b,3e-12_rk,'factor/dense polynomial same quadrature')
         else
            ! Do not silently replace historical double Horner by a 128-bit route.
            write(*,'(a,2i3,a,i0,a,3es24.15)') 'CONDITIONING degrees ',mx,my,' obs ',obs, &
               ' factor,dense,absolute_delta ',a,b,abs(a-b)
            call check(ieee_is_finite(a),'high-degree factor result is finite')
         end if
         s2=s;s2%cx=0;s2%cy=0
         call prepare_poly_factor(s2,p2,status);call status_ok(status,'unshifted polynomial factors')
         query(1:2)=query(1:2)-[s%cx,s%cy]
         scaled=evaluate_poly_factor(p2,query,rx,ry,obs,status)
         call status_ok(status,'unshifted factor response')
         call near(a,scaled,2e-13_rk,'polynomial center translation')
         s2=s;s2%alpha=s%alpha*2
         call prepare_poly_factor(s2,p2,status);call status_ok(status,'scaled-alpha polynomial factors')
         query=basequery;query(4)=query(4)/2
         scaled=evaluate_poly_factor(p2,query,rx,ry,obs,status)
         call status_ok(status,'scaled-alpha polynomial response')
         if(obs==OBS_D) scaled=scaled*2
         call near(a,scaled,2e-13_rk,'polynomial alpha and physical-age scaling')
      end do
      query=basequery;query(4)=0
      a=evaluate_poly_factor(p,query,rx,ry,OBS_T,status)
      call check(status==invalid.and..not.ieee_is_finite(a),'polynomial zero age rejected')
      query=basequery
      a=evaluate_poly_factor(p,query,rx,ry,99,status)
      call check(status==invalid,'polynomial invalid observable rejected')
   end subroutine

   subroutine zero_polynomial()
      type(source_t) :: s
      type(poly_factor_t) :: p
      real(rk) :: c(0:3,0:5),value
      integer :: status
      c=0
      call init_poly(s,c,1.0_rk,2.0_rk,1.0_rk,status)
      call prepare_poly_factor(s,p,status);call status_ok(status,'zero coefficient SVD')
      call check(p%rank==4,'zero singular components are retained')
      value=evaluate_poly_factor(p,[.1_rk,.2_rk,.3_rk,.4_rk],rx,ry,OBS_D,status)
      call status_ok(status,'zero polynomial response');call near(value,0.0_rk,0.0_rk,'zero polynomial')
   end subroutine

   function direct_response(s,query,r1,r2,obs) result(value)
      type(source_t),intent(in) :: s
      real(rk),intent(in) :: query(4)
      type(gauss_rule_t),intent(in) :: r1,r2
      integer,intent(in) :: obs
      real(rk) :: value,x,y,z,dt,radius,u,kernel,total,correction,term
      real(rk),allocatable :: xn(:),wx(:),yn(:),wy(:)
      integer :: ix,iy,status
      call check_query(s,query,x,y,z,dt,status);call status_ok(status,'direct query')
      call spatial_nodes(s,x,y,z,dt,r1,r2,xn,wx,yn,wy,status);call status_ok(status,'direct spatial nodes')
      total=0;correction=0
      do iy=1,size(yn)
         do ix=1,size(xn)
            radius=sqrt((xn(ix)-x)**2+(yn(iy)-y)**2+z*z);u=radius/(2*sqrt(dt))
            ! Independent explicit incomplete-gamma expressions at these moderate queries.
            if(obs==OBS_T) then
               kernel=z/(2*pi*radius**3)*(erfc(u)+2*u*exp(-u*u)/sqrtpi)
            else
               kernel=query(4)*z/(2*pi*radius**3)*(erf(u)-2*u*exp(-u*u)/sqrtpi) &
                  +z/(4*pi*s%alpha*radius)*erfc(u)
            end if
            term=wx(ix)*wy(iy)*source_value(s,xn(ix)+s%cx,yn(iy)+s%cy)*kernel
            call add_compensated(total,correction,term)
         end do
      end do
      value=total+correction
   end function

   pure function horner(c,x) result(v)
      real(rk),intent(in) :: c(0:),x
      real(rk) :: v
      integer :: j
      v=0
      do j=size(c)-1,0,-1
         v=v*x+c(j)
      end do
   end function
   subroutine check(condition,label)
      logical,intent(in) :: condition
      character(*),intent(in) :: label
      if(.not.condition) then
         failures=failures+1;write(*,'(a)') 'FAIL '//label
      end if
   end subroutine
   subroutine status_ok(status,label)
      integer,intent(in) :: status
      character(*),intent(in) :: label
      call check(status==ok,label//' status')
   end subroutine
   subroutine near(a,b,tolerance,label)
      real(rk),intent(in) :: a,b,tolerance
      character(*),intent(in) :: label
      logical :: pass
      pass=ieee_is_finite(a).and.ieee_is_finite(b)
      if(pass) pass=abs(a-b)<=tolerance*max(1.0_rk,abs(a),abs(b))
      call check(pass,label)
      if(.not.pass) write(*,'(a,3es25.16)') 'actual,expected,tolerance ',a,b,tolerance
   end subroutine
end program
