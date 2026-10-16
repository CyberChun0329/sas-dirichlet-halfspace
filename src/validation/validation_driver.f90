program validation_driver
   use iso_fortran_env, only: output_unit
   use validation_sas
   use mp_moments, only: mp_steady,mp_gaussian_coeff
   use validation_reference, only: ref_temperature,ref_moment,ref_deficit,ref_t90,ref_cheb,ref_operator
   use three_dimensional_analysis, only: steady_integral_semianalytic
   implicit none
   character(32) :: bank,arg
   integer :: chunk,ios
   call get_command_argument(1,bank)
   chunk=1
   call get_command_argument(2,arg)
   if(len_trim(arg)>0) read(arg,*)chunk
   select case(trim(bank))
   case('core');call core_bank()
   case('case1');call case1_bank()
   case('gaussian');call gaussian_bank()
   case('table_s2');call table_s2_bank()
   case('roots');call roots_bank(chunk)
   case default;error stop 'bank must be core/case1/gaussian/table_s2/roots'
   end select
contains
   subroutine open_result(name,unit)
      character(*),intent(in) :: name
      integer,intent(out) :: unit
      open(newunit=unit,file='results/'//name,status='new',action='write',iostat=ios)
      if(ios/=0) error stop 'output exists or cannot open; do not overwrite evidence'
   end subroutine
   subroutine checkpoint(message)
      character(*),intent(in) :: message
      write(output_unit,'(a)')trim(message)
      flush(output_unit)
   end subroutine
   function steady(x,y,z,bounds,cx,cy,shifted) result(v)
      real(rk),intent(in) :: x,y,z,bounds(4),cx(0:),cy(0:)
      logical,optional,intent(in) :: shifted
      real(rk) :: v
      logical :: sh
      sh=.false.;if(present(shifted))sh=shifted
      v=real(mp_steady(real(x,qk),real(y,qk),real(z,qk),real(bounds,qk),real(cx,qk),real(cy,qk),224,sh),rk)
   end function
   subroutine core_bank()
      real(rk),parameter :: zs(9)=[.1_rk,.5_rk,1.0_rk,2.0_rk,5.0_rk,10.0_rk,20.0_rk,50.0_rk,100.0_rk]
      integer,parameter :: mm(12)=[0,1,0,1,2,2,3,4,6,1,6,8],nn(12)=[0,0,1,1,0,3,2,4,1,6,6,8]
      real(rk),parameter :: xx(6)=[0.0_rk,.3_rk,1.0_rk,1.7_rk,0.0_rk,-.9_rk]
      real(rk),parameter :: yy(6)=[0.0_rk,-.6_rk,2.0_rk,-2.4_rk,0.0_rk,.7_rk]
      real(rk),parameter :: zz(6)=[1.0_rk,.5_rk,.7_rk,1.3_rk,.9_rk,.15_rk]
      real(rk) :: bd(4),c(0:12),d(0:12),v,q,q2,q20,v64,e,er,mx,mref,ms,mfar,A,B,Cc,gi,t,formula,R
      real(rk) :: pert,mass,maxop,scaleerror,afferr,pp(0:2),qq(0:2),vscaled,err64,epsb
      integer :: u,us,i,k,m,n,ic,s,ug,ut,uf,uo,j
      call open_result('stability.csv',u)
      write(u,'(a)')'z,m,n,binary64,mp224,GL20,GL200,GL400,error64,error224,reference_change,GL20_change'
      bd=[-1.0_rk,1.0_rk,-1.0_rk,1.0_rk];mx=0;mref=0;mfar=0;err64=0
      do i=1,9
         do k=0,4
            m=2*k;c=0;c(m)=1
            v=steady(0.0_rk,0.0_rk,zs(i),bd,c(:m),c(:m))
            v64=steady_integral_semianalytic(m,m,bd(:2),bd(3:),0.0_rk,0.0_rk,zs(i),s)
            if(s/=0)error stop 'binary64 recurrence status'
            q20=ref_moment(0.0_rk,0.0_rk,zs(i),bd,m,m,20,.false.,.false.)
            q=ref_moment(0.0_rk,0.0_rk,zs(i),bd,m,m,200,.false.,.false.)
            q2=ref_moment(0.0_rk,0.0_rk,zs(i),bd,m,m,400,.false.,.false.)
            e=abs(v-q2)/abs(q2);er=abs(q-q2)/abs(q2)
            mx=max(mx,e);mref=max(mref,er)
            if(zs(i)>=5)mfar=max(mfar,abs(q20-q2)/abs(q2))
            if(m==8.and.i==6)err64=abs(v64-q2)/abs(q2)
            write(u,'(f8.3,2(",",i0),9(",",es25.16))')zs(i),m,m,v64,v,q20,q,q2,abs(v64-q2)/abs(q2),e,er,abs(q20-q2)/abs(q2)
         end do
      end do
      close(u)
      call open_result('core_summary.txt',us)
      write(us,*)'stability_count=',45,' max_mp224_relative=',mx,' max_reference_change=',mref,' max_far_GL20_change=',mfar,' binary64_88_z10=',err64
      if(mx>=1e-12_rk.or.mref>=1e-12_rk.or.mfar>=1e-12_rk.or.err64<=.1_rk)error stop 'stability failed'
      call checkpoint('core: 45 stability comparisons complete')
      call open_result('moment_regression.csv',u)
      write(u,'(a)')'case,x,y,z,m,n,mp224,reference140,scaled_error,reference100_140'
      mx=0;mref=0
      do ic=1,6
         bd=[-1.0_rk,1.0_rk,-2.0_rk,2.0_rk]
         if(ic==1)bd=[-1.0_rk,1.0_rk,-1.0_rk,1.0_rk]
         if(ic==5)bd=[0.0_rk,1.0_rk,0.0_rk,2.0_rk]
         do k=1,12
            m=mm(k);n=nn(k);c=0;d=0;c(m)=1;d(n)=1
            v=steady(xx(ic),yy(ic),zz(ic),bd,c(:m),d(:n),.true.)
            q=ref_moment(xx(ic),yy(ic),zz(ic),bd,m,n,100,.true.,.true.)
            q2=ref_moment(xx(ic),yy(ic),zz(ic),bd,m,n,140,.true.,.true.)
            e=abs(v-q2)/(1+abs(q2));er=abs(q-q2)/(1+abs(q2));mx=max(mx,e);mref=max(mref,er)
            write(u,'(i0,3(",",es25.16),2(",",i0),4(",",es25.16))')ic,xx(ic),yy(ic),zz(ic),m,n,v,q2,e,er
         end do
      end do
      close(u);write(us,*)'moment_count=',72,' max_scaled=',mx,' max_reference_change=',mref
      if(mx>=1e-11_rk.or.mref>=1e-11_rk)error stop 'moment failed'
      call checkpoint('core: 72 general moments complete')
      call open_result('gamma_regression.csv',ug);write(ug,'(a)')'A,B,C,m,n,gamma_integral,spatial,scaled_error'
      mx=0
      do ic=1,3
         select case(ic)
         case(1);A=1.3_rk;B=.7_rk;Cc=.9_rk
         case(2);A=.2_rk;B=1.8_rk;Cc=.3_rk
         case(3);A=2;B=.4_rk;Cc=2
         end select
         bd=[0.0_rk,A,0.0_rk,B]
         do k=1,7
            if(k<=6)then
               m=mm(k);n=nn(k)
               if(k==5)then
                  m=2;n=3
               else if(k==6)then
                  m=3;n=2
               end if
            else
               m=6;n=6
            end if
            gi=sas_gamma_product(A,B,Cc,m,n,1e-13_rk,s)
            if(s/=0)error stop 'gamma status'
            q=2*pi_value()**1.5_rk/Cc*ref_moment(0.0_rk,0.0_rk,Cc,bd,m,n,100,.true.,.false.)
            e=abs(gi-q)/(1+abs(q));mx=max(mx,e)
            write(ug,'(3(es25.16,","),2(i0,","),3(es25.16,:,","))')A,B,Cc,m,n,gi,q,e
         end do
      end do
      close(ug);write(us,*)'gamma_count=',21,' max_scaled=',mx
      if(mx>=1e-11_rk)error stop 'Gamma failed'
      call open_result('signed_transient.csv',ut);write(ut,'(a)')'case,m,n,t,SAS,reference,scaled_error'
      mx=0
      do ic=1,3
         if(ic==3)then
            bd=[0.0_rk,1.0_rk,0.0_rk,2.0_rk];i=5
         else
            bd=[-1.0_rk,1.0_rk,-2.0_rk,2.0_rk];i=ic+1
         end if
         do k=1,4
            select case(k)
            case(1);m=0;n=0
            case(2);m=1;n=1
            case(3);m=2;n=3
            case(4);m=3;n=2
            end select
            c=0;d=0;c(m)=1;d(n)=1
            do j=1,3
               select case(j)
               case(1);t=.03_rk
               case(2);t=.5_rk
               case(3);t=3
               end select
               v=sas_temperature(xx(i),yy(i),zz(i),t,bd,c(:m),d(:n),2e-12_rk,s)
               if(s/=0)error stop 'transient status'
               q=ref_temperature(xx(i),yy(i),zz(i),bd,c(:m),d(:n),t,100,.true.)
               e=abs(v-q)/(1+abs(q));mx=max(mx,e)
               write(ut,'(3(i0,","),4(es25.16,:,","))')ic,m,n,t,v,q,e
            end do
         end do
      end do
      close(ut);write(us,*)'signed_transient_count=',36,' max_scaled=',mx
      if(mx>=1e-11_rk)error stop 'signed transient failed'
      bd=[-.8_rk,1.1_rk,-1.4_rk,.9_rk];scaleerror=0
      do k=1,3
         select case(k)
         case(1);m=0;n=0
         case(2);m=2;n=3
         case(3);m=4;n=4
         end select
         c=0;d=0;c(m)=1;d(n)=1
         v=steady(.2_rk,-.3_rk,.7_rk,bd,c(:m),d(:n))
         vscaled=steady(.6_rk,-.9_rk,2.1_rk,3*bd,c(:m),d(:n))/3.0_rk**(m+n)
         scaleerror=max(scaleerror,abs(vscaled-v))
      end do
      bd=[-.6_rk,1.2_rk,-1.3_rk,.7_rk]
      call affine_coeff([.4_rk,-.7_rk,.2_rk],bd(:2),pp)
      call affine_coeff([.9_rk,.3_rk,-.1_rk],bd(3:),qq)
      v=steady(.1_rk,-.2_rk,.5_rk,bd,pp,qq)
      q=ref_temperature(.1_rk,-.2_rk,.5_rk,bd,[.4_rk,-.7_rk,.2_rk],[.9_rk,.3_rk,-.1_rk],-1.0_rk,100,.true.,normalized=.true.)
      afferr=abs(v-q)
      write(us,*)'scaling_error=',scaleerror,' affine_error=',afferr
      if(scaleerror>=1e-11_rk)error stop 'dimensional scaling failed'
      if(afferr>=1e-12_rk)error stop 'affine coordinate conversion failed'
      call open_result('final_corner_checks.csv',uf);write(uf,'(a)')'A,B,C,m,n,elementary,gamma_integral,scaled_difference'
      mx=0
      do ic=1,3
         select case(ic)
         case(1);A=1.7_rk;B=.45_rk;Cc=.83_rk
         case(2);A=.03_rk;B=2.7_rk;Cc=.19_rk
         case(3);A=.55_rk;B=.82_rk;Cc=3.1_rk
         end select
         R=sqrt(A*A+B*B+Cc*Cc)
         do m=0,1
            do n=0,1
               if(m==0.and.n==0)formula=sqrt(pi_value())/Cc*atan(A*B/(Cc*R))
               if(m==1.and.n==0)formula=sqrt(pi_value())*(asinh(B/Cc)-asinh(B/sqrt(A*A+Cc*Cc)))
               if(m==0.and.n==1)formula=sqrt(pi_value())*(asinh(A/Cc)-asinh(A/sqrt(B*B+Cc*Cc)))
               if(m==1.and.n==1)formula=sqrt(pi_value())*(sqrt(A*A+Cc*Cc)+sqrt(B*B+Cc*Cc)-R-Cc)
               gi=sas_gamma_product(A,B,Cc,m,n,1e-13_rk,s)
               if(s/=0)error stop 'fresh Gamma status'
               e=abs(formula-gi)/(1+abs(gi));mx=max(mx,e)
               write(uf,'(3(es25.16,","),2(i0,","),3(es25.16,:,","))')A,B,Cc,m,n,formula,gi,e
            end do
         end do
      end do
      close(uf);write(us,*)'fresh_corner_count=',12,' max_scaled=',mx
      if(mx>=1e-12_rk)error stop 'fresh corner failed'
      call open_result('final_operator_checks.csv',uo);write(uo,'(a)')'case,t,perturbation,mass,bound'
      bd=[-.7_rk,1.4_rk,-1.1_rk,.9_rk];epsb=.013_rk;maxop=0
      do ic=1,3
         do k=1,4
            select case(k)
            case(1);t=.02_rk
            case(2);t=.4_rk
            case(3);t=3
            case(4);t=-1
            end select
            select case(ic)
            case(1);A=.1_rk;B=-.2_rk;Cc=.12_rk
            case(2);A=1.4_rk;B=.9_rk;Cc=.5_rk
            case(3);A=2.1_rk;B=-.7_rk;Cc=1.3_rk
            end select
            call ref_operator(A,B,Cc,bd,t,120,epsb,pert,mass)
            write(uo,'(i0,4(",",es25.16))')ic,t,pert,mass,epsb*mass
            maxop=max(maxop,abs(pert)-epsb*mass)
            if(abs(pert)>epsb*mass+1e-13_rk.or.mass<0.or.mass>1+1e-13_rk)error stop 'operator bound failed'
         end do
      end do
      close(uo);write(us,*)'operator_count=',12,' max_excess=',maxop
      v=sas_uniform_deficit(0.0_rk,0.0_rk,1000.0_rk,100.0_rk,1.0_rk,[-1.0_rk,1.0_rk,-2.0_rk,2.0_rk],2e-12_rk,s)
      q=800/(2*pi_value()*1e6_rk);e=abs(v/q-1)
      write(us,*)'finite_pulse_farfield_relative=',e
      if(s/=0.or.e>=4e-6_rk)error stop 'farfield failed'
      write(us,*)'status=PASS';close(us);call checkpoint('core: PASS')
   end subroutine
   subroutine affine_coeff(c,bounds,p)
      real(rk),intent(in) :: c(0:2),bounds(2)
      real(rk),intent(out) :: p(0:2)
      real(rk) :: a,b
      a=-(bounds(1)+bounds(2))/(bounds(2)-bounds(1));b=2/(bounds(2)-bounds(1))
      p=[c(0)+c(1)*a+c(2)*a*a,c(1)*b+2*c(2)*a*b,c(2)*b*b]
   end subroutine
   pure function pi_value()result(p)
      real(rk)::p
      p=acos(-1.0_rk)
   end function
   subroutine case1_bank()
      real(rk) :: x,z,v,q,q2,q20,q200,q400,mx,mref,mwhole,bd(4),e,er
      integer :: ix,i,s,u,us,uo,nodes(6),n,ip
      call open_result('case1_profiles.csv',u)
      write(u,'(a)')'x,y,z,SAS,partition80,partition120,GL20,GL200,GL400,relative_error,reference_change'
      bd=[-1.0_rk,1.0_rk,-2.0_rk,2.0_rk];mx=0;mref=0;mwhole=0
      do ix=0,2
         x=real(2*ix,rk)
         do i=0,49
            z=.1_rk+1.4_rk*i/49
            v=sas_uniform_deficit(x,0.0_rk,z,100.0_rk,1.0_rk,bd,2e-12_rk,s)
            if(s/=0)error stop 'Case I status'
            q=ref_deficit(x,0.0_rk,z,100.0_rk,1.0_rk,bd,80,.true.)
            q2=ref_deficit(x,0.0_rk,z,100.0_rk,1.0_rk,bd,120,.true.)
            q20=ref_deficit(x,0.0_rk,z,100.0_rk,1.0_rk,bd,20,.false.)
            q200=ref_deficit(x,0.0_rk,z,100.0_rk,1.0_rk,bd,200,.false.)
            q400=ref_deficit(x,0.0_rk,z,100.0_rk,1.0_rk,bd,400,.false.)
            e=abs(v-q2)/abs(q2);er=abs(q-q2)/abs(q2);mx=max(mx,e);mref=max(mref,er);mwhole=max(mwhole,abs(q200-q400)/abs(q400))
            write(u,'(11(es25.16,:,","))')x,0.0_rk,z,v,q,q2,q20,q200,q400,e,er
         end do
         call checkpoint('case1: one 50-position line complete')
      end do
      close(u);call open_result('case1_summary.txt',us)
      write(us,*)'count=',150,' max_relative=',mx,' max_reference_change=',mref,' whole200_400_change=',mwhole
      if(max(mx,mref)>=1e-12_rk)error stop 'Case I comparison failed'
      call open_result('case1_order_check.csv',uo)
      write(uo,'(a)')'x,y,z,order,SAS,whole_panel_reference,relative_error'
      nodes=[10,20,40,80,200,400]
      do ip=1,2
         x=.5_rk;z=1
         if(ip==2)then
            x=0;z=.1_rk
         end if
         v=sas_uniform_deficit(x,x,z,100.0_rk,1.0_rk,bd,2e-12_rk,s)
         do n=1,6
            q=ref_deficit(x,x,z,100.0_rk,1.0_rk,bd,nodes(n),.false.)
            write(uo,'(3(es25.16,","),i0,3(",",es25.16))')x,x,z,nodes(n),v,q,abs(v-q)/abs(v)
            if(n==5)q200=q
            if(n==6.and.ip==1)write(us,*)'representative_GL200_400_change=',abs(q-q200)/abs(q)
         end do
      end do
      close(uo)
      write(us,*)'status=PASS';close(us)
   end subroutine
   subroutine gaussian_bank()
      integer,parameter :: orders(6)=[2,4,6,8,10,12]
      real(rk) :: coeff(0:12),cr(0:12),cr2(0:12),ray,q,q2,v,pr,tv,tv2,tr,gref,bd(4),means(6),maxeval(6),maxr,ce,mtol
      real(qk) :: cq(0:12)
      integer :: k,i,M,s,u,us,uc,up
      bd=[-1.0_rk,1.0_rk,-1.0_rk,1.0_rk];means=0;maxeval=0;maxr=0
      call open_result('case2_profiles.csv',up)
      write(up,'(a)')'M,r,x,y,Gaussian_reference,reference_change,SAS,same_polynomial_reference,relative_total,relative_evaluation'
      call open_result('case2_convergence.csv',uc)
      write(uc,'(a)')'M,steady_mean,max_steady_evaluation,transient,transient_polyref,transient_Gaussian,transient_relative,transient_evaluation,transient_tolerance_change,coefficient_refinement'
      gref=ref_temperature(.5_rk,.5_rk,1.0_rk,bd,[1.0_rk],[1.0_rk],1.0_rk,120,.true.,.true.)
      mtol=0
      do k=1,6
         M=orders(k)
         call mp_gaussian_coeff(M,224,cq(:M));coeff(:M)=real(cq(:M),rk)
         call ref_cheb(M,512,cr(:M));call ref_cheb(M,1024,cr2(:M));ce=maxval(abs(cr(:M)-cr2(:M)))
         do i=0,99
            ray=1.5_rk*i/99
            q=ref_temperature(ray/sqrt(2.0_rk),ray/sqrt(2.0_rk),1.0_rk,bd,[1.0_rk],[1.0_rk],-1.0_rk,80,.true.,.true.)
            q2=ref_temperature(ray/sqrt(2.0_rk),ray/sqrt(2.0_rk),1.0_rk,bd,[1.0_rk],[1.0_rk],-1.0_rk,120,.true.,.true.)
            maxr=max(maxr,abs(q-q2)/abs(q2))
            v=real(mp_steady(real(ray/sqrt(2.0_rk),qk),real(ray/sqrt(2.0_rk),qk),1.0_qk,real(bd,qk),cq(:M),cq(:M),224),rk)
            pr=ref_temperature(ray/sqrt(2.0_rk),ray/sqrt(2.0_rk),1.0_rk,bd,cr(:M),cr(:M),-1.0_rk,120,.true.,chebyshev=.true.)
            means(k)=means(k)+abs(v-q2)/abs(q2)/100
            maxeval(k)=max(maxeval(k),abs(v-pr)/abs(pr))
            write(up,'(i0,9(",",es25.16))')M,ray,ray/sqrt(2.0_rk),ray/sqrt(2.0_rk),q2,abs(q-q2)/abs(q2),v,pr,abs(v-q2)/abs(q2),abs(v-pr)/abs(pr)
         end do
         tv=sas_temperature(.5_rk,.5_rk,1.0_rk,1.0_rk,bd,coeff(:M),coeff(:M),2e-12_rk,s)
         if(s/=0)error stop 'Gaussian transient status'
         tv2=sas_temperature(.5_rk,.5_rk,1.0_rk,1.0_rk,bd,coeff(:M),coeff(:M),2e-13_rk,s)
         if(s/=0)error stop 'Gaussian refined status'
         tr=ref_temperature(.5_rk,.5_rk,1.0_rk,bd,cr(:M),cr(:M),1.0_rk,120,.true.,chebyshev=.true.)
         mtol=max(mtol,abs(tv-tv2))
         write(uc,'(i0,9(",",es25.16))')M,means(k),maxeval(k),tv,tr,gref,abs(tv-gref)/abs(gref),abs(tv-tr)/abs(tr),abs(tv-tv2),ce
         if(max(maxeval(k),abs(tv-tr)/abs(tr))>=1e-11_rk)error stop 'Gaussian evaluation failed'
         write(output_unit,*)'gaussian: M=',M,' mean=',means(k),' eval=',maxeval(k);flush(output_unit)
      end do
      close(up);close(uc);call open_result('gaussian_summary.txt',us)
      write(us,*)'orders=',6,' ray_positions=',100,' max_reference_change=',maxr,' max_steady_eval=',maxval(maxeval),' max_tolerance_change=',mtol
      write(us,*)'status=PASS';close(us)
   end subroutine
   subroutine cosine_to_powers(M,c,p)
      integer,intent(in) :: M
      real(rk),intent(in) :: c(0:)
      real(rk),intent(out) :: p(0:)
      real(qk) :: basis(0:M,0:M),temp(0:M)
      integer :: k,j
      basis=0;basis(0,0)=1
      if(M>=1)basis(1,1)=1
      do k=2,M
         do j=1,k
            basis(j,k)=2*basis(j-1,k-1)
         end do
         basis(:,k)=basis(:,k)-basis(:,k-2)
      end do
      temp=matmul(basis,real(c(:M),qk));p=0;p(:M)=real(temp,rk)
   end subroutine
   function sas_root(z,a,b,profile,Sout) result(root)
      real(rk),intent(in) :: z,a,b
      character,intent(in) :: profile
      real(rk),optional,intent(out) :: Sout
      real(rk) :: root,S,lo,hi,mid,T,bd(4),cx(0:2),cy(0:2)
      integer :: i,status
      bd=[-a,a,-b,b];cx=0;cy=0;cx(0)=1;cy(0)=1
      if(profile=='Q')then
         cx(2)=-1/(a*a);cy(2)=-1/(b*b)
      end if
      S=steady(0.0_rk,0.0_rk,z,bd,cx,cy)
      if(present(Sout))Sout=S
      lo=1e-12_rk;hi=max(z*z,a*a,b*b)
      do
         T=sas_center_step(hi,z,a,b,profile,2e-12_rk,status)
         if(status/=0)error stop 'root integration status'
         if(T>=.9_rk*S)exit
         hi=2*hi
      end do
      do i=1,200
         mid=sqrt(lo*hi);T=sas_center_step(mid,z,a,b,profile,2e-12_rk,status)
         if(status/=0)error stop 'root integration status'
         if(T<.9_rk*S)then
            lo=mid
         else
            hi=mid
         end if
         if(hi/lo-1<2e-11_rk)exit
      end do
      root=sqrt(lo*hi)
   end function
   subroutine table_s2_bank()
      ! Original six cases, separate from the general 72-moment bank. These
      ! indices are physical-coordinate powers, as in the original Table S2.
      integer,parameter :: mm(6)=[1,0,2,3,4,5],nn(6)=[0,2,1,2,3,5]
      real(rk),parameter :: bd(4)=[.3_rk,1.5_rk,-.8_rk,1.2_rk]
      real(rk) :: c(0:5),d(0:5),v,q20,q64,q96,e,delta,mx,mref
      integer :: k,u,us,m,n
      call open_result('table_s2.csv',u)
      write(u,'(a)')'m,n,mp224,GL20,reference64,reference96,relative_error,reference_change'
      mx=0;mref=0
      do k=1,6
         m=mm(k);n=nn(k);c=0;d=0;c(m)=1;d(n)=1
         v=steady(.6_rk,.2_rk,.9_rk,bd,c(:m),d(:n))
         q20=ref_moment(.6_rk,.2_rk,.9_rk,bd,m,n,20,.false.,.false.)
         q64=ref_moment(.6_rk,.2_rk,.9_rk,bd,m,n,64,.true.,.false.)
         q96=ref_moment(.6_rk,.2_rk,.9_rk,bd,m,n,96,.true.,.false.)
         e=abs(v-q96)/abs(q96);delta=abs(q64-q96)/abs(q96)
         mx=max(mx,e);mref=max(mref,delta)
         write(u,'(i0,",",i0,6(",",es25.16))')m,n,v,q20,q64,q96,e,delta
      end do
      close(u);call open_result('table_s2_summary.txt',us)
      write(us,*)'case_count=',6,' max_relative=',mx,' max_reference_refinement=',mref
      if(mx>=1e-11_rk.or.mref>=1e-11_rk)error stop 'Table S2 independent reference failed'
      write(us,*)'status=PASS';close(us)
   end subroutine
   subroutine roots_bank(part)
      integer,intent(in) :: part
      real(rk) :: z,v,q,q2,maxsas,maxref,a,b,S,r(2,4,3)
      integer :: i,k,j,u,us,ur,computed
      character :: profile
      character(16) :: tag
      if(part<1.or.part>8)error stop 'roots chunk must be 1..8'
      write(tag,'(i0)')part
      call open_result('case3_grid_'//trim(tag)//'.csv',u)
      write(u,'(a)')'node,z,profile,SAS,reference64,reference96,SAS_relative,reference_refinement,ln_z,SAS_t90_over_z2,reference_t90_over_z2'
      maxsas=0;maxref=0;computed=0
      do i=0,400
         z=.1_rk+.00975_rk*i
         if(mod(i,8)/=part-1)cycle
         do k=1,2
            profile='U';if(k==2)profile='Q'
            v=sas_root(z,1.0_rk,.75_rk,profile)
            q=ref_t90(z,1.0_rk,.75_rk,profile,64);q2=ref_t90(z,1.0_rk,.75_rk,profile,96)
            maxsas=max(maxsas,abs(v-q2)/abs(q2));maxref=max(maxref,abs(q-q2)/abs(q2))
            write(u,'(i0,",",es25.16,",",a,8(",",es25.16))')i,z,profile,v,q,q2,abs(v-q2)/abs(q2),abs(q-q2)/abs(q2),log(z),v/z**2,q2/z**2
            computed=computed+1
         end do
         if(mod(i,40)<8)then
            write(output_unit,*)'roots chunk=',part,' input row=',i;flush(output_unit)
         end if
      end do
      close(u);call open_result('roots_'//trim(tag)//'_summary.txt',us)
      write(us,*)'input_nodes=',401,' chunk=',part,' computed_roots=',computed,' max_SAS_relative=',maxsas,' max_reference_refinement=',maxref
      if(computed/=2*((400-(part-1))/8+1).or.maxsas>=1e-9_rk.or.maxref>=1e-11_rk)error stop 'fresh grid failed'
      if(part==1)then
         call open_result('case3_table.csv',ur)
         write(ur,'(a)')'geometry,profile,z,SAS,reference64,reference96,relative_error,reference_change'
         do j=1,3
            select case(j)
            case(1);a=.2_rk;b=.25_rk
            case(2);a=.6_rk;b=.75_rk
            case(3);a=1;b=.75_rk
            end select
            do k=1,2
               profile='U';if(k==2)profile='Q'
               do i=1,4
                  select case(i)
                  case(1);z=.1_rk
                  case(2);z=.5_rk
                  case(3);z=1
                  case(4);z=2
                  end select
                  v=sas_root(z,a,b,profile,S);q=ref_t90(z,a,b,profile,64);q2=ref_t90(z,a,b,profile,96)
                  r(k,i,j)=q2
                  write(ur,'(i0,",",a,6(",",es25.16))')j,profile,z,v,q,q2,abs(v-q2)/abs(q2),abs(q-q2)/abs(q2)
                  if(abs(v-q2)/abs(q2)>=1e-9_rk)error stop 'table root failed'
               end do
            end do
         end do
         close(ur)
         do j=1,3
            write(us,*)'geometry=',j,' shallow_reduction_percent=',100*(1-r(2,1,j)/r(1,1,j)),' deep_reduction_percent=',100*(1-r(2,4,j)/r(1,4,j))
         end do
      end if
      write(us,*)'status=PASS';close(us);call checkpoint('roots chunk '//trim(tag)//': PASS')
   end subroutine
end program
