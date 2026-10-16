module replay_case_bank
   use replay_types, only: rk,ok,invalid,OBS_T,OBS_D,source_t,profile_source
   implicit none
   private
   public :: case_spec_t,build_cases,node_list,pair_list

   type case_spec_t
      character(len=48) :: id=''
      character(len=48) :: family=''
      character(len=8) :: profile=''
      type(source_t) :: source
      real(rk),allocatable :: q(:,:)
      integer :: observable=OBS_T
      real(rk) :: Fo=0.0_rk,rho=0.0_rk
      integer :: maxnode=0,bank=0
      logical :: audit=.false.
   end type

   real(rk),parameter :: GEOM_A(4)=[1.0_rk,.2_rk,.6_rk,1.0_rk]
   real(rk),parameter :: GEOM_B(4)=[2.0_rk,.25_rk,.75_rk,.75_rk]
   character(len=8),parameter :: FIVE_PROFILES(5)=[character(len=8) :: 'Q','N1','N2','B1','B2']
   real(rk),parameter :: TEMP_FO(4)=[.001_rk,.1_rk,1.0_rk,1000.0_rk]
   real(rk),parameter :: TEMP_Z(4)=[.1_rk,.5_rk,1.0_rk,2.0_rk]
   real(rk),parameter :: PAPER_FO(2)=[.75_rk,1.0_rk]
   real(rk),parameter :: PAPER_RHO(5)=[.025_rk,.1_rk,.2_rk,.4_rk,.8_rk]
   real(rk),parameter :: EXT_DEPTHS(8)=[.05_rk,.1_rk,.2_rk,.5_rk,1.0_rk,2.0_rk,5.0_rk,10.0_rk]
contains
   subroutine build_cases(family,cases,status)
      character(len=*),intent(in) :: family
      type(case_spec_t),allocatable,intent(out) :: cases(:)
      integer,intent(out) :: status
      real(rk) :: positions(2,9),a,b,L,Fo,rho,z,original_depths(50)
      character(len=8) :: profile
      integer :: i,j,k,m,idx,st
      logical :: audit

      status=invalid
      select case(trim(family))
      case('temperature_scale')
         allocate(cases(128));idx=0
         do i=1,4
            a=GEOM_A(i);b=GEOM_B(i);L=min(a,b)
            do j=1,2
               if(j==1) then
                  profile='U'
               else
                  profile='Q'
               end if
               do k=1,4
                  Fo=TEMP_FO(k)
                  do m=1,4
                     z=TEMP_Z(m);rho=z/L
                     call point_positions(a,b,positions,3)
                     idx=idx+1
                     call fill_case(cases(idx),family,idx,profile,a,b,Fo,rho,z,positions,3, &
                        OBS_T,2,32,.false.,st)
                     if(st/=ok) return
                  end do
               end do
            end do
         end do
      case('paper')
         allocate(cases(100));idx=0
         do i=1,5
            profile=FIVE_PROFILES(i)
            do j=1,2
               call geometry(j,a,b);L=min(a,b)
               call point_positions(a,b,positions,9)
               do k=1,2
                  Fo=PAPER_FO(k)
                  do m=1,5
                     rho=PAPER_RHO(m);z=rho*L;audit=(m==1)
                     idx=idx+1
                     call fill_case(cases(idx),family,idx,profile,a,b,Fo,rho,z,positions,9, &
                        OBS_T,3,32,audit,st)
                     if(st/=ok) return
                  end do
               end do
            end do
         end do
      case('deficit_original')
         allocate(cases(8));idx=0
         a=1.0_rk;b=2.0_rk;L=1.0_rk
         call single_position(0.0_rk,0.0_rk,positions)
         do i=1,50
            original_depths(i)=.1_rk+1.4_rk*real(i-1,rk)/49.0_rk
         end do
         idx=idx+1
         call fill_case(cases(idx),family,idx,'U',a,b,100.0_rk,.1_rk,.1_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=original_depths)
         if(st/=ok) return
         call single_position(.5_rk,.5_rk,positions)
         idx=idx+1
         call fill_case(cases(idx),family,idx,'U',a,b,100.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
         call single_position(.9_rk,1.8_rk,positions)
         idx=idx+1
         call fill_case(cases(idx),family,idx,'U',a,b,100.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
         call single_position(2.0_rk,0.0_rk,positions)
         idx=idx+1
         call fill_case(cases(idx),family,idx,'U',a,b,100.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
         call single_position(0.0_rk,0.0_rk,positions)
         idx=idx+1
         call fill_case(cases(idx),family,idx,'U',a,b,1.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
         idx=idx+1
         call fill_case(cases(idx),family,idx,'U',a,b,1000.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
         idx=idx+1
         call fill_case(cases(idx),family,idx,'Q',a,b,100.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
         idx=idx+1
         call fill_case(cases(idx),family,idx,'Q',a,b,1.0_rk,.05_rk,.05_rk,positions,1, &
            OBS_D,5,32,.false.,st,depths=EXT_DEPTHS)
         if(st/=ok) return
      case default
         return
      end select
      status=ok
   end subroutine

   subroutine fill_case(c,family,index,profile,a,b,Fo,rho,z,positions,npoints,observable,bank,maxnode,audit,status, &
      depths)
      type(case_spec_t),intent(out) :: c
      character(len=*),intent(in) :: family,profile
      integer,intent(in) :: index,npoints,observable,bank,maxnode
      real(rk),intent(in) :: a,b,Fo,rho,z,positions(:,:)
      logical,intent(in) :: audit
      integer,intent(out) :: status
      real(rk),intent(in),optional :: depths(:)
      integer :: i,st,nq
      real(rk) :: L,zi
      character(len=48) :: id
      c=case_spec_t();status=invalid
      write(id,'(A,"_",I4.4)') trim(family),index
      c%id=trim(id);c%family=trim(family);c%profile=trim(profile)
      c%Fo=Fo;c%rho=rho;c%observable=observable;c%bank=bank;c%maxnode=maxnode;c%audit=audit
      call profile_source(profile,a,b,c%source,st)
      if(st/=ok) return
      L=min(a,b)
      nq=npoints
      if(present(depths)) nq=size(depths)
      ! A depth-sweep batch has no single rho. Zero is the metadata sentinel;
      ! all physical positive depths remain in q and reference/value CSVs.
      if(present(depths)) c%rho=0
      allocate(c%q(4,nq))
      do i=1,nq
         zi=z
         if(present(depths)) then
            zi=depths(i)
            c%q(:,i)=[positions(1,1),positions(2,1),zi,Fo*L*L]
         else
            c%q(:,i)=[positions(1,i),positions(2,i),zi,Fo*L*L]
         end if
      end do
      status=ok
   end subroutine

   subroutine point_positions(a,b,positions,npoints)
      real(rk),intent(in) :: a,b
      real(rk),intent(out) :: positions(2,9)
      integer,intent(in) :: npoints
      positions=0.0_rk
      if(npoints==3) then
         positions(:,1)=[0.0_rk,0.0_rk]
         positions(:,2)=[.9_rk*a,.9_rk*b]
         positions(:,3)=[1.2_rk*a,0.0_rk]
      else
         positions(:,1)=[0.0_rk,0.0_rk]
         positions(:,2:5)=reshape([.9_rk*a,.9_rk*b, .9_rk*a,-.9_rk*b, &
            -.9_rk*a,.9_rk*b, -.9_rk*a,-.9_rk*b],[2,4])
         positions(:,6:9)=reshape([1.2_rk*a,0.0_rk, -1.2_rk*a,0.0_rk, &
            0.0_rk,1.2_rk*b, 0.0_rk,-1.2_rk*b],[2,4])
      end if
   end subroutine

   subroutine single_position(x,y,positions)
      real(rk),intent(in) :: x,y
      real(rk),intent(out) :: positions(2,9)
      positions=0.0_rk;positions(:,1)=[x,y]
   end subroutine

   subroutine geometry(index,a,b)
      integer,intent(in) :: index
      real(rk),intent(out) :: a,b
      if(index==1) then
         a=1.0_rk;b=2.0_rk
      else
         a=1.0_rk;b=.75_rk
      end if
   end subroutine

   subroutine node_list(c,nodes)
      type(case_spec_t),intent(in) :: c
      integer,allocatable,intent(out) :: nodes(:)
      integer :: n
      select case(c%bank)
      case(2,3)
         allocate(nodes(31));nodes=[(n,n=2,32)]
      case(4)
         allocate(nodes(47));nodes=[(n,n=2,48)]
      case(5)
         allocate(nodes(5));nodes=[2,4,8,16,32]
      case default
         allocate(nodes(0))
      end select
   end subroutine

   subroutine pair_list(c,nx,ny)
      type(case_spec_t),intent(in) :: c
      integer,allocatable,intent(out) :: nx(:),ny(:)
      integer :: grid(14),n,i,j,k,count,ng
      select case(c%bank)
      case(3)
         n=31;grid(1:10)=[4,6,8,10,12,16,20,24,28,32];ng=10
      case(4)
         n=47;grid=[4,6,8,10,12,16,20,24,28,32,36,40,44,48];ng=14
      case default
         call node_list(c,nx)
         allocate(ny(size(nx)));ny=nx
         return
      end select
      count=n+ng*ng-ng
      allocate(nx(count),ny(count));k=0
      do i=2,n+1
         k=k+1;nx(k)=i;ny(k)=i
      end do
      do i=1,ng
         do j=1,ng
            if(grid(i)==grid(j)) cycle
            k=k+1;nx(k)=grid(i);ny(k)=grid(j)
         end do
      end do
   end subroutine
end module replay_case_bank
