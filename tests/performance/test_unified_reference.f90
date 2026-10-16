program test_unified_reference
   use replay_types
   use replay_case_bank, only: case_spec_t
   use unified_poly_reference
   implicit none
   type(case_spec_t) :: c,changed
   type(reference_t) :: r,loaded
   integer :: st,u
   c%id='synthetic_reference';c%observable=OBS_T
   allocate(c%q(4,1));c%q(:,1)=[.1_rk,.2_rk,.5_rk,.75_rk]
   call profile_source('U',1.0_rk,2.0_rk,c%source,st)
   call require(st==ok,'fixture source')
   call fresh_reference(c,r,st)
   call require(st==ok.and.r%certified,'independent spatial/shared fixture refinement')
   open(newunit=u,status='scratch',action='readwrite')
   call write_reference(u,c,r);rewind(u)
   call read_reference(u,c,loaded,st)
   call require(st==ok.and.loaded%certified,'reference stream roundtrip')
   call require(all(loaded%values==r%values).and.loaded%floor==r%floor,'binary64 roundtrip exact')
   changed=c;changed%q(1,1)=changed%q(1,1)+1e-6_rk;rewind(u)
   call read_reference(u,changed,loaded,st)
   call require(st/=ok,'changed query rejected')
   changed=c;changed%source%c(0,0)=2;rewind(u)
   call read_reference(u,changed,loaded,st)
   call require(st/=ok,'changed source rejected')
   close(u)
   print '(a)','PASS fresh reference finite fixture, exact serialization and changed-input rejection'
contains
   subroutine require(test,label)
      logical,intent(in) :: test
      character(*),intent(in) :: label
      if(.not.test) then
         print *,trim(label);error stop 'reference test failed'
      end if
   end subroutine
end program
