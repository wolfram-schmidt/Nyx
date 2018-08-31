   subroutine fort_enforce_consistent_e(lo,hi,state, &
                                        state_l1,state_l2,state_l3,state_h1,state_h2,state_h3) & 
     bind(C,name="fort_enforce_consistent_e")

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT

     implicit none

     integer          :: lo(3), hi(3)
     integer          :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     real(rt) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     ! Local variables
     integer          :: i,j,k
     real(rt) :: u, v, w, rhoInv

     ! 
     ! Make sure to enforce (rho E) = (rho e) + 1/2 rho (u^2 +_ v^2 + w^2)
     !
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              rhoInv = 1.0d0 / state(i,j,k,URHO)

              u = state(i,j,k,UMX) * rhoInv
              v = state(i,j,k,UMY) * rhoInv
              w = state(i,j,k,UMZ) * rhoInv

              state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                     0.5d0 * state(i,j,k,URHO) * (u*u + v*v + w*w)

           end do
        end do
     end do

   end subroutine fort_enforce_consistent_e

   subroutine fort_mhd_enforce_consistent_e(lo,hi, state, state_l1, state_l2, state_l3, state_h1, state_h2, state_h3, &
          bx, bx_l1, bx_l2, bx_l3, bx_h1, bx_h2, bx_h3, &
          by, by_l1, by_l2, by_l3, by_h1, by_h2, by_h3, &
          bz, bz_l1, bz_l2, bz_l3, bz_h1, bz_h2, bz_h3) &
          bind(C,name="fort_mhd_enforce_consistent_e")

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT

     implicit none

     integer          :: lo(3), hi(3)
     integer          :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     integer          :: bx_l1,bx_l2,bx_l3,bx_h1,bx_h2,bx_h3
     integer          :: by_l1,by_l2,by_l3,by_h1,by_h2,by_h3
     integer          :: bz_l1,bz_l2,bz_l3,bz_h1,bz_h2,bz_h3
     real(rt) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
     real(rt) :: bx(bx_l1:bx_h1,bx_l2:bx_h2,bx_l3:bx_h3)
     real(rt) :: by(by_l1:by_h1,by_l2:by_h2,by_l3:by_h3)
     real(rt) :: bz(bz_l1:bz_h1,bz_l2:bz_h2,bz_l3:bz_h3)

     ! Local variables
     integer          :: i,j,k
     real(rt) :: u, v, w, rhoInv, bcx, bcy, bcz

     ! 
     ! Make sure to enforce (rho E) = (rho e) + 1/2 rho (u^2 + v^2 + w^2) + 1/2 (|B|^2)
     !
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rhoInv = 1.0d0 / state(i,j,k,URHO)

              u = state(i,j,k,UMX) * rhoInv
              v = state(i,j,k,UMY) * rhoInv
              w = state(i,j,k,UMZ) * rhoInv
              if(i.lt.hi(1)) then
                 bcx    = 0.5d0 * (bx(i+1,j,k) + bx(i,j,k))
              else
                 bcx = bx(i,j,k)
              endif
              if(j.lt.hi(2)) then
                 bcy    = 0.5d0 * (by(i,j+1,k) + by(i,j,k))
              else
                 bcy    = by(i,j,k)
              endif
              if(k.lt.hi(3)) then
                 bcz    = 0.5d0 * (bz(i,j,k+1) + bz(i,j,k))
              else
                 bcz    = bz(i,j,k)
              endif

              state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
              0.5d0 * state(i,j,k,URHO) * (u*u + v*v + w*w) + &
              0.5d0 *(bcx**2 + bcy**2 + bcz**2)

           end do
        end do
     end do
   end subroutine fort_mhd_enforce_consistent_e

