      subroutine fort_init_e_from_rhoe(state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns, &
                      lo,hi,a_old) bind(C, name="fort_init_e_from_rhoe")

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN
      use  eos_params_module

      implicit none

      integer  ,intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns
      real(rt) ,intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      integer  ,intent(in   ) :: lo(3), hi(3)
      real(rt) ,intent(in   ) :: a_old

      ! Local variables
      integer          :: i,j,k
      !
      ! Compute energy from the EOS
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                  0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

            enddo
         enddo
      enddo

      end subroutine fort_init_e_from_rhoe

      subroutine fort_init_e_from_rhoe_mhd(state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns, &
            bx,bx_l1,bx_l2,bx_l3,bx_h1,bx_h2,bx_h3, &
            by,by_l1,by_l2,by_l3,by_h1,by_h2,by_h3, &
            bz,bz_l1,bz_l2,bz_l3,bz_h1,bz_h2,bz_h3, &
            lo,hi,a_old) bind(C, name="fort_init_e_from_rhoe_mhd")

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN
      use  eos_params_module

      implicit none

      integer  ,intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns
      integer  ,intent(in   ) :: bx_l1, bx_l2, bx_l3, bx_h1, bx_h2, bx_h3
      integer  ,intent(in   ) :: by_l1, by_l2, by_l3, by_h1, by_h2, by_h3
      integer  ,intent(in   ) :: bz_l1, bz_l2, bz_l3, bz_h1, bz_h2, bz_h3
      real(rt) ,intent(in   ) :: bx(bx_l1:bx_h1,bx_l2:bx_h2,bx_l3:bx_h3)
      real(rt) ,intent(in   ) :: by(by_l1:by_h1,by_l2:by_h2,by_l3:by_h3)
      real(rt) ,intent(in   ) :: bz(bz_l1:bz_h1,bz_l2:bz_h2,bz_l3:bz_h3)
      real(rt) ,intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      integer  ,intent(in   ) :: lo(3), hi(3)
      real(rt) ,intent(in   ) :: a_old

      ! Local variables
      integer          :: i,j,k
      real(rt)         :: bcx, bcy, bcz
      !
      ! Compute energy from the EOS
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               bcx = 0.5d0*(bx(i+1,j,k) + bx(i,j,k))
               bcy = 0.5d0*(by(i,j+1,k) + by(i,j,k))
               bcz = 0.5d0*(bz(i,j,k+1) + bz(i,j,k))
               state(i,j,k,UEDEN) = state(i,j,k,UEINT) &
                   + 0.5d0*(state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 &
                   + state(i,j,k,UMZ)**2) / state(i,j,k,URHO) &
                   + 0.5d0*(bcx**2 + bcy**2 + bcz**2)
            enddo
         enddo
      enddo

      end subroutine fort_init_e_from_rhoe_mhd

                                                  
