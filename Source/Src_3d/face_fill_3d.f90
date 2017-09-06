! ::: -----------------------------------------------------------

      subroutine face_fillx(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3, &
                              domlo,domhi,delta,xlo,time,bc) bind(C, name="face_fillx")

      use amrex_fort_module, only : rt => amrex_real
      use fc_fill_module

      implicit none

      integer var_l1,var_l2,var_l3,var_h1,var_h2,var_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      real(rt) delta(3), xlo(3), time
      real(rt) var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)
      integer dir

      dir = 1

      call filfc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc,dir)

      end subroutine face_fillx

! ::: -----------------------------------------------------------

      subroutine face_filly(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3, &
                              domlo,domhi,delta,xlo,time,bc) bind(C, name="face_filly")

      use amrex_fort_module, only : rt => amrex_real
      use fc_fill_module

      implicit none

      integer var_l1,var_l2,var_l3,var_h1,var_h2,var_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      real(rt) delta(3), xlo(3), time
      real(rt) var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)
      integer dir

      dir = 2

      call filfc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc,dir)

      end subroutine face_filly

! ::: -----------------------------------------------------------

      subroutine face_fillz(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3, &
                              domlo,domhi,delta,xlo,time,bc) bind(C, name="face_fillz")

      use amrex_fort_module, only : rt => amrex_real
      use fc_fill_module

      implicit none

      integer var_l1,var_l2,var_l3,var_h1,var_h2,var_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      real(rt) delta(3), xlo(3), time
      real(rt) var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)
      integer dir

      dir = 3

      call filfc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc,dir)

      end subroutine face_fillz
