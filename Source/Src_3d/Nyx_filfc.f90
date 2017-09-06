module fc_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for face centered data. 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: q          <=  array to fill
! ::: q_l1,q_l2,q_l3   => index lower bounds of q
! ::: q_h1,q_h2,q_h3   => index upper bounds of q
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of q array
! ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
! ::: dir	=> direction of face centered data (i.e. x-direction)
! ::: 
! ::: -----------------------------------------------------------

      subroutine filfc(q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,domlo,domhi,dx,xlo,bc,dir)

      integer    domlo(3), domhi(3)
      integer    q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
      real(rt)   xlo(3), dx(3)
      real(rt)   q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      integer    bc(3,2)
      integer    dir

      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    i, j, k, n

      nlft = max(0,domlo(1)-q_l1)
      nrgt = max(0,q_h1-(domhi(1)+1))
      nbot = max(0,domlo(2)-q_l2)
      ntop = max(0,q_h2-(domhi(2)+1))
      ndwn = max(0,domlo(3)-q_l3)
      nup  = max(0,q_h3-(domhi(3)+1))

      if (dir.eq.1) then
!
!     ::::: X-face extending in lo-x direction
!
      if (q_l1 .lt. domlo(1)) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP .or. &
             bc(1,1) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-1,j,k) = 2.d0*q(ilo,j,k) - q(ilo+1,j,k)
            end do
	    end do

            do n = 2, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-n,j,k) = q(ilo-1,j,k)
	       end do
   	       end do
	    end do

	 else if (bc(1,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-n,j,k) = q(ilo+n,j,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(1,1)")
	 end if
      end if
!
!     ::::: X-face extending in hi-x direction
!
      if (q_h1 .gt. domhi(1)+1) then
         ihi = domhi(1)+1

	 if (bc(1,2) .eq. FOEXTRAP .or. &
             bc(1,2) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+1,j,k) = 2.d0*q(ihi,j,k) - q(ihi-1,j,k)
            end do
	    end do

            do n = 2, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+n,j,k) = q(ihi+1,j,k)
	       end do
   	       end do
	    end do

	 else if (bc(1,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ihi-n,j,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(1,2)")
	 end if
      end if
!
!     ::::: X-face extending in lo-y direction
!
      if (q_l2 .lt. domlo(2)) then

         jlo = domlo(2)

	 if (bc(2,1) .eq. FOEXTRAP .or. &
             bc(2,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-n,k) = q(i,jlo,k)
	       end do
   	       end do
	    end do

	 else if (bc(2,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jlo+n-1,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(2,1)")
	 end if
      end if

!
!     ::::: X-face extending in hi-y direction
!
      if (q_h2 .gt. domhi(2)) then

         jhi = domhi(2)

	 if (bc(2,2) .eq. FOEXTRAP .or. &
             bc(2,2) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+n,k) = q(i,jhi,k)
	       end do
   	       end do
	    end do

	 else if (bc(2,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jhi-n+1,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(2,2)")
	 end if
      end if

!
!     ::::: X-face extending in lo-z direction
!
      if (q_l3 .lt. domlo(3)) then

         klo = domlo(3)

	 if (bc(3,1) .eq. FOEXTRAP .or. &
             bc(3,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-n) = q(i,j,klo)
	       end do
   	       end do
	    end do

	 else if (bc(3,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,klo+n-1)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(3,1)")
	 end if
      end if

!
!     ::::: X-face extending in hi-z direction
!
      if (q_h3 .gt. domhi(3)) then

         khi = domhi(3)

	 if (bc(3,2) .eq. FOEXTRAP .or. &
             bc(3,2) .eq. HOEXTRAP) then

            do n = 1, nlft
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+n) = q(i,j,khi)
	       end do
   	       end do
	    end do

	 else if (bc(3,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,khi-n+1)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(3,2)")
	 end if
      end if

      else if (dir.eq.2) then
!
!     ::::: Y-face extending in lo-x direction
!
      if (q_l1 .lt. domlo(1)) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP .or. &
             bc(1,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-n,j,k) = q(ilo,j,k)
	       end do
   	       end do
	    end do

	 else if (bc(1,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-n,j,k) = q(ilo+n-1,j,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(1,1)")
	 end if
      end if
!
!     ::::: Y-face extending in hi-x direction
!
      if (q_h1 .gt. domhi(1)) then
         ihi = domhi(1)

	 if (bc(1,2) .eq. FOEXTRAP .or. &
             bc(1,2) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+n,j,k) = q(ihi,j,k)
	       end do
   	       end do
	    end do

	 else if (bc(1,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ihi-n+1,j,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(1,2)")
	 end if
      end if
!
!     ::::: Y-face extending in lo-y direction
!
      if (q_l2 .lt. domlo(2)) then

         jlo = domlo(2)

	 if (bc(2,1) .eq. FOEXTRAP .or. &
             bc(2,1) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-1,k) = 2.d0*q(i,jlo,k) - q(i,jlo+1,k)
            end do
            end do

            do n = 2, nlft
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-n,k) = q(i,jlo-1,k)
	       end do
   	       end do
	    end do

	 else if (bc(2,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jlo+n,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(2,1)")
	 end if
      end if

!
!     ::::: Y-face extending in hi-y direction
!
      if (q_h2 .gt. domhi(2)+1) then

         jhi = domhi(2)+1

	 if (bc(2,2) .eq. FOEXTRAP .or. &
             bc(2,2) .eq. HOEXTRAP) then

            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+1,k) = 2.d0*q(i,jhi,k) - q(i,jhi-1,k)
            end do
            end do

            do n = 2, nlft
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+n,k) = q(i,jhi+1,k)
	       end do
   	       end do
	    end do

	 else if (bc(2,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jhi-n,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(2,2)")
	 end if
      end if

!
!     ::::: Y-face extending in lo-z direction
!
      if (q_l3 .lt. domlo(3)) then

         klo = domlo(3)

	 if (bc(3,1) .eq. FOEXTRAP .or. &
             bc(3,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-n) = q(i,j,klo)
	       end do
   	       end do
	    end do

	 else if (bc(3,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,klo+n-1)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(3,1)")
	 end if
      end if

!
!     ::::: Y-face extending in hi-z direction
!
      if (q_h3 .gt. domhi(3)) then

         khi = domhi(3)

	 if (bc(3,2) .eq. FOEXTRAP .or. &
             bc(3,2) .eq. HOEXTRAP) then

            do n = 1, nlft
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+n) = q(i,j,khi)
	       end do
   	       end do
	    end do

	 else if (bc(3,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,khi-n+1)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(3,2)")
	 end if
      end if

      else if (dir.eq.3) then
!
!     ::::: Z-face extending in lo-x direction
!
      if (q_l1 .lt. domlo(1)) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. FOEXTRAP .or. &
             bc(1,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ilo-n,j,k) = q(ilo,j,k)
	       end do
   	       end do
	    end do

	 else if (bc(1,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ilo-n,j,k) = q(ilo+n-1,j,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(1,1)")
	 end if
      end if
!
!     ::::: Z-face extending in hi-x direction
!
      if (q_h1 .gt. domhi(1)) then
         ihi = domhi(1)

	 if (bc(1,2) .eq. FOEXTRAP .or. &
             bc(1,2) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do j = q_l2,q_h2
                  q(ihi+n,j,k) = q(ihi,j,k)
	       end do
   	       end do
	    end do

	 else if (bc(1,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do j = q_l2,q_h2
               q(ihi+n,j,k) = q(ihi-n+1,j,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(1,2)")
	 end if
      end if
!
!     ::::: Z-face extending in lo-y direction
!
      if (q_l2 .lt. domlo(2)) then

         jlo = domlo(2)

	 if (bc(2,1) .eq. FOEXTRAP .or. &
             bc(2,1) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jlo-n,k) = q(i,jlo,k)
	       end do
   	       end do
	    end do

	 else if (bc(2,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jlo-n,k) = q(i,jlo+n-1,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(2,1)")
	 end if
      end if

!
!     ::::: Z-face extending in hi-y direction
!
      if (q_h2 .gt. domhi(2)+1) then

         jhi = domhi(2)

	 if (bc(2,2) .eq. FOEXTRAP .or. &
             bc(2,2) .eq. HOEXTRAP) then

            do n = 1, nlft
               do k = q_l3,q_h3
               do i = q_l1,q_h1
                  q(i,jhi+n,k) = q(i,jhi,k)
	       end do
   	       end do
	    end do

	 else if (bc(2,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do k = q_l3,q_h3
            do i = q_l1,q_h1
               q(i,jhi+n,k) = q(i,jhi-n+1,k)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(2,2)")
	 end if
      end if

!
!     ::::: Z-face extending in lo-z direction
!
      if (q_l3 .lt. domlo(3)) then

         klo = domlo(3)

	 if (bc(3,1) .eq. FOEXTRAP .or. &
             bc(3,1) .eq. HOEXTRAP) then

            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-1) = 2.d0*q(i,j,klo) - q(i,j,klo+1)
            end do
            end do

            do n = 2, nlft
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,klo-n) = q(i,j,klo-1)
	       end do
   	       end do
	    end do

	 else if (bc(3,1) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,klo-n) = q(i,j,klo+n)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(3,1)")
	 end if
      end if

!
!     ::::: Z-face extending in hi-z direction
!
      if (q_h3 .gt. domhi(3)+1) then

         khi = domhi(3)+1

	 if (bc(3,2) .eq. FOEXTRAP .or. &
             bc(3,2) .eq. HOEXTRAP) then

            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+1) = 2.d0*q(i,j,khi) - q(i,j,khi-1)
            end do
            end do

            do n = 2, nlft
               do j = q_l2,q_h2
               do i = q_l1,q_h1
                  q(i,j,khi+n) = q(i,j,khi+1)
	       end do
   	       end do
	    end do

	 else if (bc(3,2) .eq. REFLECT_EVEN) then

	    do n = 1, nlft
            do j = q_l2,q_h2
            do i = q_l1,q_h1
               q(i,j,khi+n) = q(i,j,khi-n)
            end do
            end do
	    end do
	 else 
	    call bl_abort("Dont know how to fill this bc(3,2)")
	 end if
      end if

      end if

     end subroutine filfc

end module fc_fill_module
