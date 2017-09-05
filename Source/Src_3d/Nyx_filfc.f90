
! ::: -----------------------------------------------------------
! ::: This routine is intended to be a generic fill function
! ::: for face centered data. 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: var        <=  array to fill
! ::: var_l1,var_l2,var_l3   => index lower bounds of var
! ::: var_h1,var_h2,var_h3   => index upper bounds of var
! ::: domlo,hi  => index extent of problem domain
! ::: dx        => cell spacing
! ::: xlo       => physical location of lower left hand
! :::	           corner of q array
! ::: bc	=> array of boundary flags bc(SPACEDIM,lo:hi)
! ::: dir	=> direction of face centered data (i.e. x-direction)
! ::: 
! ::: -----------------------------------------------------------

      subroutine filfc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc,dir)

      implicit none

      integer    domlo(3), domhi(3)
      REAL_T     xlo(3), dx(3)
      REAL_T     var(var_l1:var_h1,var_l2:var_h2,var_l)
      integer    bc(3,2)
      integer    dir

      integer    nlft, nrgt, nbot, ntop, nup, ndwn
      integer    ilo, ihi, jlo, jhi, klo, khi
      integer    is,  ie,  js,  je,  ks,  ke
      integer    i, j, k

      is = max(var_l1,domlo(1))
      ie = min(var_h1,domhi(1))
      js = max(var_l2,domlo(2))
      je = min(var_h2,domhi(2))
      ks = max(var_l3,domlo(3))
      ke = min(var_h3,domhi(3))

      nlft = max(0,domlo(1)-var_l1)
      nrgt = max(0,var_h1-domhi(1))
      nbot = max(0,domlo(2)-var_l2)
      ntop = max(0,var_h2-domhi(2))
      ndwn = max(0,domlo(3)-var_l3)
      nup  = max(0,var_h3-domhi(3))

!
!     ::::: first fill sides
!
      if (nlft .gt. 0) then
         ilo = domlo(1)

	 if (bc(1,1) .eq. EXT_DIR) then
!     set all ghost cell values to a prescribed dirichlet
!     value; in this example, we have chosen 1
!	    do i = 1, nlft
!           do k = ARG_L3(q),ARG_H3(q)
!           do j = ARG_L2(q),ARG_H2(q)
!              q(ilo-i,j,k) = 1.d0
!           end do
!           end do
!	    end do
	 else if (bc(1,1) .eq. FOEXTRAP) then
	    do i = 1, nlft
            do k = ARG_L3(q),ARG_H3(q)
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-i,j,k) = q(ilo,j,k)
            end do
            end do
	    end do
	 else if (bc(1,1) .eq. HOEXTRAP) then
	    do i = 2, nlft
            do k = ARG_L3(q),ARG_H3(q)
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-i,j,k) = q(ilo,j,k)
            end do
            end do
	    end do
            if (ilo+2 .le. ie) then
               do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ilo-1,j,k) = (15*q(ilo,j,k) - 10*q(ilo+1,j,k) + &
		                 3*q(ilo+2,j,k)) * eighth
               end do
               end do
            else  
               do k = ARG_L3(q),ARG_H3(q)
               do j = ARG_L2(q),ARG_H2(q)
                  q(ilo-1,j,k) = half*(3*q(ilo,j,k) - q(ilo+1,j,k))
               end do
               end do
            end if
	 else if (bc(1,1) .eq. REFLECT_EVEN) then
	    do i = 1, nlft
            do k = ARG_L3(q),ARG_H3(q)
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-i,j,k) = q(ilo+i-1,j,k)
            end do
            end do
	    end do
	 else if (bc(1,1) .eq. REFLECT_ODD) then
	    do i = 1, nlft
            do k = ARG_L3(q),ARG_H3(q)
            do j = ARG_L2(q),ARG_H2(q)
               q(ilo-i,j,k) = -q(ilo+i-1,j,k)
            end do
            end do
	    end do
	 end if
      end if


