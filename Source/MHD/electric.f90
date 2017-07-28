module electric_field

implicit none

contains

subroutine electric(Q, E, comp) !Use ideal Ohm's Law
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : QVAR, QU,QV, QW, QMAGX, QMAGY, QMAGZ

implicit none

 real(rt), intent(in)	::Q(QVAR)
 real(rt), intent(out) 	::E
 integer, intent(in)    ::comp

		!E = -v X B
	if(comp.eq. 1) then
	E	= -Q(QV)*Q(QMAGZ) + Q(QW)*Q(QMAGY)
	elseif(comp.eq. 2) then
	E	= -Q(QW)*Q(QMAGX) + Q(QU)*Q(QMAGZ)
	elseif(comp.eq. 3) then
	E	= -Q(QU)*Q(QMAGY) + Q(QV)*Q(QMAGX)
	else
	endif
 
end subroutine electric


!=============== Interpolate the Cell Centered/Face Centered Magnetic Field Vars to Get Edge Centered Electric Field Variables ==================

subroutine elec_interp(E, q, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
			flx, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)

use amrex_fort_module, only : rt => amrex_real
use meth_params_module

implicit none
	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	integer, intent(in)		::flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	real(rt), intent(in)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	real(rt), intent(in) 	::flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)

	real(rt), intent(out)	::E(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,3,4) !12 Edges total
	
	real(rt)				::Ecen
	real(rt)				::a ,b ,d1 ,d2 ,dd1 ,dd2 
	real(rt)				::u_face ,v_face ,w_face
	
	integer					::i ,j ,k	
	E = 0.d0
!Interpolate Electric Fields to edges
	do k = flx_l3+1,flx_h3-2
		do j = flx_l2+1,flx_h2-2
			do i = flx_l1+1,flx_h1-2
!-----------------------------------Calculate Edge 1 := i + 1/2, j, k - 1/2 -------------------------------------------
		!Ey 	
		        !x-derivative
		        u_face = flx(i+1,j,k-1,QRHO,1)	
		        call electric(q(i,j,k-1,:),Ecen,2)
		        a = 2.d0*(flx(i+1,j,k-1,QMAGZ,1) - Ecen)
		        call electric(q(i+1,j,k-1,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i+1,j,k-1,QMAGZ,1))
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+2,j,k-1,QMAGZ,1) - Ecen)
                u_face = q(i+1,j,k-1,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j,k,QRHO,3)
                call electric(q(i,j,k-1,:),Ecen,3)
                a = 2.d0*(flx(i,j,k,QMAGY,1) - Ecen)
                call electric(q(i,j,k,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j,k,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i,j,k,QW)
                a = 2.d0*(Ecen - flx(i,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ey
                E(i,j,k,2,1) = 0.25d0*(flx(i+1,j,k-1,QMAGZ,1) + flx(i+1,j,k,QMAGZ,1) &
                                + flx(i,j,k,QMAGX,3) + flx(i+1,j,k,QMAGX,3)) + dd1 + dd2  
                                
!-----------------------------------Calculate Edge 2 := i, j + 1/2, k - 1/2 -------------------------------------------
         !Ex
		        !y-derivative
		        v_face = flx(i,j+1,k-1,QRHO,2)	
		        call electric(q(i,j,k-1,:),Ecen,2)
		        a = 2.d0*(flx(i,j+1,k-1,QMAGZ,2) - Ecen)
		        call electric(q(i,j+1,k-1,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j+1,k-1,QMAGZ,2))
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i,j+2,k-1,QMAGZ,2) - Ecen)
                v_face = q(i,j+1,k-1,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j,k,QRHO,3)
                call electric(q(i,j,k-1,:),Ecen,3)
                a = 2.d0*(flx(i,j,k,QMAGY,3) - Ecen)
                call electric(q(i,j,k,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j,k,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i,j,k,QW)
                a = 2.d0*(Ecen - flx(i,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ex
                E(i,j,k,1,1) = 0.25d0*(flx(i,j+1,k-1,QMAGZ,2) + flx(i,j+1,k,QMAGZ,2) &
                                + flx(i,j,k,QMAGX,3) + flx(i,j+1,k,QMAGX,3)) + dd1 + dd2  
                                

!-----------------------------------Calculate Edge 3 := i - 1/2, j, k - 1/2 -------------------------------------------
 		!Ey 	
		        !x-derivative
		        u_face = flx(i,j,k-1,QRHO,1)	
		        call electric(q(i-1,j,k-1,:),Ecen,2)
		        a = 2.d0*(flx(i,j,k-1,QMAGZ,1) - Ecen)
		        call electric(q(i,j,k-1,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j,k-1,QMAGZ,1))
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+1,j,k-1,QMAGZ,1) - Ecen)
                u_face = q(i,j,k-1,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j,k,QRHO,3)
                call electric(q(i,j,k-1,:),Ecen,3)
                a = 2.d0*(flx(i,j,k,QMAGY,1) - Ecen)
                call electric(q(i,j,k,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j,k,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i-1,j,k,QW)
                a = 2.d0*(Ecen - flx(i,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ey
                E(i,j,k,2,2) = 0.25d0*(flx(i,j,k-1,QMAGZ,1) + flx(i,j,k,QMAGZ,1) &
                                + flx(i-1,j,k,QMAGX,3) + flx(i,j,k,QMAGX,3)) + dd1 + dd2
                                
!-----------------------------------Calculate Edge 4 := i, j - 1/2, k - 1/2 -------------------------------------------
         !Ex
		        !y-derivative
		        v_face = flx(i,j,k-1,QRHO,2)	
		        call electric(q(i,j,k-1,:),Ecen,2)
		        a = 2.d0*(flx(i,j,k-1,QMAGZ,2) - Ecen)
		        call electric(q(i,j,k-1,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j,k-1,QMAGZ,2))
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i,j+1,k-1,QMAGZ,2) - Ecen)
                v_face = q(i,j,k-1,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j-1,k,QRHO,3)
                call electric(q(i,j-1,k-1,:),Ecen,3)
                a = 2.d0*(flx(i,j-1,k,QMAGY,3) - Ecen)
                call electric(q(i,j-1,k,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j-1,k,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i,j-1,k,QW)
                a = 2.d0*(Ecen - flx(i,j-1,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ex
                E(i,j,k,1,2) = 0.25d0*(flx(i,j,k-1,QMAGZ,2) + flx(i,j,k,QMAGZ,2) &
                                + flx(i,j-1,k,QMAGX,3) + flx(i,j,k,QMAGX,3)) + dd1 + dd2   

!-----------------------------------Calculate Edge 5 := i + 1/2, j - 1/2, k -------------------------------------------
         !Ez
		        !y-derivative
		        v_face = flx(i,j,k,QRHO,2)	
		        call electric(q(i,j-1,k,:),Ecen,2)
		        a = 2.d0*(flx(i,j,k,QMAGX,2) - Ecen) !-3/4
		        call electric(q(i,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j,k,QMAGX,3)) !-1/4
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i,j+1,k,QMAGX,2) - Ecen) !1/4
                v_face = q(i,j,k,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !x-derivative
		        u_face = flx(i+1,j-1,k,QRHO,1)	
		        call electric(q(i,j-1,k,:),Ecen,2)
		        a = 2.d0*(flx(i+1,j-1,k,QMAGY,1) - Ecen) ! 1/4
		        call electric(q(i+1,j-1,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i+1,j-1,k,QMAGY,1)) ! 3/4
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+2,j,k,QMAGY,1) - Ecen)
                u_face = q(i+1,j,k,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)  
                !Edge centered Ez
                E(i,j,k,3,1) = 0.25d0*(flx(i+1,j-1,k,QMAGX,2) + flx(i+1,j,k,QMAGX,2) &
                                + flx(i+1,j,k,QMAGY,1) + flx(i,j,k,QMAGY,1)) + dd1 + dd2  

!-----------------------------------Calculate Edge 6 := i + 1/2, j + 1/2, k -------------------------------------------
         !Ez
		        !y-derivative
		        v_face = flx(i,j+1,k,QRHO,2)	
		        call electric(q(i,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i,j+1,k,QMAGX,2) - Ecen) !1/4
		        call electric(q(i,j+1,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j+1,k,QMAGX,3)) !3/4
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i,j+2,k,QMAGX,2) - Ecen) !5/4
                v_face = q(i,j+1,k,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !x-derivative
		        u_face = flx(i+1,j,k,QRHO,1)	
		        call electric(q(i,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i+1,j,k,QMAGY,1) - Ecen) ! 1/4
		        call electric(q(i+1,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i+1,j,k,QMAGY,1)) ! 3/4
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+2,j,k,QMAGY,1) - Ecen)
                u_face = q(i+1,j,k,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)  
                !Edge centered Ez
                E(i,j,k,3,2) = 0.25d0*(flx(i+1,j+1,k,QMAGX,2) + flx(i,j+1,k,QMAGX,2) &
                                + flx(i+1,j+1,k,QMAGY,1) + flx(i+1,j,k,QMAGY,1)) + dd1 + dd2  
!-----------------------------------Calculate Edge 7 := i - 1/2, j + 1/2, k -------------------------------------------
         !Ez
		        !y-derivative
		        v_face = flx(i-1,j+1,k,QRHO,2)	
		        call electric(q(i-1,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i-1,j+1,k,QMAGX,2) - Ecen) !1/4
		        call electric(q(i-1,j+1,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j+1,k,QMAGX,3)) !3/4
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i-1,j+2,k,QMAGX,2) - Ecen) !5/4
                v_face = q(i-1,j+1,k,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !x-derivative
		        u_face = flx(i,j,k,QRHO,1)	
		        call electric(q(i-1,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i,j,k,QMAGY,1) - Ecen) ! -3/4
		        call electric(q(i,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j,k,QMAGY,1)) ! -1/4
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+1,j,k,QMAGY,1) - Ecen)! 1/4
                u_face = q(i,j,k,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)  
                !Edge centered Ez
                E(i,j,k,3,3) = 0.25d0*(flx(i-1,j+1,k,QMAGX,2) + flx(i,j+1,k,QMAGX,2) &
                                + flx(i-1,j+1,k,QMAGY,1) + flx(i-1,j,k,QMAGY,1)) + dd1 + dd2  			
!-----------------------------------Calculate Edge 8 := i - 1/2, j - 1/2, k -------------------------------------------
         !Ez
		        !y-derivative
		        v_face = flx(i-1,j,k,QRHO,2)	
		        call electric(q(i-1,j-1,k,:),Ecen,2)
		        a = 2.d0*(flx(i-1,j,k,QMAGX,2) - Ecen) !-3/4
		        call electric(q(i-1,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i-1,j,k,QMAGX,3)) !-1/4
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i-1,j+1,k,QMAGX,2) - Ecen) !1/4
                v_face = q(i-1,j+1,k,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !x-derivative
		        u_face = flx(i,j-1,k,QRHO,1)	
		        call electric(q(i-1,j,k-1,:),Ecen,2)
		        a = 2.d0*(flx(i,j-1,k,QMAGY,1) - Ecen) ! -3/4
		        call electric(q(i,j-1,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j-1,k,QMAGY,1)) ! -1/4
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+1,j-1,k,QMAGY,1) - Ecen)! 1/4
                u_face = q(i,j-1,k,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)  
                !Edge centered Ez
                E(i,j,k,3,4) = 0.25d0*(flx(i-1,j-1,k,QMAGX,2) + flx(i,j-1,k,QMAGX,2) &
                                + flx(i-1,j-1,k,QMAGY,1) + flx(i-1,j,k,QMAGY,1)) + dd1 + dd2  	

!-----------------------------------Calculate Edge 9 := i, j - 1/2, k + 1/2 -------------------------------------------
         !Ex
		        !y-derivative
		        v_face = flx(i,j,k,QRHO,2)	
		        call electric(q(i,j-1,k,:),Ecen,2)
		        a = 2.d0*(flx(i,j,k,QMAGZ,2) - Ecen) !-3/4
		        call electric(q(i,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j,k,QMAGZ,2)) !-1/4
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i,j+1,k,QMAGZ,2) - Ecen) !1/4
                v_face = q(i,j,k,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j-1,k+1,QRHO,3)
                call electric(q(i,j-1,k,:),Ecen,3)
                a = 2.d0*(flx(i,j-1,k+1,QMAGY,3) - Ecen) !1/4
                call electric(q(i,j-1,k+1,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j-1,k+1,QMAGX,3)) !3/4
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i,j-1,k+1,QW)
                a = 2.d0*(Ecen - flx(i,j-1,k+1,QMAGX,3)) !5/4
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ex
                E(i,j,k,1,3) = 0.25d0*(flx(i,j,k+1,QMAGZ,2) + flx(i,j,k,QMAGZ,2) &
                                + flx(i,j-1,k+1,QMAGX,3) + flx(i,j,k+1,QMAGX,3)) + dd1 + dd2  

!-----------------------------------Calculate Edge 10 := i + 1/2, j, k + 1/2 -------------------------------------------
		!Ey 	
		        !x-derivative
		        u_face = flx(i+1,j,k,QRHO,1)	
		        call electric(q(i,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i+1,j,k,QMAGZ,1) - Ecen)
		        call electric(q(i+1,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i+1,j,k,QMAGZ,1))
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+2,j,k,QMAGZ,1) - Ecen)
                u_face = q(i+1,j,k,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j,k+1,QRHO,3)
                call electric(q(i,j,k,:),Ecen,3)
                a = 2.d0*(flx(i,j,k+1,QMAGY,1) - Ecen) !1/4
                call electric(q(i,j,k+1,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j,k+1,QMAGX,3)) !3/4
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i,j,k+1,QW)
                a = 2.d0*(Ecen - flx(i,j,k+2,QMAGX,3)) !5/4
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ey
                E(i,j,k,2,3) = 0.25d0*(flx(i+1,j,k,QMAGZ,1) + flx(i+1,j,k+1,QMAGZ,1) &
                                + flx(i,j,k+1,QMAGX,3) + flx(i+1,j,k+1,QMAGX,3)) + dd1 + dd2   

!-----------------------------------Calculate Edge 11 := i, j + 1/2, k + 1/2 -------------------------------------------
         !Ex
		        !y-derivative
		        v_face = flx(i,j+1,k,QRHO,2)	
		        call electric(q(i,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i,j+1,k,QMAGZ,2) - Ecen)
		        call electric(q(i,j+1,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j+1,k,QMAGZ,2))
		        if(v_face.gt. 0.d0) then
		            d1 = a
		        elseif(v_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i,j+2,k,QMAGZ,2) - Ecen)
                v_face = q(i,j+1,k,QV)
                if(v_face.gt. 0.d0) then
                    d2 = b
                elseif(v_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i,j,k+1,QRHO,3)
                call electric(q(i,j,k,:),Ecen,3)
                a = 2.d0*(flx(i,j,k+1,QMAGY,3) - Ecen)
                call electric(q(i,j,k+1,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i,j,k+1,QW)
                a = 2.d0*(Ecen - flx(i,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ex
                E(i,j,k,1,4) = 0.25d0*(flx(i,j+1,k+1,QMAGZ,2) + flx(i,j+1,k,QMAGZ,2) &
                                + flx(i,j,k+1,QMAGX,3) + flx(i,j+1,k+1,QMAGX,3)) + dd1 + dd2  

!-----------------------------------Calculate Edge 12 := i - 1/2, j, k + 1/2 -------------------------------------------                                
 		!Ey 	
		        !x-derivative
		        u_face = flx(i,j,k,QRHO,1)	
		        call electric(q(i-1,j,k,:),Ecen,2)
		        a = 2.d0*(flx(i,j,k,QMAGZ,1) - Ecen)
		        call electric(q(i,j,k,:),Ecen,2)
		        b = 2.d0*(Ecen - flx(i,j,k,QMAGZ,1))
		        if(u_face.gt. 0.d0) then
		            d1 = a
		        elseif(u_face.lt. 0.d0) then
		            d1 = b
		        else 
		            d1 = 0.5d0*(a + b)
		        endif
                a = 2.d0*(flx(i+1,j,k,QMAGZ,1) - Ecen)
                u_face = q(i,j,k,QU)
                if(u_face.gt. 0.d0) then
                    d2 = b
                elseif(u_face.lt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd1 = 0.125d0*(d1 - d2)   
                
                !z-derivative
                w_face = flx(i-1,j,k+1,QRHO,3)
                call electric(q(i-1,j,k,:),Ecen,3)
                a = 2.d0*(flx(i-1,j,k+1,QMAGY,1) - Ecen)
                call electric(q(i-1,j,k+1,:),Ecen,3)
                b = 2.d0*(Ecen - flx(i-1,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d1 = a
                elseif(w_face.lt. 0.d0) then
                    d1 = b
                else 
                    d1 = 0.5d0*(a + b)
                endif
                
                w_face = q(i-1,j,k+1,QW)
                a = 2.d0*(Ecen - flx(i-1,j,k+1,QMAGX,3))
                if(w_face.gt. 0.d0) then
                    d2 = b
                elseif(w_face.gt. 0.d0) then
                    d2 = a
                else
                    d2 = 0.5d0*(a + b)
                endif
                !double derivative
                dd2 = 0.125d0*(d1 - d2)
                !Edge centered Ey
                E(i,j,k,2,4) = 0.25d0*(flx(i,j,k+1,QMAGZ,1) + flx(i,j,k,QMAGZ,1) &
                                + flx(i-1,j,k+1,QMAGX,3) + flx(i,j,k+1,QMAGX,3)) + dd1 + dd2
			enddo
		enddo
	enddo
end subroutine elec_interp	


end module electric_field
