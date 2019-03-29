	MODULE mod_affiche
!
	USE mod_var
!
	CONTAINS

!***********************************************************************!
	SUBROUTINE affiche ()
!***********************************************************************!
!
	IMPLICIT NONE
!            
!--Declaration des arguments
!            
!--Declaration des variables locales
	real(rp), dimension(L) :: vec_err
!            
!***********************************************************************!
!	vec_err = vec_uex -vec_uapp
!	print*, NORM2(vec_err) 

	print*, x_0 
	RETURN
	END SUBROUTINE affiche
!                  
	END MODULE mod_affiche
                  
