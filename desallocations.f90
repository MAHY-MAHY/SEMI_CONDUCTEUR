!> desallocations des vecteurs de travail present dans mod_var


MODULE mod_desallocations

  USE mod_var
  !
CONTAINS


  !***********************************************************************!
  SUBROUTINE desallocations ()

    !***********************************************************************!
    !
    IMPLICIT NONE
    !
    !--Declaration des arguments
    !
    !--Declaration des variables locales
    !
    !
!!!
    !
    ! desallocation des matrice et vecteurs globaux

    deallocate(h,h1_demi,mat_a,vec_b,x,lu_piv)
    deallocate(mat_a_eq,vec_b_eq)

    RETURN



  END SUBROUTINE desallocations


END MODULE mod_desallocations
