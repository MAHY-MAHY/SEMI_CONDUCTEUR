!>allocation des vecteurs de travail présent dans mod_var


MODULE mod_allocations

  USE mod_var
  !
CONTAINS


  !***********************************************************************!
  SUBROUTINE allocations ()

    !***********************************************************************!
    !
    IMPLICIT NONE
    !
    !--Declaration des arguments

    !--Declaration des variables locales


    ! Affichage des parametres de l'execution

    ! Allocations des matrice et vecteur globaux

    allocate(h(L),h1_demi(L+1))
    allocate(X(L+2),Y(L+1))

!!!!!******* Alloaction pour l'evolution ****!!!!!!
    allocate(mat_a(L,L),vec_b(L),lu_piv(L))
    allocate(VEC_N(L),VEC_P(L),VEC_C(L))
    allocate(VEC_PSI_APP(L))
    allocate(C(L),bD(L))
    allocate(mat_semi_imp_N(L,L),vec_semi_imp_N(L))
    allocate(mat_semi_imp_P(L,L),vec_semi_imp_P(L))

!!!!!!***** Allocaton pour l'equilibre****!!!
    allocate(mat_a_eq(L,L),vec_b_eq(L))
    allocate(VEC_PSI_APP_eq(L))
    allocate(VEC_P_eq(L))
    allocate(VEC_N_eq(L))

    RETURN
  END SUBROUTINE allocations

END MODULE mod_allocations
