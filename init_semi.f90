  !>Initialisation de variable et
  !>de fonction utile à la résolution par la méthode semi implicite

MODULE init_semi

  use mod_var

  use numerics

  use mod_initialisation_syst_lin

  IMPLICIT NONE

CONTAINS

  !@fn
  !> Choix du schéma, Scharfetter-Gummel,décentré amont, ou centré.
  !>Le choix se fait à l'execution
  Function func_B(x)
    REAL(RP) ::x
    REAL(RP) ::func_B
    select case(choice2)
    case(0)
       func_B=1-0.5*x
    case(1)
       func_B=1-min(x,0._RP)
    case(2)
       func_B=x/(exp(x)-1)
    case default
       func_B=x/(exp(x)-1)
    END select
    
    return
  END Function func_B

  !@fn
  !> Initalisation de la Matrice pour trouver N par la méthode semi implicite
  SUBROUTINE init_MAT_semi_imp_N()
    INTEGER  :: i
    REAL(RP) ::TMP_PSI1
    REAL(RP) ::TMP_PSI2
    MAT_semi_imp_N=0
    TMP_PSI1=VEC_PSI_APP(1)-VEC_PSI_APP(2)
    TMP_PSI2=VEC_PSI_APP(1)-psi_l
    MAT_semi_imp_N(1,1)=h(1)/DELTA_T+func_B(TMP_PSI1)/h1_demi(2)&
         +func_B(TMP_PSI2)/h1_demi(1)
    TMP_PSI1=VEC_PSI_APP(2)-VEC_PSI_APP(1)
    MAT_semi_imp_N(1,2)=-func_B(TMP_PSI1)/h1_demi(2)
    DO i=2,L-1
       TMP_PSI1=VEC_PSI_APP(i)-VEC_PSI_APP(i+1)
       TMP_PSI2=VEC_PSI_APP(i)-VEC_PSI_APP(i-1)
       MAT_semi_imp_N(i,i)=h(i)/DELTA_T+func_B(TMP_PSI1)/h1_demi(i+1)&
            +func_B(TMP_PSI2)/h1_demi(i)
       TMP_PSI1=VEC_PSI_APP(i+1)-VEC_PSI_APP(i)
       TMP_PSI2=VEC_PSI_APP(i-1)-VEC_PSI_APP(i)
       MAT_semi_imp_N(i,i+1)=-func_B(TMP_PSI1)/h1_demi(i+1)
       MAT_semi_imp_N(i,i-1)=-func_B(TMP_PSI2)/h1_demi(i)
    END DO
    TMP_PSI1=VEC_PSI_APP(L)-psi_r
    TMP_PSI2=VEC_PSI_APP(L)-VEC_PSI_APP(L-1)
    MAT_semi_imp_N(L,L)=h(L)/DELTA_T+func_B(TMP_PSI1)/h1_demi(L+1)&
         +func_B(TMP_PSI2)/h1_demi(L)
    TMP_PSI2=VEC_PSI_APP(L-1)-VEC_PSI_APP(L)
    MAT_semi_imp_N(L,L-1)=-func_B(TMP_PSI2)/h1_demi(L)
  END SUBROUTINE init_MAT_semi_imp_N

  !@fn
  !> Initalisation de la Matrice pour trouver P par la méthode semi implicite
    SUBROUTINE init_MAT_semi_imp_P()
    INTEGER  :: i
    REAL(RP) ::TMP_PSI1
    REAL(RP) ::TMP_PSI2
    MAT_semi_imp_P=0
    TMP_PSI1=VEC_PSI_APP(2)-VEC_PSI_APP(1)
    TMP_PSI2=psi_l-VEC_PSI_APP(1)
    MAT_semi_imp_P(1,1)=h(1)/DELTA_T+func_B(TMP_PSI1)/h1_demi(2)&
         +func_B(TMP_PSI2)/h1_demi(1)
    TMP_PSI1=VEC_PSI_APP(1)-VEC_PSI_APP(2)
    MAT_semi_imp_P(1,2)=-func_B(TMP_PSI1)/h1_demi(2)
    DO i=2,L-1
       TMP_PSI1=VEC_PSI_APP(i+1)-VEC_PSI_APP(i)
       TMP_PSI2=VEC_PSI_APP(i-1)-VEC_PSI_APP(i)
       MAT_semi_imp_P(i,i)=h(i)/DELTA_T+func_B(TMP_PSI1)/h1_demi(i+1)&
            +func_B(TMP_PSI2)/h1_demi(i)
       TMP_PSI1=VEC_PSI_APP(i)-VEC_PSI_APP(i+1)
       TMP_PSI2=VEC_PSI_APP(i)-VEC_PSI_APP(i-1)
       MAT_semi_imp_P(i,i+1)=-func_B(TMP_PSI1)/h1_demi(i+1)
       MAT_semi_imp_P(i,i-1)=-func_B(TMP_PSI2)/h1_demi(i)
    END DO
    TMP_PSI1=psi_r-VEC_PSI_APP(L)
    TMP_PSI2=VEC_PSI_APP(L-1)-VEC_PSI_APP(L)
    MAT_semi_imp_P(L,L)=h(L)/DELTA_T+func_B(TMP_PSI1)/h1_demi(L+1)&
         +func_B(TMP_PSI2)/h1_demi(L)
    TMP_PSI2=VEC_PSI_APP(L)-VEC_PSI_APP(L-1)
    MAT_semi_imp_P(L,L-1)=-func_B(TMP_PSI2)/h1_demi(L)
  END SUBROUTINE init_MAT_semi_imp_P


  !@fn
  !>Initalisation des vecteurs 2nd membres pour trouver P et N
  !>par la méthode semi implicite.
  SUBROUTINE init_2nd_semi_imp()
    INTEGER  :: i
    REAL(RP) :: TMP_PSI
    VEC_semi_imp_N=0
    VEC_semi_imp_P=0
    DO i=1,L
       VEC_semi_imp_N(i)=h(i)*VEC_N(i)/DELTA_T
       VEC_semi_imp_P(i)=h(i)*VEC_P(i)/DELTA_T
    END DO
!!!!!***** CALCUL pour N******!!!!!!
    TMP_PSI=psi_l-VEC_PSI_APP(1)
    VEC_semi_imp_N(1)=VEC_semi_imp_N(1)+n_l*func_B(TMP_PSI)/h1_demi(1)
    TMP_PSI=psi_r-VEC_PSI_APP(L)
    VEC_semi_imp_N(L)=VEC_semi_imp_N(L)+n_r*func_B(TMP_PSI)/h1_demi(L+1)
!!!!!!!******CALCUL pour P******!!!!!!!!
    TMP_PSI=VEC_PSI_APP(1)-psi_l
    VEC_semi_imp_P(1)=VEC_semi_imp_P(1)+p_l*func_B(TMP_PSI)/h1_demi(1)
    TMP_PSI=VEC_PSI_APP(L)-psi_r
    VEC_semi_imp_P(L)=VEC_semi_imp_P(L)+p_r*func_B(TMP_PSI)/h1_demi(L+1)
   
  END SUBROUTINE init_2nd_semi_imp
END MODULE init_semi

  
