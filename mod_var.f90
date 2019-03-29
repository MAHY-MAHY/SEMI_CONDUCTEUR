!>Declaration des variables globales de travail
MODULE mod_var
  !
  USE numerics
  !
  IMPLICIT NONE
  !
  !---Declaration  des  variables
  !
  ! Matrice MAT_A pour explicite
  !@enum MAT_A
  !>Matrice pour le calcul de PSI
  REAL(rp),  DIMENSION(:,:),  ALLOCATABLE  :: MAT_A
  !
  ! Vecteur second membre global VEC_B pour explicite
  !@enum VEC_B
  !> 2nd memebre pour le calcul de PSI
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_B
  !
  ! MATRice A pour schema semi implicite
  !@enum
  !>Matrice pour la calcul de N pour le shcéma semi implicite
  REAL(rp),  DIMENSION(:,:),  ALLOCATABLE  :: MAT_semi_imp_N
  !@enum
  !>Matrice pour la calcul de P pour le shcéma semi implicite
  REAL(rp),  DIMENSION(:,:),  ALLOCATABLE  :: MAT_semi_imp_P
  ! Vecteur second membre global VEC_B pour explicite
  !@enum
  !>Vecteur N pour la calcul de P pour le shcéma semi implicite
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_semi_imp_N
  !@enum
  !>Vecteur N pour la calcul de P pour le shcéma semi implicite
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_semi_imp_P
  
  ! Vecteur apparochée PSI
  !@enum
  !>Vecteur Psi calculée par volumes finies
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_PSI_APP
  ! Vecteur de N de P et C
  !@enum
  !>Vecteur N pour la calcul de n pour le shcéma explicite
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_N
  !@enum
  !>Vecteur N pour la calcul de p pour le shcéma explicite
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_P
  !@enum
  !>Vecteur C correspodant au dopage du systeme
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_C
  !@enum X
  !>VECTEUR de X correspondant au x_i 
  REAL(rp),  DIMENSION(:),  AllOCATABLE  ::  X
  !@enum y
  !>VECTEUR de Y correspondant au x_(i+1/2) 
  real(rp),dimension(:),ALLOCATABLE :: Y  ! point de discrétisation de (0,1)
  !

!!!!!****!!!!!
  ! Matrice MAT_A pour l'equilibre
  !@enum MAT_A_eq
  !>Matrice pour le calcul de psi a l'equilibre
  REAL(rp),  DIMENSION(:,:),  ALLOCATABLE  :: MAT_A_eq
  !
  ! Vecteur second membre global VEC_B pour l'equilibre
  !@enum VEC_B_eq
  !>2nd memebre pour le calcul de psi à l'equilibre
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_B_eq
  !
  ! Vecteur apparochée PSI pour l'equilibre
  !@enum VEC_PSI_APP_eq
  !>Vecteur Psi à l'equilibre
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_PSI_APP_eq
  !@enum VEC_N_eq
  !>Vecteur N à l'equilibre
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_N_eq
  !@enum VEC_p_eq
  !>Vecteur P à l'equilibr
  REAL(rp),  DIMENSION(:),  ALLOCATABLE  ::  VEC_P_eq

  !
  ! Stockage des pivots de la factorisation LU des blocs diagonaux de MAT_A (in place)
  !@enum
  !> Vecteur pibot pour la decomposition LU
  INTEGER,   DIMENSION(:),    ALLOCATABLE  ::  LU_PIV
  !
  !	MAT_A et VEC_B sont decomposes en blocs carres, de dimension NB
  !@enum
  !>Nombre de points de la ddiscretisation
  INTEGER,  SAVE  ::  L
  !@enum
  !>Initialisation pour la résolution par la méthode de Newton
  real(rp)        :: x_0
  !@enum
  !>vecteur h correspondant au pas x_(i+1/2)+1-x_(i-1/2)
  REAL(rp),  DIMENSION(:),    ALLOCATABLE  :: h
  !@enum
  !>vecteur h correspondant au pas x_(i+1)+1-x_i
  REAL(rp),  DIMENSION(:),    ALLOCATABLE  :: h1_demi
  !@enum
  !> Vecteur C de dopage pour l'equilibre = VEC_C 
  REAL(rp),  DIMENSION(:),    ALLOCATABLE  :: C
  !@enum
  !>Vecteur contenant les conditions aux bords de psi
  REAL(rp),  DIMENSION(:),    ALLOCATABLE  :: bD
  !@enum
  !>Pas de temps
  Real(rp)                  ::  DELTA_T
  !@enum
  !> temps final
  Real(rp)                  ::  Tf
  !@enum
  !>Constante correspondant à delta_n fixée à un 1
  real(rp),SAVE             :: delta_n=one
  !@enum
  !>Constante correspondant à delta_p fixée à 1
  real(rp),SAVE             :: delta_p=one
  !@enum
  !>Condition aux bords gauche de psi
  real(rp)                  :: psi_l
  !@enum
  !>Condition aux bords gauche de psi
  real(rp)                  :: psi_r
  !@enum
  !>Condition aux bords gauche de n
  real(rp)                  :: n_l
  !@enum
  !>Condition aux bords gauche de p
  real(rp)                  :: p_l
  !@enum
  !>Condition aux bords droit de n
  real(rp)                  :: n_r
  !@enum
  !>Condition aux bords droit de p
  real(rp)                  :: p_r
  !@enum
  !>Constante correspondant à l'energie pour la sortie en fonction du temps
  REAL(RP)                  :: Energie
  !@enum
  !>Constante correspondant à la dissipation pour la sortie en fonction du temps
  REAL(RP)                  :: Dissipation
  !@enum
  !>Nombre d'étapes en temps
  INTEGER :: M
  !@enum
  !>variable pour le choix du cas test précisé à l'execution
  Integer::choice
  !@enum
  !>variable pour le choix du schéma  précisé à l'execution
  Integer::choice2

END MODULE mod_var
