!>Initialisation des variables et fonctions utiles à la résolution

MODULE mod_initialisation_syst_lin
  !
  USE numerics
  !
  USE mod_var

  USE mod_tri

  implicit none
  !
CONTAINS
!!!!!!****** INITIALISATION de vecteur h h1_demi  x et x1demi
  !@fn
  !>Initialisation des vecteurs d'espace et des pas d'espace
  SUBROUTINE init_h()
    INTEGER ::i

    do i = 1, L+1
       y(i) = (i-1)*1._RP/L
    end do
    do i = 1, L
       h(i) = y(i+1) - y(i)
    end do
    
    x(1)   = zero
    x(L+2) = one
    do i= 2, L+1
       x(i) = half*(y(i) + y(i-1))
    end do
    ! vecteur h1_demi
    do i=1, L+1
       h1_demi(i) = x(i+1)-x(i)
    end do
  END SUBROUTINE init_h

!!!!**** Initialisation de la matrice et du vecteur pour Apsi=b pour l'evolution ****!!!!!

  !@fn
  !>initialisation de la matrice utile à la résolution de A*psi=B
  SUBROUTINE Init_A()

    IMPLICIT NONE

    !--Declaration des variables locales
    INTEGER  ::  i,j 

    real(rp)                 :: err = zero
    !*********************************************************    
    !Initialisation de la matrice A = mat_a
    mat_a = zero
    DO i=2,L-1
       mat_a(i,i)   =-one/h1_demi(i)-one/h1_demi(i+1)
       mat_a(i,i-1) = one/h1_demi(i)
       mat_a(i,i+1) = one/h1_demi(i+1)	
    END DO

    mat_a(1,1)   = -one/h1_demi(1)-one/h1_demi(2)	
    mat_a(1,2)   =  one/h1_demi(2)
    mat_a(L,L)   = -one/h1_demi(L) -one/h1_demi(L+1)
    mat_a(L,L-1) = one/h1_demi(L)
  END SUBROUTINE Init_A

  !@fn
  !>Calcul du 2nd memebre B utilse pour la résolution de psi
  Subroutine Init_B

    INTEGER :: i
    ! Initialisation du vecteur b=vec_b
    DO i=1, L
       vec_b(i) =  h(i)*(VEC_N(i)-VEC_P(i)-VEC_C(i))
    END DO
    VEC_B(1) = VEC_B(1)-psi_l/h1_demi(1)
    VEC_B(L) = VEC_B(L)-psi_r/h1_demi(L+1)
  END Subroutine Init_B
!!!!!!***********!!!!!!!!!!!!!!!!

!!!!!!!!!!***** Conditions initiales*********!!!!!!!!!
  !@fn
  !>Choix des conditions initiales dependant du choix 'choice'
  Subroutine init_n0_p0()
    implicit none
    INTEGER ::i
    if(choice<2)THEN
       psi_l=zero
       psi_r=zero
       n_l=one
       p_l=one
       n_r=one
       p_r=one
       DO i=1,L
          VEC_N(i)=0.5
          VEC_P(i)=1
          if(Y(i)<0.5)then
             VEC_C(i)=20
          ELSE
             VEC_C(i)=-20
          END if

          !VEC_C(i)=5
       END DO
    ELSE
       psi_l=-one
       psi_r=one
       n_l=exp(-one)
       p_l=exp(one)
       n_r=exp(one)
       p_r=exp(-one)
       DO i=1,L

          if(Y(i)<0.7)then
             VEC_C(i)=10
             VEC_P(i)=1
          ELSE
             VEC_C(i)=-10
             VEC_P(i)=0.1
          END if
          
          if(Y(i)<0.8)then
             VEC_N(i)=0.1
          ELSE
             VEC_N(i)=1
          END if

          !VEC_C(i)=5
       END DO
    END if
    
  END Subroutine init_n0_p0
!!!!!!!!!********!!!!!!!!!!!!!!!!!

!!!!!!***** Allocation pour l'equilibre *****!!!!!

  !@fn
  !>initialisation de la matrice A pour la résolution de psi à 'equilibre
  SUBROUTINE Init_A_eq()

    IMPLICIT NONE
    
    !--Declaration des variables locales
    INTEGER  ::  i,j 
 
    real(rp)                 :: err = zero
    !*********************************************************
    !!!call init_h()   !! initialisation de la subdivision
    
     !Initialisation de la matrice A = mat_a
    mat_a_eq = zero
    DO i=2,L-1
       mat_a_eq(i,i)   = one/h1_demi(i) + one/h1_demi(i+1)
       mat_a_eq(i,i-1) = -one/h1_demi(i)
       mat_a_eq(i,i+1) = -one/h1_demi(i+1)	
    END DO

    mat_a_eq(1,1)   = one/h1_demi(1) + one/h1_demi(2)	
    mat_a_eq(1,2)   = -one/h1_demi(2)
    mat_a_eq(L,L)   = one/h1_demi(L) +one/h1_demi(L+1)
    mat_a_eq(L,L-1) = -one/h1_demi(L)
  END SUBROUTINE Init_A_eq

  !@fn
  !>initialisation du vecteur 2nd membre pour trouver psi à l'equilibre  

  !@fn
  !>Definition de la fonction du 2nd membre pour la resolution de psi à l'equilibre
  function func_b_eq(psi)
    implicit none
    real(rp),dimension(L),intent(in) :: psi
    real(rp),dimension(L)            :: func_b_eq
    integer                          :: i
    ! initialisation 
    func_b_eq=zero
    !dopage C(x) = 0

    do i=1,L
       func_b_eq(i) = h(i)*(delta_n*exp(psi(i))-delta_p*exp(-psi(i)) -VEC_C(i)) 
    end do
  end function func_b_eq

!!!!!!!!!!!********fonction pour résoudre avec newton àl'equilibre****!!!!!!!

  !@fn
  !>Definition de la fonction 2nd membre pour trouver psi par la methode de Newton
  function func_Eq(psi)
    implicit none
    real(rp),dimension(L),intent(in) :: psi
    real(rp),dimension(L)            :: func_Eq

    ! Initialisation
    bD = zero
    bD(1) = -psi_l/h1_demi(1)
    bD(L) = -psi_r/h1_demi(L+1)
    func_Eq = matmul(mat_a_eq,psi) + func_b_eq(psi) + bD

  end function func_Eq 
!!!!!!!*********!!!!!!!!!!!!

  
!!!!!!******Jacobienne pour newton!!!!!!!
  !@fn
  !>Definition de la jacobienne de la fonction F pour
  !>trouver psi avec la méthode de Newton
  function Jac_f(psi)
    implicit none
    real(rp), dimension(L),intent(in) :: psi
    real(rp), dimension(L,L)          :: diag
    real(rp), dimension(L,L)          :: Jac_f
    integer                           :: i

    ! initialisation
    diag = zero
    do i=1,L
       diag(i,i) = (delta_n*exp(psi(i)) + delta_p*exp(-psi(i)))*h(i)
    end do

    Jac_f = mat_a_eq + diag

  end function Jac_F
!!!!!!!!*********!!!!!!!!!!!!


END MODULE mod_initialisation_syst_lin

