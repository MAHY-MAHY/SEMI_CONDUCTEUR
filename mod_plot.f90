!>Ecriture des vecteur des vecteurs de travail
!>dans des fichiers pour faire l'affichage.
!>Pour un affichage rapide avec python : python.plot.py
!>sinon avec gnuplot :
!> \n  plot 'Vec_psi.txt' pour psi à Tf
!> \n  plot 'Vec_N.txt'  pour N à Tf
!> \n  plot  'Vec_P.txt'  pour P à Tf
!> \n  plot 'Vec_psi_eq.txt' pour psi à l'equilibre
!> \n  plot   'Vec_n_eq.txt' pour n à l'equilibre
!> \n plot    'Vec_p_eq.txt'  pour p à l'equilibre



module mod_plot

  use mod_var

  use Resolution

  implicit none

  !
contains

  !@fn
  !> Calcul de la différence de psi, étape pour le calcul de l'energie
  function calc_D_psi(i)
    integer ,intent(in)::i
    REAL(RP) :: calc_D_psi
    REAL(RP) ::TMP_1,TMP_2
    IF(i==0)THEN
       TMP_1=VEC_PSI_APP(i+1)-psi_l
       TMP_2=VEC_PSI_APP_eq(i+1)-psi_l
       calc_D_psi=(TMP_1-TMP_2)/h1_demi(i+1)
    ELSEIF(i==L)THEN
       TMP_1=psi_r-VEC_PSI_APP(i)
       TMP_2=psi_r-VEC_PSI_APP_eq(i)
       calc_D_psi=(TMP_1-TMP_2)/h1_demi(i+1)
    ELSE
       TMP_1=VEC_PSI_APP(i+1)-VEC_PSI_APP(i)
       TMP_2=VEC_PSI_APP_eq(i+1)-VEC_PSI_APP_eq(i)
       calc_D_psi=(TMP_1-TMP_2)/h1_demi(i+1)
    END IF
  END function calc_D_psi
  !@fn
  !> Calcul de la fonction H, étape pour le calcul de l'energie
  function func_H(x)
    REAl(RP),intent(in) ::x
    REAL(RP) ::func_H
    func_H=x*log(x)-x+1
  END function func_H

  !@fn
  !> Calcul de l'energie au temps Tf
  Subroutine Calc_E()
    INTEGER ::i
    Energie=0
    DO i=1,L
       Energie =Energie+h(i)*(func_H(VEC_N(i))-func_H(VEC_N_eq(i))-log(VEC_N_eq(i))*&
            (VEC_N(i)-VEC_N_eq(i)))
       Energie=Energie+h(i)*(func_H(VEC_P(i))-func_H(VEC_P_eq(i))-log(VEC_P_eq(i))*&
            (VEC_p(i)-VEC_P_eq(i)))
    END DO
    DO i=0,L
       Energie=Energie+0.5*h1_demi(i+1)*calc_D_psi(i)
    END DO
    WRITE(6,*)Energie
  END Subroutine Calc_E

  !@fn
  !>Fonction de calcul de différence de N , une étape pour calculer la dissipation
  function calc_D_log_N(i)
    integer,intent(in) ::i
    REAL(RP) :: calc_D_log_N
    REAL(RP) ::TMP_1,TMP_2   
    IF(i==0)THEN
       TMP_1=log(n_l)-psi_l
       TMP_2=log(VEC_N(i+1))-VEC_PSI_APP(i+1)
       calc_D_log_N=(TMP_2-TMP_1)/h1_demi(i+1)
    ELSEIF(i==L)THEN
       TMP_1=log(VEC_N(i))-VEC_PSI_APP(i)
       TMP_2=log(n_r)-psi_r
       calc_D_log_N=(TMP_2-TMP_1)/h1_demi(i+1)
    ELSE
       TMP_1=log(VEC_N(i+1))-log(VEC_N(i))
       TMP_2=VEC_PSI_APP(i+1)-VEC_PSI_APP(i)
       calc_D_log_N=(TMP_1-TMP_2)/h1_demi(i+1)
    END IF
  END function calc_D_log_N

  !@fn
  !>Fonction de calcul de différence de P , une étape pour calculer la dissipation
  function calc_D_log_P(i)
    integer,intent(in) ::i
    REAL(RP) :: calc_D_log_P
    REAL(RP) ::TMP_1,TMP_2   
    IF(i==0)THEN
       TMP_1=log(p_l)+psi_l
       TMP_2=log(VEC_P(i+1))+VEC_PSI_APP(i+1)
       calc_D_log_P=(TMP_2-TMP_1)/h1_demi(i+1)
    ELSEIF(i==L)THEN
       TMP_1=log(VEC_P(i))+VEC_PSI_APP(i)
       TMP_2=log(p_r)+psi_r
       calc_D_log_P=(TMP_2-TMP_1)/h1_demi(i+1)
    ELSE
       TMP_1=log(VEC_P(i+1))-log(VEC_P(i))
       TMP_2=VEC_PSI_APP(i+1)-VEC_PSI_APP(i)
       calc_D_log_P=(TMP_1+TMP_2)/h1_demi(i+1)
    END IF
  END function calc_D_log_P



  !@fn
  !>Calcul de la dissipation à Tf 
  Subroutine Calc_I()
    INTEGER :: i
    Dissipation=0
    DO i=1,L-1
       Dissipation=Dissipation+h1_demi(i+1)*min(VEC_N(i+1),VEC_N(i))*&
            calc_D_log_N(i)**2
       Dissipation=Dissipation+h1_demi(i+1)*min(VEC_P(i+1),VEC_P(i))*&
            calc_D_log_P(i)**2
    END DO
    WRITE(6,*)Dissipation
  END Subroutine Calc_I




  !@fn
  !>écriture des vecteurs dans  des fichiers
  subroutine gnuplot()
    implicit none
    integer  :: j


    OPEN(UNIT=61,FILE='Vec_psi.txt')
    OPEN(UNIT=62,FILE='Vec_N.txt')
    OPEN(UNIT=63,FILE='Vec_P.txt')
    OPEN(UNIT=64,FILE='Vec_psi_eq.txt')
    OPEN(UNIT=65,FILE='Vec_n_eq.txt')
    OPEN(UNIT=66,FILE='Vec_p_eq.txt')
    Do j=1,L
       Write(61,*)X(j+1),VEC_PSI_APP(j)
       WRITE(62,*)X(j+1),Vec_N(j)
       WRITE(63,*)X(j+1),Vec_P(j)
       Write(64,*)X(j+1),VEC_PSI_APP_eq(j)
       Write(65,*)X(j+1),VEC_N_eq(j)
       Write(66,*)X(j+1),VEC_P_eq(j)
    end Do

  end subroutine gnuplot





end module mod_plot
