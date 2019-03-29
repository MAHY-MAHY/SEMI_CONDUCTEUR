!>Resolution de la méthode explicite et de l'equilibre

Module Resolution


  USE numerics
  !
  USE mod_var

  USE mod_initialisation_syst_lin

  USE init_semi

  IMPLICIT NONE
  
CONTAINS


  !@fn
  !>Calcul du flux du schéma pour trouver N pour la méthode explicite ,
  !>avec prise en compte des conditions de bords
  Function Flux_G(i)
    INTEGER ,intent(in)::i
    REAL(RP)::Flux_G
    REAL(RP) ::TMP_PSI
    IF(i>0 .AND. i<L)THEN
       TMP_PSI=VEC_PSI_APP(i)-VEC_PSI_APP(i+1)
       FLux_G=func_B(TMP_PSI)*VEC_N(i)
       TMP_PSI=-TMP_PSI
       FLux_G=FLux_G-func_B(TMP_PSI)*VEC_N(i+1)
       FLux_G=FLux_G/h1_demi(i+1)
    ELSEIF(i<1)THEN
       TMP_PSI=psi_l-VEC_PSI_APP(i+1)
       FLux_G=func_B(TMP_PSI)*n_l
       TMP_PSI=-TMP_PSI
       FLux_G=FLux_G-func_B(TMP_PSI)*VEC_N(i+1)
       FLux_G=FLux_G/h1_demi(1)
    ELSE
       TMP_PSI=VEC_PSI_APP(i)-psi_r
       FLux_G=func_B(TMP_PSI)*VEC_N(i)
       TMP_PSI=-TMP_PSI
       FLux_G=FLux_G-func_B(TMP_PSI)*n_r
       FLux_G=FLux_G/h1_demi(L+1)
    END IF
  END Function Flux_G

    !@fn
  !>Calcul du flux du schéma pour trouver P pour la méthode explicite ,
  !>avec prise en compte des conditions de bords
  Function Flux_H(i)
    INTEGER ,intent(in)::i
    REAL(RP)::Flux_H
    REAL(RP) ::TMP_PSI
    IF(i>0 .AND. i<L)THEN
       TMP_PSI=VEC_PSI_APP(i+1)-VEC_PSI_APP(i)
       FLux_H=func_B(TMP_PSI)*VEC_P(i)
       TMP_PSI=-TMP_PSI
       FLux_H=FLux_H-func_B(TMP_PSI)*VEC_P(i+1)
       FLux_h=FLux_H/h1_demi(i+1)
    ELSEIF(i<1)THEN
       TMP_PSI=VEC_PSI_APP(i+1)-psi_l
       FLux_H=func_B(TMP_PSI)*p_l
       TMP_PSI=-TMP_PSI
       FLux_H=FLux_H-func_B(TMP_PSI)*VEC_P(i+1)
       FLux_H=FLux_H/h1_demi(1)
    ELSE
       TMP_PSI=psi_r-VEC_PSI_APP(i)
       FLux_H=func_B(TMP_PSI)*VEC_P(i)
       TMP_PSI=-TMP_PSI
       FLux_H=FLux_H-func_B(TMP_PSI)*p_r
       FLux_H=FLux_H/h1_demi(L+1)
    END IF

  END Function Flux_H


  !@fn
  !> Resolution du système linéaire Apsi=B pour trouver psi
  SUBROUTINE RES_PSI()      !!!!!cette subroutine est deja verif
    real(rp), dimension(L)   :: TMP_PSI
    real(RP), DIMENSION(L,L)::  TMP_A
    ! DEcomposition LU de mat_a
    TMP_A=mat_A
    CALL INIT_B()
    TMP_PSI = vec_b
    call dgetrf(L,L,TMP_A,L,lu_piv,ierr)
    ! solver u = A^{-1}b
    call dgetrs('N',L,1,TMP_A,L,Lu_piv,TMP_PSI,L,ierr)
    ! solution exacte 
    VEC_PSI_APP = TMP_PSI
  END SUBROUTINE RES_PSI

  
  !@fn
  !>Mis en oeuvre du schéma de volumes finis et calcul de p,n et psi
  !>par la méthode explicite
  subroutine RES_P_N()
    implicit none
    INTEGER::i,j,k
    REAL(rp),DIMENSION(L) :: TMP_P
    REAL(rp),DIMENSION(L) :: TMP_N
    k=1
    CALL init_n0_p0
    DO i=1,M
       CALL RES_PSI()
       TMP_P=VEC_P
       TMP_N=VEC_N
       DO j=1,L
          TMP_P(j)=VEC_P(j)+DELTA_T*(Flux_H(j-1)-Flux_H(j))/h(j)
          TMP_N(j)=VEC_N(j)+DELTA_T*(Flux_G(j-1)-Flux_G(j))/h(j)
       END DO
       IF(i*1._RP/M>0.1*k)THEN
          WRITE(6,*)i*1._RP/M*100
          k=k+1
       END IF
       VEC_P=TMP_P
       VEC_N=TMP_N
       IF(mod(i,100)==0)then
          WRITE(6,*)i*1._RP/M*100
       END IF

    END DO

  end subroutine RES_P_N




!!!!!!!***** RESOLUTION SEMI_IMPLICITE****!!!!!

  !

  !@fn
  !>Resoltuion des differents systèmes A1n=B1, A2p=B2, A3psi=B3.
  !>Résolution par la méthode semi implicite
    SUBROUTINE RES_SEMi_IMP()      !!!!!cette subroutine est deja verif
    real(rp), dimension(L)   :: PIV_N,PIV_P
    real(RP), DIMENSION(L)::  TMP_N,TMP_P
    integer ::i,k
    k=1
    ! DEcomposition LU de mat_a
    CALL init_n0_p0
    !WRITE(6,*)"pour p=",ierr
    ! solver u = A^{-1}b
    DO i=1,M
       call RES_PSI()
       CALL init_MAT_semi_imp_P()
       CALL init_MAT_semi_imp_N()
       call dgetrf(L,L,MAT_semi_imp_N,L,PIV_N,ierr)
       call dgetrf(L,L,MAT_semi_imp_P,L,PIV_P,ierr)
       call init_2nd_semi_imp()
       TMP_N=VEC_semi_imp_N
       TMP_P=VEC_semi_imp_P
       call dgetrs('N',L,1,MAT_semi_imp_N,L,PIV_N,TMP_N,L,ierr)
       call dgetrs('N',L,1,MAT_semi_imp_P,L,PIV_N,TMP_P,L,ierr)
       IF(i*1._RP/M>0.1*k)THEN
          WRITE(6,*)i*1._RP/M*100
          k=k+1
       END IF
       VEC_N=TMP_N
       VEC_P=TMP_P
       PIV_P=0
       PIV_N=0
    END DO
    
  END SUBROUTINE RES_SEMI_IMP
  


!!!!!!!******* Resolution pour l'equilibre *****!!!!!!!
  !@fn
  !>Calcul de psi à l'equilibre
  subroutine eqt()
    implicit none
    real(rp),dimension(L)               :: v,w
    real(rp),dimension(L,L)             :: matr
    real(rp)                            :: tol,err
    integer                             :: iter_max,i

    ! initialisation
    err = one
    i   = 0
    print*,"tol =","iter_max ="
    read*,tol,iter_max

    call init_A_eq
    do while ( err>tol.and.i<iter_max)	
       ! DEcomposition LU de mat_a
       w = func_Eq(VEC_PSI_APP_eq)
       matr=Jac_f(VEC_PSI_APP_eq)
       !Write(6,*)matr    
       call dgetrf(L,L,matr,L,lu_piv,ierr)
       ! solver y = Jf^{-1}f(x_0)
       call dgetrs('N',L,1,matr,L,Lu_piv,w,L,ierr)
       v = VEC_PSI_APP_eq - w
       err = norm2(v-VEC_PSI_APP_eq)
       VEC_PSI_APP_eq = v
       i = i+1
    end do
    call calc_p_n_eqt()

  end subroutine eqt

  !@fn
  !>calcul de n et p à l'equilibre sachant n=delta_n*exp(psi) et p=delta_p*exp(-psi)
  subroutine calc_p_n_eqt()
    integer :: i
    Do i=1,L
       VEC_N_eq(i)=delta_n*exp(VEC_PSI_APP_eq(i))
       VEC_P_eq(i)=delta_n*exp(-VEC_PSI_APP_eq(i))
    END Do
  END subroutine calc_p_n_eqt
  

END Module Resolution
