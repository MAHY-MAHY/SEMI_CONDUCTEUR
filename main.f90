!>
!> 
!! @mainpage
  !>TITRE         : Projet effet Landau
  !! @date          : 28/03/2019
  !> @author
  !>Godfred AGYEKUM OHENEBA Vincent MAHY
!! \section Présentation
!! Le but de ce projet est de montrer la convergence à taux exponentielle d'un modèle de dérive diffusion.
!!On s'appuiera sur des schémas de volumes finies explicite et semi implicite.

 


PROGRAM semi_cond

  use mod_var

  use mod_affiche

  use mod_plot

  use mod_new

  use Resolution

  use init_semi

  use mod_allocations

  use mod_desallocations

  use mod_initialisation_syst_lin

  IMPLICIT NONE

  REAL      :: Debut, Fin,choice3 

  Call cpu_time(Debut)

  ! Initialisation des variables
  !
  machep = MAX (EPSILON(machep), 1E-15_rp)

  !  ! DESCRIPTION:
  !> Choisir L le nombre de point de points de la discretisation en x
  !>Choisir le schéma explicite ou semi implicite
  !> Choisir Delta_T le pas en temps du schéma , Attention à la condition CFL pour le schéma explicite
  !> Choisir  Tf le temps final
  !> Choisir le cas test et les conditions initales
  !> Choisir Le schéma
  ! lecture provisoire des dimensions
  WRITE(6,*)"Entrez L"
  READ(5,*) L
  WRITE(6,*)"select 0 pour le schéma explicite, et 1 pour le schéma semi implicite"
  READ(5,*)choice3
  WRITE(6,*)"Entrez DELTAT T, Attention  à la condition si schéma explicite"
  READ(5,*) DELTA_T
  WRITE(6,*)"ENTREz le temps final"
  READ(5,*) Tf
  M=aint(Tf/DELTA_T)
  WRITE(6,*)"M=",M
  WRIte(6,*)"select 0 pour cas test 1 avec dopage en creneaux 5 pour cas test 2"
  READ(5,*)choice
  WRIte(6,*)"select 0 pour schéma centré 1 pour decentré amont 2 pour Scharfetter"
  READ(5,*)choice2
  
  ! lecture de la donnéé initial
  !print*,"donnéé initiale"
  !read*,y_0	


  ! Allocation des matrices
  CALL allocations ()

  C=5   
  ! Initialisation des tableaux globaux
  call init_h()  

  CALL init_n0_p0()
  
  CALL Init_A()

  WRITE(6,*)"nl=",n_l,'nr=',n_r,'p_r=',p_r,'p_l=',p_l
  if(choice3<1)then
     CALL RES_P_N()
  ELSE
     CALL RES_SEMI_IMP()
  END if

  
  CALL eqt()

  call  gnuplot()

  !call Calc_E()

  !call Calc_I()

  ! Desallocation des matrices
  CALL desallocations ()

  CALL cpu_time(Fin)

  print '("Temps = ",f6.3," seconds.")',Fin-Debut 


  STOP
END PROGRAM semi_cond
