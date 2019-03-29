!>
!> 
!! @mainpage
  !>TITRE         : Projet effet Landau
  !! @date          : 28/03/2019
  !> @author
  !> Vincent MAHY Godfred AGYEKUM OHENEBA
!! \section Présentation
!! Le but de ce projet est de mettre en evidence l'amortissement landau au cours du temps en discretisant soit par differences finies soit par une méthode semi Lagrangienne.
PROGRAM MAIN
  
  ! DESCRIPTION:
  !> Choisir Nx le nombre de point de points de la discretisation en x
  !> Choisir Nv le nombre de point de points de la discretisation en v
  !> Choisir  Tf le temps final
  !> Choisir le choix de la méthode de resolution
  !> Choisir 0.001<epsilon<0.1 et k

  use mod_var

  use mod_affiche

  use mod_plot

  use mod_new

  use Resolution

  use mod_allocations

  use mod_desallocations

  use mod_initialisation_syst_lin

  IMPLICIT NONE
  REAL      :: Debut, Fin

  Call cpu_time(Debut)

  ! Initialisation des variables
  !
  machep = MAX (EPSILON(machep), 1E-15_rp)

  !
  ! lecture provisoire des dimensions
  WRITE(6,*)"Entrer Nx"
  READ(5,*) Nx
  WRITE(6,*)"Entrer Nv"
  READ(5,*) Nv
  WRITE(6,*)"ENTRER le temps final"
  READ(5,*) Tf
  WRITE(6,*) " taper 1 pour differences finies 2 pour semi lagrngien"
  READ(5,*) choice
  print*,"epsilon"
  read*,eps
  WRITE(6,*)"entrer k"
  read(5,*)k
  L=4._RP*PI
  Vmax=10
  hv=2._rp*Vmax/Nv
  hx=L/Nx
  WRITE(6,*)"hx=",hx,"hv=",hv
  ! Allocation des matrices
  CALL allocations ()

  call init_h()

  ! Initialisation des tableaux globaux
  call init_f()

  call init_B()
  
  CALL Init_A()

  if(choice<2)then
     call RES_F()
  ELSE
     call RES_SEMI_LAG()
  END IF
  
  call  gnuplot()

  ! Desallocation des matrices
  CALL desallocations ()

  CALL cpu_time(Fin)

  print '("Temps = ",f6.3," seconds.")',Fin-Debut 


  STOP
END PROGRAM MAIN
