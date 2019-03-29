
!>Calcul de certains paramètre utile comme 0.5, l'erreur machine etc
Module numerics
  ! ---------------------------------------------------------------------
  Integer, Parameter  ::  sp = kind(1.0  )
  Integer, Parameter  ::  dp = kind(1.0D0)
  Integer, Parameter  ::  rp = dp

  Integer, Save       ::  ierr = 0

  Real(rp), Parameter ::  one  = 1.0_rp
  Real(rp), Parameter ::  half = 0.5_rp
  Real(rp), Parameter ::  zero = 0.0_rp

  Real(rp), Save      ::  machep
  
  ! ---------------------------------------------------------------------
End Module numerics
