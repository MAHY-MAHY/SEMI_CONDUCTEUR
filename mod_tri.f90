!> Tri d'un vecteur, MODULE NON UTILISé

module mod_tri
  use mod_var

contains

  subroutine tri(t)
    implicit none

    real(rp),dimension(:),intent(inout) :: t
    integer                             :: i,k

    real(rp)                            :: temp

    do i = 2, L
       k=i-1; temp = t(i)
       do while (k>=1 .and. temp < t(k) )
          t(k+1) = t(k); k=k-1
       end do
       t(k+1) = temp
    end do

  end subroutine tri


end module mod_tri
