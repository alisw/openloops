!******************************************************************************!
! Copyright (C) 2014-2019 OpenLoops Collaboration. For authors see authors.txt !
!                                                                              !
! This file is part of OpenLoops.                                              !
!                                                                              !
! OpenLoops is free software: you can redistribute it and/or modify            !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! OpenLoops is distributed in the hope that it will be useful,                 !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                !
! GNU General Public License for more details.                                 !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with OpenLoops.  If not, see <http://www.gnu.org/licenses/>.           !
!******************************************************************************!


module ol_loop_propagators_/**/REALKIND
  implicit none
  contains

! **********************************************************************
subroutine prop_loop_A_Q(rank_in,rank_out,G_A,K,M,Gout_A)
! dressing anti-quark current with propagator
! ----------------------------------------------------------------------
! G_A(alpha,sigma,l)    = incoming anti-quark loop current
! rank_i/o              = length of tensor array up to highest incoming/outgoing rank
! K(4)                  = incoming momentum (light-cone rep)
! M                     = complex mass
! Gout_A(alpha,beta,lp) = outgoing dressed anti-quark loop current without I/(K^2-M^2) with increased rank
! HR(i,l)               = rank-raising function
! Gout_A(alpha,beta,l)  = G_A(alpha,sigma,l) * [-slash(K)+M](sigma,beta) + sum_{i=1}^4 G_A(alpha,sigma,LR(i,l))*[-gamma_i(sigma,beta)]
! gamma_i               = Dirac gamma matrix with covariant light-cone components
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  implicit none
  integer,           intent(in)  :: rank_in, rank_out
  complex(REALKIND), intent(in)  :: G_A(4,rank_in), K(5), M
  complex(REALKIND), intent(out) :: Gout_A(4,rank_out)
  integer :: l

  Gout_A = 0

  do l = 1, rank_in

    ! dressing proportional to loop momentum, raise the rank
    Gout_A(1,HR(1,l)) = Gout_A(1,HR(1,l)) + G_A(3,l)
    Gout_A(2,HR(1,l)) = Gout_A(2,HR(1,l))
    Gout_A(3,HR(1,l)) = Gout_A(3,HR(1,l))
    Gout_A(4,HR(1,l)) = Gout_A(4,HR(1,l)) + G_A(2,l)

    Gout_A(1,HR(2,l)) = Gout_A(1,HR(2,l))
    Gout_A(2,HR(2,l)) = Gout_A(2,HR(2,l)) + G_A(4,l)
    Gout_A(3,HR(2,l)) = Gout_A(3,HR(2,l)) + G_A(1,l)
    Gout_A(4,HR(2,l)) = Gout_A(4,HR(2,l))

    Gout_A(1,HR(3,l)) = Gout_A(1,HR(3,l)) + G_A(4,l)
    Gout_A(2,HR(3,l)) = Gout_A(2,HR(3,l))
    Gout_A(3,HR(3,l)) = Gout_A(3,HR(3,l)) - G_A(2,l)
    Gout_A(4,HR(3,l)) = Gout_A(4,HR(3,l))

    Gout_A(1,HR(4,l)) = Gout_A(1,HR(4,l))
    Gout_A(2,HR(4,l)) = Gout_A(2,HR(4,l)) + G_A(3,l)
    Gout_A(3,HR(4,l)) = Gout_A(3,HR(4,l))
    Gout_A(4,HR(4,l)) = Gout_A(4,HR(4,l)) - G_A(1,l)

    ! dressing proportional to tree level momentum, same rank
    Gout_A(1,l) = Gout_A(1,l) + K(1)*G_A(3,l) + K(3)*G_A(4,l) + M*G_A(1,l)
    Gout_A(2,l) = Gout_A(2,l) + K(2)*G_A(4,l) + K(4)*G_A(3,l) + M*G_A(2,l)
    Gout_A(3,l) = Gout_A(3,l) + K(2)*G_A(1,l) - K(3)*G_A(2,l) + M*G_A(3,l)
    Gout_A(4,l) = Gout_A(4,l) + K(1)*G_A(2,l) - K(4)*G_A(1,l) + M*G_A(4,l)
  end do
end subroutine prop_loop_A_Q


! **********************************************************************
subroutine prop_loop_Q_A(rank_in,rank_out,G_Q,K,M,Gout_Q)
! dressing quark current with propagator
! ----------------------------------------------------------------------
! G_Q(alpha,sigma,l)    = incoming quark loop current
! rank_i/o              = length of tensor array up to highest incoming/outgoing rank
! K(4)                  = incoming momentum (light-cone rep)
! M                     = complex mass
! Gout_Q(alpha,beta,lp) = outgoing dressed quark loop current without I/(K^2-M^2) with increased rank
! HR(i,l)               = rank-raising function
! Gout_Q(alpha,beta,l)  = [slash(K)+M](beta,sigma)*G_Q(alpha,sigma,l) + sum_{i=1}^4 [gamma_i(beta,sigma)]*G_Q(alpha,sigma,LR(i,l))
! gamma_i               = Dirac gamma matrix with covariant light-cone components
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_tensor_bookkeeping, only: HR
  implicit none
  integer,           intent(in)  :: rank_in, rank_out
  complex(REALKIND), intent(in)  :: G_Q(4,rank_in), K(5), M
  complex(REALKIND), intent(inout) :: Gout_Q(4,rank_out)
  integer :: l

  Gout_Q = 0

  do l = 1, rank_in

    ! dressing proportional to loop momentum, raise the rank
    Gout_Q(1,HR(1,l)) = Gout_Q(1,HR(1,l))
    Gout_Q(2,HR(1,l)) = Gout_Q(2,HR(1,l)) - G_Q(4,l)
    Gout_Q(3,HR(1,l)) = Gout_Q(3,HR(1,l)) - G_Q(1,l)
    Gout_Q(4,HR(1,l)) = Gout_Q(4,HR(1,l))

    Gout_Q(1,HR(2,l)) = Gout_Q(1,HR(2,l)) - G_Q(3,l)
    Gout_Q(2,HR(2,l)) = Gout_Q(2,HR(2,l))
    Gout_Q(3,HR(2,l)) = Gout_Q(3,HR(2,l))
    Gout_Q(4,HR(2,l)) = Gout_Q(4,HR(2,l)) - G_Q(2,l)

    Gout_Q(1,HR(3,l)) = Gout_Q(1,HR(3,l))
    Gout_Q(2,HR(3,l)) = Gout_Q(2,HR(3,l)) + G_Q(3,l)
    Gout_Q(3,HR(3,l)) = Gout_Q(3,HR(3,l))
    Gout_Q(4,HR(3,l)) = Gout_Q(4,HR(3,l)) - G_Q(1,l)

    Gout_Q(1,HR(4,l)) = Gout_Q(1,HR(4,l)) + G_Q(4,l)
    Gout_Q(2,HR(4,l)) = Gout_Q(2,HR(4,l))
    Gout_Q(3,HR(4,l)) = Gout_Q(3,HR(4,l)) - G_Q(2,l)
    Gout_Q(4,HR(4,l)) = Gout_Q(4,HR(4,l))

    ! dressing proportional to tree level momentum, same rank
    Gout_Q(1,l) = Gout_Q(1,l) - K(2)*G_Q(3,l) + K(4)*G_Q(4,l) + M*G_Q(1,l)
    Gout_Q(2,l) = Gout_Q(2,l) - K(1)*G_Q(4,l) + K(3)*G_Q(3,l) + M*G_Q(2,l)
    Gout_Q(3,l) = Gout_Q(3,l) - K(1)*G_Q(1,l) - K(4)*G_Q(2,l) + M*G_Q(3,l)
    Gout_Q(4,l) = Gout_Q(4,l) - K(2)*G_Q(2,l) - K(3)*G_Q(1,l) + M*G_Q(4,l)
  end do
end subroutine prop_loop_Q_A

end module ol_loop_propagators_/**/REALKIND


! **********************************************************************
!                            ALPHA INDEX INTERFACE
! ----------------------------------------------------------------------
! When is needed (e.g. fermion currents) it calls the vertex routine
! four times for the four components of the cut-open loop line
! **********************************************************************

module ol_prop_interface_/**/REALKIND
  implicit none
  contains

! **********************************************************************
subroutine loop_A_Q(Gin, K, M, Gout)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_propagators_/**/REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: Gin(:,:,:), K(5), M
  complex(REALKIND), intent(out) :: Gout(:,:,:)
  integer :: rank_in, rank_out
  rank_in  = size(Gin,2)
  rank_out = size(Gout,2)

  call prop_loop_A_Q(rank_in, rank_out, Gin(:,:,1), K, M, Gout(:,:,1))
  call prop_loop_A_Q(rank_in, rank_out, Gin(:,:,2), K, M, Gout(:,:,2))
  call prop_loop_A_Q(rank_in, rank_out, Gin(:,:,3), K, M, Gout(:,:,3))
  call prop_loop_A_Q(rank_in, rank_out, Gin(:,:,4), K, M, Gout(:,:,4))

end subroutine loop_A_Q


! **********************************************************************
subroutine loop_Q_A(G_Q, K, M, Gout_Q)
! **********************************************************************
  use KIND_TYPES, only: REALKIND
  use ol_loop_propagators_/**/REALKIND
  implicit none
  complex(REALKIND), intent(in)  :: G_Q(:,:,:), K(5), M
  complex(REALKIND), intent(out) :: Gout_Q(:,:,:)
  integer :: rank_in, rank_out
  rank_in  = size(G_Q,2)
  rank_out = size(Gout_Q,2)

  call prop_loop_Q_A(rank_in, rank_out, G_Q(:,:,1), K, M, Gout_Q(:,:,1))
  call prop_loop_Q_A(rank_in, rank_out, G_Q(:,:,2), K, M, Gout_Q(:,:,2))
  call prop_loop_Q_A(rank_in, rank_out, G_Q(:,:,3), K, M, Gout_Q(:,:,3))
  call prop_loop_Q_A(rank_in, rank_out, G_Q(:,:,4), K, M, Gout_Q(:,:,4))

end subroutine loop_Q_A

end module ol_prop_interface_/**/REALKIND
