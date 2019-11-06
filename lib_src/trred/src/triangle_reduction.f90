!******************************************************************************!
!                                                                              !
!    triangle_reduction.f90                                                    !
!    is part of trred & OpenLoops2                                             !
!    Copyright (C) 2017-2019  For authors see authors.txt.                     !
!                                                                              !
!    trred has been developed by J.-N. Lang, H. Zhang and F. Buccioni          !
!    trred is licenced under the GNU GPL version 3,                            !
!    see COPYING for details.                                                  !
!                                                                              !
!******************************************************************************!

module triangle_reduction_DP
  use triangle_expansion_DP, only: B0d_coeff,B0d_list,B0d_matrix,         &
                                B0d_0_coeff,B0d_0_list,                   &
                                C0d_0mm_coeff,C0d_0mm_matrix,             &
                                C0d_m00_coeff,C0d_000_coeff,              &
                                C0d_000_EP1_coeff,                        &
                                C0d_mmm_table,B0d_mm_opt,B0d_mm_opt_list, &
                                C0d_m0m1m1_coeff,C0d_m0m1m1_list,         &
                                B0d_m0m1_coeff,B0d_m0m1_list
  use triangle_aux_DP, only: dp,cone,cnul,A0,B0_zero,rnul
  implicit none

  real(dp), parameter :: large_z_threshold = 0.1_dp

  interface B0d
    module procedure B0d_coeff, B0d_list, B0d_matrix
  end interface

  interface B0d_0
    module procedure B0d_0_coeff, B0d_0_list
  end interface

  interface B0d_mm
    module procedure B0d_mm_opt, B0d_mm_opt_list
  end interface

  interface B0d_m0m1
    module procedure B0d_m0m1_coeff, B0d_m0m1_list
  end interface

  interface C0d_0mm
    module procedure C0d_0mm_coeff,C0d_0mm_matrix
  end interface

  interface C0d_m00
    module procedure C0d_m00_coeff
  end interface

  interface C0d_mmm
    module procedure C0d_mmm_table
  end interface

  interface C0d_000
    module procedure C0d_000_coeff
  end interface

  interface C0d_000_EP1
    module procedure C0d_000_EP1_coeff
  end interface

  interface C0d_m0m1m1
    module procedure C0d_m0m1m1_coeff, C0d_m0m1m1_list
  end interface

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     0mm                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_0mm_P12_DP(p2,m2,d,expansion) result(C_P12)
    ! computes the coefficient of the rank 1 p1[mu]-p2[mu] form factor
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    logical, optional, intent(in) :: expansion
    complex(dp) :: C_P12

    C_P12 = B0d(p2,m2,d,[cnul,cone,2*cone])/p2  + (m2/p2-cone)*C0d_0mm(p2,m2,d,1)

  end function

  function C_0mm_p1_DP(p2,m2,d) result(C_p1)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1

    C_p1 = B0d(p2,m2,d,1)/p2 - C0d_0mm(p2,m2,d)

  end function

  subroutine C_0mm_rank1(p2,m2,d,C_P12,C_p1)
    complex(dp),     intent(in)  :: m2
    real(dp),        intent(in)  :: p2,d
    complex(dp),     intent(out) :: C_P12,C_p1
    complex(dp) :: C(2)

    C = B0d_matrix(p2,m2,d,2,3,[[cnul,cnul], [cone/p2, cone/p2], [2*cone/p2,cnul]],[.true.,.false.,.false.]) &
        + C0d_0mm_matrix(p2,m2,d,2,2,[[cnul,-cone],[(m2/p2-cone),cnul]],[.false.,.false.])

    C_P12 = C(1)
    C_p1 = C(2)

  end subroutine C_0mm_rank1

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_0mm_p1p1_DP(p2,m2,d) result(C_p1p1)
    use triangle_aux_DP, only: A0mB0_0mm_p1p1
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1,b0,ab,a,c0

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_p1p1 = + C0d_0mm(p2,m2,d,0)                &
               - (A0(m2)/(2*p2**2*(1 + d)))        &
               + B0d(p2,m2,d,[+m2/(2*p2**2*(1+d)), &
                                    -m2/(2*p2**2)-3/(2*p2)])
    else
      ab = -cone/(2*p2**2*(1 + d))
      a = ab*A0mB0_0mm_p1p1(p2,m2,d)
      b0 = -d*(m2*ab*d+m2/(2*p2**2)+3/(2*p2))*B0d(p2,m2,d,2)
      c0 = C0d_0mm(p2,m2,d,0)
      C_p1p1 = a + b0 + c0
    end if
  end function

  function C_0mm_p1P12_DP(p2,m2,d) result(C_p1P12)
    use triangle_aux_DP, only: A0mB0_0mm_p1P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12,ab,a,b0,c0

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_p1P12 = + A0(m2)/(2*p2**2*(1 + d))                       &
                + B0d(p2,m2,d,[-m2/(2*p2**2*(1+d)),              &
                                          +(m2/(2*p2**2)-1/p2),        &
                                          +(m2/(2*p2**2) - 5/(2*p2))]) &
                + (1-(2*m2)/(p2))*C0d_0mm(p2,m2,d,1)
    else
      ab = cone/(2*p2**2*(1 + d))
      a = ab*A0mB0_0mm_p1P12(p2,m2,d)
      b0 = (m2+2*d*m2-(1+d)*(5+2*d)*p2)/(2*(1+d)*p2**2)*B0d(p2,m2,d,2)
      c0 = (1 - (2*m2)/(p2))*C0d_0mm(p2,m2,d,1)
      C_p1P12 = a + b0 + c0
    end if
  end function

  function C_0mm_P12P12_DP(p2,m2,d) result(C_P12P12)
    use triangle_aux_DP, only: A0mB0_0mm_P12P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12,ab,a,b0,c0

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_P12P12 = - A0(m2)/(2*p2**2*(1 + d))             &
                 + B0d(p2,m2,d, [m2/(2*p2**2*(1 + d)),  &
                                 (1/(2*p2)- m2/(2*p2**2)), &
                                 (2*m2/p2**2 - 1/p2),      &
                                 (3*m2/p2**2 - 3/p2)])     &
                 + (1 + m2**2/(p2**2) - 4*m2/p2)*C0d_0mm(p2,m2,d,2)
    else
      ab = -cone/(2*p2**2*(1 + d))
      a = ab*A0mB0_0mm_P12P12(p2,m2,d)
      b0 = B0d(p2,m2,d,[cnul,cnul,-ab*((4+3*d)*m2+(-2+d)*(1+d)*p2), &
                        (3*m2/p2**2 - 3/p2)])
      c0 = (1 + m2**2/(p2**2) - 4*m2/p2)*C0d_0mm(p2,m2,d,2)
      C_P12P12 = a + b0 + c0
    end if
  end function

  function C_0mm_g_DP(p2,m2,d) result(C_gMunu)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gMunu

    C_gMunu = B0d(p2,m2,d,[cone/4,(p2-m2)/(4*p2)]) &
              + (2*m2*C0d_0mm(p2,m2,d) + 1)/4
  end function

  subroutine C_0mm_rank2(p2,m2,d,C_P12P12,C_p1p1,C_p1P12,C_g)
    complex(dp),     intent(in)  :: m2
    real(dp),        intent(in)  :: p2,d
    complex(dp),     intent(out) :: C_P12P12,C_p1p1,C_p1P12,C_g
    complex(dp) :: C(4),a1

    C = B0d_matrix(p2,m2,d,4,4, &
    [[m2/(2*p2**2*(1+d)),       m2/(2*p2**2*(1+d)),    -m2/(2*p2**2*(1+d)), cone/4], &
     [(1/(2*p2)- m2/(2*p2**2)),-m2/(2*p2**2)-3/(2*p2), +(m2/(2*p2**2)-1/p2),(p2-m2)/(4*p2)], &
     [(2*m2/p2**2-1/p2),        cnul,                  +(m2/(2*p2**2)-5/(2*p2)), cnul], &
     [(3*m2/p2**2 - 3/p2), cnul, cnul,cnul]],[.false.,.false.,.false.,.false.]) &
    + C0d_0mm_matrix(p2,m2,d,4,3, &
    [[cnul,cone,cnul,m2/2], &
     [cnul,cnul,(1 - (2*m2)/(p2)),cnul],&
     [(1 + m2**2/(p2**2) - 4*m2/p2),cnul,cnul,cnul]],[.false.,.false.,.false.])

    a1 = A0(m2)
    C = C + [-a1/(2*p2**2*(1 + d)),-a1/(2*p2**2*(1 + d)),a1/(2*p2**2*(1 + d)),cone/4]
    C_P12P12 = C(1)
    C_p1p1 = C(2)
    C_p1P12 = C(3)
    C_g = C(4)

  end subroutine C_0mm_rank2

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_0mm_p1p1p1_DP(p2,m2,d) result(C_p1p1p1)
    use triangle_aux_DP, only: A0mB0_0mm_p1p1p1
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1p1,ab,a,b0,c0,r

    ab = m2/(3*p2**3*(1+d)**2)+(2*m2+7*p2)/(6*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,(2*m2**2+7*m2*p2+11*p2**2)/(6*p2**3)])
      c0 = -C0d_0mm(p2,m2,d,0)
      r = - m2/(6*(1+d)*p2**2)
      C_p1p1p1 = a + b0 + c0 + r
    else
      a = ab*A0mB0_0mm_p1p1p1(p2,m2,d)
      b0 = (-d**3*m2*ab+ &
            d**2*(2*m2**2+7*m2*p2+11*p2**2)/(6*(p2**3)))*B0d(p2,m2,d,3)
      c0 = -C0d_0mm(p2,m2,d,0)
      C_p1p1p1 = a + b0 + c0
    end if
  end function

  function C_0mm_p1p1P12_DP(p2,m2,d) result(C_p1p1P12)
    use triangle_aux_DP, only: A0mB0_0mm_p1p1P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1P12,c1,b,a,ab,r

    ab = -(m2/(3*p2**3*(1+d)**2))+(-m2-7*p2)/(6*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b = B0d(p2,m2,d,[-ab*m2,                             &
                       (-m2**2-7*m2*p2+6*p2**2)/(6*p2**3), &
                       -(m2**2+8*m2*p2-17*p2**2)/(6*p2**3)])
      c1 = -(p2-3*m2)*C0d_0mm(p2,m2,d,1)/p2
      r = m2*p2/(6*(1+d)*p2**3)
      C_p1p1P12 =  a + b + c1 + r
    else
      a = ab*A0mB0_0mm_p1p1P12(p2,m2,d)
      c1 = -(p2-3*m2)/p2*C0d_0mm(p2,m2,d,1)
      b = (d**3*((m2**2)/(3*(1+d)**2*p2**3)             &
           +(m2**2 + 7*m2*p2)/(6*(1+d)*p2**3))          &
           +d**2*(-m2**2 - 7*m2*p2 + 6*p2**2)/(6*p2**3) &
           -d*(m2**2 + 8*m2*p2 - 17*p2**2)/(6*p2**3))*B0d(p2,m2,d,3)

      C_p1p1P12 = a + b + c1
    end if

  end function

  function C_0mm_p1P12P12_DP(p2,m2,d) result(C_p1P12P12)
    use triangle_aux_DP, only: A0mB0_0mm_p1P12P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12P12,b0,c2,ab,a,r

    if (abs(p2/m2) .gt. large_z_threshold) then
      ab = m2/(3*p2**3*(1+d)**2)+7/(6*p2**2*(1+d))
      a = ab*A0(m2)
      b0 = B0d_list(p2,m2,d,[-ab*m2,(7*m2 - 3*p2)/(6*p2**2),         &
                             (m2**2 - 11*m2*p2 + 3*p2**2)/(3*p2**3), &
                             (m2**2 - 19*m2*p2 + 10*p2**2)/(3*p2**3) ])

      r = -m2/(6*p2**2*(1+d))
      c2 = (-3*m2**2+6*m2*p2-p2**2)/(p2**2)*C0d_0mm(p2,m2,d,2)
      C_p1P12P12 = a + b0 + c2 + r
    else
      ab = m2/(3*p2**3*(1+d)**2)+7/(6*p2**2*(1+d))
      a = ab*A0mB0_0mm_p1P12P12(p2,m2,d)
      b0 = (-d**3*ab*m2 +                          &
            d**2*(7*m2-3*p2)/(6*p2**2) +           &
            d*(m2**2-11*m2*p2+3*p2**2)/(3*p2**3) + &
            (m2**2-19*m2*p2+10*p2**2)/(3*p2**3))*B0d(p2,m2,d,3)
      c2 = (-3*m2**2+6*m2*p2-p2**2)/(p2**2)*C0d_0mm(p2,m2,d,2)
      C_p1P12P12 = a + b0 + c2

    end if
  end function

  function C_0mm_P12P12P12_DP(p2,m2,d) result(C_P12P12P12)
    use triangle_aux_DP, only: A0mB0_0mm_P12P12P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12P12,b0,c3,ab,a,r

    ab = -(m2/(3*p2**3*(1+d)**2))+(m2-7*p2)/(6*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,                                    &
                        (m2**2 - 7*m2*p2 + 2*p2**2)/(6*p2**3),     &
                        (-3*m2**2 + 14*m2*p2 - 3*p2**2)/(6*p2**3), &
                        (8*m2**2 - 19*m2*p2 + 3*p2**2)/(3*p2**3),  &
                        (11*m2**2 - 38*m2*p2 + 11*p2**2)/(3*p2**3)])
      c3 = (3*m2**4 - 27*m2**3*p2 + &
            27*m2**2*p2**2 - 3*m2*p2**3)/(3*m2*p2**3)*C0d_0mm(p2,m2,d,3)
      r = (m2*p2)/(6*(1+d)*p2**3)
      C_P12P12P12 = a + b0 + c3 + r
    else
      a = ab*A0mB0_0mm_P12P12P12(p2,m2,d)
      b0 = B0d(p2,m2,d,[cnul,cnul,cnul,-d**3*ab*m2 + &
                        d**2*(m2**2-7*m2*p2+2*p2**2)/(6*p2**3) + &
                        d*(-3*m2**2+14*m2*p2-3*p2**2)/(6*p2**3) + &
                        (8*m2**2-19*m2*p2+3*p2**2)/(3*p2**3), &
                        (11*m2**2-38*m2*p2+11*p2**2)/(3*p2**3)])
      c3 = (m2**3-9*m2**2*p2+9*m2*p2**2-p2**3)/(p2**3)*C0d_0mm(p2,m2,d,3)
      C_P12P12P12 = a + b0 + c3
    end if
  end function

  function C_0mm_gp1_DP(p2,m2,d) result(C_gp1)
    use triangle_aux_DP, only: A0mB0_0mm_gp1
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gp1,b0,c0,b,ab,a

    if (abs(p2/m2) .gt. large_z_threshold) then
      b = B0d(p2,m2,d,[-cone/6 - m2**2/(12*(1+d)*p2**2), &
                       (m2**2 + 5*m2*p2 - 2*p2**2)/(12*p2**2)])
      c0 = C0d_0mm(p2,m2,d,0)
      C_gp1 = b +(- 18*m2*c0 - 7)/36                   &
                +(m2*A0(m2))/(12*(1+d)*p2**2)
    else
      ab = m2/(12*p2**2*(1+d))
      a = ab*A0mB0_0mm_gp1(p2,m2,d)
      b0 = B0d(p2,m2,d,[-cone/6,cnul, &
               -d**2*(m2*ab)+d*(m2**2+5*m2*p2-2*p2**2)/(12*p2**2)])
      c0 = -m2*C0d_0mm(p2,m2,d,0)/2
      C_gp1 = a + b0 + c0
    end if
  end function

  function C_0mm_gP12_DP(p2,m2,d) result(C_gP12)
    use triangle_aux_DP, only: A0mB0_0mm_gP12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gP12,b0,c1,b,ab,a

    if (abs(p2/m2) .gt. large_z_threshold) then
      b = B0d(p2,m2,d,[cone/12 + m2**2/(12*(1+d)*p2**2), &
                       (-m2**2 + 5*m2*p2)/(12*p2**2) , &
                       (-m2**2 + 10*m2*p2 - p2**2)/(12*p2**2)])
      c1 = C0d_0mm(p2,m2,d,1)
      C_gP12 = (6*m2**2*p2 - 6*m2*p2**2)/(12*p2**2)*c1 &
                  +(-m2*A0(m2))/(12*(1+d)*p2**2) &
                  +cone/18 + b
    else
      ab = -m2/(12*p2**2*(1+d))
      a = ab*A0mB0_0mm_gP12(p2,m2,d)
      b0 = +B0d(p2,m2,d,[cone/12,cnul,(-m2*ab*d**2) + &
                         d*(-m2**2+5*m2*p2)/(12*p2**2) + &
                         (-m2**2+10*m2*p2-p2**2)/(12*p2**2)])
      c1 = (m2**2 - m2*p2)/(2*p2)*C0d_0mm(p2,m2,d,1)
      C_gP12 = a + b0 + c1
    end if
  end function


  subroutine C_0mm_rank3(p2,m2,d,C_P12P12P12,C_p1P12P12,C_p1p1P12,C_p1p1p1,C_gP12,C_gp1)
    complex(dp),     intent(in)  :: m2
    real(dp),        intent(in)  :: p2,d
    complex(dp),     intent(out) :: C_P12P12P12,C_p1P12P12,C_p1p1P12,C_p1p1p1,C_gP12,C_gp1
    complex(dp) :: C(6),a1,b,c0

    C = B0d_matrix(p2,m2,d,6,5, &
    [[ m2**2/(3*p2**3*(1+d)**2)-m2*(m2-7*p2)/(6*p2**3*(1+d)),&
      -m2**2/(3*p2**3*(1+d)**2)-7*m2/(6*p2**2*(1+d)), &
       m2**2/(3*p2**3*(1+d)**2) + m2*(m2+7*p2)/(6*p2**3*(1+d)), &
      -m2**2/(3*p2**3*(1+d)**2) - m2*(2*m2+7*p2)/(6*p2**3*(1+d)), &
       cone/12+m2**2/(12*(1+d)*p2**2), &
      -cone/6-m2**2/(12*(1+d)*p2**2)], &
     [ (m2**2 - 7*m2*p2 + 2*p2**2)/(6*p2**3), &
       (7*m2-3*p2)/(6*p2**2), &
       (-m2**2-7*m2*p2+6*p2**2)/(6*p2**3), &
       (2*m2**2+7*m2*p2+11*p2**2)/(6*p2**3), &
       (-m2**2 + 5*m2*p2)/(12*p2**2), &
       (m2**2+5*m2*p2-2*p2**2)/(12*p2**2)], &
     [ (-3*m2**2 + 14*m2*p2 - 3*p2**2)/(6*p2**3), &
       (m2**2 - 11*m2*p2 + 3*p2**2)/(3*p2**3), &
      -(m2**2+8*m2*p2-17*p2**2)/(6*p2**3), &
       cnul,(-m2**2 + 10*m2*p2 - p2**2)/(12*p2**2),cnul], &
     [ (8*m2**2 - 19*m2*p2 + 3*p2**2)/(3*p2**3), &
       (m2**2 - 19*m2*p2 + 10*p2**2)/(3*p2**3), &
       cnul,cnul,cnul,cnul], &
     [(11*m2**2 - 38*m2*p2 + 11*p2**2)/(3*p2**3),cnul,cnul,cnul,cnul,cnul]],&
     [.false.,.false.,.false.,.false.,.false.])


    C = C + C0d_0mm_matrix(p2,m2,d,6,4, &
   [[cnul,cnul,cnul,-cone,cnul,-m2/2], &
    [cnul,cnul,-cone+3*m2/p2,cnul,(m2**2 - m2*p2)/(2*p2),cnul], &
    [cnul,(-3*m2**2+6*m2*p2-p2**2)/(p2**2),cnul,cnul,cnul,cnul], &
    [(m2**3-9*m2**2*p2+9*m2*p2**2-p2**3)/(p2**3),cnul,cnul,cnul,cnul,cnul]],&
   [.false.,.false.,.false.,.false.])

    a1 = A0(m2)
    C = C + [-(m2/(3*p2**3*(1+d)**2))+(m2-7*p2)/(6*p2**3*(1+d)),  &
               m2/(3*p2**3*(1+d)**2)+7/(6*p2**2*(1+d)),           &
             -(m2/(3*p2**3*(1+d)**2))+(-m2-7*p2)/(6*p2**3*(1+d)), &
               m2/(3*p2**3*(1+d)**2)+(2*m2+7*p2)/(6*p2**3*(1+d)), &
             -m2/(12*(1+d)*p2**2),                                &
             +m2/(12*(1+d)*p2**2)] * a1

    C = C + [(m2*p2)/(6*(1+d)*p2**3), &
             -m2/(6*p2**2*(1+d)), &
              m2*p2/(6*(1+d)*p2**3), &
            - m2/(6*(1+d)*p2**2), &
            + cone/18, &
            - 7*cone/36]

    C_P12P12P12 = C(1)
    C_p1P12P12 = C(2)
    C_p1p1P12 = C(3)
    C_p1p1p1 = C(4)
    C_gP12 = C(5)
    C_gp1 = C(6)

  end subroutine C_0mm_rank3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     m00                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_m00_p1_DP(p2,m2,d,expansion) result(C_p1)
    ! computes the coefficient of the rank 1 p1[mu] form factor
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    logical, optional, intent(in) :: expansion
    complex(dp) :: C_p1

    if (present(expansion) .and. .not. expansion) then
      ! tensor reduction result with no \delta expansion
      C_p1 = - C0d_m00(p2,m2,d)                 &
             + B0d(p2,m2,d,0)/(p2*d)            &
             - B0d(p2,m2,rnul,0)/(p2*d)
    else
      C_p1 = B0d(p2,m2,d,1)/p2 - C0d_m00(p2,m2,d)
    end if

  end function

  function C_m00_P12_DP(p2,m2,d,expansion) result(C_P12)
    ! computes the coefficient of the rank 1 p1[mu]-p2[mu] form factor
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    logical, optional, intent(in) :: expansion
    complex(dp) :: C_P12

    if (present(expansion) .and. .not. expansion) then
      ! tensor reduction result with no \delta expansion
      C_P12 = -(m2 + p2)/(p2*d)*C0d_m00(p2,m2,d,0) &
              + (2 + d)/(p2*d**2)*B0d(p2,m2,d,0)   &
              - B0_zero()/(p2*d)                   &
              - 2/(p2*d**2)*B0d(p2,m2,rnul,0)
    else
      C_P12 = B0d(p2,m2,d,[cnul,cone,2*cone])/p2 &
              + (-m2/p2-cone)*C0d_m00(p2,m2,d,1)
    end if

  end function

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_m00_P12P12_DP(p2,m2,d) result(C_P12P12)
    use triangle_aux_DP, only: A0mB0
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12,c1,c2,c3

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_P12P12 = + A0(m2)/(2*(cone + d))          &
                 + B0d(p2,m2,d,[-m2/(2*(cone+d)), &
                                +(p2+m2)/2,       &
                                -(2*m2 + p2),     &
                                -3*(m2 + p2)])    &
                 + (p2**2 + m2**2 + 2*m2*p2)*C0d_m00(p2,m2,d,2)
      C_P12P12 = C_P12P12/p2**2
    else
      ! For small p2/m2 the formula above is numerically unstable due to large
      ! cancellations between A0 and B0, DB0.

      ! Observation:
      ! The numerical instability resides in A0 and the first two (truncated)
      ! B0 functions:
      ! A0 - m2 B0^0 + (1+d) (p2+m2) B0^1

      ! Used identities:
      ! -B0^0 + d B0^1 = B0^{0,(0)}    (only the first coefficients of B0^0)
      ! B0^1 = B0^{1,(0)} + d B0^2
      ! A0mB0 = A0(m2)-m2 B0(p2,0,m2) + (p2+m2) B0^{1,(0)}(p2,0,m2)
      c1 = + A0mB0(p2,m2)/(2*p2**2*(cone + d))
      c2 = + B0d(p2,m2,d,[cnul,d*p2/(2*(cone + d)),               &
                          -(((4+3*d)*m2+(2+d)*p2)/(2*(1+d))), &
                          -3*(m2+p2)])/p2**2
      c3 = + (cone+m2**2/(p2**2)+2*m2/p2)*C0d_m00(p2,m2,d,2)
      C_P12P12 = c1+c2+c3
    end if

  end function

  function C_m00_p1p1_DP(p2,m2,d) result(C_p1p1)
    use triangle_aux_DP, only: A0mB0
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1,c1,c2,c3

    if (abs(p2/m2) .gt. large_z_threshold) then
      c1 = A0(m2)/(2*(cone + d))
      c2 = B0d(p2,m2,d,[-m2/(2*(cone+d)),(m2-3*p2)/2])
      c3 = C0d_m00(p2,m2,d,0)
      C_p1p1 = (c1+c2)/p2**2 + c3
    else
      c1 = A0mB0(p2,m2)
      c2 = B0d(p2,m2,d,[cnul,cone*p2*(-4-3*d),d*(m2+p2)])
      c3 = C0d_m00(p2,m2,d,0)
      C_p1p1 = (c1+c2)/(2*p2**2*(cone + d)) + c3
    end if
  end function

  function C_m00_p1P12_DP(p2,m2,d) result(C_p1P12)
    use triangle_aux_DP, only: A0mB0
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12,c1,c2,c3

    if (abs(p2/m2) .gt. large_z_threshold) then
      c1 = -A0(m2)/(cone + d)
      c2 = B0d(p2,m2,d,[m2/(cone+d),-(m2+2*p2),-(m2+5*p2)])
      c3 = C0d_m00(p2,m2,d,1)*(m2/p2+cone)
      C_p1P12 = (c1+c2)/(2*p2**2) + c3
    else
      c1 = -A0mB0(p2,m2)/(cone + d)
      c2 = B0d(p2,m2,d,[cnul,-p2*(cone+2*d)/(cone+d), &
                        -(m2+5*p2)-d*(m2+p2)/(cone+d)])
      c3 = C0d_m00(p2,m2,d,1)*(m2/p2+cone)
      C_p1P12 = (c1+c2)/(2*p2**2) + c3
    end if

  end function

  function C_m00_g_DP(p2,m2,d) result(C_g)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_g

    C_g = (cone+B0d(p2,m2,d,[cone,cone+m2/p2]))/4

  end function

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_m00_p1p1p1_DP(p2,m2,d) result(C_p1p1p1)
    use triangle_aux_DP, only: A0mB0_p1p1p1
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1p1,a,b0,ab,c,r

    if (abs(p2/m2) .gt. large_z_threshold) then
      ab = m2/(3*p2**3*(1+d))*(cone/(1+d) + cone-5*p2/(2*m2))
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,(2*m2**2-5*m2*p2+11*p2**2)/(6*p2**3)])
      c = -C0d_m00(p2,m2,d)
      r = -m2/(6*p2**2*(1+d))
      C_p1p1p1 = a + b0 + c + r
    else
      ab = m2/(3*p2**3*(1+d))*(cone/(1+d)+cone-5*p2/(2*m2))
      a = ab*A0mB0_p1p1p1(p2,m2,d)
      b0 = (-ab*d**3*m2 + &
            d**2*(2*m2**2-5*m2*p2+11*p2**2)/(6*p2**3))*B0d(p2,m2,d,3)
      c =  -C0d_m00(p2,m2,d)
      C_p1p1p1 = a + b0 + c
    end if

  end function

  function C_m00_p1p1P12_DP(p2,m2,d) result(C_p1p1P12)
    use triangle_aux_DP, only: A0mB0_p1p1P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1P12,a,b0,ab,c,r

    ab = m2/(3*p2**3*(1+d))*((5*p2/m2-cone)/2 - cone/(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,                            &
                        (6*p2**2+5*m2*p2-m2**2)/(6*p2**3), &
                        (17*p2**2+4*m2*p2-m2**2)/(6*p2**3)])
      c =  -(m2/p2+cone)*C0d_m00(p2,m2,d,1)
      r = m2/(6*p2**2*(1+d))
      C_p1p1P12 = a + b0 + c + r
    else
      a = ab*A0mB0_p1p1P12(p2,m2,d)
      b0 = (-ab*d**3*m2 + d**2*(6*p2**2+5*m2*p2-m2**2)/(6*p2**3) + &
            d*(17*p2**2+4*m2*p2-m2**2)/(6*p2**3))*B0d(p2,m2,d,3)
      c =  -(m2/p2+cone)*C0d_m00(p2,m2,d,1)
      C_p1p1P12 = a + b0 + c
    end if

  end function

  function C_m00_p1P12P12_DP(p2,m2,d) result(C_p1P12P12)
    use triangle_aux_DP, only: A0mB0_p1P12P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12P12,a,b0,ab,c,r

    ab = (m2/(1+d) - 5*p2/2)/(3*p2**3*(1+d))

    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d, [-ab*m2,                            &
                         -(5*m2+3*p2)/(6*p2**2),            &
                         (m2**2+7*p2*m2+3*p2**2)/(3*p2**3), &
                         (m2**2+11*p2*m2+10*p2**2)/(3*p2**3)])
      c = -(m2**2+2*m2*p2+p2**2)/p2**2*C0d_m00(p2,m2,d,2)
      r = -m2/(6*p2**2*(1+d))
      C_p1P12P12 = a + b0 + c + r
    else
      a = ab*A0mB0_p1P12P12(p2,m2,d)
      b0 = (-ab*d**3*m2                          &
            -d**2*(5*m2+3*p2)/(6*p2**2)          &
            +d*(m2**2+7*p2*m2+3*p2**2)/(3*p2**3) &
            +(m2**2+11*p2*m2+10*p2**2)/(3*p2**3))*B0d(p2,m2,d,3)
      c = -(m2**2+2*m2*p2+p2**2)/p2**2*C0d_m00(p2,m2,d,2)
      C_p1P12P12 = a + b0 + c
    end if

  end function

  function C_m00_P12P12P12_DP(p2,m2,d) result(C_P12P12P12)
    use triangle_aux_DP, only: A0mB0_P12P12P12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12P12,a,b0,ab,c,r

    ab = ((m2+5*p2)/2-m2/(1+d))/(3*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,                                &
                        (m2**2+5*m2*p2+2*p2**2)/(6*p2**3),     &
                        (-3*m2**2-10*m2*p2-3*p2**2)/(6*p2**3), &
                        (8*m2**2+11*m2*p2+3*p2**2)/(3*p2**3),  &
                        (11*(m2**2+2*m2*p2+p2**2))/(3*p2**3)])
      c = (-m2**3-3*m2**2*p2-3*m2*p2**2-p2**3)/(p2**3)*C0d_m00(p2,m2,d,3)
      r = m2/(6*p2**2*(1+d))
      C_P12P12P12 = a + b0 + c + r
    else
      a = ab*A0mB0_P12P12P12(p2,m2,d)
      b0 = B0d(p2,m2,d,[cnul,cnul,cnul,-d**3*ab*m2                      &
                        + d**2*(m2**2+5*m2*p2+2*p2**2)/(6*p2**3)  &
                        + d*(-3*m2**2-10*m2*p2-3*p2**2)/(6*p2**3) &
                        + (8*m2**2+11*m2*p2+3*p2**2)/(3*p2**3),   &
                        (11*(m2**2+2*m2*p2+p2**2))/(3*p2**3)])
      c = (-m2**3-3*m2**2*p2-3*m2*p2**2-p2**3)/(p2**3)*C0d_m00(p2,m2,d,3)
      C_P12P12P12 = a + b0 + c
    end if

  end function

  function C_m00_gP12_DP(p2,m2,d) result(C_gP12)
    use triangle_aux_DP, only: A0mB0_gP12
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gP12,a,b0,ab,r

    ab = -m2/(12*p2**2*(1+d))
    r = cone/18
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d, [(cone/12-ab*m2),           &
                         (-m2**2-m2*p2)/(12*p2**2), &
                         (-m2**2-2*m2*p2-p2**2)/(12*p2**2)])
    else
      a = ab*A0mB0_gP12(p2,m2,d)
      b0 = cone*B0d(p2,m2,d)/12 + &
           (-d**3*ab*m2 + d**2*(-m2**2-m2*p2)/(12*p2**2) + &
            d*(-m2**2-2*m2*p2-p2**2)/(12*p2**2))*B0d(p2,m2,d,3)
    end if
    C_gP12 = a + b0 + r
  end function

  function C_m00_gp1_DP(p2,m2,d) result(C_gp1)
    use triangle_aux_DP, only: A0mB0_gp1
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gp1,a,b0,ab,r

    ab = m2/(12*p2**2*(1+d))
    r = -7._dp/36
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-(cone/6+ab*m2),(m2**2-m2*p2-2*p2**2)/(12*p2**2)])
    else
      a = ab*A0mB0_gp1(p2,m2,d)
      b0 = B0d(p2,m2,d,[-cone/6,cnul,d**2*(m2**2-m2*p2-2*p2**2)/(12*p2**2),&
                        - d**3*ab*m2])
    end if
    C_gp1 = a + b0 + r
  end function

! --------------------------------------------------------
!               massless triangle
! --------------------------------------------------------

! ---------- rank 1 ----------------
  function C_000_p1_DP(p2,d) result(C_p1)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_p1

    C_p1 = B0d_0(p2,d,1)/p2 - C0d_000(p2,d,0)
  end function

  function C_000_P12_DP(p2,d) result(C_P12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_P12

    C_P12 = B0d_0(p2,d,[cnul,cone,2*cone])/p2 - C0d_000(p2,d,1)
  end function

! ---------- rank 2 ----------------
  function C_000_p1p1_DP(p2,d) result(C_p1p1)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_p1p1

    C_p1p1 = -3*B0d_0(p2,d,1)/(2*p2) + C0d_000(p2,d)
  end function

  function C_000_P12P12_DP(p2,d) result(C_P12P12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_P12P12

    C_P12P12 = B0d_0(p2,d,[cnul,cone/2,-cone,-3*cone])/p2 &
               + C0d_000(p2,d,2)
  end function

  function C_000_p1P12_DP(p2,d) result(C_p1P12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_p1P12

    C_p1P12 = B0d_0(p2,d,[cnul,-cone, -5*cone/2])/p2 + C0d_000(p2,d,1)
  end function

  function C_000_g_DP(p2,d) result(C_g)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_g

    C_g = cone/4 + B0d_0(p2,d,[cone/4,cone/4])
  end function

! ---------- rank 3 ----------------
  function C_000_p1p1p1_DP(p2,d) result(C_p1p1p1)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_p1p1p1

    C_p1p1p1 = 11*B0d_0(p2,d,1)/(6*p2) - C0d_000(p2,d)
  end function

  function C_000_p1p1P12_DP(p2,d) result(C_p1p1P12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_p1p1P12

    C_p1p1P12 = B0d_0(p2,d,[cnul, cone,17*cone/6])/p2 &
                - C0d_000(p2,d,1)
  end function

  function C_000_p1P12P12_DP(p2,d) result(C_p1P12P12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_p1P12P12

    C_p1P12P12 = B0d_0(p2,d,[cnul, -cone/2, cone, 10*cone/3])/p2 &
                 - C0d_000(p2,d,2)
  end function

  function C_000_P12P12P12_DP(p2,d) result(C_P12P12P12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_P12P12P12

    C_P12P12P12 = B0d_0(p2,d,[cnul, cone/3, -cone/2, cone, 11*cone/3])/p2 &
                  - C0d_000(p2,d,3)
  end function

  function C_000_gp1_DP(p2,d) result(C_gp1)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_gp1

    C_gp1 = -7._dp/36 + B0d_0(p2,d,[-cone/6, -cone/6])
  end function

  function C_000_gP12_DP(p2,d) result(C_gP12)
    real(dp), intent(in) :: p2,d
    complex(dp) :: C_gP12

    C_gP12 = B0d_0(p2,d,[cone/12, cnul, -cone/12]) + cone/18
  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     MMM                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_mmm_p1_DP(p2,m2,d) result(C_p1)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1

    C_p1 = B0d_mm(p2,m2,d,1)/p2 - C0d_mmm(p2,m2,d)

  end function

  function C_mmm_P12_DP(p2,m2,d) result(C_P12)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12, b

    b = B0d_mm(p2,m2,d,[cnul, cone/p2, 2*cone/p2])

    ! C_P12 = B0d_mm(p2,m2,muUV2,d,1)/p2 - C0d_mmm(p2,m2,d,1) + 2*B0d_mm(p2,m2,muUV2,d,2)/p2
    C_P12 = b - C0d_mmm(p2,m2,d,1)

  end function

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_mmm_p1p1_DP(p2,m2,d) result(C_p1p1)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1

    C_p1p1 = - 3._dp*B0d_mm(p2,m2,d,1)/(2._dp*p2) + C0d_mmm(p2,m2,d)

  end function

  function C_mmm_p1P12_DP(p2,m2,d) result(C_p1P12)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12,b

    b = B0d_mm(p2,m2,d,[cnul, -cone/p2, - 5*cone/(2*p2)])

    ! C_p1P12 = - B0d_mm(p2,m2,muUV2,d,1)/p2 + (cone - m2/p2)*C0d_mmm(p2,m2,d,1) &
    !           - 5._dp/(2*p2)*B0d_mm(p2,m2,muUV2,d,2)
    C_p1P12 = b + (cone - m2/p2)*C0d_mmm(p2,m2,d,1)

  end function

  function C_mmm_P12P12_DP(p2,m2,d) result(C_P12P12)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12, b

    b = B0d_mm(p2,m2,d,[cnul, cone/(2*p2), -cone/p2, -3*cone/p2])

    ! C_P12P12 = - 3._dp*B0d_mm(p2,m2,muUV2,d,3)/p2 &
    !            - B0d_mm(p2,m2,muUV2,d,2)/p2 + (cone - 2*m2/p2)*C0d_mmm(p2,m2,d,2) &
    !            + B0d_mm(p2,m2,muUV2,d,1)/(2*p2)
    C_P12P12 = b + (cone - 2*m2/p2)*C0d_mmm(p2,m2,d,2)

  end function

  function C_mmm_g_DP(p2,m2,d) result(C_g)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_g, b1, b0, c0, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    b = B0d_mm(p2,m2,d,[cone/4._dp, cone/4._dp])
    c0 = C0d_mmm(p2,m2,d,0)

    ! C_g = b1/4._dp + (b0+2._dp*m2*c0+cone)/4._dp
    C_g = b + (2._dp*m2*c0+cone)/4._dp

  end function

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_mmm_p1p1p1_DP(p2,m2,d) result(C_p1p1p1)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1p1, b1, b0, c0, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    b = B0d_mm(p2,m2,d,[ -m2/(3*(1+d)*p2**2), (2*m2 + 11*p2)/(6*p2**2)])
    c0 = C0d_mmm(p2,m2,d,0)

    ! C_p1p1p1 = (2*m2 + 11*p2)/(6*p2**2)*b1 &
    !            +(- m2*b0 + A0(m2,muUV2) -m2)/(3*(1+d)*p2**2) - c0
    C_p1p1p1 = b + ( A0(m2) -m2)/(3*(1+d)*p2**2) - c0

  end function

  function C_mmm_p1p1P12_DP(p2,m2,d) result(C_p1p1P12)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1P12, b3,b2,b1, b0, c1, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    b = B0d_mm(p2,m2,d,[m2/(3*(1+d)*p2**2), (-m2 + 3*p2)/(3*p2**2), (-4*m2+17*p2)/(6*p2**2) ])
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    c1 = C0d_mmm(p2,m2,d,1)

    ! C_p1p1P12 = (-4*m2+17*p2)/(6*p2**2)*b2 &
    !            + (-m2 + 3*p2)/(3*p2**2)*b1 &
    !            + (2*m2/p2 - cone)*c1 &
    !            + (m2*b0 - A0(m2,muUV2) + m2)/(3*(1+d)*p2**2)
    C_p1p1P12 = b + (2*m2/p2 - cone)*c1 + (- A0(m2) + m2)/(3*(1+d)*p2**2)


  end function

  function C_mmm_p1P12P12_DP(p2,m2,d) result(C_p1P12P12)
  !!!!!!!!!!!!!!!!!!!!
  !! for this function, using B0_mm_opt_list will cause meld_DP numerical instability
  !!!!!!!!!!!!!!!!!!!!
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12P12, b3,b2,b1, b0, c2, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    b = B0d_mm(p2,m2,d,[-m2/(3*(1+d)*p2**2), (2*m2 - 3*p2)/(6*p2**2), (-4*m2+3*p2)/(3*p2**2), 2*(-4*m2+5*p2)/(3*p2**2) ])
    c2 = C0d_mmm(p2,m2,d,2)

    ! C_p1P12P12 = 2*(-4*m2+5*p2)/(3*p2**2)*b3 &
    !            + (-4*m2+3*p2)/(3*p2**2)*b2 &
    !            + (4*m2/p2-cone)*c2 &
    !            + (2*m2 - 3*p2)/(6*p2**2)*b1 &
    !            + (-m2*b0 + A0(m2,muUV2) -m2)/(3*(1+d)*p2**2)
    C_p1P12P12 = b + (4*m2/p2-cone)*c2 + (A0(m2) -m2)/(3*(1+d)*p2**2)


  end function

  function C_mmm_P12P12P12_DP(p2,m2,d) result(C_P12P12P12)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12P12, b4,b3,b2,b1, b0, c3, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    ! b4 = B0d_mm(p2,m2,muUV2,d,4)
    b = B0d_mm(p2,m2,d,[ m2/(3*(1+d)*p2**2), (-m2+p2)/(3*p2**2), (4*m2-3*p2)/(6*p2**2), &
                              (-8*m2+3*p2)/(3*p2**2), (-16*m2+11*p2)/(3*p2**2) ])
    c3 = C0d_mmm(p2,m2,d,3)

    ! C_P12P12P12 = (-16*m2+11*p2)/(3*p2**2)*b4 &
    !             + (-8*m2+3*p2)/(3*p2**2)*b3 &
    !             + (18*m2-3*p2)/(3*p2)*c3 &
    !             + (4*m2-3*p2)/(6*p2**2)*b2 &
    !             + (-m2+p2)/(3*p2**2)*b1 &
    !             + (m2*b0 - A0(m2,muUV2)+m2)/(3*(1+d)*p2**2)

    C_P12P12P12 = b + (18*m2-3*p2)/(3*p2)*c3 + (- A0(m2)+m2)/(3*(1+d)*p2**2)

  end function

  function C_mmm_gp1_DP(p2,m2,d) result(C_gp1)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gp1, b4,b3,b2,b1, b0, c0, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    b = B0d_mm(p2,m2,d,[ -cone/6._dp, (2*m2-p2)/(6*p2) ])
    c0 = C0d_mmm(p2,m2,d,0)

    C_gp1 = b + (-18*m2*c0-7._dp)/36._dp

  end function

    function C_mmm_gP12_DP(p2,m2,d) result(C_gP12)
    complex(dp),     intent(in) :: m2
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gP12, b4,b3,b2,b1, b0, c1, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    b = B0d_mm(p2,m2,d,[cone/12._dp, 4*m2/(12*p2), (8*m2-p2)/(12*p2)])
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    ! b4 = B0d_mm(p2,m2,muUV2,d,4)
    c1 = C0d_mmm(p2,m2,d,1)

    C_gP12 = b +(-m2*p2*c1)/(2*p2) + cone/18._dp

  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   m0m1m1                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_m0m1m1_p1_DP(p2,m02,m12,d) result(C_p1)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1

    C_p1 = B0d_m0m1(p2,m02,m12,d,1)/p2 - C0d_m0m1m1(p2,m02,m12,d)

  end function

  function C_m0m1m1_P12_DP(p2,m02,m12,d) result(C_P12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12

    C_P12 = B0d_m0m1(p2,m02,m12,d,[cnul, cone/p2, 2*cone/p2])  &
          - (cone + (m02-m12)/p2)*C0d_m0m1m1(p2,m02,m12,d,1)

  end function

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_m0m1m1_p1p1_DP(p2,m02,m12,d) result(C_p1p1)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1,z0,z1,b1,a1

    z0 = m02/p2
    z1 = m12/p2
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [(z1-z0)/(2*p2*(cone+d)),-(3*cone-z0+z1)/(2*p2)])
    a1 = (A0(m02)-A0(m12))/(2*p2**2*(cone+d))
    C_p1p1 =  a1 + b1 + C0d_m0m1m1(p2,m02,m12,d)

  end function

  function C_m0m1m1_p1P12_DP(p2,m02,m12,d) result(C_p1P12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12,z0,z1,a1,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    a1 = -(A0(m02)-A0(m12))/(2*p2**2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d,       &
                  [(z0-z1)/(2*p2*(cone+d)), &
                  -((2*cone+z0-z1)/(2*p2)), &
                  -((5*cone+z0-z1)/(2*p2))])
    c1 = (cone+z0-2*z1)*C0d_m0m1m1(p2,m02,m12,d,1)
    C_p1P12 = a1 + b1 + c1
  end function

  function C_m0m1m1_P12P12_DP(p2,m02,m12,d) result(C_P12P12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12,z0,z1,a1,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    a1 = (A0(m02)-A0(m12))/(2*p2**2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d,        &
                  [(z1-z0)/(2*p2*(1+d)),   &
                   (cone+z0-z1)/(2*p2),    &
                   (-cone-2*z0+2*z1)/(p2), &
                   -(3*(cone+z0-z1))/(p2)])

    c1 = ((cone + z0)**2 - 2*(2*cone+z0)*z1+z1**2)*C0d_m0m1m1(p2,m02,m12,d,2)
    C_P12P12 = a1 + b1 + c1
  end function

  function C_m0m1m1_g_DP(p2,m02,m12,d) result(C_g)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_g,z0,z1,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    b1 = B0d_m0m1(p2,m02,m12,d,[cone,cone+z0-z1])/4
    c1 = p2*z1*C0d_m0m1m1(p2,m02,m12,d,0)/2
    C_g = cone/4 + b1 + c1
  end function

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_m0m1m1_p1p1p1_DP(p2,m02,m12,d) result(C_p1p1p1)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1p1,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2
    a11 = ((z0-z1)/(cone+d) + (-5*cone+2*z0-2*z1)/2)*A0(m02)/(3*p2**2*(cone+d))
    a12 = -((z0-z1)/(cone+d) + (-7*cone+2*z0-2*z1)/2)*A0(m12)/(3*p2**2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d, &
                [(-z0**2+2*z0*z1-z1**2)/(3*p2*(cone+d)**2) + &
                (5*z0-2*z0**2-7*z1+4*z0*z1-2*z1**2)/(6*p2*(cone+d)), &
                (11*cone-5*z0+2*z0**2+7*z1-4*z0*z1+2*z1**2)/(6*p2)])
    c1 = -C0d_m0m1m1(p2,m02,m12,d)
    C_p1p1p1 =  a11 + a12 + b1 + c1 -(z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_p1p1P12_DP(p2,m02,m12,d) result(C_p1p1P12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1p1P12,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2
    a11 = ((z1-z0)/(cone+d) + (5*cone-z0+z1)/2)*A0(m02)/(3*p2**2*(cone+d))
    a12 = -((z1-z0)/(cone+d) + (7*cone-z0+z1)/2)*A0(m12)/(3*p2**2*(cone+d))

    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [(z0-z1)**2/(3*p2*(cone+d)**2) + &
                  (-5*z0+z0**2+7*z1-2*z0*z1+z1**2) /(6*p2*(cone+d)), &
                   -(-6*cone+z0**2+z1*(7*cone+z1)-z0*(5*cone+2*z1))/(6*p2), &
                   -(-17*cone+z0**2-2*z0*(2*cone+z1)+z1*(8*cone+z1))/(6*p2)])

    c1 = -(cone+z0-3*z1)*C0d_m0m1m1(p2,m02,m12,d,1)
    C_p1p1P12 = a11 + a12 + b1 + c1 + (z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_p1P12P12_DP(p2,m02,m12,d) result(C_p1P12P12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_p1P12P12,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    a11 = ((z0-z1)/(3*p2**2*(cone+d)**2)-5/(6*p2**2*(cone+d)))*A0(m02)
    a12 = (-(z0-z1)/(3*p2**2*(cone+d)**2)+7/(6*p2**2*(cone+d)))*A0(m12)

    b1 = B0d_m0m1(p2,m02,m12,d, &
          [-((z0 - z1)**2/(3*p2*(cone+d)**2))+(5*z0-7*z1)/( 6*p2*(cone+d)), &
           (-3*cone-5*z0+7*z1)/(6*p2), &
           (3*cone+7*z0+z0**2-11*z1-2*z0*z1+z1**2)/(3*p2), &
           (10*cone+11*z0+z0**2-19*z1-2*z0*z1+z1**2)/(3*p2)])

    c1 = (-cone-2*z0-z0**2+6*z1+4*z0*z1-3*z1**2)*C0d_m0m1m1(p2,m02,m12,d,2)
    C_p1P12P12 = a11 + a12 + b1 + c1 - (z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_P12P12P12_DP(p2,m02,m12,d) result(C_P12P12P12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_P12P12P12,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2
    a11 = ((-z0+z1)/(3*p2**2*(cone+d)**2) + (5*cone+z0-z1)/(6*p2**2*(cone+d)))*A0(m02)
    a12 = -((-z0+z1)/(3*p2**2*(cone+d)**2) + (7*cone+z0-z1)/(6*p2**2*(cone+d)))*A0(m12)
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [(z0-z1)**2/(3*p2*(cone+d)**2)-(z0**2+z0*(5*cone-2*z1)+ &
                   (-7*cone+z1)*z1)/(6*p2*(cone+d)), &
                   (2*cone+z0**2+z0*(5*cone-2*z1)+(-7*cone+z1)*z1)/(6*p2), &
                   -((3*cone+3*z0**2+z0*(10*cone-6*z1)+z1*(-14*cone+3*z1))/(6*p2)), &
                   (3*cone+8*z0**2+z0*(11*cone-16*z1)+z1*(-19*cone+8*z1))/(3*p2), &
                   (11*cone-38*z1+11*(z0*(2*cone+z0)-2*z0*z1+z1**2))/(3*p2)])

    c1 = -(cone+z0-z1)*((cone+z0)**2-2*(4*cone+z0)*z1+z1**2)*C0d_m0m1m1(p2,m02,m12,d,3)
    C_P12P12P12 = a11 + a12 + b1 + c1 + (z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_gp1_DP(p2,m02,m12,d) result(C_gp1)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gp1,z0,z1,a11,a12,b1,b2,c1,r1

    z0 = m02/p2
    z1 = m12/p2
    a11 = (z0-z1)*A0(m02)/(12*p2*(cone+d))
    a12 = -(z0-z1)*A0(m12)/(12*p2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [-(z0-z1)**2/(12*(cone+d))-cone/6, &
                   (-2*cone-z0+z0**2+5*z1-2*z0*z1+z1**2)/(12)])
    c1 = -(p2*z1)*C0d_m0m1m1(p2,m02,m12,d)/2
    r1 = -(7*cone)/36
    C_gp1 = a11 + a12 + b1 + b2 + c1 + r1

  end function

  function C_m0m1m1_gP12_DP(p2,m02,m12,d) result(C_gP12)
    complex(dp),     intent(in) :: m02,m12
    real(dp),        intent(in) :: p2,d
    complex(dp) :: C_gP12,z0,z1,a1,b1,c1,r1

    z0 = m02/p2
    z1 = m12/p2
    a1 = (-z0+z1)/(12*p2*(1+d))*(A0(m02)-A0(m12))
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [cone/12 + &
                   (z0**2-2*z0*z1+z1**2)/(12*(cone+d)), &
                  -(z0+z0**2-2*z0*z1+(-5*cone+z1)*z1)/12, &
                  -((cone+z0)**2-2*(5*cone+z0)*z1+z1**2)/12])

    c1 = p2*z1*(-cone-z0+z1)*C0d_m0m1m1(p2,m02,m12,d,1)/2
    C_gP12 = a1 + b1 + c1 + cone/18

  end function

end module triangle_reduction_DP

module triangle_reduction_QP
  use triangle_expansion_QP, only: B0d_coeff,B0d_list,B0d_matrix,         &
                                B0d_0_coeff,B0d_0_list,                   &
                                C0d_0mm_coeff,C0d_0mm_matrix,             &
                                C0d_m00_coeff,C0d_000_coeff,              &
                                C0d_000_EP1_coeff,                        &
                                C0d_mmm_table,B0d_mm_opt,B0d_mm_opt_list, &
                                C0d_m0m1m1_coeff,C0d_m0m1m1_list,         &
                                B0d_m0m1_coeff,B0d_m0m1_list
  use triangle_aux_QP, only: qp,cone,cnul,A0,B0_zero,rnul
  implicit none

  real(qp), parameter :: large_z_threshold = 0.1_qp

  interface B0d
    module procedure B0d_coeff, B0d_list, B0d_matrix
  end interface

  interface B0d_0
    module procedure B0d_0_coeff, B0d_0_list
  end interface

  interface B0d_mm
    module procedure B0d_mm_opt, B0d_mm_opt_list
  end interface

  interface B0d_m0m1
    module procedure B0d_m0m1_coeff, B0d_m0m1_list
  end interface

  interface C0d_0mm
    module procedure C0d_0mm_coeff,C0d_0mm_matrix
  end interface

  interface C0d_m00
    module procedure C0d_m00_coeff
  end interface

  interface C0d_mmm
    module procedure C0d_mmm_table
  end interface

  interface C0d_000
    module procedure C0d_000_coeff
  end interface

  interface C0d_000_EP1
    module procedure C0d_000_EP1_coeff
  end interface

  interface C0d_m0m1m1
    module procedure C0d_m0m1m1_coeff, C0d_m0m1m1_list
  end interface

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     0mm                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_0mm_P12_QP(p2,m2,d,expansion) result(C_P12)
    ! computes the coefficient of the rank 1 p1[mu]-p2[mu] form factor
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    logical, optional, intent(in) :: expansion
    complex(qp) :: C_P12

    C_P12 = B0d(p2,m2,d,[cnul,cone,2*cone])/p2  + (m2/p2-cone)*C0d_0mm(p2,m2,d,1)

  end function

  function C_0mm_p1_QP(p2,m2,d) result(C_p1)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1

    C_p1 = B0d(p2,m2,d,1)/p2 - C0d_0mm(p2,m2,d)

  end function

  subroutine C_0mm_rank1(p2,m2,d,C_P12,C_p1)
    complex(qp),     intent(in)  :: m2
    real(qp),        intent(in)  :: p2,d
    complex(qp),     intent(out) :: C_P12,C_p1
    complex(qp) :: C(2)

    C = B0d_matrix(p2,m2,d,2,3,[[cnul,cnul], [cone/p2, cone/p2], [2*cone/p2,cnul]],[.true.,.false.,.false.]) &
        + C0d_0mm_matrix(p2,m2,d,2,2,[[cnul,-cone],[(m2/p2-cone),cnul]],[.false.,.false.])

    C_P12 = C(1)
    C_p1 = C(2)

  end subroutine C_0mm_rank1

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_0mm_p1p1_QP(p2,m2,d) result(C_p1p1)
    use triangle_aux_QP, only: A0mB0_0mm_p1p1
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1,b0,ab,a,c0

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_p1p1 = + C0d_0mm(p2,m2,d,0)                      &
               - (A0(m2)/(2*p2**2*(1 + d)))        &
               + B0d(p2,m2,d,[+m2/(2*p2**2*(1+d)), &
                                    -m2/(2*p2**2)-3/(2*p2)])
    else
      ab = -cone/(2*p2**2*(1 + d))
      a = ab*A0mB0_0mm_p1p1(p2,m2,d)
      b0 = -d*(m2*ab*d+m2/(2*p2**2)+3/(2*p2))*B0d(p2,m2,d,2)
      c0 = C0d_0mm(p2,m2,d,0)
      C_p1p1 = a + b0 + c0
    end if
  end function

  function C_0mm_p1P12_QP(p2,m2,d) result(C_p1P12)
    use triangle_aux_QP, only: A0mB0_0mm_p1P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12,ab,a,b0,c0

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_p1P12 = + A0(m2)/(2*p2**2*(1 + d))                       &
                + B0d(p2,m2,d,[-m2/(2*p2**2*(1+d)),              &
                                          +(m2/(2*p2**2)-1/p2),        &
                                          +(m2/(2*p2**2) - 5/(2*p2))]) &
                + (1-(2*m2)/(p2))*C0d_0mm(p2,m2,d,1)
    else
      ab = cone/(2*p2**2*(1 + d))
      a = ab*A0mB0_0mm_p1P12(p2,m2,d)
      b0 = (m2+2*d*m2-(1+d)*(5+2*d)*p2)/(2*(1+d)*p2**2)*B0d(p2,m2,d,2)
      c0 = (1 - (2*m2)/(p2))*C0d_0mm(p2,m2,d,1)
      C_p1P12 = a + b0 + c0
    end if
  end function

  function C_0mm_P12P12_QP(p2,m2,d) result(C_P12P12)
    use triangle_aux_QP, only: A0mB0_0mm_P12P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12,ab,a,b0,c0

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_P12P12 = - A0(m2)/(2*p2**2*(1 + d))             &
                 + B0d(p2,m2,d, [m2/(2*p2**2*(1 + d)),  &
                                 (1/(2*p2)- m2/(2*p2**2)), &
                                 (2*m2/p2**2 - 1/p2),      &
                                 (3*m2/p2**2 - 3/p2)])     &
                 + (1 + m2**2/(p2**2) - 4*m2/p2)*C0d_0mm(p2,m2,d,2)
    else
      ab = -cone/(2*p2**2*(1 + d))
      a = ab*A0mB0_0mm_P12P12(p2,m2,d)
      b0 = B0d(p2,m2,d,[cnul,cnul,-ab*((4+3*d)*m2+(-2+d)*(1+d)*p2), &
                        (3*m2/p2**2 - 3/p2)])
      c0 = (1 + m2**2/(p2**2) - 4*m2/p2)*C0d_0mm(p2,m2,d,2)
      C_P12P12 = a + b0 + c0
    end if
  end function

  function C_0mm_g_QP(p2,m2,d) result(C_gMunu)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gMunu

    C_gMunu = B0d(p2,m2,d,[cone/4,(p2-m2)/(4*p2)]) &
              + (2*m2*C0d_0mm(p2,m2,d) + 1)/4
  end function

  subroutine C_0mm_rank2(p2,m2,d,C_P12P12,C_p1p1,C_p1P12,C_g)
    complex(qp),     intent(in)  :: m2
    real(qp),        intent(in)  :: p2,d
    complex(qp),     intent(out) :: C_P12P12,C_p1p1,C_p1P12,C_g
    complex(qp) :: C(4),a1

    C = B0d_matrix(p2,m2,d,4,4, &
    [[m2/(2*p2**2*(1+d)),       m2/(2*p2**2*(1+d)),    -m2/(2*p2**2*(1+d)), cone/4], &
     [(1/(2*p2)- m2/(2*p2**2)),-m2/(2*p2**2)-3/(2*p2), +(m2/(2*p2**2)-1/p2),(p2-m2)/(4*p2)], &
     [(2*m2/p2**2-1/p2),        cnul,                  +(m2/(2*p2**2)-5/(2*p2)), cnul], &
     [(3*m2/p2**2 - 3/p2), cnul, cnul,cnul]],[.false.,.false.,.false.,.false.]) &
    + C0d_0mm_matrix(p2,m2,d,4,3, &
    [[cnul,cone,cnul,m2/2], &
     [cnul,cnul,(1 - (2*m2)/(p2)),cnul],&
     [(1 + m2**2/(p2**2) - 4*m2/p2),cnul,cnul,cnul]],[.false.,.false.,.false.])

    a1 = A0(m2)
    C = C + [-a1/(2*p2**2*(1 + d)),-a1/(2*p2**2*(1 + d)),a1/(2*p2**2*(1 + d)),cone/4]
    C_P12P12 = C(1)
    C_p1p1 = C(2)
    C_p1P12 = C(3)
    C_g = C(4)

  end subroutine C_0mm_rank2

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_0mm_p1p1p1_QP(p2,m2,d) result(C_p1p1p1)
    use triangle_aux_QP, only: A0mB0_0mm_p1p1p1
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1p1,ab,a,b0,c0,r

    ab = m2/(3*p2**3*(1+d)**2)+(2*m2+7*p2)/(6*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,(2*m2**2+7*m2*p2+11*p2**2)/(6*p2**3)])
      c0 = -C0d_0mm(p2,m2,d,0)
      r = - m2/(6*(1+d)*p2**2)
      C_p1p1p1 = a + b0 + c0 + r
    else
      a = ab*A0mB0_0mm_p1p1p1(p2,m2,d)
      b0 = (-d**3*m2*ab+ &
            d**2*(2*m2**2+7*m2*p2+11*p2**2)/(6*(p2**3)))*B0d(p2,m2,d,3)
      c0 = -C0d_0mm(p2,m2,d,0)
      C_p1p1p1 = a + b0 + c0
    end if
  end function

  function C_0mm_p1p1P12_QP(p2,m2,d) result(C_p1p1P12)
    use triangle_aux_QP, only: A0mB0_0mm_p1p1P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1P12,c1,b,a,ab,r

    ab = -(m2/(3*p2**3*(1+d)**2))+(-m2-7*p2)/(6*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b = B0d(p2,m2,d,[-ab*m2,                             &
                       (-m2**2-7*m2*p2+6*p2**2)/(6*p2**3), &
                       -(m2**2+8*m2*p2-17*p2**2)/(6*p2**3)])
      c1 = -(p2-3*m2)*C0d_0mm(p2,m2,d,1)/p2
      r = m2*p2/(6*(1+d)*p2**3)
      C_p1p1P12 =  a + b + c1 + r
    else
      a = ab*A0mB0_0mm_p1p1P12(p2,m2,d)
      c1 = -(p2-3*m2)/p2*C0d_0mm(p2,m2,d,1)
      b = (d**3*((m2**2)/(3*(1+d)**2*p2**3)             &
           +(m2**2 + 7*m2*p2)/(6*(1+d)*p2**3))          &
           +d**2*(-m2**2 - 7*m2*p2 + 6*p2**2)/(6*p2**3) &
           -d*(m2**2 + 8*m2*p2 - 17*p2**2)/(6*p2**3))*B0d(p2,m2,d,3)

      C_p1p1P12 = a + b + c1
    end if

  end function

  function C_0mm_p1P12P12_QP(p2,m2,d) result(C_p1P12P12)
    use triangle_aux_QP, only: A0mB0_0mm_p1P12P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12P12,b0,c2,ab,a,r

    if (abs(p2/m2) .gt. large_z_threshold) then
      ab = m2/(3*p2**3*(1+d)**2)+7/(6*p2**2*(1+d))
      a = ab*A0(m2)
      b0 = B0d_list(p2,m2,d,[-ab*m2,(7*m2 - 3*p2)/(6*p2**2),         &
                             (m2**2 - 11*m2*p2 + 3*p2**2)/(3*p2**3), &
                             (m2**2 - 19*m2*p2 + 10*p2**2)/(3*p2**3) ])

      r = -m2/(6*p2**2*(1+d))
      c2 = (-3*m2**2+6*m2*p2-p2**2)/(p2**2)*C0d_0mm(p2,m2,d,2)
      C_p1P12P12 = a + b0 + c2 + r
    else
      ab = m2/(3*p2**3*(1+d)**2)+7/(6*p2**2*(1+d))
      a = ab*A0mB0_0mm_p1P12P12(p2,m2,d)
      b0 = (-d**3*ab*m2 +                          &
            d**2*(7*m2-3*p2)/(6*p2**2) +           &
            d*(m2**2-11*m2*p2+3*p2**2)/(3*p2**3) + &
            (m2**2-19*m2*p2+10*p2**2)/(3*p2**3))*B0d(p2,m2,d,3)
      c2 = (-3*m2**2+6*m2*p2-p2**2)/(p2**2)*C0d_0mm(p2,m2,d,2)
      C_p1P12P12 = a + b0 + c2

    end if
  end function

  function C_0mm_P12P12P12_QP(p2,m2,d) result(C_P12P12P12)
    use triangle_aux_QP, only: A0mB0_0mm_P12P12P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12P12,b0,c3,ab,a,r

    ab = -(m2/(3*p2**3*(1+d)**2))+(m2-7*p2)/(6*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,                                    &
                        (m2**2 - 7*m2*p2 + 2*p2**2)/(6*p2**3),     &
                        (-3*m2**2 + 14*m2*p2 - 3*p2**2)/(6*p2**3), &
                        (8*m2**2 - 19*m2*p2 + 3*p2**2)/(3*p2**3),  &
                        (11*m2**2 - 38*m2*p2 + 11*p2**2)/(3*p2**3)])
      c3 = (3*m2**4 - 27*m2**3*p2 + &
            27*m2**2*p2**2 - 3*m2*p2**3)/(3*m2*p2**3)*C0d_0mm(p2,m2,d,3)
      r = (m2*p2)/(6*(1+d)*p2**3)
      C_P12P12P12 = a + b0 + c3 + r
    else
      a = ab*A0mB0_0mm_P12P12P12(p2,m2,d)
      b0 = B0d(p2,m2,d,[cnul,cnul,cnul,-d**3*ab*m2 + &
                        d**2*(m2**2-7*m2*p2+2*p2**2)/(6*p2**3) + &
                        d*(-3*m2**2+14*m2*p2-3*p2**2)/(6*p2**3) + &
                        (8*m2**2-19*m2*p2+3*p2**2)/(3*p2**3), &
                        (11*m2**2-38*m2*p2+11*p2**2)/(3*p2**3)])
      c3 = (m2**3-9*m2**2*p2+9*m2*p2**2-p2**3)/(p2**3)*C0d_0mm(p2,m2,d,3)
      C_P12P12P12 = a + b0 + c3
    end if
  end function

  function C_0mm_gp1_QP(p2,m2,d) result(C_gp1)
    use triangle_aux_QP, only: A0mB0_0mm_gp1
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gp1,b0,c0,b,ab,a

    if (abs(p2/m2) .gt. large_z_threshold) then
      b = B0d(p2,m2,d,[-cone/6 - m2**2/(12*(1+d)*p2**2), &
                       (m2**2 + 5*m2*p2 - 2*p2**2)/(12*p2**2)])
      c0 = C0d_0mm(p2,m2,d,0)
      C_gp1 = b +(- 18*m2*c0 - 7)/36                   &
                +(m2*A0(m2))/(12*(1+d)*p2**2)
    else
      ab = m2/(12*p2**2*(1+d))
      a = ab*A0mB0_0mm_gp1(p2,m2,d)
      b0 = B0d(p2,m2,d,[-cone/6,cnul, &
               -d**2*(m2*ab)+d*(m2**2+5*m2*p2-2*p2**2)/(12*p2**2)])
      c0 = -m2*C0d_0mm(p2,m2,d,0)/2
      C_gp1 = a + b0 + c0
    end if
  end function

  function C_0mm_gP12_QP(p2,m2,d) result(C_gP12)
    use triangle_aux_QP, only: A0mB0_0mm_gP12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gP12,b0,c1,b,ab,a

    if (abs(p2/m2) .gt. large_z_threshold) then
      b = B0d(p2,m2,d,[cone/12 + m2**2/(12*(1+d)*p2**2), &
                       (-m2**2 + 5*m2*p2)/(12*p2**2) , &
                       (-m2**2 + 10*m2*p2 - p2**2)/(12*p2**2)])
      c1 = C0d_0mm(p2,m2,d,1)
      C_gP12 = (6*m2**2*p2 - 6*m2*p2**2)/(12*p2**2)*c1 &
                  +(-m2*A0(m2))/(12*(1+d)*p2**2) &
                  +cone/18 + b
    else
      ab = -m2/(12*p2**2*(1+d))
      a = ab*A0mB0_0mm_gP12(p2,m2,d)
      b0 = +B0d(p2,m2,d,[cone/12,cnul,(-m2*ab*d**2) + &
                         d*(-m2**2+5*m2*p2)/(12*p2**2) + &
                         (-m2**2+10*m2*p2-p2**2)/(12*p2**2)])
      c1 = (m2**2 - m2*p2)/(2*p2)*C0d_0mm(p2,m2,d,1)
      C_gP12 = a + b0 + c1
    end if
  end function


  subroutine C_0mm_rank3(p2,m2,d,C_P12P12P12,C_p1P12P12,C_p1p1P12,C_p1p1p1,C_gP12,C_gp1)
    complex(qp),     intent(in)  :: m2
    real(qp),        intent(in)  :: p2,d
    complex(qp),     intent(out) :: C_P12P12P12,C_p1P12P12,C_p1p1P12,C_p1p1p1,C_gP12,C_gp1
    complex(qp) :: C(6),a1,b,c0

    C = B0d_matrix(p2,m2,d,6,5, &
    [[ m2**2/(3*p2**3*(1+d)**2)-m2*(m2-7*p2)/(6*p2**3*(1+d)),&
      -m2**2/(3*p2**3*(1+d)**2)-7*m2/(6*p2**2*(1+d)), &
       m2**2/(3*p2**3*(1+d)**2) + m2*(m2+7*p2)/(6*p2**3*(1+d)), &
      -m2**2/(3*p2**3*(1+d)**2) - m2*(2*m2+7*p2)/(6*p2**3*(1+d)), &
       cone/12+m2**2/(12*(1+d)*p2**2), &
      -cone/6-m2**2/(12*(1+d)*p2**2)], &
     [ (m2**2 - 7*m2*p2 + 2*p2**2)/(6*p2**3), &
       (7*m2-3*p2)/(6*p2**2), &
       (-m2**2-7*m2*p2+6*p2**2)/(6*p2**3), &
       (2*m2**2+7*m2*p2+11*p2**2)/(6*p2**3), &
       (-m2**2 + 5*m2*p2)/(12*p2**2), &
       (m2**2+5*m2*p2-2*p2**2)/(12*p2**2)], &
     [ (-3*m2**2 + 14*m2*p2 - 3*p2**2)/(6*p2**3), &
       (m2**2 - 11*m2*p2 + 3*p2**2)/(3*p2**3), &
      -(m2**2+8*m2*p2-17*p2**2)/(6*p2**3), &
       cnul,(-m2**2 + 10*m2*p2 - p2**2)/(12*p2**2),cnul], &
     [ (8*m2**2 - 19*m2*p2 + 3*p2**2)/(3*p2**3), &
       (m2**2 - 19*m2*p2 + 10*p2**2)/(3*p2**3), &
       cnul,cnul,cnul,cnul], &
     [(11*m2**2 - 38*m2*p2 + 11*p2**2)/(3*p2**3),cnul,cnul,cnul,cnul,cnul]],&
     [.false.,.false.,.false.,.false.,.false.])


    C = C + C0d_0mm_matrix(p2,m2,d,6,4, &
   [[cnul,cnul,cnul,-cone,cnul,-m2/2], &
    [cnul,cnul,-cone+3*m2/p2,cnul,(m2**2 - m2*p2)/(2*p2),cnul], &
    [cnul,(-3*m2**2+6*m2*p2-p2**2)/(p2**2),cnul,cnul,cnul,cnul], &
    [(m2**3-9*m2**2*p2+9*m2*p2**2-p2**3)/(p2**3),cnul,cnul,cnul,cnul,cnul]],&
   [.false.,.false.,.false.,.false.])

    a1 = A0(m2)
    C = C + [-(m2/(3*p2**3*(1+d)**2))+(m2-7*p2)/(6*p2**3*(1+d)),  &
               m2/(3*p2**3*(1+d)**2)+7/(6*p2**2*(1+d)),           &
             -(m2/(3*p2**3*(1+d)**2))+(-m2-7*p2)/(6*p2**3*(1+d)), &
               m2/(3*p2**3*(1+d)**2)+(2*m2+7*p2)/(6*p2**3*(1+d)), &
             -m2/(12*(1+d)*p2**2),                                &
             +m2/(12*(1+d)*p2**2)] * a1

    C = C + [(m2*p2)/(6*(1+d)*p2**3), &
             -m2/(6*p2**2*(1+d)), &
              m2*p2/(6*(1+d)*p2**3), &
            - m2/(6*(1+d)*p2**2), &
            + cone/18, &
            - 7*cone/36]

    C_P12P12P12 = C(1)
    C_p1P12P12 = C(2)
    C_p1p1P12 = C(3)
    C_p1p1p1 = C(4)
    C_gP12 = C(5)
    C_gp1 = C(6)

  end subroutine C_0mm_rank3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     m00                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_m00_p1_QP(p2,m2,d,expansion) result(C_p1)
    ! computes the coefficient of the rank 1 p1[mu] form factor
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    logical, optional, intent(in) :: expansion
    complex(qp) :: C_p1

    if (present(expansion) .and. .not. expansion) then
      ! tensor reduction result with no \delta expansion
      C_p1 = - C0d_m00(p2,m2,d)                 &
             + B0d(p2,m2,d,0)/(p2*d)            &
             - B0d(p2,m2,rnul,0)/(p2*d)
    else
      C_p1 = B0d(p2,m2,d,1)/p2 - C0d_m00(p2,m2,d)
    end if

  end function

  function C_m00_P12_QP(p2,m2,d,expansion) result(C_P12)
    ! computes the coefficient of the rank 1 p1[mu]-p2[mu] form factor
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    logical, optional, intent(in) :: expansion
    complex(qp) :: C_P12

    if (present(expansion) .and. .not. expansion) then
      ! tensor reduction result with no \delta expansion
      C_P12 = -(m2 + p2)/(p2*d)*C0d_m00(p2,m2,d,0) &
              + (2 + d)/(p2*d**2)*B0d(p2,m2,d,0)   &
              - B0_zero()/(p2*d)                   &
              - 2/(p2*d**2)*B0d(p2,m2,rnul,0)
    else
      C_P12 = B0d(p2,m2,d,[cnul,cone,2*cone])/p2 &
              + (-m2/p2-cone)*C0d_m00(p2,m2,d,1)
    end if

  end function

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_m00_P12P12_QP(p2,m2,d) result(C_P12P12)
    use triangle_aux_QP, only: A0mB0
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12,c1,c2,c3

    if (abs(p2/m2) .gt. large_z_threshold) then
      C_P12P12 = + A0(m2)/(2*(cone + d))          &
                 + B0d(p2,m2,d,[-m2/(2*(cone+d)), &
                                +(p2+m2)/2,       &
                                -(2*m2 + p2),     &
                                -3*(m2 + p2)])    &
                 + (p2**2 + m2**2 + 2*m2*p2)*C0d_m00(p2,m2,d,2)
      C_P12P12 = C_P12P12/p2**2
    else
      ! For small p2/m2 the formula above is numerically unstable due to large
      ! cancellations between A0 and B0, DB0.

      ! Observation:
      ! The numerical instability resides in A0 and the first two (truncated)
      ! B0 functions:
      ! A0 - m2 B0^0 + (1+d) (p2+m2) B0^1

      ! Used identities:
      ! -B0^0 + d B0^1 = B0^{0,(0)}    (only the first coefficients of B0^0)
      ! B0^1 = B0^{1,(0)} + d B0^2
      ! A0mB0 = A0(m2)-m2 B0(p2,0,m2) + (p2+m2) B0^{1,(0)}(p2,0,m2)
      c1 = + A0mB0(p2,m2)/(2*p2**2*(cone + d))
      c2 = + B0d(p2,m2,d,[cnul,d*p2/(2*(cone + d)),               &
                          -(((4+3*d)*m2+(2+d)*p2)/(2*(1+d))), &
                          -3*(m2+p2)])/p2**2
      c3 = + (cone+m2**2/(p2**2)+2*m2/p2)*C0d_m00(p2,m2,d,2)
      C_P12P12 = c1+c2+c3
    end if

  end function

  function C_m00_p1p1_QP(p2,m2,d) result(C_p1p1)
    use triangle_aux_QP, only: A0mB0
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1,c1,c2,c3

    if (abs(p2/m2) .gt. large_z_threshold) then
      c1 = A0(m2)/(2*(cone + d))
      c2 = B0d(p2,m2,d,[-m2/(2*(cone+d)),(m2-3*p2)/2])
      c3 = C0d_m00(p2,m2,d,0)
      C_p1p1 = (c1+c2)/p2**2 + c3
    else
      c1 = A0mB0(p2,m2)
      c2 = B0d(p2,m2,d,[cnul,cone*p2*(-4-3*d),d*(m2+p2)])
      c3 = C0d_m00(p2,m2,d,0)
      C_p1p1 = (c1+c2)/(2*p2**2*(cone + d)) + c3
    end if
  end function

  function C_m00_p1P12_QP(p2,m2,d) result(C_p1P12)
    use triangle_aux_QP, only: A0mB0
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12,c1,c2,c3

    if (abs(p2/m2) .gt. large_z_threshold) then
      c1 = -A0(m2)/(cone + d)
      c2 = B0d(p2,m2,d,[m2/(cone+d),-(m2+2*p2),-(m2+5*p2)])
      c3 = C0d_m00(p2,m2,d,1)*(m2/p2+cone)
      C_p1P12 = (c1+c2)/(2*p2**2) + c3
    else
      c1 = -A0mB0(p2,m2)/(cone + d)
      c2 = B0d(p2,m2,d,[cnul,-p2*(cone+2*d)/(cone+d), &
                        -(m2+5*p2)-d*(m2+p2)/(cone+d)])
      c3 = C0d_m00(p2,m2,d,1)*(m2/p2+cone)
      C_p1P12 = (c1+c2)/(2*p2**2) + c3
    end if

  end function

  function C_m00_g_QP(p2,m2,d) result(C_g)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_g

    C_g = (cone+B0d(p2,m2,d,[cone,cone+m2/p2]))/4

  end function

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_m00_p1p1p1_QP(p2,m2,d) result(C_p1p1p1)
    use triangle_aux_QP, only: A0mB0_p1p1p1
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1p1,a,b0,ab,c,r

    if (abs(p2/m2) .gt. large_z_threshold) then
      ab = m2/(3*p2**3*(1+d))*(cone/(1+d) + cone-5*p2/(2*m2))
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,(2*m2**2-5*m2*p2+11*p2**2)/(6*p2**3)])
      c = -C0d_m00(p2,m2,d)
      r = -m2/(6*p2**2*(1+d))
      C_p1p1p1 = a + b0 + c + r
    else
      ab = m2/(3*p2**3*(1+d))*(cone/(1+d)+cone-5*p2/(2*m2))
      a = ab*A0mB0_p1p1p1(p2,m2,d)
      b0 = (-ab*d**3*m2 + &
            d**2*(2*m2**2-5*m2*p2+11*p2**2)/(6*p2**3))*B0d(p2,m2,d,3)
      c =  -C0d_m00(p2,m2,d)
      C_p1p1p1 = a + b0 + c
    end if

  end function

  function C_m00_p1p1P12_QP(p2,m2,d) result(C_p1p1P12)
    use triangle_aux_QP, only: A0mB0_p1p1P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1P12,a,b0,ab,c,r

    ab = m2/(3*p2**3*(1+d))*((5*p2/m2-cone)/2 - cone/(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,                            &
                        (6*p2**2+5*m2*p2-m2**2)/(6*p2**3), &
                        (17*p2**2+4*m2*p2-m2**2)/(6*p2**3)])
      c =  -(m2/p2+cone)*C0d_m00(p2,m2,d,1)
      r = m2/(6*p2**2*(1+d))
      C_p1p1P12 = a + b0 + c + r
    else
      a = ab*A0mB0_p1p1P12(p2,m2,d)
      b0 = (-ab*d**3*m2 + d**2*(6*p2**2+5*m2*p2-m2**2)/(6*p2**3) + &
            d*(17*p2**2+4*m2*p2-m2**2)/(6*p2**3))*B0d(p2,m2,d,3)
      c =  -(m2/p2+cone)*C0d_m00(p2,m2,d,1)
      C_p1p1P12 = a + b0 + c
    end if

  end function

  function C_m00_p1P12P12_QP(p2,m2,d) result(C_p1P12P12)
    use triangle_aux_QP, only: A0mB0_p1P12P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12P12,a,b0,ab,c,r

    ab = (m2/(1+d) - 5*p2/2)/(3*p2**3*(1+d))

    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d, [-ab*m2,                            &
                         -(5*m2+3*p2)/(6*p2**2),            &
                         (m2**2+7*p2*m2+3*p2**2)/(3*p2**3), &
                         (m2**2+11*p2*m2+10*p2**2)/(3*p2**3)])
      c = -(m2**2+2*m2*p2+p2**2)/p2**2*C0d_m00(p2,m2,d,2)
      r = -m2/(6*p2**2*(1+d))
      C_p1P12P12 = a + b0 + c + r
    else
      a = ab*A0mB0_p1P12P12(p2,m2,d)
      b0 = (-ab*d**3*m2                          &
            -d**2*(5*m2+3*p2)/(6*p2**2)          &
            +d*(m2**2+7*p2*m2+3*p2**2)/(3*p2**3) &
            +(m2**2+11*p2*m2+10*p2**2)/(3*p2**3))*B0d(p2,m2,d,3)
      c = -(m2**2+2*m2*p2+p2**2)/p2**2*C0d_m00(p2,m2,d,2)
      C_p1P12P12 = a + b0 + c
    end if

  end function

  function C_m00_P12P12P12_QP(p2,m2,d) result(C_P12P12P12)
    use triangle_aux_QP, only: A0mB0_P12P12P12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12P12,a,b0,ab,c,r

    ab = ((m2+5*p2)/2-m2/(1+d))/(3*p2**3*(1+d))
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-ab*m2,                                &
                        (m2**2+5*m2*p2+2*p2**2)/(6*p2**3),     &
                        (-3*m2**2-10*m2*p2-3*p2**2)/(6*p2**3), &
                        (8*m2**2+11*m2*p2+3*p2**2)/(3*p2**3),  &
                        (11*(m2**2+2*m2*p2+p2**2))/(3*p2**3)])
      c = (-m2**3-3*m2**2*p2-3*m2*p2**2-p2**3)/(p2**3)*C0d_m00(p2,m2,d,3)
      r = m2/(6*p2**2*(1+d))
      C_P12P12P12 = a + b0 + c + r
    else
      a = ab*A0mB0_P12P12P12(p2,m2,d)
      b0 = B0d(p2,m2,d,[cnul,cnul,cnul,-d**3*ab*m2                      &
                        + d**2*(m2**2+5*m2*p2+2*p2**2)/(6*p2**3)  &
                        + d*(-3*m2**2-10*m2*p2-3*p2**2)/(6*p2**3) &
                        + (8*m2**2+11*m2*p2+3*p2**2)/(3*p2**3),   &
                        (11*(m2**2+2*m2*p2+p2**2))/(3*p2**3)])
      c = (-m2**3-3*m2**2*p2-3*m2*p2**2-p2**3)/(p2**3)*C0d_m00(p2,m2,d,3)
      C_P12P12P12 = a + b0 + c
    end if

  end function

  function C_m00_gP12_QP(p2,m2,d) result(C_gP12)
    use triangle_aux_QP, only: A0mB0_gP12
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gP12,a,b0,ab,r

    ab = -m2/(12*p2**2*(1+d))
    r = cone/18
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d, [(cone/12-ab*m2),           &
                         (-m2**2-m2*p2)/(12*p2**2), &
                         (-m2**2-2*m2*p2-p2**2)/(12*p2**2)])
    else
      a = ab*A0mB0_gP12(p2,m2,d)
      b0 = cone*B0d(p2,m2,d)/12 + &
           (-d**3*ab*m2 + d**2*(-m2**2-m2*p2)/(12*p2**2) + &
            d*(-m2**2-2*m2*p2-p2**2)/(12*p2**2))*B0d(p2,m2,d,3)
    end if
    C_gP12 = a + b0 + r
  end function

  function C_m00_gp1_QP(p2,m2,d) result(C_gp1)
    use triangle_aux_QP, only: A0mB0_gp1
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gp1,a,b0,ab,r

    ab = m2/(12*p2**2*(1+d))
    r = -7._qp/36
    if (abs(p2/m2) .gt. large_z_threshold) then
      a = ab*A0(m2)
      b0 = B0d(p2,m2,d,[-(cone/6+ab*m2),(m2**2-m2*p2-2*p2**2)/(12*p2**2)])
    else
      a = ab*A0mB0_gp1(p2,m2,d)
      b0 = B0d(p2,m2,d,[-cone/6,cnul,d**2*(m2**2-m2*p2-2*p2**2)/(12*p2**2),&
                        - d**3*ab*m2])
    end if
    C_gp1 = a + b0 + r
  end function

! --------------------------------------------------------
!               massless triangle
! --------------------------------------------------------

! ---------- rank 1 ----------------
  function C_000_p1_QP(p2,d) result(C_p1)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_p1

    C_p1 = B0d_0(p2,d,1)/p2 - C0d_000(p2,d,0)
  end function

  function C_000_P12_QP(p2,d) result(C_P12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_P12

    C_P12 = B0d_0(p2,d,[cnul,cone,2*cone])/p2 - C0d_000(p2,d,1)
  end function

! ---------- rank 2 ----------------
  function C_000_p1p1_QP(p2,d) result(C_p1p1)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_p1p1

    C_p1p1 = -3*B0d_0(p2,d,1)/(2*p2) + C0d_000(p2,d)
  end function

  function C_000_P12P12_QP(p2,d) result(C_P12P12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_P12P12

    C_P12P12 = B0d_0(p2,d,[cnul,cone/2,-cone,-3*cone])/p2 &
               + C0d_000(p2,d,2)
  end function

  function C_000_p1P12_QP(p2,d) result(C_p1P12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_p1P12

    C_p1P12 = B0d_0(p2,d,[cnul,-cone, -5*cone/2])/p2 + C0d_000(p2,d,1)
  end function

  function C_000_g_QP(p2,d) result(C_g)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_g

    C_g = cone/4 + B0d_0(p2,d,[cone/4,cone/4])
  end function

! ---------- rank 3 ----------------
  function C_000_p1p1p1_QP(p2,d) result(C_p1p1p1)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_p1p1p1

    C_p1p1p1 = 11*B0d_0(p2,d,1)/(6*p2) - C0d_000(p2,d)
  end function

  function C_000_p1p1P12_QP(p2,d) result(C_p1p1P12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_p1p1P12

    C_p1p1P12 = B0d_0(p2,d,[cnul, cone,17*cone/6])/p2 &
                - C0d_000(p2,d,1)
  end function

  function C_000_p1P12P12_QP(p2,d) result(C_p1P12P12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_p1P12P12

    C_p1P12P12 = B0d_0(p2,d,[cnul, -cone/2, cone, 10*cone/3])/p2 &
                 - C0d_000(p2,d,2)
  end function

  function C_000_P12P12P12_QP(p2,d) result(C_P12P12P12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_P12P12P12

    C_P12P12P12 = B0d_0(p2,d,[cnul, cone/3, -cone/2, cone, 11*cone/3])/p2 &
                     - C0d_000(p2,d,3)
  end function

  function C_000_gp1_QP(p2,d) result(C_gp1)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_gp1

    C_gp1 = -7._qp/36 + B0d_0(p2,d,[-cone/6, -cone/6])
  end function

  function C_000_gP12_QP(p2,d) result(C_gP12)
    real(qp), intent(in) :: p2,d
    complex(qp) :: C_gP12

    C_gP12 = B0d_0(p2,d,[cone/12, cnul, -cone/12]) + cone/18
  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     MMM                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_mmm_p1_QP(p2,m2,d) result(C_p1)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1

    C_p1 = B0d_mm(p2,m2,d,1)/p2 - C0d_mmm(p2,m2,d)

  end function

  function C_mmm_P12_QP(p2,m2,d) result(C_P12)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12, b

    b = B0d_mm(p2,m2,d,[cnul, cone/p2, 2*cone/p2])

    ! C_P12 = B0d_mm(p2,m2,muUV2,d,1)/p2 - C0d_mmm(p2,m2,d,1) + 2*B0d_mm(p2,m2,muUV2,d,2)/p2
    C_P12 = b - C0d_mmm(p2,m2,d,1)

  end function

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_mmm_p1p1_QP(p2,m2,d) result(C_p1p1)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1

    C_p1p1 = - 3._qp*B0d_mm(p2,m2,d,1)/(2._qp*p2) + C0d_mmm(p2,m2,d)

  end function

  function C_mmm_p1P12_QP(p2,m2,d) result(C_p1P12)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12,b

    b = B0d_mm(p2,m2,d,[cnul, -cone/p2, - 5*cone/(2*p2)])

    ! C_p1P12 = - B0d_mm(p2,m2,muUV2,d,1)/p2 + (cone - m2/p2)*C0d_mmm(p2,m2,d,1) &
    !           - 5._qp/(2*p2)*B0d_mm(p2,m2,muUV2,d,2)
    C_p1P12 = b + (cone - m2/p2)*C0d_mmm(p2,m2,d,1)

  end function

  function C_mmm_P12P12_QP(p2,m2,d) result(C_P12P12)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12, b

    b = B0d_mm(p2,m2,d,[cnul, cone/(2*p2), -cone/p2, -3*cone/p2])

    ! C_P12P12 = - 3._qp*B0d_mm(p2,m2,muUV2,d,3)/p2 &
    !            - B0d_mm(p2,m2,muUV2,d,2)/p2 + (cone - 2*m2/p2)*C0d_mmm(p2,m2,d,2) &
    !            + B0d_mm(p2,m2,muUV2,d,1)/(2*p2)
    C_P12P12 = b + (cone - 2*m2/p2)*C0d_mmm(p2,m2,d,2)

  end function

  function C_mmm_g_QP(p2,m2,d) result(C_g)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_g, b1, b0, c0, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    b = B0d_mm(p2,m2,d,[cone/4._qp, cone/4._qp])
    c0 = C0d_mmm(p2,m2,d,0)

    ! C_g = b1/4._qp + (b0+2._qp*m2*c0+cone)/4._qp
    C_g = b + (2._qp*m2*c0+cone)/4._qp

  end function

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_mmm_p1p1p1_QP(p2,m2,d) result(C_p1p1p1)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1p1, b1, b0, c0, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    b = B0d_mm(p2,m2,d,[ -m2/(3*(1+d)*p2**2), (2*m2 + 11*p2)/(6*p2**2)])
    c0 = C0d_mmm(p2,m2,d,0)

    ! C_p1p1p1 = (2*m2 + 11*p2)/(6*p2**2)*b1 &
    !            +(- m2*b0 + A0(m2,muUV2) -m2)/(3*(1+d)*p2**2) - c0
    C_p1p1p1 = b + ( A0(m2) -m2)/(3*(1+d)*p2**2) - c0

  end function

  function C_mmm_p1p1P12_QP(p2,m2,d) result(C_p1p1P12)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1P12, b3,b2,b1, b0, c1, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    b = B0d_mm(p2,m2,d,[m2/(3*(1+d)*p2**2), (-m2 + 3*p2)/(3*p2**2), (-4*m2+17*p2)/(6*p2**2) ])
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    c1 = C0d_mmm(p2,m2,d,1)

    ! C_p1p1P12 = (-4*m2+17*p2)/(6*p2**2)*b2 &
    !            + (-m2 + 3*p2)/(3*p2**2)*b1 &
    !            + (2*m2/p2 - cone)*c1 &
    !            + (m2*b0 - A0(m2,muUV2) + m2)/(3*(1+d)*p2**2)
    C_p1p1P12 = b + (2*m2/p2 - cone)*c1 + (- A0(m2) + m2)/(3*(1+d)*p2**2)


  end function

  function C_mmm_p1P12P12_QP(p2,m2,d) result(C_p1P12P12)
  !!!!!!!!!!!!!!!!!!!!
  !! for this function, using B0_mm_opt_list will cause meld_QP numerical instability
  !!!!!!!!!!!!!!!!!!!!
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12P12, b3,b2,b1, b0, c2, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    b = B0d_mm(p2,m2,d,[-m2/(3*(1+d)*p2**2), (2*m2 - 3*p2)/(6*p2**2), (-4*m2+3*p2)/(3*p2**2), 2*(-4*m2+5*p2)/(3*p2**2) ])
    c2 = C0d_mmm(p2,m2,d,2)

    ! C_p1P12P12 = 2*(-4*m2+5*p2)/(3*p2**2)*b3 &
    !            + (-4*m2+3*p2)/(3*p2**2)*b2 &
    !            + (4*m2/p2-cone)*c2 &
    !            + (2*m2 - 3*p2)/(6*p2**2)*b1 &
    !            + (-m2*b0 + A0(m2,muUV2) -m2)/(3*(1+d)*p2**2)
    C_p1P12P12 = b + (4*m2/p2-cone)*c2 + (A0(m2) -m2)/(3*(1+d)*p2**2)


  end function

  function C_mmm_P12P12P12_QP(p2,m2,d) result(C_P12P12P12)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12P12, b4,b3,b2,b1, b0, c3, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    ! b4 = B0d_mm(p2,m2,muUV2,d,4)
    b = B0d_mm(p2,m2,d,[ m2/(3*(1+d)*p2**2), (-m2+p2)/(3*p2**2), (4*m2-3*p2)/(6*p2**2), &
                              (-8*m2+3*p2)/(3*p2**2), (-16*m2+11*p2)/(3*p2**2) ])
    c3 = C0d_mmm(p2,m2,d,3)

    ! C_P12P12P12 = (-16*m2+11*p2)/(3*p2**2)*b4 &
    !             + (-8*m2+3*p2)/(3*p2**2)*b3 &
    !             + (18*m2-3*p2)/(3*p2)*c3 &
    !             + (4*m2-3*p2)/(6*p2**2)*b2 &
    !             + (-m2+p2)/(3*p2**2)*b1 &
    !             + (m2*b0 - A0(m2,muUV2)+m2)/(3*(1+d)*p2**2)

    C_P12P12P12 = b + (18*m2-3*p2)/(3*p2)*c3 + (- A0(m2)+m2)/(3*(1+d)*p2**2)

  end function

  function C_mmm_gp1_QP(p2,m2,d) result(C_gp1)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gp1, b4,b3,b2,b1, b0, c0, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    b = B0d_mm(p2,m2,d,[ -cone/6._qp, (2*m2-p2)/(6*p2) ])
    c0 = C0d_mmm(p2,m2,d,0)

    C_gp1 = b + (-18*m2*c0-7._qp)/36._qp

  end function

    function C_mmm_gP12_QP(p2,m2,d) result(C_gP12)
    complex(qp),     intent(in) :: m2
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gP12, b4,b3,b2,b1, b0, c1, b

    ! b0 = B0d_mm(p2,m2,muUV2,d,0)
    ! b1 = B0d_mm(p2,m2,muUV2,d,1)
    ! b2 = B0d_mm(p2,m2,muUV2,d,2)
    b = B0d_mm(p2,m2,d,[cone/12._qp, 4*m2/(12*p2), (8*m2-p2)/(12*p2)])
    ! b3 = B0d_mm(p2,m2,muUV2,d,3)
    ! b4 = B0d_mm(p2,m2,muUV2,d,4)
    c1 = C0d_mmm(p2,m2,d,1)

    C_gP12 = b +(-m2*p2*c1)/(2*p2) + cone/18._qp

  end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                   m0m1m1                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!
!  rank1  !
!!!!!!!!!!!

  function C_m0m1m1_p1_QP(p2,m02,m12,d) result(C_p1)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1

    C_p1 = B0d_m0m1(p2,m02,m12,d,1)/p2 - C0d_m0m1m1(p2,m02,m12,d)

  end function

  function C_m0m1m1_P12_QP(p2,m02,m12,d) result(C_P12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12

    C_P12 = B0d_m0m1(p2,m02,m12,d,[cnul, cone/p2, 2*cone/p2])  &
          - (cone + (m02-m12)/p2)*C0d_m0m1m1(p2,m02,m12,d,1)

  end function

!!!!!!!!!!!
!  rank2  !
!!!!!!!!!!!

  function C_m0m1m1_p1p1_QP(p2,m02,m12,d) result(C_p1p1)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1,z0,z1,b1,a1

    z0 = m02/p2
    z1 = m12/p2
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [(z1-z0)/(2*p2*(cone+d)),-(3*cone-z0+z1)/(2*p2)])
    a1 = (A0(m02)-A0(m12))/(2*p2**2*(cone+d))
    C_p1p1 =  a1 + b1 + C0d_m0m1m1(p2,m02,m12,d)

  end function

  function C_m0m1m1_p1P12_QP(p2,m02,m12,d) result(C_p1P12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12,z0,z1,a1,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    a1 = -(A0(m02)-A0(m12))/(2*p2**2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d,       &
                  [(z0-z1)/(2*p2*(cone+d)), &
                  -((2*cone+z0-z1)/(2*p2)), &
                  -((5*cone+z0-z1)/(2*p2))])
    c1 = (cone+z0-2*z1)*C0d_m0m1m1(p2,m02,m12,d,1)
    C_p1P12 = a1 + b1 + c1
  end function

  function C_m0m1m1_P12P12_QP(p2,m02,m12,d) result(C_P12P12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12,z0,z1,a1,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    a1 = (A0(m02)-A0(m12))/(2*p2**2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d,        &
                  [(z1-z0)/(2*p2*(1+d)),   &
                   (cone+z0-z1)/(2*p2),    &
                   (-cone-2*z0+2*z1)/(p2), &
                   -(3*(cone+z0-z1))/(p2)])

    c1 = ((cone + z0)**2 - 2*(2*cone+z0)*z1+z1**2)*C0d_m0m1m1(p2,m02,m12,d,2)
    C_P12P12 = a1 + b1 + c1
  end function

  function C_m0m1m1_g_QP(p2,m02,m12,d) result(C_g)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_g,z0,z1,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    b1 = B0d_m0m1(p2,m02,m12,d,[cone,cone+z0-z1])/4
    c1 = p2*z1*C0d_m0m1m1(p2,m02,m12,d,0)/2
    C_g = cone/4 + b1 + c1
  end function

!!!!!!!!!!!
!  rank3  !
!!!!!!!!!!!

  function C_m0m1m1_p1p1p1_QP(p2,m02,m12,d) result(C_p1p1p1)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1p1,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2
    a11 = ((z0-z1)/(cone+d) + (-5*cone+2*z0-2*z1)/2)*A0(m02)/(3*p2**2*(cone+d))
    a12 = -((z0-z1)/(cone+d) + (-7*cone+2*z0-2*z1)/2)*A0(m12)/(3*p2**2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d, &
                [(-z0**2+2*z0*z1-z1**2)/(3*p2*(cone+d)**2) + &
                (5*z0-2*z0**2-7*z1+4*z0*z1-2*z1**2)/(6*p2*(cone+d)), &
                (11*cone-5*z0+2*z0**2+7*z1-4*z0*z1+2*z1**2)/(6*p2)])
    c1 = -C0d_m0m1m1(p2,m02,m12,d)
    C_p1p1p1 =  a11 + a12 + b1 + c1 -(z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_p1p1P12_QP(p2,m02,m12,d) result(C_p1p1P12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1p1P12,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2
    a11 = ((z1-z0)/(cone+d) + (5*cone-z0+z1)/2)*A0(m02)/(3*p2**2*(cone+d))
    a12 = -((z1-z0)/(cone+d) + (7*cone-z0+z1)/2)*A0(m12)/(3*p2**2*(cone+d))

    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [(z0-z1)**2/(3*p2*(cone+d)**2) + &
                  (-5*z0+z0**2+7*z1-2*z0*z1+z1**2) /(6*p2*(cone+d)), &
                   -(-6*cone+z0**2+z1*(7*cone+z1)-z0*(5*cone+2*z1))/(6*p2), &
                   -(-17*cone+z0**2-2*z0*(2*cone+z1)+z1*(8*cone+z1))/(6*p2)])

    c1 = -(cone+z0-3*z1)*C0d_m0m1m1(p2,m02,m12,d,1)
    C_p1p1P12 = a11 + a12 + b1 + c1 + (z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_p1P12P12_QP(p2,m02,m12,d) result(C_p1P12P12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_p1P12P12,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2

    a11 = ((z0-z1)/(3*p2**2*(cone+d)**2)-5/(6*p2**2*(cone+d)))*A0(m02)
    a12 = (-(z0-z1)/(3*p2**2*(cone+d)**2)+7/(6*p2**2*(cone+d)))*A0(m12)

    b1 = B0d_m0m1(p2,m02,m12,d, &
          [-((z0 - z1)**2/(3*p2*(cone+d)**2))+(5*z0-7*z1)/( 6*p2*(cone+d)), &
           (-3*cone-5*z0+7*z1)/(6*p2), &
           (3*cone+7*z0+z0**2-11*z1-2*z0*z1+z1**2)/(3*p2), &
           (10*cone+11*z0+z0**2-19*z1-2*z0*z1+z1**2)/(3*p2)])

    c1 = (-cone-2*z0-z0**2+6*z1+4*z0*z1-3*z1**2)*C0d_m0m1m1(p2,m02,m12,d,2)
    C_p1P12P12 = a11 + a12 + b1 + c1 - (z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_P12P12P12_QP(p2,m02,m12,d) result(C_P12P12P12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_P12P12P12,z0,z1,a11,a12,b1,c1

    z0 = m02/p2
    z1 = m12/p2
    a11 = ((-z0+z1)/(3*p2**2*(cone+d)**2) + (5*cone+z0-z1)/(6*p2**2*(cone+d)))*A0(m02)
    a12 = -((-z0+z1)/(3*p2**2*(cone+d)**2) + (7*cone+z0-z1)/(6*p2**2*(cone+d)))*A0(m12)
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [(z0-z1)**2/(3*p2*(cone+d)**2)-(z0**2+z0*(5*cone-2*z1)+ &
                   (-7*cone+z1)*z1)/(6*p2*(cone+d)), &
                   (2*cone+z0**2+z0*(5*cone-2*z1)+(-7*cone+z1)*z1)/(6*p2), &
                   -((3*cone+3*z0**2+z0*(10*cone-6*z1)+z1*(-14*cone+3*z1))/(6*p2)), &
                   (3*cone+8*z0**2+z0*(11*cone-16*z1)+z1*(-19*cone+8*z1))/(3*p2), &
                   (11*cone-38*z1+11*(z0*(2*cone+z0)-2*z0*z1+z1**2))/(3*p2)])

    c1 = -(cone+z0-z1)*((cone+z0)**2-2*(4*cone+z0)*z1+z1**2)*C0d_m0m1m1(p2,m02,m12,d,3)
    C_P12P12P12 = a11 + a12 + b1 + c1 + (z0+z1)/(6*p2*(cone+d))

  end function

  function C_m0m1m1_gp1_QP(p2,m02,m12,d) result(C_gp1)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gp1,z0,z1,a11,a12,b1,b2,c1,r1

    z0 = m02/p2
    z1 = m12/p2
    a11 = (z0-z1)*A0(m02)/(12*p2*(cone+d))
    a12 = -(z0-z1)*A0(m12)/(12*p2*(cone+d))
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [-(z0-z1)**2/(12*(cone+d))-cone/6, &
                   (-2*cone-z0+z0**2+5*z1-2*z0*z1+z1**2)/(12)])
    c1 = -(p2*z1)*C0d_m0m1m1(p2,m02,m12,d)/2
    r1 = -(7*cone)/36
    C_gp1 = a11 + a12 + b1 + b2 + c1 + r1

  end function

  function C_m0m1m1_gP12_QP(p2,m02,m12,d) result(C_gP12)
    complex(qp),     intent(in) :: m02,m12
    real(qp),        intent(in) :: p2,d
    complex(qp) :: C_gP12,z0,z1,a1,b1,c1,r1

    z0 = m02/p2
    z1 = m12/p2
    a1 = (-z0+z1)/(12*p2*(1+d))*(A0(m02)-A0(m12))
    b1 = B0d_m0m1(p2,m02,m12,d, &
                  [cone/12 + &
                   (z0**2-2*z0*z1+z1**2)/(12*(cone+d)), &
                  -(z0+z0**2-2*z0*z1+(-5*cone+z1)*z1)/12, &
                  -((cone+z0)**2-2*(5*cone+z0)*z1+z1**2)/12])

    c1 = p2*z1*(-cone-z0+z1)*C0d_m0m1m1(p2,m02,m12,d,1)/2
    C_gP12 = a1 + b1 + c1 + cone/18

  end function

end module triangle_reduction_QP
