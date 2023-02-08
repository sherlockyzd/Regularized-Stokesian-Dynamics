subroutine trajectory
! ------------------------------------------------------------ !
! The routine advances the particle trajectories using the
! displacements and velocities computed earlier in the code.
! 
! Inputs: none
! Outputs: none
! Changes: x - this is the updated position of the particles
!          x0 - this is the position of the particles in the
!              current time step.
! ------------------------------------------------------------ !
use global

implicit none

    integer*4 :: ii, jj, kk
    real*8 :: rot

    ! Advance the particles using the velocity due to the swimming gait 
    ! and an external force/ torque.  Note the components of u have 
    ! already been multiplied by delta t.  
    if ( ( istep .eq. 1 ) .or. ( iloop .eq. 1 ) ) then

        do ii = 1, n6

            x( ii ) = x0( ii ) + u( ii )

            uold( ii, 1 ) = u( ii )

        end do

    else if ( ( istep .eq. 2 ) .or. ( iloop .eq. 2 ) ) then

        do ii = 1, n6

            x( ii ) = x0( ii ) + 1.d0 / 2.d0 * ( 3.d0 * u( ii ) - uold( ii, 1 ) )

            uold( ii, 2 ) = u( ii )

        end do

    else if ( ( istep .eq. 3 ) .or. ( iloop .eq. 3 ) ) then

        do ii = 1, n6

            x( ii ) = x0( ii ) + 1.d0 / 12.d0 * ( 23.d0 * u( ii ) - 16.d0 * uold( ii, 2 ) &
                    + 5.d0 * uold( ii, 1 ) )

            uold( ii, 3 ) = u( ii )

        end do

    else

        do ii = 1, n6

            x( ii ) = x0( ii ) + 1.d0 / 24.d0 * ( 55.d0 * u( ii ) - 59.d0 * uold( ii, 3 ) &
                    + 37.d0 * uold( ii, 2 ) - 9.d0 * uold( ii, 1 ) )

            uold( ii, 1 ) = uold( ii, 2 )
            uold( ii, 2 ) = uold( ii, 3 )
            uold( ii, 3 ) = u( ii )

        end do

    end if


    ! Update the rotational displacement of the particles
    do ii = 1, np

        jj = 6 * ( ii - 1 ) + 3

        rot = dsqrt( x( jj + 3 ) * x( jj + 3 ) + x( jj + 1 ) * x( jj + 1 ) &
            + x( jj + 2 ) * x( jj + 2 ) )

        if ( rot .gt. twopi ) then

            do kk = 1, 3

                x( jj + kk ) = x( jj + kk ) * ( rot - twopi ) / rot

            end do

        end if

    end do

end subroutine trajectory
