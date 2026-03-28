! Copyright (C) 2005 Barbara Ercolano
!
! Version 2.02
module atom_He
    use constants_mod
    use common_mod
    use continuum_mod
    use grid_mod
    use xSec_mod
    use elements_mod
    use atom_H

    implicit none

    contains

    ! Calculate He fb and ff emission coefficients
    subroutine atom_He_fb_ff(TeUsed, sqrTeUsed, gammaHeI, gammaHeII, ffCoeff2, logGammaHeIloc, logGammaHeIIloc, ffCoeff1)
        implicit none

        real, intent(in) :: TeUsed, sqrTeUsed
        double precision, dimension(nbins), intent(out) :: gammaHeI, gammaHeII, ffCoeff2
        double precision, dimension(nbins), intent(in) :: ffCoeff1
        real, dimension(nlimGammaHeI), intent(in) :: logGammaHeIloc
        real, dimension(nlimGammaHeII), intent(in) :: logGammaHeIIloc

        integer :: i, j, iup, ilow, HeIPnuP, HeIIPnuP
        double precision :: aFit, factor, expFactor, constant, phXSecHeI, phXSecHeII
        real, dimension(3) :: statW

        gammaHeI = 0.
        gammaHeII = 0.
        ffCoeff2 = 0.

        ! calculate gammaHeI (nu<0.799Ryd)
        do i = 1, nlimGammaHeI/2
            iup  = HeINuEdgeP(2*i)
            ilow = HeINuEdgeP(2*i - 1)

            if (iup == ilow) then
                print*, "! fb_ff (HeI): frequency grid too course - divide&
                     & by zero                          [nu]",&
                     & nuArray(iup)
                stop
            end if
            aFit = (logGammaHeIloc(2*i)-logGammaHeIloc(2*i - 1)) /&
                 & (nuArray(iup) - nuArray(ilow))
            ! carry out interpolation between edges
            do j = ilow, iup-1
                gammaHeI(j) = logGammaHeIloc(2*i)+aFit* (nuArray(j)&
                     &-nuArray(iup))
                gammaHeI(j) = 10.**gammaHeI(j)
            end do
        end do

        ! calculate gammaHeII (nu<4.0Ryd)
        do i = 1, nlimGammaHeII/2
            iup  = HeIINuEdgeP(2*i)
            ilow = HeIINuEdgeP(2*i - 1)

            if (iup == ilow) then
                print*, "! fb_ff (HeII): frequency grid too course - divide&
                     & by zero                          [nu]",&
                     & nuArray(iup)
                stop
            end if

            aFit = (logGammaHeIIloc(2*i)-logGammaHeIIloc(2*i - 1)) /&
                 & (nuArray(iup) - nuArray(ilow))
            ! carry out interpolation between edges
            do j = ilow, iup-1
                gammaHeII(j) = logGammaHeIIloc(2*i)+aFit* (nuArray(j)&
                     &-nuArray(iup))
                gammaHeII(j) = 10.**gammaHeII(j)
            end do
        end do

        constant = 4.9874105e-6 ! Ryd*Ryd*(h^2/(2*Pi*Me*K))**(3/2) [cm*K^(3/2)]
        HeIPnuP = HeIlevNuP(1)
        HeIIPnuP = HeIIlevNuP(1)
        statW = (/2., 0.5, 2./)

        ! Set up free-free coefficient for Z=2
        do i = 1, nbins
            ffCoeff2(i) = ffCoeff(nuArray(i), 2, gauntFFHeII(i), TeUsed, sqrTeUsed)
        end do

        do i = 1, nbins
            factor = (nuArray(i)*nuArray(i)*nuArray(i)) / (TeUsed*sqrTeUsed)
            ! calculate gammaHeI
            if (i >= HeINuEdgeP(nlimGammaHeI) ) then
                expFactor = exp(dble( (-nuArray(i) + nuArray(HeIPnuP)) * hcRyd_k / TeUsed))
                phXSecHeI = xSecArray(i-HeIPnuP+1+HeISingXSecP(1)-1)

                gammaHeI(i) = fourPi * phXSecHeI * statW(2) * hcRyd * constant*&
&                          factor * expFactor * 1.e20 *1.e20
                ! add the free free contribution
                gammaHeI(i) = gammaHeI(i) + ffCoeff1(i)
                if (i >= HeIINuEdgeP(nlimGammaHeII)) then
                    expFactor  = exp(dble( (-nuArray(i) + nuArray(HeIIPnuP)) * hcRyd_k / TeUsed))
                    phXSecHeII = xSecArray(i-HeIIPnuP+1+HeIIXSecP(1)-1)
                    gammaHeII(i) = fourPi * phXSecHeII * statW(3) * hcRyd * constant*&
&                              factor * expFactor * 1.e20 *1.e20
                    ! add the free free contribution
                    gammaHeII(i) = gammaHeII(i) + ffCoeff2(i)
                end if
            end if
        end do

    end subroutine atom_He_fb_ff

    subroutine atom_He_twoPhoton(twoPhotHeI, TeUsed)
        implicit none

        double precision, dimension(nbins), intent(out) :: twoPhotHeI
        real, intent(in) :: TeUsed

        ! local variables
        integer             :: i, j          ! counters
        real, parameter     :: A21S11S=51.3  ! HeI 2q [1/s]
        real, parameter     :: nu0=1.514     ! 21s -> 11s frequency [Ryd]
        real                :: alphaEff21SHeI! effective rec coeff to HeI 21S [e-14 cm^3/s]
        real                :: Ay            ! interpolated rate [1/s]
        real                :: fit1, fit2    ! interpolation coefficients
        real                :: y             ! nu/nu0

        twoPhotHeI = 0.

        alphaEff21SHeI = 6.23*((TeUsed/10000.)**(-0.827))

        j = 1
        do i = 1, nbins
            y  = nuArray(i)/nu0

            if ( y > y_dat(j) ) j=j+1
            if ( j >= 41 ) exit

            fit1 = (Ay_dat(j+1)-Ay_dat(j)) / &
&                                  (y_dat(j+1)-y_dat(j))
            fit2 = Ay_dat(j) - y_dat(j)*fit1

            Ay = fit1*y + fit2

            if ( Ay < 9.2163086E-03 ) Ay = 0.

            twoPhotHeI(i) = alphaEff21SHeI*Ay*0.66262*y/A21S11S
        end do

    end subroutine atom_He_twoPhoton

    subroutine atom_He_RecLinesEmission(TeUsed, NeUsed, log10Te, grids, iG, ix, iy, iz, HeIRecLines)
        implicit none

        real, intent(in) :: TeUsed, NeUsed, log10Te
        type(grid_type), intent(in) :: grids(*)
        integer, intent(in) :: iG, ix, iy, iz
        double precision, dimension(34), intent(out) :: HeIRecLines

        integer :: i, ilow, iup, denint
        real :: HeII4686, T4, x1, x2, coeff

        T4 = TeUsed / 10000.

        ! copy array. this is where the file used to be read in
        HeIIRecLines = HeIIRecLinedata

        ! calculate HeII 4686 [E-25 ergs*cm^3/s]
        HeII4686 = 10.**(-.997*log10Te+5.16)
        HeII4686 = HeII4686*NeUsed*grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),2)*ionDenUsed(elementXref(2),3)

        ! calculate emission due to HeII recombination lines [e-25 ergs/s/cm^3]
        do iup = 30, 3, -1
            do ilow = 2, min(16, iup-1)
                HeIIRecLines(iup, ilow) = HeIIRecLines(iup, ilow)*HeII4686
             end do
        end do

        ! now do HeI

        if (NeUsed <= 100.) then
           denint=0
        elseif (NeUsed > 100. .and. NeUsed <= 1.e4) then
           denint=1
        elseif (NeUsed > 1.e4 .and. NeUsed <= 1.e6) then
            denint=2
        elseif (NeUsed > 1.e6) then
           denint=3
        end if

        if (TeUsed > 5000.) then
           if (denint>0.and.denint<3) then
              do i = 1, 34
                  x1=HeIrecLineCoeff(i,denint,1)*(T4**(HeIrecLineCoeff(i,denint,2)))*exp(HeIrecLineCoeff(i,denint,3)/T4)
                  x2=HeIrecLineCoeff(i,denint+1,1)*(T4**(HeIrecLineCoeff(i,denint+1,2)))*exp(HeIrecLineCoeff(i,denint+1,3)/T4)
                  HeIRecLines(i) = x1+((x2-x1)*(NeUsed-100.**denint)/(100.**(denint+1)-100.**(denint)))
              end do
          elseif(denint==0) then
              do i = 1, 34
                 HeIRecLines(i) = HeIrecLineCoeff(i,1,1)*(T4**(HeIrecLineCoeff(i,1,2)))*exp(HeIrecLineCoeff(i,1,3)/T4)
              end do
          elseif(denint==3) then
              do i = 1, 34
                 HeIRecLines(i) = HeIrecLineCoeff(i,3,1)*(T4**(HeIrecLineCoeff(i,3,2)))*exp(HeIrecLineCoeff(i,3,3)/T4)
              end do
          end if
       else
           if (denint>0.and.denint<3) then
              do i = 1, 34
                 x1 = HeIrecLineCoeff(i,1,1)*((0.5)**(HeIrecLineCoeff(i,1,2)))*exp(HeIrecLineCoeff(i,1,3)/(0.5))
                 x2 = HeIrecLineCoeff(i,1,1)*((0.6)**(HeIrecLineCoeff(i,1,2)))*exp(HeIrecLineCoeff(i,1,3)/(0.6))
                 coeff = (LOG10(x1)-LOG10(x2))/(-0.079181246)
                 HeIRecLines(i) = (x1/(0.5**coeff))*T4**coeff
              end do
          elseif(denint==0) then
              do i = 1, 34
                 x1 = HeIrecLineCoeff(i,1,1)*((0.5)**(HeIrecLineCoeff(i,1,2)))*exp(HeIrecLineCoeff(i,1,3)/(0.5))
                 x2 = HeIrecLineCoeff(i,1,1)*((0.6)**(HeIrecLineCoeff(i,1,2)))*exp(HeIrecLineCoeff(i,1,3)/(0.6))
                 coeff = (LOG10(x1)-LOG10(x2))/(-0.079181246)
                 HeIRecLines(i) = (x1/(0.5**coeff))*T4**coeff
              end do
          elseif(denint==3) then
              do i = 1, 34
                 x1 = HeIrecLineCoeff(i,3,1)*((0.5)**(HeIrecLineCoeff(i,3,2)))*exp(HeIrecLineCoeff(i,3,3)/(0.5))
                 x2 = HeIrecLineCoeff(i,3,1)*((0.6)**(HeIrecLineCoeff(i,3,2)))*exp(HeIrecLineCoeff(i,3,3)/(0.6))
                 coeff = (LOG10(x1)-LOG10(x2))/(-0.079181246)
                 HeIRecLines(i) = (x1/(0.5**coeff))*T4**coeff
              end do
          end if
       endif

       HeIRecLines=HeIRecLines*NeUsed*grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),2)*ionDenUsed(elementXref(2),2)

    end subroutine atom_He_RecLinesEmission

    subroutine atom_He_inOpacity(xSecP, nuLowP, nuHighP, den, b, contBoltz, opacity)
        implicit none

        integer, intent(in)               :: nuLowP, nuHighP  ! pointers to lower and higher limits in nuArray
        integer, intent(in)               :: xSecP            ! x section pointer
        real, intent(in)                  :: b                ! departure coefficient
        real, intent(in)                  :: den              ! density of the lower level [cm^-3]
        real, dimension(nbins), intent(in):: contBoltz        ! Boltzmann factors
        real, dimension(:), intent(inout) :: opacity          ! Opacity array

        call atom_H_inOpacity(xSecP, nuLowP, nuHighP, den, b, contBoltz, opacity)
    end subroutine atom_He_inOpacity

end module atom_He
