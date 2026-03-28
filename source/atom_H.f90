! Copyright (C) 2005 Barbara Ercolano
!
! Version 2.02
module atom_H
    use constants_mod
    use common_mod
    use continuum_mod
    use grid_mod
    use xSec_mod
    use elements_mod

    implicit none

    contains

    ! Calculate H fb and ff emission coefficients
    subroutine atom_H_fb_ff(TeUsed, sqrTeUsed, gammaHI, ffCoeff1, logGammaHIloc)
        implicit none
        real, intent(in) :: TeUsed, sqrTeUsed
        double precision, dimension(nbins), intent(out) :: gammaHI, ffCoeff1
        real, dimension(nlimGammaHI), intent(in) :: logGammaHIloc

        integer :: i, j, iup, ilow, HIPnuP
        double precision :: aFit, factor, expFactor, constant, phXSecHI
        real, dimension(3) :: statW

        gammaHI = 0.
        ffCoeff1 = 0.

        ! calculate gammaHI (nu<1.Ryd)
        do i = 1, nlimGammaHI/2
            iup  = HINuEdgeP(2*i)
            ilow = HINuEdgeP(2*i - 1)

            if (iup == ilow) then
                print*, "! fb_ff (H): frequency grid too course - divide&
                     & by zero                          [nu]",&
                     & nuArray(iup)
                stop
            end if
            aFit = (logGammaHIloc(2*i)-logGammaHIloc(2*i - 1)) /&
                 & (nuArray(iup) - nuArray(ilow))
            ! carry out interpolation between edges
            do j = ilow, iup-1
                gammaHI(j) = logGammaHIloc(2*i)+aFit* (nuArray(j)&
                     &-nuArray(iup))
                gammaHI(j) = 10.**gammaHI(j)
            end do
        end do

        ! calculate f-b emission for HI  (nu > 1.) by means of Milne relation
        constant = 4.9874105e-6 ! Ryd*Ryd*(h^2/(2*Pi*Me*K))**(3/2) [cm*K^(3/2)]
        HIPnuP   = HlevNuP(1)
        statW = (/2., 0.5, 2./)

        ! Set up free-free coefficient for Z=1
        do i = 1, nbins
            ffCoeff1(i) = ffCoeff(nuArray(i), 1, gauntFF(i), TeUsed, sqrTeUsed)
        end do

        do i = HIPnuP, nbins
            factor = (nuArray(i)*nuArray(i)*nuArray(i)) / (TeUsed*sqrTeUsed)
            expFactor = exp(dble( (-nuArray(i) + nuArray(HIPnuP)) * hcRyd_k / TeUsed))

            phXSecHI =  xSecArray(i-HIPnuP+1+HlevXSecP(1)-1)
            gammaHI(i) = fourPi * phXSecHI * statW(1) * hcRyd * constant*&
&                          factor * expFactor * 1.e20 *1.e20
            ! add the free free contribution
            gammaHI(i) = gammaHI(i)+ffCoeff1(i)
        end do
    end subroutine atom_H_fb_ff

    ! this function calculates the ff emission coeficient for interactions between
    ! ions of nuclear charge Z and electrons.
    ! see Allen pg 103.
    function ffCoeff(nu, Z, g, TeUsed, sqrTeUsed)
        implicit none

        double precision     :: ffCoeff     ! ff emission coefficient [e-40 erg*cm^3/s.Hz]
        real, intent(in)     :: g           ! ff gaunt
        real, intent(in)     :: nu          ! frequency [Ryd]
        real, intent(in)     :: TeUsed, sqrTeUsed

        integer, intent(in)  :: Z           ! nuclear charge

        ! local variables
        double precision     :: expFactor   ! exponential factor

        expFactor = exp(dble(-nu*hcRyd_k/TeUsed))

        ffCoeff = fourPi*54.43*Z*Z*g*expFactor/sqrTeUsed

    end function ffCoeff

    ! this function calculates the two photon emission coefficient
    ! for Hydrogenic ions. See Nussbaumer and Schmutz, A&A 138(1984)495
    function hydro2phot(nu, Z, TeUsed, NeUsed)
        implicit none

        double precision                :: hydro2phot   ! emission coeff for 2-phot cont [e-40 erg/cm^3/s/Hz]
                                            ! for the hydrogenic ion [e-40 erg*cm^3/s/Hz]
        real, intent(in)    :: nu           ! frequency [Ryd]
        real, intent(in)    :: TeUsed, NeUsed

        integer, intent(in) :: Z            ! nuclear charge

        ! local variables
        real                :: alphaEffH22S ! effective recombination coefficient
                                            ! to the 22S state of HI [e-13 cm^3/s]
        real                :: alphaEffHeII21S ! effective recombination coefficient
                                            ! to the 21S state of HeII [e-13 cm^3/s]
        real, parameter     :: A22S12S=8.23 ! Einstein A of forbibben LyAlpha [1/s]
        real, parameter     :: A2qH = 8.2249! total hydrogen 2s->1s 2 photon
                                            ! transition probability [1/s]
        real                :: A2qZ         ! total Z-ion 2s->1s 2 photon
                                            ! transition probability [1/s]
        real                :: Ay           ! transition probability (nu-dependent)
        real                :: factor       ! general interpolation factor
        real                :: gNu          ! spec distribution of 2phot emission
                                            ! [e-27 erg/Hz]
        real                :: nu0          ! frequencyof the 2s->1s transition
        real                :: Q22S22P      ! collisional transitional rates
                                            ! for HI 22S-22P
        real                :: y            ! nu/nu0 where nu0is the frequency
                                            ! of the 2s->1s transition
        ! select the hydrogenic ion
        select case(Z)
        ! for HI
        case (1)
            ! nu [Ryd] of 2s->1s transition (=1215.7 A)
            nu0 = 0.7496
            y = nu/nu0

            ! check for nu >= nu0
            if (nu >= nu0 ) then
                hydro2phot = 0.
                return
            end if

            ! for alphaEffH22S see Pengelly, MNRAS, 127(1964)145
            alphaEffH22S = 0.8368*((TeUsed*1.e-4)**(-0.723))

            Ay = 202.0*(y*(1-y)*(1-(4*y*(1-y))**0.8) + &
&                 0.88*((y*(1-y))**1.53)*(4*y*(1-y))**0.8)

            gNu = hPlanck*y*Ay/A2qH
            gNu = gNu*1e27

            ! calculate the collisional transitional rates for HI 22S-22P
            ! interpolate over temperature the value given in table 4.10
            ! from Osterbrock
            if (TeUsed <= 10000.) then
                Q22S22P = 5.31e-4
            else if (TeUsed >= 20000.) then
                Q22S22P = 4.71e-4
            else
                factor = log10(4.71e-4/5.31e-4) / log10(2.)
                Q22S22P = 10**( log10(5.31e-4) + factor*log10(TeUsed/10000.) )
            end if

            hydro2phot = alphaEffH22s*gNu / (1.+(NeUsed*Q22S22P)/A22S12S)

        case (2) ! HeII
            !  nu [Ryd] of 2s->1s transition (=303.8 A)
            nu0 = 3.00
            y = nu/nu0

            ! check for nu >= nu0
            if (nu >= nu0 ) then
                hydro2phot = 0.
                return
            end if

            ! for alphaEffHeII21S see Storey & Hummer, MNRAS, 272(1995)41
            ! the values used are for Ne = 100 cm^-3; however alphaEffHeII21S is not
            ! very sensitive to density.. only interpolate over temperature
            if (TeUsed <= 5000.) then
                alphaEffHeII21S = 6.161
            else if (TeUsed >= 30000.) then
                alphaEffHeII21S = 2.035
            else if (TeUsed > 5000. .and. TeUsed <= 10000. ) then
                factor = log10(4.091/6.161) / log10(2.)
                alphaEffHeII21S = 10**( log10(6.161) + factor*log10(TeUsed/5000.) )
            else if (TeUsed > 10000. .and. TeUsed <= 15000. ) then
                factor = log10(3.189/4.091) / log10(1.5)
                alphaEffHeII21S = 10**( log10(4.091) + factor*log10(TeUsed/10000.) )
            else if (TeUsed > 15000. .and. TeUsed < 30000. ) then
                factor = log10(2.035/3.189) / log10(2.)
                alphaEffHeII21S = 10**( log10(3.189) + factor*log10(TeUsed/15000.) )
            end if

            Ay = (Z**6)*0.9994667*202.0*(y*(1-y)*(1-(4*y*(1-y))**0.8) + &
&                 0.88*((y*(1-y))**1.53)*(4*y*(1-y))**0.8)

            A2qZ = 8.226*Z**6

            gNu = hPlanck*y*Ay/A2qZ
            gNu = gNu*1e27

            ! collisional de-excitation of the 22S of HeII is negligible

            ! caalculate the 2 photon emission
            hydro2phot = alphaEffHeII21s*gNu


        end select

    end function hydro2phot

    subroutine atom_H_RecLinesEmission(TeUsed, NeUsed, log10Te, grids, iG, ix, iy, iz, Hbeta)
        implicit none

        real, intent(in) :: TeUsed, NeUsed, log10Te
        type(grid_type), intent(in) :: grids(*)
        integer, intent(in) :: iG, ix, iy, iz
        real, intent(out) :: Hbeta

        integer :: iup, ilow
        real :: Lalpha

        ! copy HI data into array. this is where reading from file used to happen
        HIRecLines = HIRecLineData

        ! calculate Hbeta
        ! fits to Storey and Hummer MNRAS 272(1995)41
        Hbeta = 10**(-0.870*log10Te + 3.57)
        Hbeta = Hbeta*NeUsed*ionDenUsed(elementXref(1),2)*grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),1)

        ! calculate emission due to HI recombination lines [e-25 ergs/s/cm^3]
        do iup = 15, 3, -1
            do ilow = 2, min(8, iup-1)
                HIRecLines(iup, ilow) = HIRecLines(iup, ilow)*Hbeta
            end do
        end do

        ! add contribution of Lyman alpha
        ! fits to Storey and Hummer MNRAS 272(1995)41
        Lalpha = 10**(-0.897*log10Te + 5.05)
        HIRecLines(15, 8) =HIRecLines(15, 8) + grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),1)*&
             & ionDenUsed(elementXref(1),2)*&
             & NeUsed*Lalpha

    end subroutine atom_H_RecLinesEmission

    subroutine atom_H_inOpacity(xSecP, nuLowP, nuHighP, den, b, contBoltz, opacity)
        implicit none

        integer, intent(in)               :: nuLowP, nuHighP  ! pointers to lower and higher limits in nuArray
        integer, intent(in)               :: xSecP            ! x section pointer
        real, intent(in)                  :: b                ! departure coefficient
        real, intent(in)                  :: den              ! density of the lower level [cm^-3]
        real, dimension(nbins), intent(in):: contBoltz        ! Boltzmann factors
        real, dimension(:), intent(inout) :: opacity          ! Opacity array

        ! local variables
        integer             :: i, iup, k
        real                :: bInv

        k = xSecP - nuLowP
        iup = min(nuHighP, nbins)
        iup = max(nuLowP, iup)
        if (b > 1e-35) then
            bInv = 1./b
            do i = nuLowP, iup
                opacity(i) = opacity(i) + xSecArray(i+k)*den*&
&                                  max(0., 1.-contBoltz(i)*bInv)
            end do
        else
            do i = nuLowP, iup
                opacity(i) = opacity(i) + xSecArray(i+k)*den
           end do
        end if
    end subroutine atom_H_inOpacity

end module atom_H
