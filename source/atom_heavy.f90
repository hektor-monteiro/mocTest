! Copyright (C) 2005 Barbara Ercolano
!
! Version 2.02
module atom_heavy
    use constants_mod
    use common_mod
    use continuum_mod
    use grid_mod
    use xSec_mod
    use elements_mod
    use interpolation_mod

    implicit none

    contains

    subroutine atom_heavy_fb_ff(TeUsed, sqrTeUsed, gammaHeavies, grids, iG, ix, iy, iz)
        implicit none

        real, intent(in) :: TeUsed, sqrTeUsed
        double precision, dimension(nbins), intent(out) :: gammaHeavies
        type(grid_type), intent(in) :: grids(*)
        integer, intent(in) :: iG, ix, iy, iz

        integer :: elem, ion, nElec, outShell, g0, g1, highNuP, IPnuP, xSecP, i
        double precision :: expFactor, phXSecM, constant, factor

        gammaHeavies = 0.
        constant = 4.9874105e-6 ! Ryd*Ryd*(h^2/(2*Pi*Me*K))**(3/2) [cm*K^(3/2)]

        do i = 1, nbins
           ! only elemnts up to z=19 calculated.. for now
           do elem = 3, 19
                do ion = 1, nstages-1
                    ! check if this element is present in the nebula
                    if (.not.lgElementOn(elem)) exit
                    if (ion > elem) exit

                    ! find the number of electron in this ion
                    nElec = elem - ion +1

                    ! find the outer shell number and stat weights
                    call getOuterShell(elem, nElec, outShell, g0, g1)

                    ! get pointer to this ion's high energy limit in nuArray
                    highNuP = elementP(elem, ion, outShell, 2)
                    if ( highNuP > nbins ) then
                        print*, "! fb_ff: high frequency limit beyond grid limit. Increase&
&                             frequency grid, or switch element off [nuMax,elem,&
&                             ion]", nuMax,elem,ion
                        stop
                    end if

                    ! get pointer to this ion's IP in nuArray
                    IPnuP = elementP(elem, ion, outShell, 1)
                    if ( IPnuP > nbins ) then
                        print*, "! fb_ff: IP frequency beyond grid limit. Increase&
&                             frequency grid, or switch element off [nuMax,elem,&
&                             ion]", nuMax,elem,ion
                        stop
                    end if

                    ! check IP of atom against current energy
                    if ( (i >= IPnuP) .and. (i < highNuP) ) then
                        xSecP = elementP(elem, ion, outShell, 3)
                        phXSecM = xSecArray(i+xSecP-IPnuP+1-1)

                        expFactor = exp( dble((-nuArray(i)+nuArray(IPnuP)) * hcRyd_k / TeUsed))

                        gammaHeavies(i) = gammaHeavies(i) + &
&                                  expFactor * phXSecM * (real(g0)/real(g1)) * &
&                                  ionDenUsed(elementXref(elem), ion+1) * grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),elem)
                    end if
                end do
            end do

            factor = nuArray(i)*nuArray(i)*nuArray(i) / (TeUsed*sqrTeUsed)
            gammaHeavies(i) = gammaHeavies(i) * fourPi * constant *&
&                                hcRyd * factor * 1.e20 * 1.e20
        end do
    end subroutine atom_heavy_fb_ff

    subroutine atom_heavy_forLines(TeUsed, NeUsed, grids, iG, ix, iy, iz)
        implicit none

        real, intent(in) :: TeUsed, NeUsed
        type(grid_type), intent(in) :: grids(*)
        integer, intent(in) :: iG, ix, iy, iz

        integer :: elem, ion

        forbiddenLines = 0.

        do elem = 3, nElements
           do ion = 1, min(elem+1, nstages)
              if (.not.lgElementOn(elem)) exit

              if (lgDataAvailable(elem, ion)) then
                 if (elem == 26 .and. ion == 2) then
                    if (nstages > 2) then
                       call equilibrium(file_name=dataFile(elem, ion), &
                            &ionDenUp=ionDenUsed(elementXref(elem),ion+1)&
                            &/ionDenUsed(elementXref(elem),ion), Te=TeUsed,&
                            &Ne=NeUsed, FlineEm=forbiddenLinesLarge(&
                            &1:nForLevelsLarge, 1:nForLevelsLarge))
                    else
                       call equilibrium(file_name=dataFile(elem, ion), &
                            &ionDenUp=0., Te=TeUsed,&
                            &Ne=NeUsed, FlineEm=forbiddenLinesLarge(&
                            &1:nForLevelsLarge, 1:nForLevelsLarge))
                    end if

                    forbiddenLinesLarge(:, :) =&
                         &forbiddenLinesLarge(:, :)*grids(iG)%elemAbun(&
                         &grids(iG)%abFileIndex(ix,iy,iz),elem)*&
                         & ionDenUsed(elementXref(elem), ion)

                 else
                    if (ion<nstages) then
                       call equilibrium(file_name=dataFile(elem, ion), &
                            &ionDenUp=ionDenUsed(elementXref(elem),ion+1)/&
                            &ionDenUsed(elementXref(elem),ion), Te=TeUsed,&
                            &Ne=NeUsed, FlineEm=forbiddenLines(elem, ion,&
                            &1:nForLevels, 1:nForLevels))
                    else
                       call equilibrium(file_name=dataFile(elem, ion), &
                            &ionDenUp=0., Te=TeUsed, Ne=NeUsed, &
                            &FlineEm=forbiddenLines(elem, ion,1:nForLevels,&
                            &1:nForLevels))
                    end if

                    forbiddenLines(elem, ion, :, :) =&
                         &forbiddenLines(elem, ion, :, :)*grids(iG)%elemAbun(&
                         &grids(iG)%abFileIndex(ix,iy,iz),elem)*&
                         & ionDenUsed(elementXref(elem), ion)

                 end if
              end if

           end do
        end do

        forbiddenLines = forbiddenLines*1.9865e9
    end subroutine atom_heavy_forLines


    subroutine equilibrium(file_name, ionDenUp, Te, Ne, fLineEm, wav)
        implicit none

        double precision, dimension(:,:), &
            & intent(inout) :: fLineEm     ! forbidden line emissivity
        double precision, dimension(:,:), &
            & intent(inout), optional :: wav         ! wavelength of transition [A]

        double precision     :: constant                  ! calculations constant
        double precision     :: delTeK                    ! Boltzmann exponent
        double precision     :: expFac                    ! calculations factor
        double precision     :: Eji                       ! energy between levels j and i
        double precision     :: sumN                      ! normalization factor for populations
        double precision     :: sqrTe                     ! sqrt(Te)

        double precision, allocatable          :: x(:,:)      ! matrix arrays
        double precision, allocatable          :: y(:)        !
        double precision, allocatable          :: n(:)        ! level population arrays

        real                 :: a_r(4),a_d(5),z,br       !
        real                 :: qomInt                    ! interpolated value of qom
        real                 :: log10Te

        real, intent(in) &
            & :: Te, &    ! electron temperature [K]
            & Ne, &       ! electron density [cm^-3]
            & ionDenUp    ! ion density of the upper ion stage

        integer  :: i, j, iT         ! counters/indeces
        integer  :: elem, ion        ! location of atomic data in the array
        integer  :: nLev             ! number of levels in atomic data file
        integer  :: nTemp            ! number of temperature points in atomic data file
        character(len = *), intent(in)  :: file_name   ! ionic data file name

        if (ionDenUp .ne. ionDenUp) then !NaN test, results from division by zero in the call. to be fixed upstream from here
            fLineEm = 0.d0
            return
        endif

        do i=3,nelements
            do j=1,10
                if (trim(atomic_data_array(i,j)%ion) .eq. trim(file_name)) then
                    goto 888
                endif
            enddo
        enddo

888     elem = i
        ion = j

        sqrTe   = sqrt(Te)
        log10Te = log10(Te)

        if (atomic_data_array(elem,ion)%nlevs > size(fLineEm(1,:))) then
           print*, '! equilibrium: model ion has more levels than &
                &allowed by nForLevels - please enlarge', atomic_data_array(elem,ion)%nlevs, &
                & size(fLineEm(1,:)),  nForLevels, file_name
           stop
        end if
        br = atomic_data_array(elem,ion)%br
        z = atomic_data_array(elem,ion)%z
        a_r = atomic_data_array(elem,ion)%a_r
        a_d = atomic_data_array(elem,ion)%a_d
        nlev = atomic_data_array(elem,ion)%nlevs
        ntemp = atomic_data_array(elem,ion)%ntemps

        atomic_data_array(elem,ion)%alphaTotal = 0.

        do j = 2, nLev
           if (atomic_data_array(elem,ion)%br .le. 0.d0) then
              exit
           else
             atomic_data_array(elem,ion)%alphaTotal(1) = 1.
           end if
           atomic_data_array(elem,ion)%alphaTotal(j)=1.e-13*br*z*a_r(1)*(Te*1.e-4/z**2)**a_r(2)/(1.+a_r(3)*&
                & (Te*1.e-4/Z**2)**a_r(4))+1.0e-12*br*(a_d(1)&
                & /(Te*1.e-4)+a_d(2)+a_d(3)*Te*1.e-4+a_d(4)*(Te*1.e-4)**2)*&
                & (Te*1.e-4)**(-1.5)*exp(-a_d(5)*1.e4/Te)
        end do

        ! form matrices
        ! set up qeff
        do i = 2, nLev
           do j = i, nLev
              do iT = 1, nTemp
                 atomic_data_array(elem,ion)%qq(iT) = real(atomic_data_array(elem,ion)%qom(iT, i-1, j))
              end do

              if (nTemp == 1) then
                 qomInt = atomic_data_array(elem,ion)%qq(1)
              else if (nTemp == 2) then
                 qomInt = atomic_data_array(elem,ion)%qq(1) + &
                      (atomic_data_array(elem,ion)%qq(2)-atomic_data_array(elem,ion)%qq(1)) / (atomic_data_array(elem,ion)%logTemp(2)-atomic_data_array(elem,ion)%logTemp(1)) * (log10Te - atomic_data_array(elem,ion)%logTemp(1))
              else
                 call linear_interpolation(atomic_data_array(elem,ion)%logTemp, atomic_data_array(elem,ion)%qq, log10Te, qomInt)
              end if

              atomic_data_array(elem,ion)%cs(i-1, j) = qomInt
              constant = 1.4388463d0
              delTeK = (atomic_data_array(elem,ion)%e(i-1)-atomic_data_array(elem,ion)%e(j))*constant
              expFac = exp( delTeK/Te )

              atomic_data_array(elem,ion)%qeff(i-1, j) = 8.63d-6 * atomic_data_array(elem,ion)%cs(i-1, j) * expFac /&
                   &(atomic_data_array(elem,ion)%g(i-1)*sqrTe)
              atomic_data_array(elem,ion)%qeff(j, i-1) = 8.63d-6 * atomic_data_array(elem,ion)%cs(i-1, j) / (atomic_data_array(elem,ion)%g(j)*sqrTe)

           end do
        end do

        allocate(x(nlev,nlev))
        allocate(y(nlev))
        allocate(n(nlev))

        x=0.
        y=0.
        n=0.

        do i= 2, nLev
           do j = 1, nLev
              x(1,:) = 1.
              y(1)   = 1.

              if (j /= i) then
                 x(i, j) = x(i, j) + Ne*atomic_data_array(elem,ion)%qeff(j, i)
                 x(i, i) = x(i, i) - Ne*atomic_data_array(elem,ion)%qeff(i, j)
                 if (j > i) then
                    x(i, j) = x(i, j) + atomic_data_array(elem,ion)%a(j, i)
                 else
                    x(i, i) = x(i, i) - atomic_data_array(elem,ion)%a(i, j)
                 end if
              end if
           end do
           y(i) = -Ne*ionDenUp*atomic_data_array(elem,ion)%alphaTotal(i)
        end do

        call luSlv(x, y, nLev)

        n = y

        sumN = 0.d0
        do i = 1, nLev
           sumN = sumN+n(i)
        end do
        do i = 1, nLev
           n(i) = n(i)/sumN
        end do

        do i = 1, nLev-1
           do j = i+1, nLev
              if (atomic_data_array(elem,ion)%a(j,i) /= 0.d0) then
                 Eji = (atomic_data_array(elem,ion)%e(j)-atomic_data_array(elem,ion)%e(i))
                 fLineEm(i,j) = atomic_data_array(elem,ion)%a(j,i) * Eji * n(j)
                 if (Eji>0.) then
                    if (present(wav)) wav(i,j) = 1.e8/Eji
                 else
                    if (present(wav)) wav(i,j) = 0.d0
                 end if
              end if
           end do
        end do

        deallocate(x)
        deallocate(y)
        deallocate(n)

    end subroutine equilibrium

    subroutine luSlv(a, b, n)
        implicit none
        integer, intent(in)                  :: n
        double precision,&
             & intent(inout), dimension(:,:) :: a
        double precision,&
             & intent(inout), dimension(:)   :: b

        call lured(a,n)
        call reslv(a,b,n)
    end subroutine luSlv

    subroutine lured(a,n)
        implicit none
        integer, intent(in)                  :: n
        double precision,&
             & intent(inout), dimension(:,:)    :: a
        integer          :: i, j, k
        double precision :: factor

        if (n == 1) return

        do i = 1, n-1
           do k = i+1, n
              factor = a(k,i)/a(i,i)
              do j = i+1, n
                 a(k, j) = a(k, j) - a(i, j) * factor
              end do
           end do
        end do
    end subroutine lured

    subroutine reslv(a,b,n)
        implicit none
        integer, intent(in)              :: n
        double precision,&
             & intent(inout), dimension(:,:) :: a
        double precision,&
             & intent(inout), dimension(:)   :: b
        integer    :: i, j, k, l

        if (n == 1) then
           b(n) = b(n) / a(n,n)
           return
        end if

        do i = 1, n-1
           do j = i+1, n
              b(j) = b(j) - b(i)*a(j, i)/ a(i, i)
           end do
        end do

        b(n) = b(n) / a(n,n)
        do i = 1, n-1
           k = n-i
           l = k+1
           do j = l, n
              b(k) = b(k) - b(j)*a(k, j)
           end do
           b(k) = b(k) / a(k,k)
        end do
    end subroutine reslv

    subroutine atom_heavy_putOpacity(nElem, densityArray, contBoltz, opacity)
        implicit none
        integer, intent(in) :: nElem
        real, dimension(:,:), intent(in) :: densityArray
        real, dimension(nbins), intent(in):: contBoltz
        real, dimension(:), intent(inout) :: opacity

        integer :: nIon, nuLowP, nuHighP, nShell, xSecP

        do nIon = 1, min(nElem, nstages)
            if ( densityArray(nElem, nIon) > 0. ) then
                do nShell = 1, nShells(nElem, nIon)
                    nuLowP = elementP(nElem, nIon, nShell, 1)
                    nuHighP = elementP(nElem, nIon, nShell, 2)
                    xSecP = elementP(nElem, nIon, nShell, 3)
                    call atom_heavy_inOpacity(xSecP, nuLowP, nuHighP, densityArray(nElem, nIon), 0., contBoltz, opacity)
                end do
            end if
        end do
    end subroutine atom_heavy_putOpacity

    subroutine atom_heavy_inOpacity(xSecP, nuLowP, nuHighP, den, b, contBoltz, opacity)
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
    end subroutine atom_heavy_inOpacity

end module atom_heavy
