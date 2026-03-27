! Copyright (C) 2005 Barbara Ercolano
!
! Version 2.02
module emission_mod
    use constants_mod
    use common_mod
    use continuum_mod
    use grid_mod
    use xSec_mod
    use elements_mod
    use atom_H
    use atom_He
    use atom_heavy

    type(resLine_type), allocatable         :: resLine(:)  ! resonant lines

    ! units are [e-25erg/s/N_gas]

    double precision, dimension(34)        :: HeIRecLines          ! emissivity from HeI rec lines

    real                                :: BjumpTemp=0   !
    real                                :: HbetaProb   ! probability of Hbeta
    real                                :: log10Ne     ! log10(Ne) at this cell
    real                                :: log10Te     ! log10(Te) at this cell
    real                                :: sqrTeUsed   ! sqr(Te) at this cell

    integer                             :: abFileUsed  ! abundance file index used
    integer                             :: cellPUsed   !  cell index
    integer                             :: nspE
        logical,save :: lgFloc=.true.                !
    contains

    ! this is the subroutine to drive the calculation
    ! of the emission spectrum at the given cell
    subroutine emissionDriver(grids, ix, iy, iz, ig)
        implicit none

        integer, intent(in) :: ix, iy, iz, ig  ! pointers to this cell


        type(grid_type), intent(inout) :: grids(*) ! the cartesian grid

        ! local variables

        integer :: err                       ! allocation error status
        integer :: i                         ! counters

        real    :: dV                        ! volume element

        double precision :: gammaHI(nbins)          ! HI fb+ff emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision :: gammaHeI(nbins)         ! HeI fb+ff emission coefficient [e-40erg*cm^3/s/Hz]
        double precision :: gammaHeII(nbins)        ! HeII fb+ff emission coefficient [e-40erg*cm^3/s/Hz]
        double precision :: gammaHeavies(nbins)     ! Heavy ions f-b emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision :: ffCoeff1(nbins)         ! ff emission coefficient for Z=1 ions [e-40 erg*cm^3/s/Hz]
        double precision :: ffCoeff2(nbins)         ! ff emission coefficient for Z=2 ions [e-40 erg*cm^3/s/Hz]
        double precision :: twoPhotHI(nbins)        ! HI 2photon emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision :: twoPhotHeI(nbins)       ! HeI 2photon emission coefficient [e-40 erg*cm^3/s/Hz]
        double precision :: twoPhotHeII(nbins)      ! HeII 2photon emission coefficient [e-40 erg*cm^3/s/Hz]

    ! units are [e-40erg/cm^3/s/Hz]
        double precision :: emissionHI(nbins)                     ! HI continuum emission coefficient
        double precision :: emissionHeI(nbins)                    ! HeI continuum emission coefficient
        double precision :: emissionHeII(nbins)                   ! HeII continuum emission coefficient

        ! check whether this cell is outside the nebula
        if (grids(iG)%active(ix, iy, iz)<=0) return

        cellPUsed = grids(iG)%active(ix, iy, iz)
        if (lgMultiDustChemistry) then
           nspE = grids(iG)%dustAbunIndex(cellPUsed)
        else
           nspE = 1
        end if

        ! set the dust emission PDF
        if (lgDust .and. .not.lgGas) call setDustPDF()

        if (.not.lgGas) return

        ! find the physical properties of this cell
        ionDenUsed= grids(iG)%ionDen(cellPUsed, :, :)

        abFileUsed= grids(iG)%abFileIndex(ix,iy,iz)
        NeUsed    = grids(iG)%Ne(cellPUsed)
        TeUsed    = grids(iG)%Te(cellPUsed)
        sqrTeUsed = sqrt(TeUsed)

        log10Te = log10(TeUsed)
        log10Ne = log10(NeUsed)

        ! check the electron density is not too low
        if (NeUsed < 1.e-5) then
             NeUsed = 1.e-5
        end if

        ! zero out arrays
        continuum      = 0.
        emissionHI     = 0.
        emissionHeI    = 0.
        emissionHeII   = 0.
        gammaHI        = 0.
        gammaHeI       = 0.
        gammaHeII      = 0.
        gammaHeavies   = 0.
        twoPhotHI      = 0.
        twoPhotHeI     = 0.
        twoPhotHeII    = 0.
        ffCoeff1       = 0.
        ffCoeff2       = 0.

        ! calculates the H, He and heavy ions f-b and f-f emission
        ! coefficients
        call fb_ff()

        ! calculate two photon emission coefficients for HI, HeI and
        ! HeII
        call twoPhoton()

        ! calculate the total emission coefficients for H and He

        emissionHI = (gammaHI + twoPhotHI) * grids(iG)&
                 &%elemAbun(grids(iG)%abFileIndex(ix,iy,iz), 1) * ionDenUsed(elementXref(1),2)*NeUsed &
                 &+ gammaHeavies*NeUsed
        emissionHeI = (gammaHeI + twoPhotHeI ) * grids(iG)&
                 &%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),2) * ionDenUsed(elementXref(2),2)*NeUsed
        emissionHeII= (gammaHeII + twoPhotHeII ) * grids(iG)&
                 &%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),2) * ionDenUsed(elementXref(2),3)*NeUsed

        ! zero out near zero contributions

        where (emissionHI < 1.e-15 ) emissionHI = 0.
        where (emissionHeI < 1.e-15 ) emissionHeI = 0.
        where (emissionHeII < 1.e-15 ) emissionHeII = 0.

        continuum = (emissionHI + emissionHeI + emissionHeII)

        dV = getVolume(grids(iG),ix,iy,iz)

        ! add contribution of this cell to the integrated balmer jump
        BjumpTemp = BjumpTemp + real(((emissionHI(BjumpP)&
             &+emissionHeI(BjumpP)+emissionHeII(BjumpP)) *&
             & cRyd*cRyd*nuArray(BjumpP)*nuArray(BjumpP)*1.e-8/c -&
             & (emissionHI(BjumpP-1)+emissionHeI(BjumpP-1)&
             &+emissionHeII(BjumpP-1)) *&
             & cRyd*cRyd*nuArray(BjumpP-1)*nuArray(BjumpP-1)*1.e-8/c)&
             & *grids(iG)%Hden(grids(iG)%active(ix,iy,iz))*dV)

        ! calculate the emission due to HI and HeI recombination lines
        call RecLinesEmission()

        ! calculate the emission due to the heavy elements forbidden lines
        call forLines()

        ! set the PDFs for diffuse radiation
        call setDiffusePDF()

        ! set the resonance line escape probabilities
        if (lgDust.and.convPercent>=resLinesTransfer .and. lgResLinesFirst &
             & .and. (.not.nIterateMC==1)) then
           call setResLineEscapeProb(grids,ig,ix,iy,iz)
        end if

        ! set the extra packets due to dust absorption of resonance lines
        if (lgDust.and.convPercent>=resLinesTransfer &
             & .and. (.not.nIterateMC==1)) then
           call initResLinePackets()
        end if

    contains

      subroutine initResLinePackets()
        implicit none

        real                :: lineIntensity

        integer             :: iRes, imul, n, ai

        species: do n = 1, nSpeciesPart(nspE)
           do ai = 1, nSizes
           
              if (grainRadius(ai) >= componentMinRadius(nspE, n) .and. grainRadius(ai) <= componentMaxRadius(nspE, n)) then

                 if (grids(iG)%Tdust(n,ai,cellPUsed)>0. .and. &
                      & grids(iG)%Tdust(n,ai,cellPUsed)<TdustSublime(dustComPoint(nspE)-1+n)) exit species
              endif

           end do
        end do species

        if (n>nSpeciesPart(nspE)) then
           ! all grains have sublimed
           grids(iG)%resLinePackets(cellPUsed) = 0
           return
        end if

        ! calculate total energy deposited into grains at location icell in grid igrid
        lineIntensity=0.
        do iRes = 1, nResLines
           do imul = 1, resLine(iRes)%nmul

              if (resLine(iRes)%elem==1) then
                 if ( resLine(iRes)%ion == 1 .and. resLine(iRes)%moclow(imul)==1 &
                      & .and. resLine(iRes)%mochigh(imul)==2 ) then

                    ! fits to Storey and Hummer MNRAS 272(1995)41
                    lineIntensity = lineIntensity+10**(-0.897*log10(TeUsed) + 5.05)* &
                         & grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),1)*ionDenUsed(elementXref(1),2)*&
                         &  NeUsed*grids(iG)%Hden(grids(iG)%active(ix,iy,iz))*dV*&
                         & (1.-grids(iG)%fEscapeResPhotons(grids(iG)%active(ix,iy,iz), iRes))

                 else

                    print*, "! initResLinePackets: [warning] only dust heating &
                         &from H Lyman alpha and resonance lines from heavy &
                         &elements is implemented in this version. &
                         &Please contact B. Ercolano -1-", ires

                 end if
              else if (resLine(iRes)%elem>2) then

                 lineIntensity = lineIntensity+&
                      &real(forbiddenLines(resLine(iRes)%elem,resLine(iRes)%ion,&
                      &resLine(iRes)%moclow(imul),resLine(iRes)%mochigh(imul))*&
                      & grids(iG)%Hden(grids(iG)%active(ix,iy,iz))*dV*&
                      & (1.-grids(iG)%fEscapeResPhotons(grids(iG)%active(ix,iy,iz), iRes)))

              else

                 print*, "! initResLinePackets: [warning] only dust heating from &
                      &H Lyman alpha and resonance lines from heavy &
                      &elements is implemented in this version. &
                      &Please contact B. Ercolano -2-", ires


              end if

           end do
        end do

        ! calculate number of dust emitted extra packets to be transfered from this location
        ! NOTE :1.e-16 from 1.e45 (from dV) * 1.e-36 (from Lstar) * 1.e-25 (from local forlines)
        grids(iG)%resLinePackets(grids(iG)%active(ix,iy,iz)) = nint(lineIntensity*1.e-16/(Lstar(1)/Nphotons(1)))

      end subroutine initResLinePackets



      ! this subroutine calculates the H, He and heavy ions f-b and f-f
      ! emission coefficients
      ! - it can only be called by emission driver -
      subroutine fb_ff()
        implicit none

        logical            :: lgTeOut         ! is Te out of range of
        integer            :: nTemp           ! =1 if 5000K<Te<10000K, =2 if 10000K<Te<20000K
        integer            :: i               ! counters

        double precision                :: aFit, factor

        real, allocatable :: logGammaHIloc(:), logGammaHeIloc(:), logGammaHeIIloc(:)
                                                      ! gamma [e-40 erg/Hz] for HI, HeI and HeII

        allocate (logGammaHIloc(nlimGammaHI), stat=err)
        if (err /= 0) then
            print*, "! fb_ff: can't allocate array memory"
            stop
        end if
        allocate (logGammaHeIloc(nlimGammaHeI), stat=err)
        if (err /= 0) then
            print*, "! fb_ff: can't allocate array memory"
            stop
        end if
        allocate (logGammaHeIIloc(nlimGammaHeII), stat=err)
        if (err /= 0) then
            print*, "! fb_ff: can't allocate array memory"
            stop
        end if

        logGammaHIloc = 0.
        logGammaHeIloc = 0.
        logGammaHeIIloc = 0.

        ! check what coefficients are needed for the local electron Temperature
        lgTeOut = .false. ! (re)initialize out of range Te flag
        if (TeUsed <= tkGamma(1)) then
            lgTeOut = .true.
            logGammaHIloc = logGammaHI(1,:)
            logGammaHeIloc = logGammaHeI(1,:)
            logGammaHeIIloc = logGammaHeII(1,:)
        else if (TeUsed >= tkGamma(ntkGamma)) then
            lgTeOut = .true.
            logGammaHIloc = logGammaHI(ntkGamma,:)
            logGammaHeIloc = logGammaHeI(ntkGamma,:)
            logGammaHeIIloc = logGammaHeII(ntkGamma,:)
         else
            call locate (tkGamma, TeUsed, nTemp)

            factor = log10(TeUsed/tkGamma(nTemp))
            do i = 1, nlimGammaHI
               aFit = (logGammaHI(nTemp+1,i)-logGammaHI(nTemp,i))/ &
                    & log10(tkGamma(ntemp+1)/tkGamma(ntemp))
               logGammaHIloc(i) = logGammaHI(nTemp,i) + real(aFit * factor)
            end do

            do i = 1, nlimGammaHeI
               aFit = (logGammaHeI(nTemp+1,i)-logGammaHeI(nTemp,i))/ &
                    & log10(tkGamma(ntemp+1)/tkGamma(ntemp))
               logGammaHeIloc(i) = logGammaHeI(nTemp,i) + real(aFit * factor)

            end do

            do i = 1, nlimGammaHeII
               aFit = (logGammaHeII(nTemp+1,i)-logGammaHeII(nTemp,i))/ &
                    & log10(tkGamma(ntemp+1)/tkGamma(ntemp))
               logGammaHeIIloc(i) = logGammaHeII(nTemp,i) + real(aFit * factor)

            end do
        end if

        call atom_H_fb_ff(TeUsed, sqrTeUsed, gammaHI, ffCoeff1, logGammaHIloc)
        call atom_He_fb_ff(TeUsed, sqrTeUsed, gammaHeI, gammaHeII, ffCoeff2, logGammaHeIloc, logGammaHeIIloc, ffCoeff1)
        call atom_heavy_fb_ff(TeUsed, sqrTeUsed, gammaHeavies, grids, iG, ix, iy, iz)

    end subroutine fb_ff

    subroutine twoPhoton()
        implicit none
        integer :: i

        do i = 1, nbins
            twoPhotHI(i) = hydro2phot(nuArray(i), 1, TeUsed, NeUsed)
            twoPhotHeII(i) = hydro2phot(nuArray(i), 2, TeUsed, NeUsed)
        end do
        call atom_He_twoPhoton(twoPhotHeI, TeUsed)
    end subroutine twoPhoton

    subroutine RecLinesEmission()
        implicit none
        real :: Hbeta
        call atom_H_RecLinesEmission(TeUsed, NeUsed, log10Te, grids, iG, ix, iy, iz, Hbeta)
        call atom_He_RecLinesEmission(TeUsed, NeUsed, log10Te, grids, iG, ix, iy, iz, HeIRecLines)
    end subroutine RecLinesEmission

    subroutine forLines()
        implicit none
        call atom_heavy_forLines(TeUsed, NeUsed, grids, iG, ix, iy, iz)
    end subroutine forLines

    subroutine setDiffusePDF()
        implicit none

        ! local variables

        real (kind=8)     :: tg
        real(kind=8),dimension(nTbins) :: Tspike,Pspike

        real              :: alpha2tS          ! effective rec coeff to the 2s trip HeI
        real              :: bb                ! blackbody
        real              :: const             ! general calculation constant
        real              :: correction        ! general correction term
        real              :: normalize         ! normHI+normHeI+normHeII
        real              :: normDust           ! normalization constant for dust
        real              :: normFor           ! normalization constant for lines
        real              :: normHI            ! normalization constant for HI
        real              :: normHeI           !                            HeI
        real              :: normHeII          !                            HeII
        real              :: normRec           !                            rec. lines
        real              :: treal
        real              :: T4                ! TeUsed/10000.

        real, allocatable     :: sumDiffuseDust(:) ! summation terms for dust emission
        real, allocatable     :: sumDiffuseHI(:)   ! summation terms for HI emission
        real, allocatable     :: sumDiffuseHeI(:)  ! summation terms for HeI emission
        real, allocatable     :: sumDiffuseHeII(:) ! summation terms for HeII emission

        integer           :: dcp                ! local dustComPointer
        integer           :: elem, ion          ! counters
        integer           :: err                ! allocation error status
        integer           :: j2TsP              ! pointer to 2s triplet state in nuArray
        integer           :: i,iup,ilow,j       ! counters
        integer           :: nuP                ! frequency pointer in nuArray
        integer           :: nS, ai, freq, iT   ! dust counters

        character(len=50) :: cShapeLoc

        real,dimension(NHeIILyman) :: HeIILyman ! HeII Lyman lines em. [e-25ergs*cm^3/s]

        cShapeLoc = 'blackbody'

        ! allocate space for the summation arrays
        allocate(sumDiffuseHI(1:nbins), stat = err)
        if (err /= 0) then
            print*, "! setDiffusePDF: can't allocate recPDF array memory"
            stop
        end if
        allocate(sumDiffuseHeI(1:nbins), stat = err)
        if (err /= 0) then
            print*, "! setDiffusePDF: can't allocate recPDF array memory"
            stop
        end if
        allocate(sumDiffuseHeII(1:nbins), stat = err)
        if (err /= 0) then
            print*, "! setDiffusePDF: can't allocate recPDF array memory"
            stop
        end if



        ! assign pointers for the He line photons
        call locate(nuArray, 1.45673, j2TsP)

        ! calculate normalization constants for the HI, HeI and HeII diffuse probs

        ! initialize constants
        normHI   = 0.
        normHeI  = 0.
        normHeII = 0.

        sumDiffuseHI = 0.
        sumDiffuseHeI = 0.
        sumDiffuseHeII = 0.

        ! perform summations
        do i = 1, nbins

           sumDiffuseHI(i) = real(cRyd*widFlx(i)*emissionHI(i)/1.e15)

           normHI = normHI + sumDiffuseHI(i)

           sumDiffuseHeI(i) =  real(cRyd*widFlx(i)*emissionHeI(i)/1.e15)

           normHeI = normHeI + sumDiffuseHeI(i)

           sumDiffuseHeII(i) = real(cRyd*widFlx(i)*emissionHeII(i)/1.e15)

           normHeII = normHeII + sumDiffuseHeII(i)

        end do


        ! Add contribution of HeI lines able to ionize HI

        ! triplet states:
        ! calculate effective recombination coefficients to triplet states[e-14 cm^3/s]
        ! Benjamin, SKillman and mits, 1999, ApJ 514, 307
        
        ! limit T4 to the data limit used in the paper fits
        T4 = MAX(5000.,TeUsed)*1.e-4
        
        alpha2tS = 27.2*(T4**(-0.678))

        ! calculate correction for partial sums and normalization constant
        ! triplet 2s
        ! NOTE: multiply by 1.e11 to make units of [e-25 erg/cm^3/s] as e-14 from
        ! rec coeff (see above).
        correction = alpha2tS * &
             & NeUsed*grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),2)*&
             &ionDenUsed(elementXref(2),2)*1.45673*hcRyd*1e11

        sumDiffuseHeI(j2TsP) = sumDiffuseHeI(j2TsP) + correction
        normHeI = normHeI + correction

        ! Add contribution of HeII lines able to ionize HI and HeII

        ! copy HeIILymanData into new array.  this is where file used to be read in

        HeIILyman = HeIILymanData

        ! calculate Lyman alpha first
        HeIILyman(4) = 10.**(-0.792*log10(TeUsed)+6.01)
        HeIILyman(4) = HeIILyman(4)*NeUsed*ionDenUsed(elementXref(2),3)*grids(iG)%elemAbun(grids(iG)%abFileIndex(ix,iy,iz),2)

        ! calculate emission due to HeII Lyman lines [e-25 ergs/s/cm^3]
        do i = 1, NHeIILyman-1
            HeIILyman(i) = HeIILyman(i)*HeIILyman(4)
        end do

        ! add contributions in respective bins
        do i = 1, NHeIILyman
            ! locate the frequency bin in the array
            call locate(NuArray, HeIILymanNu(i), NuP)
            ! add contribution to HeII partial sums ...
            sumDiffuseHeII(NuP) = sumDiffuseHeII(NuP) + HeIILyman(i)
            ! ... and to HeII normilazation constant
            normHeII = normHeII + HeIILyman(i)
        end do


        ! calculate dust emission
        if (lgDust) then

           dcp = dustComPoint(nspE)

           allocate(sumDiffuseDust(1:nbins), stat = err)
           if (err /= 0) then
              print*, "! setDiffusePDF: can't allocate sumDiffuseDust array memory"
              stop
           end if

           sumDiffuseDust = 0.
           normDust       = 0.
           do ai = 1, nSizes
              do nS = 1, nSpeciesPart(nspE)
                 
                 if (grainRadius(ai) >= componentMinRadius(nspE, nS) .and. grainRadius(ai) <= componentMaxRadius(nspE, nS)) then

                    if (grids(iG)%Tdust(nS,ai,cellPUsed)<TdustSublime(dcp-1+nS)) then

                       if (lgQHeat .and. grainRadius(ai)<minaQHeat .and. convPercent>minConvQheat .and. nIterateMC>1) then

                          tg =  grids(iG)%Tdust(nS,ai,cellPused)
                          Tspike=0.
                          Pspike=0.

                          call qHeat(nS, ai,tg,Tspike,Pspike)

                          if (lgWritePss .and. taskid==0) then
                             write(89, *) '              ix iy iz cellp nSpecies aRadius Teq'
                             write(89, *)  ix, iy, iz, cellpUsed, ns, grainRadius(ai), tg
                             write(89, *) '              i Tspike(i) Pspike(i)'
                             do iT = 1, nTbins
                                write(89, *) it, Tspike(it), Pspike(it)
                             end do
                             write(89, *) 'Tmean: ', tg
                             write(89, *) ' '
                          end if

                          do freq = 1, nbins
                             do iT=1,nTbins
                                treal = real(Tspike(iT))
                                bb = getFlux(nuArray(freq), treal, cShapeLoc)
                                sumDiffuseDust(freq) = sumDiffuseDust(freq) + &
                                     &  real(xSecArray(dustAbsXsecP(dcp-1+nS,ai)+freq-1)*bb*widFlx(freq)*&
                                     & grainWeight(ai)*grainAbun(nspE, nS)*Pspike(iT))
                                normDust = normDust+real(xSecArray(dustAbsXsecP(dcp+nS-1,ai)+freq-1)*bb*widFlx(freq)*&
                                     & grainWeight(ai)*grainAbun(nspE, nS)*Pspike(iT))
                             end do
                          end do

                       else

                          do freq = 1, nbins
                             bb = getFlux(nuArray(freq), grids(iG)%Tdust(nS,ai,cellPUsed), cShapeLoc)
                             sumDiffuseDust(freq) = sumDiffuseDust(freq) + &
                                  &  xSecArray(dustAbsXsecP(dcp-1+nS,ai)+freq-1)*bb*widFlx(freq)*&
                                  & grainWeight(ai)*grainAbun(nspE,nS)
                             normDust = normDust+xSecArray(dustAbsXsecP(dcp-1+nS,ai)+freq-1)*bb*widFlx(freq)*&
                                  & grainWeight(ai)*grainAbun(nspE, nS)
                          end do
                       end if
                    end if 
                    
                 endif !!! grain size interval if block

              end do
           end do

           ! the hPlanck is re-introduced here (was excluded in the bb calcs)
           const = hPlanck*4.*Pi*1.e25*cRyd
           sumDiffuseDust = sumDiffuseDust*const*grids(iG)%Ndust(cellPUsed)/grids(iG)%Hden(cellPUsed)
           normDust = normDust*const*grids(iG)%Ndust(cellPUsed)/grids(iG)%Hden(cellPUsed)

        end if

        ! Total normalization constant
        normalize = normHI + normHeI + normHeII

        if (lgDust) normalize = normalize+normDust

        ! Sum  up energy in recombination lines

        normRec = 0.
        ! HI
        do iup = 3, 15
            do ilow = 2, min0(8, iup-1)
                normRec = normRec+real(HIRecLines(iup, ilow))
            end do
        end do
        ! HeII
        do iup = 3, 30
            do ilow = 2, min0(16, iup-1)
                normRec = normRec+real(HeIIRecLines(iup, ilow))
            end do
        end do
        ! HeI singlets
        do i = 1, 34
            normRec = normRec+real(HeIRecLines(i))
        end do

        ! Sum  up energy in forbidden lines

        normFor = 0.
        do elem = 3, nElements
           do ion = 1, min(elem+1, nstages)
              if (elem ==26 .and. ion==2) then
                 do ilow = 1, nForLevelsLarge
                    do iup = 1, nForLevelsLarge
                       if (forbiddenLinesLarge(ilow, iup) > 1.e-35)  then
                          normFor = normFor + real(forbiddenLinesLarge(ilow, iup))
                       end if
                    end do
                 end do
              else
                 do ilow = 1, nForLevels
                    do iup = 1, nForLevels
                       if (forbiddenLines(elem, ion, ilow, iup) > 1.e-35)  then
                          normFor = normFor + real(forbiddenLines(elem, ion, ilow, iup))
                       end if
                    end do
                 end do
              end if
           end do
        end do

        ! total energy in lines (note: this variable will later be turned into the fraction of
        ! non-ionizing line photons as in declaration)

        grids(iG)%totalLines(cellPUsed) = normRec + normFor

        if (lgDust) then

           do i = 1, nbins
              if (i == 1) then
                 grids(iG)%recPDF(cellPUsed, i) = sumDiffuseHI(i) + sumDiffuseHeI(i) +&
                      & sumDiffuseHeII(i) + sumDiffuseDust(i)
              else
                 grids(iG)%recPDF(cellPUsed, i) = grids(iG)%recPDF(cellPUsed, i-1) +&
                      & sumDiffuseHI(i) + sumDiffuseHeI(i) +&
                      & sumDiffuseHeII(i) + sumDiffuseDust(i)
              end if

           end do

        else

           do i = 1, nbins
              if (i == 1) then
                 grids(iG)%recPDF(cellPUsed, i) = (sumDiffuseHI(i) + sumDiffuseHeI(i) +&
                      & sumDiffuseHeII(i) )
              else
                 grids(iG)%recPDF(cellPUsed, i) = grids(iG)%recPDF(cellPUsed, i-1) +&
                      & (sumDiffuseHI(i) + sumDiffuseHeI(i) +&
                      & sumDiffuseHeII(i) )
              end if
           end do

        end if

        grids(iG)%recPDF(cellPUsed, :) = grids(iG)%recPDF(cellPUsed, :)/grids(iG)%recPDF(cellPUsed, nbins)

        grids(iG)%recPDF(cellPUsed, nbins) = 1.

        do i = 1, nbins
           if (grids(iG)%recPDF(cellPUsed, i) > 0.999998) grids(iG)%recPDF(cellPUsed, i) = 1.
        enddo

        ! calculate the PDF for recombination and forbidden lines
        if (lgDebug) then
           i = 1
           ! HI recombination lines
           do iup = 3, 15
              do ilow = 2, min(8, iup-1)
                 if (i == 1) then

                    grids(iG)%linePDF(cellPUsed, i) = real(HIRecLines(iup, ilow) / &
&                                                  grids(iG)%totalLines(grids(iG)%active(ix,iy,iz)))



                 else

                    grids(iG)%linePDF(cellPUsed, i) = grids(iG)%linePDF(cellPUsed, i-1) + &
&                        real(HIRecLines(iup, ilow) / grids(iG)%totalLines(cellPUsed))



                 end if
                 i = i+1

              end do
           end do

           ! HeI singlet recombination lines
           do j = 1, 34
              grids(iG)%linePDF(cellPUsed, i) = grids(iG)%linePDF(cellPUsed, i-1) + &
&                           real(HeIRecLines(j) / grids(iG)%totalLines(cellPUsed))
              i = i+1
           end do

           ! HeII recombination lines
           do iup = 3, 30
              do ilow = 2, min(16, iup -1)
                 grids(iG)%linePDF(cellPUsed, i) = grids(iG)%linePDF(grids(iG)%active(ix, iy,&
                      & iz), i-1) + real(HeIIRecLines(iup, ilow) / grids(iG)&
                      &%totalLines(grids(iG)%active(ix,iy,iz)))
                 i = i+1
              end do
           end do

           ! heavy elements forbidden lines
           do elem = 3, nElements
              do ion = 1, min(elem+1, nstages)

                 if (.not.lgElementOn(elem)) exit
                 if (lgDataAvailable(elem, ion)) then
                    if (elem == 26 .and. ion == 2) then

                       do iup = 1, nForLevelsLarge
                          do ilow = 1, nForLevelsLarge
                             grids(iG)%linePDF(cellPUsed, i) = grids(iG)%linePDF(grids(iG)%active(ix, iy&
                                  &, iz), i-1) + real(forbiddenLinesLarge(iup, ilow)&
                                  & / grids(iG)%totalLines(grids(iG)%active(ix,iy,iz)))

                             if (grids(iG)%linePDF(cellPUsed, i) > 1. ) grids(iG)&
                                  &%linePDF(cellPUsed, i) = 1.
                             i = i+1

                          end do
                       end do

                    else

                       do iup = 1, nForLevels
                          do ilow = 1, nForLevels
                             grids(iG)%linePDF(cellPUsed, i) = grids(iG)%linePDF(grids(iG)%active(ix, iy&
                                  &, iz), i-1) + real(forbiddenLines(elem, ion, iup, ilow)&
                                  & / grids(iG)%totalLines(grids(iG)%active(ix,iy,iz)))

                             if (grids(iG)%linePDF(cellPUsed, i) > 1. ) grids(iG)&
                                  &%linePDF(cellPUsed, i) = 1.
                             i = i+1
                          end do
                       end do

                    end if
                 end if

              end do
           end do
           grids(iG)%linePDF(cellPUsed, nLines) = 1.

           ! calculate the probability of Hbeta
           HbetaProb = real(HIRecLines(4,2) / (grids(iG)%totalLines(grids(iG)%active(ix,iy,iz))&
                &+normalize))

        end if ! debug condition

        ! now compute the fraction of non-ionizing photons

        normalize = normalize + grids(iG)%totalLines(grids(iG)%active(ix,iy,iz))

        grids(iG)%totalLines(grids(iG)%active(ix,iy,iz)) = grids(iG)%totalLines(grids(iG)%active(ix,iy,iz)) / normalize
        ! deallocate arrays
        if ( allocated(sumDiffuseHI) ) deallocate(sumDiffuseHI)
        if ( allocated(sumDiffuseHeI) ) deallocate(sumDiffuseHeI)
        if ( allocated(sumDiffuseHeII) ) deallocate(sumDiffuseHeII)
        if (lgDust) then
           if ( allocated(sumDiffuseDust) ) deallocate(sumDiffuseDust)
        end if

    end subroutine setDiffusePDF

    subroutine setDustPDF()
      implicit none

      real :: bb,treal
      real (kind=8) :: tg
      real(kind=8),dimension(nTbins) :: Tspike,Pspike

      integer :: i, n, ai, iT, dcp ! counters

      character(len=50) :: cShapeLoc

      Tspike=0.
      Pspike=0.

      cShapeLoc = 'blackbody'

      dcp = dustComPoint(nspE)

      grids(iG)%dustPDF(grids(iG)%active(ix,iy,iz),:)=0.
      do n = 1, nSpeciesPart(nspE)
         do ai = 1, nSizes
         
            if (grainRadius(ai) >= componentMinRadius(nspE, n) .and. grainRadius(ai) <= componentMaxRadius(nspE, n)) then

               if (grids(iG)%Tdust(n,ai,cellPUsed) >0. .and. &
                    & grids(iG)%Tdust(n,ai,cellPused)<TdustSublime(dcp-1+n)) then


                  if (lgQHeat .and. grainRadius(ai)<minaQHeat .and. &
                       & convPercent>minConvQheat.and. nIterateMC>1) then

                     tg =  grids(iG)%Tdust(n,ai,cellPused)

                     Tspike=0.
                     Pspike=0.

                     call qHeat(n, ai,tg,Tspike,Pspike)

                     do i = 1, nbins
                        do iT = 1, nTbins
                           treal = real(Tspike(iT))
                           bb = getFlux(nuArray(i), treal, cShapeLoc)
                           grids(iG)%dustPDF(cellPUsed, i) = grids(iG)%dustPDF(cellPUsed, i)+ &
                                & real(xSecArray(dustAbsXsecP(dcp-1+n,ai)+i-1)*bb*widFlx(i)*&
                                & grainWeight(ai)*&
                                & grainAbun(nspE, n)*Pspike(iT))
                        end do
                     end do

                  else

                     do i = 1, nbins
                        treal = grids(iG)%Tdust(n,ai,cellPused)
                        bb = getFlux(nuArray(i), treal, cShapeLoc)
                        grids(iG)%dustPDF(cellPused, i) = grids(iG)%dustPDF(cellPused, i)+&
                             & xSecArray(dustAbsXsecP(n+dcp-1,ai)+i-1)*bb*widFlx(i)*&
                             & grainWeight(ai)*grainAbun(nspE, n)
                     end do
                  end if

               end if 
               
            endif   ! end of grain size interval block           
            
         end do
      end do

      ! normalise
      do i = 2, nbins

         grids(iG)%dustPDF(cellPused,i) = &
              & grids(iG)%dustPDF(cellPUsed,i-1)+grids(iG)%dustPDF(cellPUsed,i)

      end do
      grids(iG)%dustPDF(cellPused,:) = &
           grids(iG)%dustPDF(cellPused,:)/grids(iG)%dustPDF(cellPused,nbins)

      grids(iG)%dustPDF(cellPused,nbins) = 1.

    end subroutine setDustPDF


    ! calculates the quantum heating according to
    ! guhathakurta & draine (G&D, 1989) ApJ 345, 230
    ! adapted from a program originally written by
    ! R. Sylvester
    subroutine qHeat(ns,na,tg,temp,pss)
      implicit none

      real ::  tbase
      integer, intent(in) :: ns, na
      integer             :: i,j,ii,if, dcp
      integer             :: wlp,wllp,wlsp,natom

      character(len=1)    :: sorc

      real(kind=8), dimension(nTbins,nTbins) :: A, B
      real(kind=8), intent(inout) :: tg
      real(kind=8), intent(inout) :: temp(nTbins),pss(nTbins)
      real(kind=8)                :: lambda(nbins), sQabs(nbins),&
           & temp1(nbins),temp2(nbins)
      real(kind=8)                :: U(nTbins),dU(nTbins)
      real(kind=8)                :: mgrain,mdUdt,ch
      real(kind=8)                :: const1, radField(nbins),hc,wlim
      real(kind=8)                :: qHeatTemp,X(nTbins)
      real(kind=8)                :: UiTrun,sumX, sumBx
      real(kind=8)                :: wl,wll,wls
      real(kind=8)                :: qint,qint1,qint2,delwav

      integer :: ifreq

      dcp = dustComPoint(nspE)

      ! radiation field at this location in Flambdas
      if (lgDebug) then
         radField = grids(iG)%Jste(cellPUsed,:) + grids(iG)%Jdif(cellPUsed,:)
      else
         radField = grids(iG)%Jste(cellPUsed,:)
      end if

      do ifreq = 1, nbins
         radField(ifreq) = radField(ifreq)/(widflx(ifreq)*fr1ryd)

        temp1(ifreq) = (radField(ifreq)* (nuArray(ifreq)*fr1Ryd)**2)/(c)
         temp2(ifreq) = xSecArray(dustAbsXSecP(ns+dcp-1,na)+ifreq-1)
      end do
      do ifreq=1, nbins
         radField(ifreq) = temp1(nbins-ifreq+1)
         lambda(ifreq) = c/(nuArray(nbins-ifreq+1)*fr1Ryd)
         sQabs(ifreq) = temp2(nbins-ifreq+1)
      end do


      const1 = (hPlanck*fr1Ryd)
      hc = hPlanck*c

      if (grainRadius(na) > minaQHeat) then
         temp(1)=tg
         pss(1) = 1.
         do i = 2, nTbins
            temp(i)=0.
            pss(i) = 0.
         end do
         return
      end if

      temp = 0.
      pss  = 0.

      sorc = grainLabel(ns)(1:1)

      select case(sorc)
      case ('S') ! silicates
         ! factor of 1.e12 to change from um to A
         natom = nint(0.44*1.e12*grainRadius(na)**3)
         if (natom<=0) then
            print*, '! qHeat: invalid value for natom', sorc, na,natom, grainRadius(na)
         end if

         mgrain = (4.*Pi/3.0)*(grainRadius(na)**3)*rho(dcp-1+ns)*1.e-12

      case ('C')  ! carbonaceous
         natom = nint(0.454*3.*1.e12*grainRadius(na)**3)
         if (natom<=0) then
            print*, '! qHeat: invalid value for natom', sorc, na,natom, grainRadius(na)
         end if
         mgrain = (4.*Pi/3.0)*(grainRadius(na)**3)*rho(dcp-1+ns)*1.e-12
      case default
         qheatTemp = tg
         return
      end select

      ! find enthalpy bins
      call getTmin(ns,na,temp(1))

      tbase=real(tg)
      U(1) = enthalpy(temp(1),natom,na,sorc)

      do i=2,10
         temp(i)= temp(1) + real(i-1)/real(9)*(tbase-temp(1))  !big bins
         U(i)=enthalpy(temp(i),natom,na,sorc)
      end do

      do i = 11, nTbins
         temp(i) = tbase + real(i-10)/real(nTbins-10)*(tmax-tbase)
         U(i) = enthalpy(temp(i), natom, na, sorc)
      end do

      dU(1) = 0.5*(U(2)-U(1))
      do i = 2, nTbins-1
         dU(i) = 0.5*(U(i+1)-U(i-1))
      end do
      dU(nTbins) = 0.5*(U(nTbins)-U(nTbins-1))

      ! discrete heating term:
      ! grain heated from state ii to higher state if
      ! by a photon of wavelength hc/( U(if)-U(ii) )
      do ii = 1, nTbins-1
         do if = ii+1, nTbins

            ! calculate wav of transition ii->if
            wl = hc/ (U(if)-U(ii))
            do wlp = nbins, 1, -1
               if (lambda(wlp)<wl) exit
            end do
            wlp = wlp+1

            if ( (dU(if) >= (0.1* (U(if)-U(ii)))) .and. &
                 & (if < nTbins) ) then

               wls = hc/( (U(if)+U(if+1))/2.0 - U(ii)) ! shorter wl
               wll = hc/( (U(if)+U(if-1))/2.0 - U(ii)) ! longer wl

               if ( wl>= 0.1 .or. wl < lambda(1)) then
                  A(if,ii) = 0.
               else

                  if (wll>0.1) wll =.1
                  if (wls<lambda(1)) wls = lambda(1)

                  do wllp = nbins, 1, -1
                     if (lambda(wllp)<wll) exit
                  end do
                  wllp = wllp+1
                  do wlsp = nbins, 1, -1
                     if (lambda(wlsp)<wls) exit
                  end do
                  wlsp = wlsp+1

                  ! here R. Sylvester program interpolates
                  ! our nu grid is finer and maybe this is
                  ! not necessary - test this and possibly change
                  ! if not too expensive

                  ! 1st trapezium
                  A(if,ii) = 0.5*(sQabs(wlsP)*radField(wlsP)&
                       &*wls**3+sQabs(wlP)*radField(wlP)&
                       &*wl**3)*(U(if+1)-U(if))/2.

                  ! 2nd trapezium
                  A(if,ii) = A(if,ii)+0.5*(sQabs(wlP)*radField(wlP)&
                       &*wl**3+sQabs(wllP)*radField(wllP)&
                       &*wll**3)*(U(if)-U(if-1))/2.

                  A(if,ii) = A(if,ii) /(hc**2)

               end if
            else
               if ( (wl>=0.1) .or. (wl<lambda(1))) then
                  A(if,ii) = 0.
               else
                  A(if,ii) = sQabs(wlP)*&
                       & radField(wlP)*dU(if)*wl**3/(hc**2)
               end if

            end if
         end do

         ! heating of enthalpies beyond the highest bin
         ! G&D p232
         UiTrun = U(nTbins)+dU(nTbins)/2.0 - U(ii)
         wlim = hc/UiTrun

         qint = 0.
         do j = 1,nbins-1

            qint1 = lambda(j)*sQabs(j)*radField(j)
            qint2 = lambda(j+1)*sQabs(j+1)*radField(j+1)
            delwav = lambda(j+1) - lambda(j)

            if (lambda(j+1)<wlim) then

               qint = qint + 0.5*(qint1+qint2)*delwav

            else if (lambda(j) < wlim) then
               qint = qint+(wlim-lambda(j))*(qint1+0.5*(wlim-lambda(j))*&
                    & (qint2-qint1)/delwav)
            end if

         end do

         qint = qint/hc

         A(nTbins,ii) = A(nTbins,ii)+qint
      end do

      ! Work out the cooling terms
      do if = 1, nTbins-1

         call clrate(temp(if+1), lambda,radField,sQabs,mdUdt)
         call cheat(if+1,nTbins,lambda,radField,sQabs,U(1:nTbins),ch)


         A(if,if+1) = mdUdt/dU(if+1)

         if (mdUdt <= ch) then
            A(if,if+1) = A(if,if+1)*1.d-20
            A(if+2,if+1) = A(if+2,if+1) + (ch-mdUdt)/dU(if+1)
         else
            A(if,if+1)=A(if,if+1) - ch/dU(if+1)
         end if

      end do


      if (A(2,1) == 0.) A(2,1) = 2.*A(1,2)


      ! find the B(f,j) terms

      B=0.
      do ii = 1, nTbins

         B(nTbins,ii) = A(nTbins,ii)
         do if = nTbins-1, ii+1,-1

            B(if,ii) = B(if+1,ii)+A(if,ii) ! G&D p233
         end do
      end do
      sumX = 0.
      x(1) = 1.d-100
      do if = 2, nTbins
         sumBx = 0.
         do j = 1, if-1

            sumBx = sumBx + B(if,j)*x(j)
         end do

         x(if) = sumBx/A(if-1,if) ! G&D eqn 2.16

         if (x(if) > 1.d250) then
            do i = 1, if
               if (x(i)< 1.) then
                  x(i) = 0.
               else
                  x(i) = x(i)*1.d-250
               end if
            end do
         end if
      end do

      do j=1,nTbins
         sumx = sumx+x(j)
      end do

      do i = 1, nTbins-1
         if (x(i) > 0.) then
            if (dlog10(x(i))-dlog10(sumx) <= -250.) x(i) = 0.
         end if
         Pss(i) = x(i)/sumx
      end do

      qHeatTemp = 0.
      do i = 1, nTbins

         qHeatTemp = qHeatTemp+Pss(i)*temp(i)
      end do

      tg=qHeatTemp

    end  subroutine qHeat


    subroutine getTmin(ns, ai, tmin)
      implicit none

      real(kind=8), intent(out)      :: tmin  ! continuous heating/cooling equilibrium temp
      real                   :: dustAbsIntegral   ! dust absorption integral
      real                   :: cutoff            ! lambda=1000um (guhathakurta & draine 89)
      real, dimension(nbins) :: radField          ! radiation field

      integer, intent(in) :: ns, ai ! grain species and size pointers
      integer :: iT,ifreq ! pointer to dust temp in dust temp array
      integer :: cutoffP ! pointer to energy cutoff in nuarray
      integer :: dcp

      dcp = dustComPoint(nspE)

      cutoff = 9.11e-5

      if (cutoff<nuArray(1)) then
         print*, "! getTmin: [warning] nuMin larger than cutoff for continuous heating &
              & in quantum heating routine"
         cutoff = nuArray(1)
         cutoffP = 1
      else if (cutoff>=nuArray(nbins)) then
         print*, "! getTmin: nuMax smaller than cutoff for continuous heating &
              & in quantum heating routine"
         stop
      else
         call locate(nuArray, cutoff, cutoffP)
      end if

      ! radiation field at this location
      if (lgDebug) then
         radField = grids(iG)%Jste(cellPUsed,:) + grids(iG)%Jdif(cellPUsed,:)
      else
         radField = grids(iG)%Jste(cellPUsed,:)
      end if

      do ifreq = 1, cutoffP
         radField(ifreq)=0.
      end do

      ! zero out dust temperature
      tmin = 0.

      ! calculate absorption integral
      dustAbsIntegral = 0.
      do i = 1, cutoffP
         dustAbsIntegral =  dustAbsIntegral + xSecArray(dustAbsXsecP(dcp-1+nS,ai)+i-1)*radField(i)/Pi
      end do


      call locate(dustEmIntegral(nS,ai,:), dustAbsIntegral, iT)
      if (iT<=0) then
         iT=1
      end if

      tmin = real(iT)

    end subroutine getTmin


          ! calculates the enthalpy for a spherical silicate
          ! or graphite grain of radius a, temperature tg,
          ! containing natom atoms
          ! based on a program originally written by
          ! R. Sylvester using G&D heat capacity
          function enthalpy(tg,natom,na,sorc)
            implicit none

            real(kind=8) :: enthalpy
            real(kind=8), intent(in) :: tg

            real(kind=8) :: enth1,  enth2,  enth3,  enth4,&
                 & enth1_50,  enth2_150,  enth3_500, &
                 & vg

            integer, intent(in) :: natom, na
            character (len=1), intent(in) :: sorc

            enth1  = (1.4e3/3.) * tg**3
            enth1_50 = (1.4e3/3.) * 50.**3
            enth2 = (2.1647e4/2.3) * (tg**2.3 - 50.**2.3)
            enth2_150 = (2.1647e4/2.3) * (150.**2.3 - 50.**2.3)
            enth3 = (4.8375e5/1.68) * (tg**1.68 - 150.**1.68)
            enth3_500 = (4.8375e5/1.68) * (500.**1.68 - 150.**1.68)
            enth4 = 3.3107e7*(tg-500.)

            vg = (4.18878*1.e-12*grainRadius(na)**3)

            select case (sorc)
            case ('S')

               if (tg<=50.) then
                  enthalpy = enth1
               else if (tg <= 150.) then
                  enthalpy = enth1_50 + enth2
               else if (tg <= 500.) then
                  enthalpy = enth1_50 + enth2_150 + enth3
               else
                  enthalpy = enth1_50 + enth2_150 + enth3_500 +&
                       & enth4
               end if
               enthalpy = (1.-(2./real(natom))) * vg *enthalpy

            case('C')
               enthalpy = real(natom)*(4.15e-22*tg**3.3)
               enthalpy = enthalpy/(1.+6.51e-3*tg + 1.5e-6*tg**2 + 8.3e-7*tg**2.3)
               enthalpy = (1.-(2./real(natom))) * enthalpy
            case default
               print*, '! enthalpy: calculations only available for C & S &
                    & identifiers'
            end select

          end function enthalpy


          ! Calculation of the Plank black body equation
          ! lam in cm; output in erg/cm2/s/cm
          ! C1 = 2.Pi.h.c^2
          ! C2 = h.c/k
          function bbody(lam,tIn)

            real(kind=8)            :: bbody
            real(kind=8), intent(in):: lam, tIn
            real(kind=8), parameter :: C1=3.74185e-5
            real(kind=8), parameter :: C2=1.43883
            real(kind=8)            :: exptest

            exptest=C2/(lam*tIn)

            if(exptest > 400.) then
               bbody=0.
            else
               bbody = C1*(1./lam**5)/(exp(exptest)-1.)
            end if
          end function bbody

          subroutine clrate(tem,lam,fl,sQa,mdudt)
            implicit none

            real(kind=8), intent(out) :: mdudt

            real(kind=8), intent(in) :: tem,lam(nbins),fl(nbins),sQa(nbins)

            real(kind=8) :: cutoff, eout(nbins), engin, engout,eomid, &
                 & delwav

            integer :: k,cutoffP

            cutoff = .1 ! cm

            ! calculate -dU/dt
            engin = 0.
            engout = 0.

            do k = 1, nbins

               eout(k) = sQa(k)*bbody(lam(k),tem)
               if (lam(k)<cutoff) cutoffP=k

               if (k>=2) then
                  eomid=(eout(k)+eout(k-1))/2.
                  delwav = lam(k)-lam(k-1)
                  engout = engout+4.*eomid*delwav

               end if
            end do

            engin  = 0.25*fl(cutoffP)*sQa(cutoffP)*cutoff
            if (engin<0.) print*, "! clrate [warning]: engin < 0."

            mdUdt = engout-engin

          end subroutine clrate

          ! continuous heating within a temoerature bin: photon not energetic
          ! enough to raise the grain to the next temp bin
          subroutine cheat(it,nTbins,lam,fl,sQa,U,ch)
            implicit none

            real(kind=8), intent(out)   :: ch
            real(kind=8), intent(in)    :: U(*),lam(nbins),fl(nbins),sQa(nbins)
            real(kind=8) :: einmax,shwav,cutoff, ein(nbins),engin,&
                 & eimid,hc,delwav

            integer, intent(in) :: it,nTbins
            integer :: nshort,nlong,k

            if (it>nTbins) then
               ch=0.
               return
            end if

            hc = hPlanck*c

            engin=0.
            nshort=0
            nlong=0

            cutoff = .1

            if (it>=nTbins) then
               ch = 0.
               return
            end if

            einmax = 0.5*(U(it+1)-(U(it)))
            shwav=hc/einmax

            if(shwav>=cutoff) then
               ch=0.
               return
            end if

            do k = 1, nbins
               if (lam(k)<shwav) nshort = k
               if (lam(k)<cutoff) nlong = k
               ein(k) = sQa(k)*fl(k)
            end do
            nshort = nshort+1
            if (nshort>=nlong) then
               ch = 0.
               return
            end if

            do k = nshort,nlong
               if (k>=2) then
                  eimid = (ein(k)+ein(k-1))/2.
                  delwav = lam(k)-lam(k-1)
                  engin=engin+eimid*delwav
               end if
            end do

            ch = engin

          end subroutine cheat


  end subroutine emissionDriver



  subroutine initResLines(grid)
    implicit none

    type(grid_type), intent(inout) :: grid(*)

      character(len=2)   :: label
      character(len=120) :: reader ! file reader

      integer :: el  ! element number
      integer :: err ! allocation error status
      integer :: i,ii,j ! counters
      integer :: iGrid ! counter
      integer :: ios ! I/O error status
      integer :: nL  ! line counter
      integer :: nmul ! number of multiplets
      integer :: nResLinesFile
      integer :: safeLimit = 1000000 !

      open(unit=19, action="read", file=PREFIX//"/share/mocassin/data/resLines.dat", status="old", position="rewind", iostat=ios)
      if (ios /= 0) then
         print*, "! initResLines: can't open file: ",PREFIX,"/share/mocassin/data/resLines.dat"
         stop
      end if

      nResLines =0
      do i = 1, safeLimit
         read(unit=19,fmt='(A2,A111,I4)',iostat=ios) label,reader, nmul
         if (ios/=0) exit
         do j = 1, nmul
            read(unit=19,fmt=*, iostat=ios)
            if (ios/=0) then
               print*, "! initResLines: error reading ",PREFIX,"/share/mocassin/data/resLines.dat file"
               stop
            end if
         end do

         ! check what element and if it's included
         el=0
         do j = 1, nElements
            if (trim(label) == trim(elemLabel(j)) ) then
               el=j
               exit
            end if
         end do
         if (el==0) then
            print*, "! initResLines: illegal label ", label, el
            stop
         end if
         if (lgElementOn(el)) nResLines=nResLines+1

      end do
      if (nResLines <= 0) then
         print*, '! initResLines: number of resonant lines is zero or negative'
         stop
      end if

      if (i>=safeLimit) then
         print*, "! initResLines : safe limit exceeded in resonant line list read"
         stop
      end if
      nResLinesFile=i-1
      rewind(19)

      ! allocate space for resLineLabel and resLineNuP arrays
      allocate(resLine(1:nResLines), stat = err)
      if (err /= 0) then
         print*, "! initResLines: can't allocate array memory, resLine"
         stop
      end if

      nL=1
      do ii = 1, nResLinesFile
         if (nL>nResLines) exit

         read(unit=19,fmt='(A2,I2,1x,A23,11x,A7,1x,F9.4,7x,F9.6,4x,F12.6,1x,I2,1x,I2,1x,E8.2,1x,E8.2,1x,I2)', &
              &iostat=ios) resLine(nL)%species, resLine(nL)%ion,&
              & resLine(nL)%transition, resLine(nL)%multiplet, resLine(nL)%wav, resLine(nL)%Elow, resLine(nL)%Ehigh,&
              & resLine(nL)%glow, &
              & resLine(nL)%ghigh, resLine(nL)%Aik, resLine(nL)%fik, resLine(nL)%nmul
         if (ios/=0) then
            print*, '! initResline: error reading file -1- '
            stop
         end if

         allocate(resLine(nL)%moclow(1:resLine(nL)%nmul), stat=err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, resLine(nL)%moclow"
            stop
         end if
         allocate(resLine(nL)%mochigh(1:resLine(nL)%nmul), stat=err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, resLine(nL)%mochigh"
            stop
         end if

         do j = 1, resLine(nL)%nmul
            read(unit=19,fmt='(A112,1x,I2,1x,I2)', iostat=ios) reader, &
                 & resLine(nL)%moclow(j), resLine(nL)%mochigh(j)
            if (ios/=0) then
               print*, "! initResLines: error reading ",PREFIX,"/share/mocassin/data/resLines.dat file - 2"
               stop
            end if
         end do

         resLine(nL)%elem=0
         do j = 1, nElements

            if (trim(resLine(nL)%species) == trim(elemLabel(j)))  then
               resLine(nL)%elem=j
               exit
            end if
         end do

         if (resLine(nL)%elem==0) then
            print*, "! initResLines: illegal label -2 -", resLine(nl)%elem, nL
            stop
         end if

         if (lgElementOn(resLine(nL)%elem)) then

            ! convert A to Ryd
            resLine(nL)%freq = c*1.e8/(resLine(nL)%wav*fr1Ryd)

            ! find the mass of the ion in [g]
            resLine(nL)%m_ion = aWeight(resLine(nL)%elem)*amu

            ! locate on nuArary
            call locate(nuArray, resLine(nL)%freq, resLine(nL)%nuP)
            if (resLine(nL)%nuP<=0 .or. resLine(nL)%nuP>=nbins) then
               print*, "! initResLines: wavelength of resonant line is outside grid range", &
                    & ii, resLine(nL)%nuP, resLine(nL)%freq, nuArray(1), nuArray(nbins)
               stop
            end if

            nL = nL +1
         end if

      end do
      close(19)

      if (taskid == 0) then
         print*, "! initResLines: The following resonance lines have been initialised and will &
              &be included in the radiative transfer"
         do i = 1, nResLines
            do j = 1, resLine(i)%nmul
               print*, i, resLine(i)%elem, resLine(i)%ion, resLine(i)%nmul, resLine(i)%wav
               print*, resLine(i)%moclow(j), resLine(i)%mochigh(j)
            end do
         end do
      end if

      do iGrid = 1, nGrids
         ! allocate fEscapeResPhotons
         allocate(grid(iGrid)%fEscapeResPhotons(0:grid(iGrid)%nCells,1:nResLines), stat = err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, fEscapeResPhotons"
            stop
         end if
         grid(iGrid)%fEscapeResPhotons(0:grid(iGrid)%nCells,1:nResLines) = 0.

         ! and the extra packets to be transferred from each location
         allocate(grid(iGrid)%resLinePackets(0:grid(iGrid)%nCells), stat = err)
         if (err /= 0) then
            print*, "! initResLines: can't allocate array memory, resLinePackets"
            stop
         end if
         grid(iGrid)%resLinePackets(0:grid(iGrid)%nCells) = 0

      end do

      allocate(dustHeatingBudget(0:nAbComponents, 0:nResLines), stat=err)
      if (err /= 0) then
         print*, "! initResLines: can't allocate array memory, dustHeatingBudget"
         stop
      end if
      dustHeatingBudget=0.

    end subroutine initResLines

    subroutine setResLineEscapeProb(grids, gPin, xPin,yPin,zPin)
      implicit none

      type(grid_type), intent(inout) :: grids(*) ! incoming grids

      type(vector) :: uHat, rVec ! direction and position vectors

      real    :: a0,DeltaNuD ! line centre xSec [cm^2] and doppler width
      real    :: bFac
      real    :: deltax,deltaxloc    !
      real    :: k_d, k_l  ! dust and line opacities
      real    :: dSx, dSy, dSz, dS ! distances from walls
      real    :: I1, I2    ! calculation integrals
      real    :: Pline     ! prob line
      real    :: radius    ! radial distance from the centre of the grid
      real    :: tau_mu    ! optical depth in direction mu at line centre
      real    :: gasdensity,gastemperature ! properties of the cells along the line of travel
      real    :: gasdensity0,gastemperature0 ! local cell properties
      real :: k_dTest=0.5e-20, HdenTest=0.252, ionFracTest=0.5, TeTest=1.e4


      integer, intent(in) :: xPin,yPin,zPin,gPin ! incoming cell pointers
      integer             :: xp,yp,zp ! cell pointers
      integer             :: gPmother,xPmother,yPmother,zPmother
      integer             :: icount  ! counter
      integer             :: idir, jdir, kdir ! freq counter
      integer, parameter  :: safeLimit=100000 ! loop limit
      integer             :: gP ! subgrid counters
      integer             :: cellP   ! cell pointer on active grid array
      integer             :: compoP  ! pointer to gas component index
      integer             :: iLine   ! res line counter
      integer             :: ndirx, ndiry,ndirz ! number of directions to sample tau

      logical             :: lgTest


      lgTest = .false.

      xp=xPin
      yp=yPin
      zp=zPin
      gp=gPin

      cellP = grids(gPin)%active(xP,yP,zP)

      if (cellP<=0) return

      ndirx = 4
      ndiry = 4
      ndirz = 4

      ! get local physical properties of the gas
      if (lgTest) then
         gasdensity0 = HdenTest
         gastemperature0 = TeTest
         nResLines=1
      else
         gasdensity0 = grids(gPin)%Hden(grids(gPin)%active(xp,yp,zp))
         gastemperature0 = grids(gPin)%Te(grids(gPin)%active(xP,yP,zP))
      end if

      do iLine = 1, nResLines
         if (.not.lgElementOn(resLine(iLine)%elem)) then
            print*, "!setResLineEscapeProb: attempted to calculate escape probability of resonant line &
                 & from an absent ion"
            stop
         end if

         xp=xPin
         yp=yPin
         zp=zPin
         gp=gPin

         cellP = grids(gPin)%active(xP,yP,zP)
         compoP = grids(gPin)%abFileIndex(xP,yP,zP)

         Pline=0.

         ! calculate mean line opacity for this cell
         DeltaNuD = (resLine(iLine)%freq*fr1Ryd/c)*(2.76124e-16*gastemperature0/resLine(iLine)%m_ion)**0.5
!         call meanLineOpacity(resLine(iLine), DeltaNuD, a0)
         call lineCentreXSec(resLine(iLine), DeltaNuD, a0)
         a0 = a0/1.e-14

         ! find k_d and k_l
         if (lgTest) then
            k_d = k_dTest
            k_l = a0* HdenTest * grids(gPin)%elemAbun(compoP,resLine(iLine)%elem)* &
                 & ionFracTest
         else
            k_d = grids(gPin)%absOpac(cellP,resLine(iLine)%nuP)
            k_l = a0* grids(gPin)%Hden(cellP) * grids(gPin)%elemAbun(compoP,resLine(iLine)%elem)* &
                 &  grids(gPin)%ionDen(cellP,elementXref(resLine(iLine)%elem),resLine(iLine)%ion)
         end if

         ! calculate line profile and calculation integrals
         I1=0.
         I2=0.
         deltax=0.

         ! profiles are too narrow for the following

         if (resLine(iLine)%nuP>=nbins) then
            print*, '! setResLineEscapeProb: res line freq pointer >= nbins', iline, &
                 & resLine(iLine)%nuP, nbins
            stop
         end if

         uHat%x = 1./sqrt(3.)
         uHat%y = uHat%x
         uHat%z = uHat%x

         deltax=0.
         do idir = 0,ndirx
            uHat = rotateX(uHat,real(idir)*Pi/2.)
            do jdir = 0,ndiry
               uHat = rotateY(uHat,real(jdir)*Pi/2.)
               do kdir = 0,ndirz
                  uHat = rotateZ(uHat,real(kdir)*Pi/2.)

                  tau_mu = 0.

                  rVec%x = grids(gPin)%xAxis(xPin)
                  rVec%y = grids(gPin)%yAxis(yPin)
                  rVec%z = grids(gPin)%zAxis(zPin)

                  xp = xpin
                  yp = ypin
                  zp = zpin
                  gp = gpin


                  do icount = 1, safeLimit



                     ! check if we are in a subgrid
                     if (grids(gPin)%active(xP,yP,zP)<0) then


                        gPmother = gP
                        xPmother = xP
                        yPmother = yP
                        zPmother = zP

                        gP = abs(grids(gPin)%active(xP,yP,zP))

                        ! where are we in the sub-grid?
                        call locate(grids(gP)%xAxis, rVec%x, xP)
                        if (xP< grids(gP)%nx) then
                           if (rVec%x >=  (grids(gP)%xAxis(xP+1)+grids(gP)%xAxis(xP))/2.) &
                                & xP = xP + 1
                        end if

                        call locate(grids(gP)%yAxis, rVec%y, yP)
                        if (yP< grids(gP)%ny) then
                           if (rVec%y >=  (grids(gP)%yAxis(yP+1)+grids(gP)%yAxis(yP))/2.) &
                                & yP = yP + 1
                        end if
                        call locate(grids(gP)%zAxis, rVec%z, zP)
                        if (zP< grids(gP)%nz) then
                           if (rVec%z >=  (grids(gP)%zAxis(zP+1)+grids(gP)%zAxis(zP))/2.) &
                                & zP = zP + 1
                        end if

                     end if

                     cellP = grids(gPin)%active(xP,yP,zP)

                     if (cellP > 0 ) then
                        compoP = grids(gPin)%abFileIndex(xP,yP,zP)
                        if (xp == xpin .and. yp == ypin .and. zp == zpin .and. gP == gpin) then
                           gasdensity = gasdensity0
                           gastemperature=gastemperature0
                        else
                           if (lgTest) then
                              gasdensity = HdenTest
                              gastemperature = TeTest
                           else
                              gasdensity = grids(gP)%Hden(grids(gP)%active(xp,yp,zp))
                              gastemperature = grids(gP)%Te(grids(gP)%active(xP,yP,zP))
                           end if
                        end if

                     else

                        gasdensity=0.
                        gastemperature=1.

                     end if

                     ! find distances from all walls

                     if (uHat%x>0.) then
                        if (xP<grids(gP)%nx) then
                           dSx = ( (grids(gP)%xAxis(xP+1)+grids(gP)%xAxis(xP))/2.-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then
                              rVec%x=(grids(gP)%xAxis(xP+1)+grids(gP)%xAxis(xP))/2.
                              xP = xP+1
                           end if
                        else
                           dSx = ( grids(gP)%xAxis(grids(gP)%nx)-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then
                              rVec%x=grids(gP)%xAxis(grids(gP)%nx)
                           end if
                        end if
                     else if (uHat%x<0.) then
                        if (xP>1) then
                           dSx = ( (grids(gP)%xAxis(xP)+grids(gP)%xAxis(xP-1))/2.-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then
                              rVec%x=(grids(gP)%xAxis(xP)+grids(gP)%xAxis(xP-1))/2.
                              xP = xP-1
                           end if
                        else
                           dSx = (grids(gP)%xAxis(1)-rVec%x)/uHat%x
                           if (abs(dSx)<1.e-10) then
                              rVec%x=grids(gP)%xAxis(1)
                           end if
                        end if
                     else if (uHat%x==0.) then
                        dSx = grids(gP)%xAxis(grids(gP)%nx)
                     end if

                     if (.not.lg1D) then
                        if (uHat%y>0.) then
                           if (yP<grids(gP)%ny) then
                              dSy = ( (grids(gP)%yAxis(yP+1)+grids(gP)%yAxis(yP))/2.-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then
                                 rVec%y=(grids(gP)%yAxis(yP+1)+grids(gP)%yAxis(yP))/2.
                                 yP = yP+1
                              end if
                           else
                              dSy = (  grids(gP)%yAxis(grids(gP)%ny)-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then
                                 rVec%y=grids(gP)%yAxis(grids(gP)%ny)
                              end if
                           end if
                        else if (uHat%y<0.) then
                           if (yP>1) then
                              dSy = ( (grids(gP)%yAxis(yP)+grids(gP)%yAxis(yP-1))/2.-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then
                                 rVec%y=(grids(gP)%yAxis(yP)+grids(gP)%yAxis(yP-1))/2.
                                 yP = yP-1
                              end if
                           else
                              dSy = ( grids(gP)%yAxis(1)-rVec%y)/uHat%y
                              if (abs(dSy)<1.e-10) then
                                 rVec%y=grids(gP)%yAxis(1)
                              end if
                           end if
                        else if (uHat%y==0.) then
                           dSy = grids(gP)%yAxis(grids(gP)%ny)
                        end if

                        if (uHat%z>0.) then
                           if (zP<grids(gP)%nz) then
                              dSz = ( (grids(gP)%zAxis(zP+1)+grids(gP)%zAxis(zP))/2.-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then
                                 rVec%z=(grids(gP)%zAxis(zP+1)+grids(gP)%zAxis(zP))/2.
                                 zP = zP+1
                              end if
                           else
                              dSz = ( grids(gP)%zAxis(grids(gP)%nz)-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then
                                 rVec%z=grids(gP)%zAxis(grids(gP)%nz)
                              end if
                           end if
                        else if (uHat%z<0.) then
                           if (zP>1) then
                              dSz = ( (grids(gP)%zAxis(zP)+grids(gP)%zAxis(zP-1))/2.-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then
                                 rVec%z=(grids(gP)%zAxis(zP)+grids(gP)%zAxis(zP-1))/2.
                                 zP = zP-1
                              end if
                           else
                              dSz = ( grids(gP)%zAxis(1)-rVec%z)/uHat%z
                              if (abs(dSz)<1.e-10) then
                                 rVec%z=grids(gP)%zAxis(1)
                              end if
                           end if
                        else if (uHat%z==0.) then
                           dSz = grids(gP)%zAxis(grids(gP)%nz)
                        end if

                     else

                        print*, '! setResLineEscapeProb: oneD option not yet implemented for res line transfer routine'
                        stop

                     end if


                     ! cater for cells on cell wall
                     if ( abs(dSx)<1.e-10 ) dSx = grids(gP)%xAxis(grids(gP)%nx)
                     if ( abs(dSy)<1.e-10 ) dSy = grids(gP)%yAxis(grids(gP)%ny)
                     if ( abs(dSz)<1.e-10 ) dSz = grids(gP)%zAxis(grids(gP)%nz)

                     ! find the nearest wall
                     dSx = abs(dSx)
                     dSy = abs(dSy)
                     dSz = abs(dSz)

                     if (dSx<=0.) then
                        print*, '! setResLineEscapeProb: [warning] dSx <= 0.'
                        dS = min(dSy, dSz)
                     else if (dSy<=0.) then
                        print*, '! setResLineEscapeProb: [warning] dSy <= 0.'
                        dS = min(dSx, dSz)
                     else if (dSz<=0.) then
                        print*, '! setResLineEscapeProb: [warning] dSz <= 0.'
                        dS = min(dSx, dSy)
                     else
                        dS = min(dSx,dSy,dSz)
                     end if

                     ! this should now never ever happen
                     if (dS <= 0.) then
                        print*, 'setResLineEscapeProb: dS <= 0'
                        stop
                     end if

                     ! update optical depth
                     if (cellP>0) then
                        if (lgTest) then
                           tau_mu  = tau_mu + a0*1.e-14*dS*gasdensity*grids(gP)%elemAbun(compoP,resLine(iLine)%elem)*&
                                &ionFracTest
                        else
                           tau_mu  = tau_mu + a0*1.e-14*dS*gasdensity*grids(gP)%elemAbun(compoP,resLine(iLine)%elem)*&
                                &grids(gP)%ionDen(cellP,elementXref(resLine(iLine)%elem),resLine(iLine)%ion)
                        end if
                     end if

                     ! update position
                     rVec = rVec + dS*uHat

                     ! keep track of where you are on mother grid
                     if (gP>1) then
                        if (uHat%x>0.) then
                           if ( xPmother < grids(gPmother)%nx ) then
                              if ( rVec%x >= (grids(gPmother)%xAxis(xPmother)+&
                                   & grids(gPmother)%xAxis(xPmother+1))/2. ) then
                                 xPmother = xPmother+1
                              end if
                           else
                              if ( rVec%x > grids(gPmother)%xAxis(xPmother)) then
                                 print*, '! pathSegment: insanity occurred at mother grid transfer (x axis +)', &
                                      & rVec%x, gP, gPmother
                                 stop
                              end if
                           end if
                        else
                           if ( xPmother > 1 ) then
                              if ( rVec%x <= (grids(gPmother)%xAxis(xPmother-1)+&
                                   & grids(gPmother)%xAxis(xPmother))/2. ) then
                                 xPmother = xPmother-1
                              end if
                           else
                              if (rVec%x < grids(gPmother)%xAxis(1)) then
                                 print*, '! pathSegment: insanity occurred at mother grids transfer (x axis -)', &
                                      & rVec%x, gP, gPmother
                                 stop
                              end if
                           end if
                        end if
                        if (uHat%y>0.) then
                           if (  yPmother < grids(gPmother)%ny ) then
                              if ( rVec%y >= (grids(gPmother)%yAxis( yPmother)+&
                                   & grids(gPmother)%yAxis( yPmother+1))/2. ) then
                                 yPmother =  yPmother+1
                              end if
                           else
                              if ( rVec%y > grids(gPmother)%yAxis( yPmother)) then
                                 print*, '! pathSegment: insanity occurred at mother grid transfer (y axis +)', &
                                      & rVec%y, gP, gPmother
                                 stop
                              end if
                           end if
                        else
                           if (  yPmother > 1 ) then
                              if ( rVec%y <= (grids(gPmother)%yAxis( yPmother-1)+&
                                   & grids(gPmother)%yAxis( yPmother))/2. ) then
                                 yPmother =  yPmother-1
                              end if
                           else
                              if (rVec%y < grids(gPmother)%yAxis(1)) then
                                 print*, '! pathSegment: insanity occurred at mother grid transfer (y axis -)', &
                                      & rVec%y, gP, gPmother
                                 stop
                              end if
                           end if
                        end if
                        if (uHat%z>0.) then
                           if (  zPmother < grids(gPmother)%nz ) then
                              if ( rVec%z >= (grids(gPmother)%zAxis( zPmother)+&
                                   & grids(gPmother)%zAxis( zPmother+1))/2. ) then
                                 zPmother =  zPmother+1
                              end if
                           else
                              if ( rVec%z > grids(gPmother)%zAxis( zPmother)) then
                                 print*, '! pathSegment: insanity occurred at mother grid transfer (z axis +)', &
                                      & rVec%z, gP, gPmother
                                 stop
                              end if
                           end if
                        else
                           if (  zPmother > 1 ) then
                              if ( rVec%z <= (grids(gPmother)%zAxis( zPmother-1)+&
                                   & grids(gPmother)%zAxis( zPmother))/2. ) then
                                 zPmother =  zPmother-1
                              end if
                           else
                              if (rVec%z < grids(gPmother)%zAxis(1)) then
                                 print*, '! pathSegment: insanity occurred at mother grid transfer (z axis -)', &
                                      & rVec%z, gP, gPmother
                                 stop
                              end if
                           end if
                        end if
                     end if

                     ! and update current grid indeces
                     if (.not.lg1D) then
                        if ( (dS == dSx) .and. (uHat%x > 0.)  ) then
                           xP = xP+1
                        else if ( (dS == dSx) .and. (uHat%x < 0.) ) then
                           xP = xP-1
                        else if ( (dS == dSy) .and. (uHat%y > 0.) ) then
                           yP = yP+1
                        else if ( (dS == dSy) .and. (uHat%y < 0.) ) then
                           yP = yP-1
                        else if ( (dS == dSz) .and. (uHat%z > 0.) ) then
                           zP = zP+1
                        else if ( (dS == dSz) .and. (uHat%z < 0.) ) then
                           zP = zP-1
                        end if
                     else

                        print*, '! setResLineEscapeProb: oneD option not yet implemented for res line transfer routine'
                        stop

                        radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                             & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                             & (rVec%z/1.e10)*(rVec%z/1.e10))
                        call locate(grids(gP)%xAxis, radius , xP)

                     end if

                     if(lgPlaneIonization) then

                        if ( (rVec%y <= grids(gP)%yAxis(1) .or. yP<1) .and. uHat%y<0.) then


                           if (gP==1) then
                              ! exit the loop if on mother or return to mother if on sub
                              yP = nint(grids(gP)%yAxis(1))
                              exit
!                              yP=1
!                              rVec%y = grids(gP)%yAxis(1)
!                              uHat%y = -uHat%y
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
!                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP', gP
                              stop
                           end if
                        end if


                        if (rVec%y >= grids(gP)%yAxis(grids(gP)%ny)+grids(gP)%geoCorrY .or. yP>grids(gP)%ny) then

                           if (gP==1) then
                              ! exit the loop if on mother or return to mother if on sub
                              yP = grids(gP)%ny
                              exit
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP', gP
                              stop
                           end if

                        end if

                        if (rVec%x <= grids(gP)%xAxis(1) .or. xP<1) then

                           if (gP==1) then
!                         rVec%x = -rVec%x
                              xP=1
                              rVec%x = grids(gP)%xAxis(1)
                              uHat%x = -uHat%x
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP', gP
                              stop
                           end if
                        end if

                        if (rVec%x >=  grids(gP)%xAxis(grids(gP)%nx) .or. xP>grids(gP)%nx)then

                           if (gP==1) then
!                         rVec%x = 2.*grid(gP)%xAxis(grid(gP)%nx)- rVec%x
                              xP = grids(gP)%nx
                              rVec%x = grids(gP)%xAxis(grids(gP)%nx)
                              uHat%x = -uHat%x
                           else if (gP>1) then
                              xP = xPmother
                              yP =  yPmother
                              zP =  zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP', gP
                              stop
                           end if

                        end if

                        if (rVec%z <= grids(gP)%zAxis(1) .or.zP<1) then

                           if (gP==1) then
!                         rVec%z = -rVec%z
                              zP=1
                              rVec%z = grids(gP)%zAxis(1)
                              uHat%z = -uHat%z
                           else if (gP>1) then
                              xP = xPmother
                              yP =  yPmother
                              zP =  zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP', gP
                              stop
                           end if

                        end if

                        if (rVec%z >=  grids(gP)%zAxis(grids(gP)%nz) .or. zP>grids(gP)%nz)then

                           if (gP==1) then
                              zP = grids(gP)%nz
!                         rVec%z = 2.*grid(gP)%zAxis(grid(gP)%nz)- rVec%z
                              rVec%z = grids(gP)%zAxis(grids(gP)%nz)
                              uHat%z = -uHat%z
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP', gP
                              stop
                           end if

                        end if

                     end if

                     ! check if the path is still within the ionized region
                     radius = 1.e10*sqrt((rVec%x/1.e10)*(rVec%x/1.e10) + &
                          & (rVec%y/1.e10)*(rVec%y/1.e10) + &
                          & (rVec%z/1.e10)*(rVec%z/1.e10))

                     if ( radius >= R_out .and. R_out >= 0. ) exit

                     if (lgSymmetricXYZ .and. gP == 1) then
                        if (lgPlaneIonization) then
                           print*, '! setResLineEscapeProb: lgSymmetric and lgPlaneionization flags both raised'
                           stop
                        end if
                        if ( rVec%x <= grids(gP)%xAxis(1) .or. xP<1) then
                           if (uHat%x<0.) uHat%x = -uHat%x
                           xP = 1
                           rVec%x = grids(gP)%xAxis(1)
                        end if
                        if ( rVec%y <= grids(gP)%yAxis(1) .or. yP<1) then
                           if (uHat%y<0.) uHat%y = -uHat%y
                           yP = 1
                           rVec%y = grids(gP)%yAxis(1)
                        end if
                        if ( rVec%z <= grids(gP)%zAxis(1) .or. zP<1) then
                           if (uHat%z<0.) uHat%z = -uHat%z
                           zP=1
                           rVec%z = grids(gP)%zAxis(1)
                        end if
                     end if


                     if (.not.lgPlaneIonization) then

                        if ( (abs(rVec%x) >= grids(gP)%xAxis(grids(gP)%nx)+grids(gP)%geoCorrX) .or.&
                             &(abs(rVec%y) >= grids(gP)%yAxis(grids(gP)%ny)+grids(gP)%geoCorrY) .or.&
                             &(abs(rVec%z) >= grids(gP)%zAxis(grids(gP)%nz)+grids(gP)%geoCorrZ) .or. &
                             & xP>grids(gP)%nx .or. yP>grids(gP)%ny .or. zP>grids(gP)%nz  ) then


                           if ((gP==1) .or.  (radius >= R_out .and. R_out >= 0.)) then
                              if (xP > grids(gP)%nx) xP = grids(gP)%nx
                              if (yP > grids(gP)%ny) yP = grids(gP)%ny
                              if (zP > grids(gP)%nz) zP = grids(gP)%nz
                              exit
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP - ', gP
                              stop
                           end if

                        end if
                        if ( .not. lgSymmetricXYZ .and. ( rVec%x < grids(gP)%xAxis(1) .or.&
                             & rVec%y < grids(gP)%yAxis(1) .or.&
                             & rVec%z < grids(gP)%zAxis(1) .or. &
                             & xP<1.or. yP<1 .or. zP<1 )) then

                           if (gP==1) then
                              if (xP < 1) xP = 1
                              if (yP < 1) yP = 1
                              if (zP < 1) zP = 1
                              exit
                           else if (gP>1) then
                              xP = xPmother
                              yP = yPmother
                              zP = zPmother
                              gP = gPmother
                           else
                              print*, '! setResLineEscapeProb: insanity occurred - invalid gP - ', gP
                              stop
                           end if

                        end if
                     end if


                     if (gP>1) then

                        if ( ( (rVec%x <= grids(gP)%xAxis(1) .or. xP<1) .and. uHat%x <=0.) .or. &
                             & ( (rVec%y <= grids(gP)%yAxis(1) .or. yP<1) .and. uHat%y <=0.) .or. &
                             & ( (rVec%z <= grids(gP)%zAxis(1) .or. zP<1) .and. uHat%z <=0.) ) then

                           ! go back to mother grid
                           xP = xPmother
                           yP = yPmother
                           zP = zPmother
                           gP = gPmother

                        end if

                     end if
                  end do

                  if (iCount>=10000) then
                     print*, '! setResLineEscapeProb: safeLimit exceeded in tau loop', icount, safelimit
                  end if

                  if (tau_mu <= 0.) then
                     tau_mu = 0.
                     if (lgTalk) &
                          & print*, '! setResLineEscapeProb: [warning] tau_mu <= 0. ', &
                          & tau_mu, uHat, xPin,yPin,zPin, cellP, xp,yp,zp
                  end if


                  ! calculate direction-dependent escape probability according to Kwan & Krolik 1981
                  if (tau_mu>=1.) then
                     bFac = sqrt(log(tau_mu))/(1.+tau_mu/1.e5)
                     Pline = 1./(tau_mu*sqrt(Pi)*(1.2+bFac))
                  else
                     Pline=(1.-exp(-2.*tau_mu))/(2.*tau_mu)
                  end if

                  if (tau_mu <= 0.) then
                     deltaxloc=0.
                  else
                     deltaxloc = log(tau_mu/sqrt(Pi))
                     if (deltaxloc >0.) then
                        deltaxloc=sqrt(deltaxloc)
                     else
                        deltaxloc=0.
                     end if
                  end if

                  ! calculate fEscapeResPhotons according to Harrington, Cohen and Hess (1984) eqn A5
                  grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) =&
                       &grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine)+ &
                       & Pline
                  deltax = deltax+deltaxloc

               end do
            end do
         end do

         grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) = &
              & grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine)/&
              & real((ndirx+1)*(ndiry+1)*(ndirz+1))

         deltax=deltax/real((ndirx+1)*(ndiry+1)*(ndirz+1))

         if (grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine)>0.) then
            grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) = &
                 & 1./(1. + k_d*1.e+14*2.*deltax/(k_l*&
                 & grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine)))
         else
            grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),iLine) = 0.
         end if

      end do


      if (lgTest) then
         print*, xPin, yPin, zPin,grids(gPin)%fEscapeResPhotons(grids(gPin)%active(xPin,yPin,zPin),1)
      end if

    end subroutine setResLineEscapeProb

    ! calculate direction-dependent escape probability according to Harrington, Cohen and Hess (1984) eqn A1
    function getM2_tau(width,lineIn,tauIn, profilein)
      implicit none

      type(resLine_type), intent(in) :: lineIn

      real, intent(in)              :: width, tauIn
      real, intent(in), optional    :: profileIn(nbins)
      real                          :: getM2_tau
      real                          :: dx

      integer                       :: i

      if (tauIn<0.) then
         print*, "! getM2_tau: negative tau", tauIn
         stop
      end if

      if (present(profilein)) then
         getM2_tau=0.
         do i = 1, nbins
            dx = widFlx(lineIn%nuP)*fr1ryd/(real(nbins-1)*width)
            getM2_tau = getM2_tau+profileIn(i)*exp(-tauIn*profileIn(i))*dx
         end do
      else
         getM2_tau = exp(-tauIn)
      end if

    end function getM2_tau

    function profile(profileID, widthIn, lineIn, nuIn)
      implicit none

      type(resLine_type), intent(in) :: lineIn
      real              :: profile ! phi(nu)
      real, intent(in)  :: nuIn ! frequency [Ryd]
      real, intent(in)  :: widthIn !

      character(len=7), intent(in) :: profileID ! what type of broadening?

      select case (profileID)
         case('doppler')
            profile = (1./(sqrt(Pi)*widthIn) ) * Exp(- ((nuIn-lineIn%freq)*fr1Ryd/widthIn)**2)
         case default
            print*, '! profile: unknown profile requested ', profileID
            stop
      end select

    end function profile

    ! use line centre xsec formula
    ! a0 = (Pi e^2 / (m_e c)) f /( (Pi^1/2 DeltaNuD))
    ! DeltaNuD = Doppler width = Nu0 /c * (2 k T /m_ion)^1/2
    !
    subroutine lineCentreXSec(line, width, line_xSec)
      implicit none

      type(resLine_type), intent(inout) :: line

      real, intent(in)         :: width
      real, intent(out)        :: line_xSec


      if (width<=0.) then
         print*, '! lineCentreXSec: width less or equal to zero'
         stop
      end if
      line_xSec = 14.969e-3*line%fik/width

    end subroutine lineCentreXSec


    ! use hummer & kunasz eqn 2.3 : k_l = N_k B_kl h v_kl/(4 Pi DeltaNuP)
    ! and B_kl = e^2 f_kl / (4 E_0 m_e c h n_kl)
    ! hence k_l = N_k e^2 f_kl / (4 E_0 m_e c 4 Pi DeltaNuP)
    ! E_0 is permissivity of free space
    ! a0 = (Pi e^2 / (m_e c)) f /( (Pi^1/2 DeltaNuD))
    ! DeltaNuD = Doppler width = Nu0 /c * (2 k T /m_ion)^1/2
    !
    subroutine meanLineOpacity(line, width, line_xSec)
      implicit none

      type(resLine_type), intent(inout) :: line

      real, intent(in)         :: width
      real, intent(out)        :: line_xSec


      if (width<=0.) then
         print*, '! lineCentreXSec: width less or equal to zero'
         stop
      end if
      line_xSec = 1.681e-4*line%fik/width

    end subroutine meanLineOpacity

  end module emission_mod


