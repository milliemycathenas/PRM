subroutine ValenceEnergy

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    real :: NJKp1, NJKn1, NJK0
    real :: IJKp1, IJKn1, IJK0

    VHamiltonian = 0

    VHamiltonianloop: &
    do bra=1, DimensionNumber
    do ket=bra, DimensionNumber

        call VInitialise

        NJKp1 = 0
        NJKn1 = 0
        NJK0 = 0

        IJKp1 = 0
        IJKn1 = 0
        IJK0 = 0

        call KEquiv
        call KMinus
        call KPlus

        KBra=NSCVector(1,bra)
        KKet=NSCVector(1,ket)

        if (KBra==KKet+1) NJKp1=sqrt(I*(I+1)-KKet*(KKet+1))
        if (KBra==KKet-1) NJKn1=sqrt(I*(I+1)-KKet*(KKet-1))
        if (KBra==KKet) NJK0=KKet

        if (KBra==-KKet+1) IJKp1=sqrt(I*(I+1)+KKet*(-KKet+1))
        if (KBra==-KKet-1) IJKn1=sqrt(I*(I+1)+KKet*(-KKet-1))
        if (KBra==-KKet) IJK0=-KKet

        CoeffNtoI = 0
        do tau3=1, tau3_max
            CoeffNtoI = CoeffNtoI + PNumI(tau3,NSCVector(2,ket))&
            + PNumN(tau3,NSCVector(2,ket))*PNumI(tau3,NSCVector(2,ket))
        end do
        CoeffNtoI = CoeffNtoI + ( I-NSCVector(1,ket) )

        VHamiltonian(bra,ket) = 0.25*A(1)*(NJSMono+NJSPMono+NJSTB) &
        + 0.25*A(2)*(NJuJdMono+NJuJdPMono+NJuJdTB+NJdJuMono+NJdJuPMono+NJdJuTB) &
        + A(3)*(NJ3SMono+NJ3SPMono+NJ3STB) - 2*CoriolisF*A(3)*NJK0*NJ3Mono &
        - 0.5*CoriolisF*A(1)*NJKp1*NJdMono - 0.5*CoriolisF*A(2)*NJKp1*NJuMono &
        - 0.5*CoriolisF*A(1)*NJKn1*NJuMono - 0.5*CoriolisF*A(2)*NJKn1*NJdMono &
        + (-1)**CoeffNtoI &
        *(0.25*A(1)*(IJSMono+IJSPMono+IJSTB) &
        + 0.25*A(2)*(IJuJdMono+IJuJdPMono+IJuJdTB+IJdJuMono+IJdJuPMono+IJdJuTB) &
        + A(3)*(IJ3SMono+IJ3SPMono+IJ3STB) - 2*CoriolisF*A(3)*IJK0*IJ3Mono &
        - 0.5*CoriolisF*A(1)*IJKp1*IJdMono - 0.5*CoriolisF*A(2)*IJKp1*IJuMono &
        - 0.5*CoriolisF*A(1)*IJKn1*IJuMono - 0.5*CoriolisF*A(2)*IJKn1*IJdMono)

        ! SignPass=- 2*CoriolisF*A(3)*NJK0*NJ3Mono &
        ! - 0.5*CoriolisF*A(1)*NJKp1*NJdMono - 0.5*CoriolisF*A(2)*NJKp1*NJuMono &
        ! - 0.5*CoriolisF*A(1)*NJKn1*NJuMono - 0.5*CoriolisF*A(2)*NJKn1*NJdMono 
        
        ! PassMark= (-1)**CoeffNtoI &
        ! *(- 2*CoriolisF*A(3)*IJK0*IJ3Mono &
        ! - 0.5*CoriolisF*A(1)*IJKp1*IJdMono - 0.5*CoriolisF*A(2)*IJKp1*IJuMono &
        ! - 0.5*CoriolisF*A(1)*IJKn1*IJuMono - 0.5*CoriolisF*A(2)*IJKn1*IJdMono)

        ! SignPass=0.25*A(1)*(NJSMono+NJSPMono+NJSTB) &
        ! + 0.25*A(2)*(NJuJdMono+NJuJdPMono+NJuJdTB+NJdJuMono+NJdJuPMono+NJdJuTB) &
        ! + A(3)*(NJ3SMono+NJ3SPMono+NJ3STB)

        ! PassMark= (-1)**CoeffNtoI &
        ! *(0.25*A(1)*(IJSMono+IJSPMono+IJSTB) &
        ! + 0.25*A(2)*(IJuJdMono+IJuJdPMono+IJuJdTB+IJdJuMono+IJdJuPMono+IJdJuTB) &
        ! + A(3)*(IJ3SMono+IJ3SPMono+IJ3STB))

        ! if (SignPass/=0) write(*,*) SignPass
        ! if (PassMark/=0) write(*,*) PassMark

    end do
    end do VHamiltonianloop

    open(unit=Voutput, file='VHamiltonian.dat')
        ! do bra=1, DimensionNumber
        !     write(Voutput,'(100000f12.6)') VHamiltonian(bra,:)
        ! end do
        do bra=1, DimensionNumber
        do ket=1, DimensionNumber
            if (VHamiltonian(bra,ket)/=0) write(Voutput,*) VHamiltonian(bra,ket)
        end do
        end do
    close(unit=Voutput)

end subroutine ValenceEnergy

