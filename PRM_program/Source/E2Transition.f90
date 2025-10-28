subroutine ETransition

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters
    use CGParameters

    implicit none

    TFileName = 'WaveFunction/TWaveFunction'//achar(int(EIniState+64))//'.dat'
    open(unit=Toutput,file=TFileName)
        do bra=1, DimensionArray(int(EIniState))
            read(Toutput,'(100000f12.8)') EIniWF(bra,:)
        end do
    close(unit=Toutput)

    TFileName = 'WaveFunction/TWaveFunction'//achar(int(EFinState+64))//'.dat'
    open(unit=Toutput,file=TFileName)
        do bra=1, DimensionArray(int(EFinState))
            read(Toutput,'(1000000f12.8)') EFinWF(bra,:)
        end do
    close(unit=Toutput)

    SCFileName = 'SCSerialNum/SCSerialNum'//achar(int(EIniState+64))//'.dat'
    open(unit=Toutput,file=SCFileName)
        do bra=1, DimensionArray(int(EIniState))
            read(Toutput,*) EIniNSCVector(:,bra)
        end do
    close(unit=Toutput)

    SCFileName = 'SCSerialNum/SCSerialNum'//achar(int(EFinState+64))//'.dat'
    open(unit=Toutput,file=SCFileName) 
        do bra=1, DimensionArray(int(EFinState))  
            read(Toutput,*) EFinNSCVector(:,bra)
        end do
    close(unit=Toutput)

    SCFileName = 'ISCSerialNum/ISCSerialNum'//achar(int(EIniState+64))//'.dat'
    open(unit=Toutput,file=SCFileName)
        do bra=1, DimensionArray(int(EIniState))
            read(Toutput,*) EIniISCVector(:,bra)
        end do
    close(unit=Toutput)

    SCFileName = 'ISCSerialNum/ISCSerialNum'//achar(int(EFinState+64))//'.dat'
    open(unit=Toutput,file=SCFileName)
        do bra=1, DimensionArray(int(EFinState))
            read(Toutput,*) EFinISCVector(:,bra)
        end do
    close(unit=Toutput)

    BEResArray = 0
    do bra=1, DimensionArray(int(EFinState))
    do ket=1, DimensionArray(int(EIniState))

        if (EFinNSCVector(2,bra)/=EIniNSCVector(2,ket)) cycle

        FinK = EFinNSCVector(1,bra)
        IniK = EIniNSCVector(1,ket)
        
        call CGFactor(EIniState,dble(2),EFinState,IniK,dble(0),FinK, CG20)

        call CGFactor(EIniState,dble(2),EFinState,IniK,dble(2),FinK, CG2P2)

        call CGFactor(EIniState,dble(2),EFinState,IniK,dble(-2),FinK, CG2N2)

        BEResArray(bra,ket) = EFinWF(bra,1)*EIniWF(ket,1) &
        *(cos(gamma_rad)*CG20 + (sin(gamma_rad)*(CG2P2+CG2N2)/(sqrt(2.))))
    
    end do
    end do

    do bra=1, DimensionArray(int(EFinState))
    do ket=1, DimensionArray(int(EIniState))

        if (EFinNSCVector(2,bra)/=EIniISCVector(2,ket)) cycle

        FinK = EFinNSCVector(1,bra)
        IniK = EIniNSCVector(1,ket)

        call CGFactor(EIniState,dble(2),EFinState,-IniK,dble(0),FinK, ICG20)

        call CGFactor(EIniState,dble(2),EFinState,-IniK,dble(2),FinK, ICG2P2)

        call CGFactor(EIniState,dble(2),EFinState,-IniK,dble(-2),FinK, ICG2N2)

        BEResArray(bra,ket) = EFinWF(bra,1)*EIniWF(ket,1) &
        *(((-1)**(EIniState-IniK))*(cos(gamma_rad)*ICG20 + sin(gamma_rad)*(ICG2P2+ICG2N2)/(sqrt(2.))))&
        + BEResArray(bra,ket)

    end do
    end do

    BERes = 5*((Q_intr*sum(BEResArray))**2)/(16*pi)

    write(*,*) BERes

end subroutine ETransition