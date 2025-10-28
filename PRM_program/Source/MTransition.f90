subroutine MTransition

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters
    use CGParameters

    implicit none

    MFinNSCVector=0
    MIniNSCVector=0

    TFileName = 'WaveFunction/TWaveFunction'//achar(int(MIniState+64))//'.dat'
    open(unit=Toutput,file=TFileName)
        do bra=1, DimensionArray(int(MIniState))
            read(Toutput,'(100000f12.8)') MIniWF(bra,:)
        end do
    close(unit=Toutput)

    TFileName = 'WaveFunction/TWaveFunction'//achar(int(MFinState+64))//'.dat'
    open(unit=Toutput,file=TFileName)
        do bra=1, DimensionArray(int(MFinState))
            read(Toutput,'(1000000f12.8)') MFinWF(bra,:)
        end do
    close(unit=Toutput)

    SCFileName = 'SCSerialNum/SCSerialNum'//achar(int(MIniState+64))//'.dat'
    open(unit=Toutput,file=SCFileName)
        do bra=1, DimensionArray(int(MIniState))
            read(Toutput,*) MIniNSCVector(:,bra)
        end do
    close(unit=Toutput)

    SCFileName = 'SCSerialNum/SCSerialNum'//achar(int(MFinState+64))//'.dat'
    open(unit=Toutput,file=SCFileName) 
        do bra=1, DimensionArray(int(MFinState))  
            read(Toutput,*) MFinNSCVector(:,bra)
        end do
    close(unit=Toutput)

    ! do bra=1, DimensionArray(int(MFinState))
    !     if (PNumI(2,MFinNSCVector(2,bra))==1) then
    !         write (*,'(10f12.1)') MFinNSCVector(1,bra)
    !         write (*,'(100I2)') SPVector(1,1:countsp_max(1),MFinNSCVector(2,bra)), SPVector(2,1:countsp_max(2),MFinNSCVector(2,bra))
    !     end if
    ! end do

    ! do bra=1, DimensionArray(int(MIniState))
    !     if (PNumI(2,MIniNSCVector(2,bra))==1) then
    !         write (*,'(10f12.1)') MIniNSCVector(1,bra)
    !         ! write (*,'(100I2)') SPVector(1,1:countsp_max(1),MIniNSCVector(2,bra)), SPVector(2,1:countsp_max(2),MIniNSCVector(2,bra))
    !     end if
    ! end do

    ! stop

    BMResArray = 0
    do bra=1, DimensionArray(int(MFinState))
    do ket=1, DimensionArray(int(MIniState))

        FinK = MFinNSCVector(1,bra)
        IniK = MIniNSCVector(1,ket)

        AJplus = 0
        AJz = 0
        AJminus = 0

        BJplus = 0
        BJz = 0
        BJminus = 0

        CJplus = 0
        CJz = 0
        CJminus = 0

        DJplus = 0
        DJz = 0
        DJminus = 0

        CGAplus = 0
        CGAz = 0
        CGAminus = 0

        CGBplus = 0
        CGBz = 0
        CGBminus = 0

        CGCplus = 0
        CGCz = 0
        CGCminus = 0

        CGDplus = 0
        CGDz = 0
        CGDminus = 0

        !!!part A begin!!!
            do tau3=1 ,tau3_max
                do count=1, countsp_max(tau3)
                BSPVector(tau3,count)=SPVector(tau3,count,MFinNSCVector(2,bra))
                KSPVector(tau3,count)=SPVector(tau3,count,MIniNSCVector(2,ket))
                end do
            end do 

            do tau3=1, tau3_max

                AnnihMuAloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle AnnihMuAloop
                AnnihNuAloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle AnnihNuAloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence=sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence=sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                CoeffAnnih=0
                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    AJplus(tau3)=AJplus(tau3)+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    AJz(tau3)=AJz(tau3)+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    AJminus(tau3)=AJminus(tau3)+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,MFinNSCVector(2,bra))
                        KSPVector(orbit,count)=SPVector(orbit,count,MIniNSCVector(2,ket))
                    end do
                end do

                end do AnnihNuAloop
                end do AnnihMuAloop
            end do

            do tau3=1, tau3_max
                if (PorN(tau3)==1) g_factor = g_proton-g_rotor
                if (PorN(tau3)==0) g_factor = g_neutron-g_rotor
                
                AJplus(tau3)=AJplus(tau3)*g_factor
                AJz(tau3)=AJz(tau3)*g_factor
                AJminus(tau3)=AJminus(tau3)*g_factor                    
            end do
        !!!part A end!!!

        !!!part B begin!!!
            do tau3=1, tau3_max
                do count=1, countsp_max(tau3)
                    BSPVector(tau3,count)=ISPVector(tau3,count,MFinNSCVector(2,bra))
                    KSPVector(tau3,count)=SPVector(tau3,count,MIniNSCVector(2,ket))
                end do
            end do

            do tau3=1, tau3_max

                AnnihMuBloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle AnnihMuBloop
                AnnihNuBloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle AnnihNuBloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence=sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence=sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    BJplus(tau3)=BJplus(tau3)+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    BJz(tau3)=BJz(tau3)+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    BJminus(tau3)=BJminus(tau3)+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=ISPVector(orbit,count,MFinNSCVector(2,bra))
                        KSPVector(orbit,count)=SPVector(orbit,count,MIniNSCVector(2,ket))
                    end do
                end do

                end do AnnihNuBloop
                end do AnnihMuBloop
            end do

            do tau3=1, tau3_max
                if (PorN(tau3)==1) g_factor = g_proton-g_rotor
                if (PorN(tau3)==0) g_factor = g_neutron-g_rotor
                
                BJplus(tau3)=BJplus(tau3)*g_factor
                BJz(tau3)=BJz(tau3)*g_factor
                BJminus(tau3)=BJminus(tau3)*g_factor                    
            end do
        !!!part B end!!!

        !!!part C begin!!!
            do tau3=1, tau3_max
                do count=1, countsp_max(tau3)
                    BSPVector(tau3,count)=SPVector(tau3,count,MFinNSCVector(2,bra))
                    KSPVector(tau3,count)=ISPVector(tau3,count,MIniNSCVector(2,ket))
                end do
            end do

            do tau3=1, tau3_max

                AnnihMuCloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle AnnihMuCloop
                AnnihNuCloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle AnnihNuCloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence=sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence=sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    CJplus(tau3)=CJplus(tau3)+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    CJz(tau3)=CJz(tau3)+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    CJminus(tau3)=CJminus(tau3)+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,MFinNSCVector(2,bra))
                        KSPVector(orbit,count)=ISPVector(orbit,count,MIniNSCVector(2,ket))
                    end do
                end do

                end do AnnihNuCloop
                end do AnnihMuCloop
            end do
            
            do tau3=1, tau3_max
                if (CJz(tau3)/=0) write(*,*) CJz(tau3)
            end do

            do tau3=1, tau3_max
                if (PorN(tau3)==1) g_factor = g_proton-g_rotor
                if (PorN(tau3)==0) g_factor = g_neutron-g_rotor

                CJplus(tau3)=CJplus(tau3)*g_factor
                CJz(tau3)=CJz(tau3)*g_factor
                CJminus(tau3)=CJminus(tau3)*g_factor
            end do
        !!!part C end!!!

        !!!part D begin!!!
            do tau3=1, tau3_max
                do count=1, countsp_max(tau3)
                    BSPVector(tau3,count)=ISPVector(tau3,count,MFinNSCVector(2,bra))
                    KSPVector(tau3,count)=ISPVector(tau3,count,MIniNSCVector(2,ket))
                end do
            end do

            do tau3=1, tau3_max

                AnnihMuDloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle AnnihMuDloop
                AnnihNuDloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle AnnihNuDloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence=sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence=sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    DJplus(tau3)=DJplus(tau3)+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    DJz(tau3)=DJz(tau3)+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    DJminus(tau3)=DJminus(tau3)+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=ISPVector(orbit,count,MFinNSCVector(2,bra))
                        KSPVector(orbit,count)=ISPVector(orbit,count,MIniNSCVector(2,ket))
                    end do
                end do

                end do AnnihNuDloop
                end do AnnihMuDloop
            end do

            do tau3=1, tau3_max
                if (PorN(tau3)==1) g_factor = g_proton-g_rotor
                if (PorN(tau3)==0) g_factor = g_neutron-g_rotor

                DJplus(tau3)=DJplus(tau3)*g_factor
                DJz(tau3)=DJz(tau3)*g_factor
                DJminus(tau3)=DJminus(tau3)*g_factor
            end do
        !!!part D end!!!

        ! if (sum(BJplus)+sum(BJz)+sum(BJminus)/=0) then
        !     write(*,*) 'B', bra, ket
        !     CoeffNtoI=0
        !     do tau3=1, tau3_max
        !         CoeffNtoI=CoeffNtoI+PNumI(tau3,MFinNSCVector(2,bra)) &
        !         +PNumN(tau3,MFinNSCVector(2,bra))*PNumI(tau3,MFinNSCVector(2,bra))
        !     end do
        !     CoeffNtoI=CoeffNtoI+MFinState-FinK
        ! end if

        ! if (sum(CJplus)+sum(CJz)+sum(CJminus)/=0) then
            CoeffNtoI=0
            do tau3=1, tau3_max
                CoeffNtoI=CoeffNtoI+PNumI(tau3,MIniNSCVector(2,ket)) &
                +PNumN(tau3,MIniNSCVector(2,ket))*PNumI(tau3,MIniNSCVector(2,ket))
                ! write(*,*) tau3,PNumN(tau3,MIniNSCVector(2,ket)),PNumI(tau3,MIniNSCVector(2,ket))
            end do
            CoeffNtoI=CoeffNtoI+MIniState-MIniNSCVector(1,ket)
        ! end if

        ! if (sum(DJplus)+sum(DJz)+sum(DJminus)/=0) then
        !     write(*,*) 'D', bra, ket
        !     CoeffNtoI=0
        !     do tau3=1, tau3_max
        !         CoeffNtoI=CoeffNtoI+PNumI(tau3,MFinNSCVector(2,bra)) &
        !         +PNumN(tau3,MFinNSCVector(2,bra))*PNumI(tau3,MFinNSCVector(2,bra)) &
        !         +PNumI(tau3,MIniNSCVector(2,ket)) &
        !         +PNumN(tau3,MIniNSCVector(2,ket))*PNumI(tau3,MIniNSCVector(2,ket))
        !     end do
        !     CoeffNtoI=CoeffNtoI+MFinState-FinK+MIniState-IniK
        ! end if

        if (sum(AJplus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),IniK,dble(1),FinK, CGAplus)
        ! if (sum(BJplus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),IniK,dble(1),-FinK, CGBplus)
        if (sum(CJplus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),-IniK,dble(1),FinK, CGCplus)
        ! if (sum(DJplus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),-IniK,dble(1),-FinK, CGDplus)

        if (sum(AJz)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),IniK,dble(0),FinK, CGAz)
        ! if (sum(BJz)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),IniK,dble(0),-FinK, CGBz)
        if (sum(CJz)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),-IniK,dble(0),FinK, CGCz)
        ! if (sum(DJz)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),-IniK,dble(0),-FinK, CGDz)

        if (sum(AJminus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),IniK,dble(-1),FinK, CGAminus)
        ! if (sum(BJminus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),IniK,dble(-1),-FinK, CGBminus)
        if (sum(CJminus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),-IniK,dble(-1),FinK, CGCminus)
        ! if (sum(DJminus)/=0) call CGFactor(dble(MIniState),dble(1),dble(MFinState),-IniK,dble(-1),-FinK, CGDminus)

        BMResArray(bra,ket)= (CGAminus*sum(AJminus)/sqrt(2.)+CGAz*sum(AJz)-CGAplus*sum(AJplus)/sqrt(2.)) &
        ! + ((-1)**(MFinState-FinK)) &
        ! *(CGBminus*sum(BJminus)/sqrt(2.)+CGBz*sum(BJz)-CGBplus*sum(BJplus)/sqrt(2.)) &
        + ((-1)**CoeffNtoI) &
        *(CGCminus*sum(CJminus)/sqrt(2.)+CGCz*sum(CJz)-CGCplus*sum(CJplus)/sqrt(2.)) 
        ! + ((-1)**(MFinState-FinK+MIniState-IniK)) &
        ! *(CGDminus*sum(DJminus)/sqrt(2.)+CGDz*sum(DJz)-CGDplus*sum(DJplus)/sqrt(2.))

        BMResArray(bra,ket)=BMResArray(bra,ket)*MFinWF(bra,2)*MIniWF(ket,2)

    end do
    end do

    BMRes = 3*(sum(BMResArray)*sum(BMResArray))/(4*pi)

    write(*,*) BMRes

end subroutine MTransition