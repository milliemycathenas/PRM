subroutine AMomentumR

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters

    implicit none

    real :: NRotorK0, NRotorKp1, NRotorKn1, NRotorKp2, NRotorKn2, NRotorKz
    real :: IRotorK0, IRotorKp1, IRotorKn1, IRotorKp2, IRotorKn2, IRotorKz

    R1Array = 0
    R2Array = 0
    R3Array = 0

    Rloop: &
    do bra=1, DimensionArray(int(AMState))
    do ket=1, DimensionArray(int(AMState))

        call RInitialise

        NRotorK0=0
        NRotorKp1=0
        NRotorKn1=0
        NRotorKp2=0
        NRotorKn2=0
        NRotorKz=0

        IRotorK0=0
        IRotorKp1=0
        IRotorKn1=0
        IRotorKp2=0
        IRotorKn2=0
        IRotorKz=0

        KBra=AMNSCVector(1,bra)
        KKet=AMNSCVector(1,ket)

        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count)=SPVector(tau3,count,AMNSCVector(2,bra))
                KSPVector(tau3,count)=SPVector(tau3,count,AMNSCVector(2,ket))
            end do
        end do

        if (KBra==KKet+1 .or. KBra==KKet-1) then
            NRUDloop: &
            do tau3=1, tau3_max

                NRUDMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle NRUDMuloop
                NRUDNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle NRUDNuloop

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
                        CoeffAnnih=(-1)**(MuSequence+NuSequence)

                        RNJUp=RNJUp+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        RNJDown=RNJDown+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    end if
                    
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do NRUDNuloop
                end do NRUDMuloop
            end do NRUDloop
        end if
 
        if (KBra==KKet) then
            NRMloop: &
            do tau3=1, tau3_max

                NRMMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle NRMMuloop
                NRMNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle NRMNuloop

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
                        CoeffAnnih=(-1)**(MuSequence+NuSequence)

                        RNJ3=RNJ3+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih

                        RNJ3SMono=RNJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        RNJUDSMono=RNJUDSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        RNJxySMono=RNJxySMono+DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        end if
                    
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do NRMNuloop
                end do NRMMuloop
            end do NRMloop
            

            NRPloop: &
            do tau3=1, tau3_max
            do tau3Apo=1, tau3_max
                if (tau3==tau3Apo) cycle

                NRPMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle NRPMuloop

                NRPNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle NRPNuloop

                NRPMuApoloop: &
                do AnnihMuApo=1, countsp_max(tau3)
                    if (BSPVector(tau3Apo,AnnihMuApo)==0) cycle NRPMuApoloop

                NRPNuApoloop: &
                do AnnihNuApo=1, countsp_max(tau3)
                    if (KSPVector(tau3Apo,AnnihNuApo)==0) cycle NRPNuApoloop

                    BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                    MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))
                    KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                    NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

                    BSPVector(tau3Apo,AnnihMuApo)=BSPVector(tau3Apo,AnnihMuApo)-1
                    MuSeApo=sum(BSPVector(tau3Apo,1:AnnihMuApo-1))
                    KSPVector(tau3Apo,AnnihNuApo)=KSPVector(tau3Apo,AnnihNuApo)-1
                    NuSeApo=sum(KSPVector(tau3Apo,1:AnnihNuApo-1))

                    SignPass=1
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                        end do
                    end do

                    if (SignPass==1) then
                        CoeffAnnih=(-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                        RNJ3SPMono=RNJ3SPMono+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RNJUDSPMono=RNJUDSPMono &
                        + DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RNJxySPMono=RNJxySPMono &
                        + 0.5*DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + 0.5*DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                    end if
                
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do

                end do NRPNuApoloop
                end do NRPMuApoloop
                end do NRPNuloop
                end do NRPMuloop
            end do
            end do NRPloop

            NRTBloop: &
            do tau3=1, tau3_max
                NRTBMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle NRTBMuloop
                NRTBNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle NRTBNuloop

                    BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                    MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))

                    KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                    NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

                NRTBMuApoloop: &
                do AnnihMuApo=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMuApo)==0) cycle NRTBMuApoloop
                NRTBNuApoloop: &
                do AnnihNuApo=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNuApo)==0) cycle NRTBNuApoloop

                    BSPVector(tau3,AnnihMuApo)=BSPVector(tau3,AnnihMuApo)-1
                    MuSeApo=sum(BSPVector(tau3,1:AnnihMuApo-1))

                    KSPVector(tau3,AnnihNuApo)=KSPVector(tau3,AnnihNuApo)-1
                    NuSeApo=sum(KSPVector(tau3,1:AnnihNuApo-1))

                    SignPass=1
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                        end do
                    end do

                    if (SignPass==1) then
                        CoeffAnnih=(-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                        RNJ3STB=RNJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RNJUDSTB=RNJUDSTB &
                        + DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RNJxySTB=RNJxySTB &
                        + 0.5*DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + 0.5*DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                    end if

                end do NRTBNuApoloop
                end do NRTBMuApoloop

                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do

                end do NRTBNuloop
                end do NRTBMuloop
            end do NRTBloop
        end if

        !!!normal rotor coefficient begin!!!
            if (KBra==KKet .and. AMNSCVector(2,bra)==AMNSCVector(2,ket)) &
            NRotorK0=AMState*(AMState+1)-KKet**2

            if (KBra==KKet+1) &
            NRotorKp1=sqrt(AMState*(AMState+1)-KKet*(KKet+1))

            if (KBra==KKet-1) &
            NRotorKn1=sqrt(AMState*(AMState+1)-KKet*(KKet-1))

            if (KBra==KKet+2 .and. AMNSCVector(2,bra)==AMNSCVector(2,ket)) &
            NRotorKp2=sqrt(AMState*(AMState+1)-KKet*(KKet+1)) &
            *sqrt(AMState*(AMState+1)-(KKet+1)*(KKet+2))

            if (KBra==KKet-2 .and. AMNSCVector(2,bra)==AMNSCVector(2,ket)) &
            NRotorKn2=sqrt(AMState*(AMState+1)-KKet*(KKet-1)) &
            *sqrt(AMState*(AMState+1)-(KKet-1)*(KKet-2))

            if (KBra==KKet .and. AMNSCVector(2,bra)==AMNSCVector(2,ket)) &
            NRotorKz=KKet**2
        !!!normal rotor coefficient end!!!

        KKet=-KKet

        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count)=SPVector(tau3,count,AMNSCVector(2,bra))
                KSPVector(tau3,count)=ISPVector(tau3,count,AMNSCVector(2,ket))
            end do
        end do

        if (KBra==KKet+1 .or. KBra==KKet-1) then
            IRUDloop: &
            do tau3=1, tau3_max

                IRUDMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle IRUDMuloop
                IRUDNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle IRUDNuloop

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
                        CoeffAnnih=(-1)**(MuSequence+NuSequence)

                        RIJUp=RIJUp+DecoupleABSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        RIJDown=RIJDown+DecoupleCSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    end if
                    
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do IRUDNuloop
                end do IRUDMuloop
            end do IRUDloop
        end if

        if (KBra==KKet) then
            IRMloop: &
            do tau3=1, tau3_max

                IRMMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle IRMMuloop
                IRMNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle IRMNuloop

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
                        CoeffAnnih=(-1)**(MuSequence+NuSequence)

                        RIJ3=RIJ3+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih

                        RIJ3SMono=RIJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        RIJUDSMono=RIJUDSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        RIJxySMono=RIJxySMono+DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    end if
                    
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do IRMNuloop
                end do IRMMuloop
            end do IRMloop

            IRPloop: &
            do tau3=1, tau3_max
            do tau3Apo=1, tau3_max
                if (tau3==tau3Apo) cycle

                IRPMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle IRPMuloop
                IRPNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle IRPNuloop
                IRPMuApoloop: &
                do AnnihMuApo=1, countsp_max(tau3)
                    if (BSPVector(tau3Apo,AnnihMuApo)==0) cycle IRPMuApoloop
                IRPNuApoloop: &
                do AnnihNuApo=1, countsp_max(tau3)
                    if (KSPVector(tau3Apo,AnnihNuApo)==0) cycle IRPNuApoloop

                    BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                    MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))
                    KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                    NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

                    BSPVector(tau3Apo,AnnihMuApo)=BSPVector(tau3Apo,AnnihMuApo)-1
                    MuSeApo=sum(BSPVector(tau3Apo,1:AnnihMuApo-1))
                    KSPVector(tau3Apo,AnnihNuApo)=KSPVector(tau3Apo,AnnihNuApo)-1
                    NuSeApo=sum(KSPVector(tau3Apo,1:AnnihNuApo-1))

                    SignPass=1
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                        end do
                    end do

                    if (SignPass==1) then
                        CoeffAnnih=(-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                        RIJ3SPMono=RIJ3SPMono + DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RIJUDSPMono=RIJUDSPMono &
                        + DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RIJxySPMono=RIJxySPMono &
                        + 0.5*DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + 0.5*DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                    end if
                
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do

                end do IRPNuApoloop
                end do IRPMuApoloop
                end do IRPNuloop
                end do IRPMuloop
            end do
            end do IRPloop

            IRTBloop: &
            do tau3=1, tau3_max
                IRTBMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle IRTBMuloop
                IRTBNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle IRTBNuloop

                    BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                    MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))
                    KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                    NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

                IRTBMuApoloop: &
                do AnnihMuApo=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMuApo)==0) cycle IRTBMuApoloop
                IRTBNuApoloop: &
                do AnnihNuApo=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNuApo)==0) cycle IRTBNuApoloop

                    BSPVector(tau3,AnnihMuApo)=BSPVector(tau3,AnnihMuApo)-1
                    MuSeApo=sum(BSPVector(tau3,1:AnnihMuApo-1))
                    KSPVector(tau3,AnnihNuApo)=KSPVector(tau3,AnnihNuApo)-1
                    NuSeApo=sum(KSPVector(tau3,1:AnnihNuApo-1))

                    SignPass=1
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                        end do
                    end do

                    if (SignPass==1) then
                        CoeffAnnih=(-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                        RIJ3STB=RIJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RIJUDSTB=RIJUDSTB &
                        + DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                        RIJxySTB=RIJxySTB &
                        + 0.5*DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                        + 0.5*DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                        *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                    end if

                end do IRTBNuApoloop
                end do IRTBMuApoloop

                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do IRTBNuloop
                end do IRTBMuloop
            end do IRTBloop
        end if

        !!!inverse rotor coefficient begin!!!
            ! if (KBra==KKet .and. AMNSCVector(2,bra)==AMISCVector(2,ket)) &
            ! IRotorK0=AMState*(AMState+1)-KKet**2

            if (KBra==KKet+1) &
            IRotorKp1=sqrt(AMState*(AMState+1)-KKet*(KKet+1))

            if (KBra==KKet-1) &
            IRotorKn1=sqrt(AMState*(AMState+1)-KKet*(KKet-1))

            ! if (KBra==KKet+2 .and. AMNSCVector(2,bra)==AMISCVector(2,ket)) &
            ! IRotorKp2=sqrt(AMState*(AMState+1)-KKet*(KKet+1)) &
            ! *sqrt(AMState*(AMState+1)-(KKet+1)*(KKet+2))

            ! if (KBra==KKet-2 .and. AMNSCVector(2,bra)==AMISCVector(2,ket)) &
            ! IRotorKn2=sqrt(AMState*(AMState+1)-KKet*(KKet-1)) &
            ! *sqrt(AMState*(AMState+1)-(KKet-1)*(KKet-2))

            ! if (KBra==KKet .and. AMNSCVector(2,bra)==AMISCVector(2,ket)) IRotorKz=KKet**2
        !!!inverse rotor coefficient end!!!

        R1Array(bra,ket) = 0.25*(NRotorKp2+NRotorKn2) + 0.5*NRotorK0 &
        + 0.25*(RNJUDSMono+RNJUDSPMono+RNJUDSTB) &
        + 0.5*(RNJxySMono+RNJxySPMono+RNJxySTB) &
        - 0.5*NRotorKp1*(RNJUp+RNJDown) - 0.5*NRotorKn1*(RNJUp+RNJDown) &
        + ((-1)**(AMState-AMNSCVector(1,ket))) &
        *(0.25*(IRotorKp2+IRotorKn2) + 0.5*IRotorK0 &
        + 0.25*(RIJUDSMono+RIJUDSPMono+RIJUDSTB) &
        + 0.5*(RIJxySMono+RIJxySPMono+RIJxySTB) &
        - 0.5*IRotorKp1*(RIJUp+RIJDown) - 0.5*IRotorKn1*(RIJUp+RIJDown))

        R2Array(bra,ket) = 0.5*NRotorK0 - 0.25*(NRotorKp2+NRotorKn2) &
        + 0.5*(RNJxySMono+RNJxySPMono+RNJxySTB) &
        - 0.25*(RNJUDSMono+RNJUDSPMono+RNJUDSTB) &
        - 0.5*(NRotorKp1*RNJUp+NRotorKn1*RNJDown) &
        + 0.5*(NRotorKp1*RNJDown+NRotorKn1*RNJUp) &
        + ((-1)**(AMState-AMNSCVector(1,ket))) &
        *(0.5*IRotorK0 - 0.25*(IRotorKp2+IRotorKn2) &
        + 0.5*(RIJxySMono+RIJxySPMono+RIJxySTB) &
        - 0.25*(RIJUDSMono+RIJUDSPMono+RIJUDSTB) &
        - 0.5*(IRotorKp1*RIJUp+IRotorKn1*RIJDown) &
        + 0.5*(IRotorKp1*RIJDown+IRotorKn1*RIJUp))

        R3Array(bra,ket) = NRotorKz - 2*AMNSCVector(1,ket)*RNJ3 + (RNJ3SMono+RNJ3SPMono+RNJ3STB) &
        + ((-1)**(AMState-AMNSCVector(1,ket))) &
        *(IRotorKz + 2*AMNSCVector(1,ket)*RIJ3 + (RIJ3SMono+RIJ3SPMono+RIJ3STB))

        R1Array(bra,ket) = R1Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)
        R2Array(bra,ket) = R2Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)
        R3Array(bra,ket) = R3Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)

    end do
    end do Rloop

    R1Res=sqrt(sum(R1Array))
    R2Res=sqrt(sum(R2Array))
    R3Res=sqrt(sum(R3Array))

    write(*,*) R1Res, R2Res, R3Res

end subroutine AMomentumR