subroutine AMomentumJ

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use FunctionParameters

    implicit none

    PJ1Array=0
    PJ2Array=0
    PJ3Array=0

    NJ1Array=0
    NJ2Array=0
    NJ3Array=0

    do bra=1, DimensionArray(int(AMState))
    do ket=1, DimensionArray(int(AMState))
       
        KBra=AMNSCVector(1,bra)
        KKet=AMNSCVector(1,ket)

        PNJUDSMono=0
        PNJUDSPMono=0
        PNJUDSTB=0

        PNJxySMono=0
        PNJxySPMono=0
        PNJxySTB=0

        PNJ3SMono=0
        PNJ3SPMono=0
        PNJ3STB=0

        NNJUDSMono=0
        NNJUDSPMono=0
        NNJUDSTB=0

        NNJxySMono=0
        NNJxySPMono=0
        NNJxySTB=0

        NNJ3SMono=0
        NNJ3SPMono=0
        NNJ3STB=0

        if (KBra==KKet) then

            do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count)=SPVector(tau3,count,AMNSCVector(2,bra))
                KSPVector(tau3,count)=SPVector(tau3,count,AMNSCVector(2,ket))
            end do
            end do

            NJMloop: &
            do tau3=1, tau3_max

                NJMMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle NJMMuloop
                NJMNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle NJMNuloop

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

                        if (PorN(tau3)==1) then
                            PNJUDSMono=PNJUDSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            PNJxySMono=PNJxySMono+2*DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            PNJ3SMono=PNJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        end if

                        if (PorN(tau3)==0) then
                            NNJUDSMono=NNJUDSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            NNJxySMono=NNJxySMono+2*DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            NNJ3SMono=NNJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        end if
                    end if
                    
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do NJMNuloop
                end do NJMMuloop
            end do NJMloop

            ! NJPloop: &
            ! do tau3=1, tau3_max
            ! do tau3Apo=1, tau3_max
            !     if (tau3==tau3Apo) cycle

            !     NJPMuloop: &
            !     do AnnihMu=1, countsp_max(tau3)
            !         if (BSPVector(tau3,AnnihMu)==0) cycle NJPMuloop

            !     NJPNuloop: &
            !     do AnnihNu=1, countsp_max(tau3)
            !         if (KSPVector(tau3,AnnihNu)==0) cycle NJPNuloop

            !     NJPMuApoloop: &
            !     do AnnihMuApo=1, countsp_max(tau3)
            !         if (BSPVector(tau3Apo,AnnihMuApo)==0) cycle NJPMuApoloop

            !     NJPNuApoloop: &
            !     do AnnihNuApo=1, countsp_max(tau3)
            !         if (KSPVector(tau3Apo,AnnihNuApo)==0) cycle NJPNuApoloop

            !         BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
            !         MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))
            !         KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
            !         NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

            !         BSPVector(tau3Apo,AnnihMuApo)=BSPVector(tau3Apo,AnnihMuApo)-1
            !         MuSeApo=sum(BSPVector(tau3Apo,1:AnnihMuApo-1))
            !         KSPVector(tau3Apo,AnnihNuApo)=KSPVector(tau3Apo,AnnihNuApo)-1
            !         NuSeApo=sum(KSPVector(tau3Apo,1:AnnihNuApo-1))

            !         SignPass=1
            !         do orbit=1, tau3_max
            !             do count=1, countsp_max(orbit)
            !                 if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
            !             end do
            !         end do

            !         if (SignPass==1) then
            !             CoeffAnnih=(-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

            !             if (PorN(tau3)==1) then
            !                 PNJUDSPMono=PNJUDSPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 PNJxySPMono=PNJxySPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 PNJ3SPMono=PNJ3SPMono+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
            !                 +DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
            !             end if

            !             if (PorN(tau3)==0) then
            !                 NNJUDSPMono=NNJUDSPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 NNJxySPMono=NNJxySPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 NNJ3SPMono=NNJ3SPMono+DecoupleESum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleESum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
            !             end if
            !         end if
                
            !         do orbit=1, tau3_max
            !             do count=1, countsp_max(orbit)
            !                 BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
            !                 KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
            !             end do
            !         end do

            !     end do NJPNuApoloop
            !     end do NJPMuApoloop
            !     end do NJPNuloop
            !     end do NJPMuloop
            ! end do
            ! end do NJPloop

            NJTBloop: &
            do tau3=1, tau3_max
                NJTBMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle NJTBMuloop
                NJTBNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle NJTBNuloop

                    BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                    MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))

                    KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                    NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

                NJTBMuApoloop: &
                do AnnihMuApo=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMuApo)==0) cycle NJTBMuApoloop
                NJTBNuApoloop: &
                do AnnihNuApo=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNuApo)==0) cycle NJTBNuApoloop

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

                        if (PorN(tau3)==1) then
                            PNJUDSTB=PNJUDSTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            PNJxySTB=PNJxySTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            PNJ3STB=PNJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                        end if

                        if (PorN(tau3)==0) then
                            NNJUDSTB=NNJUDSTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            NNJxySTB=NNJxySTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            NNJ3STB=NNJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                        end if
                    end if

                end do NJTBNuApoloop
                end do NJTBMuApoloop

                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do

                end do NJTBNuloop
                end do NJTBMuloop
            end do NJTBloop
        end if

        PIJUDSMono=0
        PIJUDSPMono=0
        PIJUDSTB=0

        PIJxySMono=0
        PIJxySPMono=0
        PIJxySTB=0

        PIJ3SMono=0
        PIJ3SPMono=0
        PIJ3STB=0

        NIJUDSMono=0
        NIJUDSPMono=0
        NIJUDSTB=0

        NIJxySMono=0
        NIJxySPMono=0
        NIJxySTB=0

        NIJ3SMono=0
        NIJ3SPMono=0
        NIJ3STB=0

        if (KBra==-KKet) then

            do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count)=SPVector(tau3,count,AMNSCVector(2,bra))
                KSPVector(tau3,count)=ISPVector(tau3,count,AMNSCVector(2,ket))
            end do
            end do

            IJMloop: &
            do tau3=1, tau3_max

                IJMMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle IJMMuloop
                IJMNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle IJMNuloop

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

                        if (PorN(tau3)==1) then
                            PIJUDSMono=PIJUDSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            PIJxySMono=PIJxySMono+2*DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            PIJ3SMono=PIJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        end if

                        if (PorN(tau3)==0) then
                            NIJUDSMono=NIJUDSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            NIJxySMono=NIJxySMono+2*DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                            NIJ3SMono=NIJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                        end if
                    end if
                    
                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do
                
                end do IJMNuloop
                end do IJMMuloop
            end do IJMloop

            ! IJPloop: &
            ! do tau3=1, tau3_max
            ! do tau3Apo=1, tau3_max
            !     if (tau3==tau3Apo) cycle

            !     IJPMuloop: &
            !     do AnnihMu=1, countsp_max(tau3)
            !         if (BSPVector(tau3,AnnihMu)==0) cycle IJPMuloop

            !     IJPNuloop: &
            !     do AnnihNu=1, countsp_max(tau3)
            !         if (KSPVector(tau3,AnnihNu)==0) cycle IJPNuloop

            !     IJPMuApoloop: &
            !     do AnnihMuApo=1, countsp_max(tau3)
            !         if (BSPVector(tau3Apo,AnnihMuApo)==0) cycle IJPMuApoloop

            !     IJPNuApoloop: &
            !     do AnnihNuApo=1, countsp_max(tau3)
            !         if (KSPVector(tau3Apo,AnnihNuApo)==0) cycle IJPNuApoloop

            !         BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
            !         MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))
            !         KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
            !         NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

            !         BSPVector(tau3Apo,AnnihMuApo)=BSPVector(tau3Apo,AnnihMuApo)-1
            !         MuSeApo=sum(BSPVector(tau3Apo,1:AnnihMuApo-1))
            !         KSPVector(tau3Apo,AnnihNuApo)=KSPVector(tau3Apo,AnnihNuApo)-1
            !         NuSeApo=sum(KSPVector(tau3Apo,1:AnnihNuApo-1))

            !         SignPass=1
            !         do orbit=1, tau3_max
            !             do count=1, countsp_max(orbit)
            !                 if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
            !             end do
            !         end do

            !         if (SignPass==1) then
            !             CoeffAnnih=(-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

            !             if (PorN(tau3)==1) then
            !                 PIJUDSPMono=PIJUDSPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 PIJxySPMono=PIJxySPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 PIJ3SPMono=PIJ3SPMono+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
            !             end if

            !             if (PorN(tau3)==0) then
            !                 NIJUDSPMono=NIJUDSPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 NIJxySPMono=NIJxySPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
            !                 + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

            !                 NIJ3SPMono=NIJ3SPMono+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
            !                 *DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
            !             end if
            !         end if
                
            !         do orbit=1, tau3_max
            !             do count=1, countsp_max(orbit)
            !                 BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
            !                 KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
            !             end do
            !         end do

            !     end do IJPNuApoloop
            !     end do IJPMuApoloop
            !     end do IJPNuloop
            !     end do IJPMuloop
            ! end do
            ! end do IJPloop

            IJTBloop: &
            do tau3=1, tau3_max
                IJTBMuloop: &
                do AnnihMu=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMu)==0) cycle IJTBMuloop
                IJTBNuloop: &
                do AnnihNu=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNu)==0) cycle IJTBNuloop

                    BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                    MuSerial=sum(BSPVector(tau3,1:AnnihMu-1))

                    KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                    NuSerial=sum(KSPVector(tau3,1:AnnihNu-1))

                IJTBMuApoloop: &
                do AnnihMuApo=1, countsp_max(tau3)
                    if (BSPVector(tau3,AnnihMuApo)==0) cycle IJTBMuApoloop
                IJTBNuApoloop: &
                do AnnihNuApo=1, countsp_max(tau3)
                    if (KSPVector(tau3,AnnihNuApo)==0) cycle IJTBNuApoloop

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

                        if (PorN(tau3)==1) then
                            PIJUDSTB=PIJUDSTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            PIJxySTB=PIJxySTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            PIJ3STB=PIJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                        end if

                        if (PorN(tau3)==0) then
                            NIJUDSTB=NIJUDSTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            NIJxySTB=NIJxySTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                            + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                            NIJ3STB=NIJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                            *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                        end if
                    end if

                end do IJTBNuApoloop
                end do IJTBMuApoloop

                    do orbit=1, tau3_max
                        do count=1, countsp_max(orbit)
                            BSPVector(orbit,count)=SPVector(orbit,count,AMNSCVector(2,bra))
                            KSPVector(orbit,count)=ISPVector(orbit,count,AMNSCVector(2,ket))
                        end do
                    end do

                end do IJTBNuloop
                end do IJTBMuloop
            end do IJTBloop
        end if

        PJ1Array(bra,ket)=0.25*(PNJUDSMono+PNJUDSPMono+PNJUDSTB) &
        + 0.25*(PNJxySMono+PNJxySPMono+PNJxySTB) &
        + ((-1)**(AMState-AMNSCVector(1,ket))) &
        *(0.25*(PIJUDSMono+PIJUDSPMono+PIJUDSTB) &
        + 0.25*(PIJxySMono+PIJxySPMono+PIJxySTB))

        PJ2Array(bra,ket)=-0.25*(PNJUDSMono+PNJUDSPMono+PNJUDSTB) &
        + 0.25*(PNJxySMono+PNJxySPMono+PNJxySTB) &
        + ((-1)**(AMState-AMNSCVector(1,ket))) &
        *(-0.25*(PIJUDSMono+PIJUDSPMono+PIJUDSTB) &
        + 0.25*(PIJxySMono+PIJxySPMono+PIJxySTB))

        PJ3Array(bra,ket)=(PNJ3SMono+PNJ3SPMono+PNJ3STB) &
        + (-1)**(AMState-KKet) &
        *(PIJ3SMono+PIJ3SPMono+PIJ3STB)

        PJ1Array(bra,ket)=PJ1Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)
        PJ2Array(bra,ket)=PJ2Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)
        PJ3Array(bra,ket)=PJ3Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)

        NJ1Array(bra,ket)=0.25*(NNJUDSMono+NNJUDSPMono+NNJUDSTB) &
        + 0.25*(NNJxySMono+NNJxySPMono+NNJxySTB) &
        + ((-1)**(AMState-KKet)) &
        *(0.25*(NIJUDSMono+NIJUDSPMono+NIJUDSTB) &
        + 0.25*(NIJxySMono+NIJxySPMono+NIJxySTB))

        NJ2Array(bra,ket)=-0.25*(NNJUDSMono+NNJUDSPMono+NNJUDSTB) &
        + 0.25*(NNJxySMono+NNJxySPMono+NNJxySTB) &
        + ((-1)**(AMState-KKet)) &
        *(-0.25*(NIJUDSMono+NIJUDSPMono+NIJUDSTB) &
        + 0.25*(NIJxySMono+NIJxySPMono+NIJ3STB))

        NJ3Array(bra,ket)=(NNJ3SMono+NNJ3SPMono+NNJ3STB) &
        + ((-1)**(AMState-KKet)) &
        *(NIJ3SMono+NIJ3SPMono+NIJ3STB)

        NJ1Array(bra,ket)=NJ1Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)
        NJ2Array(bra,ket)=NJ2Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)
        NJ3Array(bra,ket)=NJ3Array(bra,ket)*AMWF(bra,2)*AMWF(ket,2)

    end do  
    end do

    PJ1Res=sqrt(sum(PJ1Array))
    PJ2Res=sqrt(sum(PJ2Array))
    PJ3Res=sqrt(sum(PJ3Array))

    NJ1Res=sqrt(sum(NJ1Array))
    NJ2Res=sqrt(sum(NJ2Array))
    NJ3Res=sqrt(sum(NJ3Array))

    write(*,*) NJ1Res, NJ2Res, NJ3Res

end subroutine AMomentumJ