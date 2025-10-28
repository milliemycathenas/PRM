subroutine KEquiv

    use ProcessVariables
    use Constants
    use PhysicalQuantities
    use DiagonalVariables
    use FunctionParameters

    implicit none

    if (NSCVector(1,bra)==NSCVector(1,ket)) then

        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count)=SPVector(tau3,count,NSCVector(2,bra))
                KSPVector(tau3,count)=SPVector(tau3,count,NSCVector(2,ket))
            end do
        end do

        NEqMloop: &
        do tau3=1, tau3_max
            NMMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle NMMuloop
            NMNuloop: &
            do AnnihNu=1, countsp_max(tau3) 
                if (KSPVector(tau3,AnnihNu)==0) cycle NMNuloop
                
                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence = sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence = sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    NJ3Mono = NJ3Mono+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    NJ3SMono = NJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    NJSMono = NJSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    NJdJuMono = NJdJuMono+DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    NJuJdMono = NJuJdMono+DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do NMNuloop
            end do NMMuloop
        end do NEqMloop

        NEqPloop: &
        do tau3=1, tau3_max
        do tau3Apo=1, tau3_max
            if (tau3==tau3Apo) cycle

            NPMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle NPMuloop
                
            NPNuloop: &
            do AnnihNu=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNu)==0) cycle NPNuloop    
                
            NPMuApoloop: &
            do AnnihMuApo=1, countsp_max(tau3Apo)
                if (BSPVector(tau3Apo,AnnihMuApo)==0) cycle NPMuApoloop
                
            NPNuApoloop: &
            do AnnihNuApo=1, countsp_max(tau3Apo)
                if (KSPVector(tau3Apo,AnnihNuApo)==0) cycle NPNuApoloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSerial = sum(BSPVector(tau3,1:AnnihMu-1))
                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSerial = sum(KSPVector(tau3,1:AnnihNu-1))
                BSPVector(tau3Apo,AnnihMuApo)=BSPVector(tau3Apo,AnnihMuApo)-1
                MuSeApo = sum(BSPVector(tau3Apo,1:AnnihMuApo-1))
                KSPVector(tau3Apo,AnnihNuApo)=KSPVector(tau3Apo,AnnihNuApo)-1
                NuSeApo = sum(KSPVector(tau3Apo,1:AnnihNuApo-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                    NJ3SPMono = NJ3SPMono+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    NJSPMono = NJSPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                    + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    NJuJdPMono = NJuJdPMono+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    NJdJuPMono = NJdJuPMono+DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do NPNuApoloop
            end do NPMuApoloop
            end do NPNuloop
            end do NPMuloop
        end do
        end do NEqPloop

        NEqTBloop: &
        do tau3=1, tau3_max
            NTBMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle NTBMuloop

            NTBNuloop: &
            do AnnihNu=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNu)==0) cycle NTBNuloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSerial = sum(BSPVector(tau3,1:AnnihMu-1))
                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSerial = sum(KSPVector(tau3,1:AnnihNu-1))

            NTBMuApoloop: &
            do AnnihMuApo=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMuApo)==0) cycle NTBMuApoloop    
                
            NTBNuApoloop: &
            do AnnihNuApo=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNuApo)==0) cycle NTBNuApoloop 

                BSPVector(tau3,AnnihMuApo)=BSPVector(tau3,AnnihMuApo)-1
                MuSeApo = sum(BSPVector(tau3,1:AnnihMuApo-1))
                KSPVector(tau3,AnnihNuApo)=KSPVector(tau3,AnnihNuApo)-1
                NuSeApo = sum(KSPVector(tau3,1:AnnihNuApo-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do
                
                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                    NJ3STB = NJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    NJSTB = NJSTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                    + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    NJuJdTB = NJuJdTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    NJdJuTB = NJdJuTB+DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                end if

            end do NTBNuApoloop
            end do NTBMuApoloop

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do NTBNuloop
            end do NTBMuloop
        end do NEqTBloop

    end if

    if (NSCVector(1,bra)==ISCVector(1,ket)) then

        do tau3=1, tau3_max
            do count=1, countsp_max(tau3)
                BSPVector(tau3,count) = SPVector(tau3,count,NSCVector(2,bra))
                KSPVector(tau3,count) = ISPVector(tau3,count,NSCVector(2,ket))
            end do
        end do

        IEqMloop: &
        do tau3=1, tau3_max
            IMMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle IMMuloop
            IMNuloop: &
            do AnnihNu=1, countsp_max(tau3) 
                if (KSPVector(tau3,AnnihNu)==0) cycle IMNuloop
                
                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSequence = sum(BSPVector(tau3,1:AnnihMu-1))

                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSequence = sum(KSPVector(tau3,1:AnnihNu-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSequence+NuSequence)

                    IJ3Mono = IJ3Mono+DecoupleDSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih

                    IJ3SMono = IJ3SMono+DecoupleESum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    IJSMono = IJSMono+DecoupleGSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    IJdJuMono = IJdJuMono+DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                    IJuJdMono = IJuJdMono+DecoupleFSum(tau3,AnnihMu,AnnihNu)*CoeffAnnih
                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=ISPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do IMNuloop
            end do IMMuloop

        end do IEqMloop

        IEqPloop: &
        do tau3=1, tau3_max
        do tau3Apo=1, tau3_max
            if (tau3==tau3Apo) cycle 

            IPMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle IPMuloop
                
            IPNuloop: &
            do AnnihNu=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNu)==0) cycle IPNuloop
                
            IPMuApoloop: &
            do AnnihMuApo=1, countsp_max(tau3Apo)
                if (BSPVector(tau3Apo,AnnihMuApo)==0) cycle IPMuApoloop
                
            IPNuApoloop: &
            do AnnihNuApo=1, countsp_max(tau3Apo)
                if (KSPVector(tau3Apo,AnnihNuApo)==0) cycle IPNuApoloop
                
                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSerial = sum(BSPVector(tau3,1:AnnihMu-1))
                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSerial = sum(KSPVector(tau3,1:AnnihNu-1))
                BSPVector(tau3Apo,AnnihMuApo)=BSPVector(tau3Apo,AnnihMuApo)-1
                MuSeApo = sum(BSPVector(tau3Apo,1:AnnihMuApo-1))
                KSPVector(tau3Apo,AnnihNuApo)=KSPVector(tau3Apo,AnnihNuApo)-1
                NuSeApo = sum(KSPVector(tau3Apo,1:AnnihNuApo-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then

                    CoeffAnnih = (-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                    IJ3SPMono = IJ3SPMono + DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleDSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    IJSPMono = IJSPMono + DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                    + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                    
                    IJuJdPMono = IJuJdPMono + DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    IJdJuPMono = IJdJuPMono + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3Apo,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                end if

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=ISPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do IPNuApoloop
            end do IPMuApoloop            
            end do IPNuloop
            end do IPMuloop
        end do
        end do IEqPloop

        IEqTBloop: &
        do tau3=1, tau3_max
            ITBMuloop: &
            do AnnihMu=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMu)==0) cycle ITBMuloop
                
            ITBNuloop: &
            do AnnihNu=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNu)==0) cycle ITBNuloop

                BSPVector(tau3,AnnihMu)=BSPVector(tau3,AnnihMu)-1
                MuSerial = sum(BSPVector(tau3,1:AnnihMu-1))
                KSPVector(tau3,AnnihNu)=KSPVector(tau3,AnnihNu)-1
                NuSerial = sum(KSPVector(tau3,1:AnnihNu-1))
                
            ITBMuApoloop: &
            do AnnihMuApo=1, countsp_max(tau3)
                if (BSPVector(tau3,AnnihMuApo)==0) cycle ITBMuApoloop
                
            ITBNuApoloop: &
            do AnnihNuApo=1, countsp_max(tau3)
                if (KSPVector(tau3,AnnihNuApo)==0) cycle ITBNuApoloop

                BSPVector(tau3,AnnihMuApo)=BSPVector(tau3,AnnihMuApo)-1
                MuSeApo = sum(BSPVector(tau3,1:AnnihMuApo-1))
                KSPVector(tau3,AnnihNuApo)=KSPVector(tau3,AnnihNuApo)-1
                NuSeApo = sum(KSPVector(tau3,1:AnnihNuApo-1))

                SignPass=1
                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        if (BSPVector(orbit,count)/=KSPVector(orbit,count)) SignPass=0
                    end do
                end do

                if (SignPass==1) then
                    CoeffAnnih = (-1)**(MuSerial+NuSerial+MuSeApo+NuSeApo)

                    IJ3STB = IJ3STB+DecoupleDSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleDSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    IJSTB = IJSTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih &
                    + DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    IJuJdTB = IJuJdTB+DecoupleABSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleCSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih

                    IJdJuTB = IJdJuTB+DecoupleCSum(tau3,AnnihMu,AnnihNu) &
                    *DecoupleABSum(tau3,AnnihMuApo,AnnihNuApo)*CoeffAnnih
                end if

            end do ITBNuApoloop
            end do ITBMuApoloop

                do orbit=1, tau3_max
                    do count=1, countsp_max(orbit)
                        BSPVector(orbit,count)=SPVector(orbit,count,NSCVector(2,bra))
                        KSPVector(orbit,count)=ISPVector(orbit,count,NSCVector(2,ket))
                    end do
                end do

            end do ITBNuloop
            end do ITBMuloop
        end do IEqTBloop

    end if

end subroutine KEquiv