%SOURCE=[1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;1,2,3,4,5,6,7,8,9,10,11,12;]
SOURCE=meshgrid(1:10,1:10)
%SOURCE=[1,1,1,1,1,1,1,1,1,1;2,2,2,2,2,2,2,2,2,2;3,3,3,3,3,3,3,3,3,3;4,4,4,4,4,4,4,4,4,4;5,5,5,5,5,5,5,5,5,5;6,6,6,6,6,6,6,6,6,6;7,7,7,7,7,7,7,7,7,7;8,8,8,8,8,8,8,8,8,8;9,9,9,9,9,9,9,9,9,9;10,10,10,10,10,10,10,10,10,10]
mysize=size(SOURCE)
ly=mysize(1)
lx=mysize(2)
maximum=sqrt((lx/2).^2+(ly/2).^2);
        %vi4islyaem ro-theta 4etverti
        for i=1:(lx/2)
        for j=1:(ly/2)
            x=i
            y=j
            Iro(j,i)=sqrt(x.^2+y.^2)/maximum
            Itheta(j,i)=atan2(y,x)
        end
        end
        for i=1:(lx/2)
        for j=1:(ly/2)
            x=i
            y=-j
            IIro(j,i)=sqrt(x.^2+y.^2)/maximum
            IItheta(j,i)=atan2(y,x)
        end
        end
        for i=1:(lx/2)
        for j=1:(ly/2)
            x=-i
            y=j
            IIIro(j,i)=sqrt(x.^2+y.^2)/maximum
            IIItheta(j,i)=atan2(y,x)
        end
        end
        for i=1:(lx/2)
        for j=1:(ly/2)
            x=-i
            y=-j
            IVro(j,i)=sqrt(x.^2+y.^2)/maximum
            IVtheta(j,i)=atan2(y,x)
        end
        end
        %perevora4ivaem 4etverti
        tempr=IIIro
        tempt=IIItheta
        for i=1:lx/2
        for j=1:ly/2
            IIIro(j,i)=tempr(round(ly/2)+1-j,round(lx/2)+1-i)
            IIItheta(j,i)=tempt(round(ly/2)+1-j,round(lx/2)+1-i)
        end
        end
        tempr=Iro
        tempt=Itheta
        
        for i=1:lx/2
        for j=1:ly/2
            Iro(j,i)=tempr(round(ly/2)+1-j,i)
            Itheta(j,i)=tempt(round(ly/2)+1-j,i)
        end
        end
        tempr=IVro
        tempt=IVtheta
        for i=1:lx/2
        for j=1:ly/2
            IVro(j,i)=tempr(j,round(lx/2)+1-i)
            IVtheta(j,i)=tempt(j,round(lx/2)+1-i)
        end
        end
        clear tempr tempt 
        %forming ro-theta matrix together
        for i=1:(lx/2)
        for j=1:(ly/2)
            ro(j+round(ly/2),i+round(lx/2))=IIro(j,i)
            ro(j,i+round(lx/2))=Iro(j,i)
            ro(j,i)=IIIro(j,i)
            ro(j+round(ly/2),i)=IVro(j,i)
        end
        end
        for i=1:(lx/2)
        for j=1:(ly/2)
            theta(j+round(ly/2),i+round(lx/2))=IItheta(j,i)
            theta(j,i+round(lx/2))=Itheta(j,i)
            theta(j,i)=IIItheta(j,i)
            theta(j+round(ly/2),i)=IVtheta(j,i)
        end
        end
        clear Iro IIro IIIro IVro Itheta IVtheta IIItheta IItheta
        %vi4islyaem polynomials
        for i=1:(lx)
        for j=1:(ly)
            %prep
            tro=ro(j,i)
            ttheta=theta(j,i)
            ro2=tro*tro
            ro3=ro2*tro
            ro4=ro2*ro2
            ro5=ro4*tro
            ro6=ro4*ro2
            ro8=ro6*ro2
            %polynomials
            P0(j,i)=0
            P1(j,i)=tro*cos(ttheta)
            P2(j,i)=tro*sin(ttheta)
            P3(j,i)=2*ro2-1
            P4(j,i)=ro2*cos(2*ttheta)
            P5(j,i)=ro2*sin(2*ttheta)
            P6(j,i)=(3*ro2-2)*tro*cos(ttheta)
            P7(j,i)=(3*ro2-2)*tro*sin(ttheta)
            P8(j,i)=6*ro4-6*ro2+1
            P9(j,i)=ro3*cos(3*ttheta)
            P10(j,i)=ro3*sin(3*ttheta)
            P11(j,i)=(4*ro2-3)*ro2*cos(2*ttheta)
            P12(j,i)=(4*ro2-3)*ro2*sin(2*ttheta)
            P13(j,i)=(10*ro4-12*ro2+3)*tro*cos(ttheta)
            P14(j,i)=(10*ro4-12*ro2+3)*tro*sin(ttheta)
            P15(j,i)=20*ro6-30*ro4+12*ro2-1
            P16(j,i)=ro4*cos(4*ttheta)
            P17(j,i)=ro4*sin(4*ttheta)
            P18(j,i)=(5*ro2-4)*ro3*cos(3*ttheta)
            P19(j,i)=(5*ro2-4)*ro3*sin(3*ttheta)
            P20(j,i)=(15*ro4-20*ro2+6)*ro2*cos(2*ttheta)
            P21(j,i)=(15*ro4-20*ro2+6)*ro2*sin(2*ttheta)
            P22(j,i)=(35*ro6-60*ro4+30*ro2-4)*tro*cos(ttheta)
            P23(j,i)=(35*ro6-60*ro4+30*ro2-4)*tro*sin(ttheta)
            P24(j,i)=70*ro8-140*ro6+90*ro4-20*ro2+1
            P25(j,i)=ro5*cos(5*ttheta)
            P26(j,i)=ro5*sin(5*ttheta)
            P27(j,i)=(6*ro2-5)*ro4*cos(4*ttheta)
            P28(j,i)=(6*ro2-5)*ro4*sin(4*ttheta)
            P29(j,i)=(21*ro4-30*ro2+10)*ro3*cos(3*ttheta)
            P30(j,i)=(21*ro4-30*ro2+10)*ro3*sin(3*ttheta)
            P31(j,i)=(56*ro6-105*ro4+60*ro2-10)*ro2*cos(2*ttheta)
            P32(j,i)=(56*ro6-105*ro4+60*ro2-10)*ro2*sin(2*ttheta)
            P33(j,i)=(126*ro8-280*ro6+210*ro4-60*ro2+5)*tro*cos(ttheta)
            P34(j,i)=(126*ro8-280*ro6+210*ro4-60*ro2+5)*tro*sin(ttheta)
            P35(j,i)=252*ro8*ro2-630*ro8+560*ro6-210*ro4+30*ro2-1
            P36(j,i)=924*ro8*ro4-2772*ro8*ro2+3150*ro8-1680*ro6+420*ro4-42*ro2+1
        end
        end
        %LSES
        amo=1
        for i=1:lx
        for j=1:ly
            tempmatp(amo,:)=[P0(j,i), P1(j,i), P2(j,i), P3(j,i), P4(j,i), P5(j,i), P6(j,i), P7(j,i),...
                P8(j,i), P9(j,i), P10(j,i), P11(j,i), P12(j,i), P13(j,i), P14(j,i), P15(j,i), P16(j,i),...
                P17(j,i), P18(j,i), P19(j,i), P20(j,i), P21(j,i), P22(j,i), P23(j,i), P24(j,i), P25(j,i),...
                P26(j,i), P27(j,i), P28(j,i), P29(j,i), P30(j,i), P31(j,i), P32(j,i), P33(j,i), P34(j,i), P35(j,i), P36(j,i)]
            tempmats(amo,:)=SOURCE(j,i)
            amo=amo+1
            if (amo==(lx*ly) ) break 
            end
        end
            if (amo==(lx*ly) ) break
            end
        end  
        C=linsolve(tempmatp,tempmats)
        clear tempmatp tempmats
        %changing coeffs - decoment and put number instead of X to change coef
        
        C(2,1)=0
        
        %forming result matrix
        for i=1:(lx)
        for j=1:(ly)
            RESULT(j,i)=C(1,1)+C(2,1)*P1(j,i)+C(3,1)*P2(j,i)...
                +C(4,1)*P3(j,i)+C(5,1)*P4(j,i)+C(6,1)*P5(j,i)+C(7,1)*P6(j,i)+...
                C(8,1)*P7(j,i)+C(9,1)*P8(j,i)+C(10,1)*P9(j,i)+C(11,1)*P10(j,i)...
                +C(12,1)*P11(j,i)+C(13,1)*P12(j,i)+C(14,1)*P13(j,i)+C(15,1)*P14(j,i)...
                +C(16,1)*P15(j,i)+C(17,1)*P16(j,i)+C(18,1)*P17(j,i)+C(19,1)*P18(j,i)...
                +C(20,1)*P19(j,i)+C(21,1)*P20(j,i)+C(22,1)*P21(j,i)+C(23,1)*P22(j,i)...
                +C(24,1)*P23(j,i)+C(25,1)*P24(j,i)+C(26,1)*P25(j,i)+C(27,1)*P26(j,i)...
                +C(28,1)*P27(j,i)+C(29,1)*P28(j,i)+C(30,1)*P29(j,i)+C(31,1)*P30(j,i)...
                +C(32,1)*P31(j,i)+C(33,1)*P32(j,i)+C(34,1)*P33(j,i)+C(35,1)*P34(j,i)...
                +C(36,1)*P35(j,i)+C(37,1)*P36(j,i)
        end
        end
        %surf(SOURCE)
        surf(RESULT)