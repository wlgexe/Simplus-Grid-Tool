% This class defines the model of 8-th order synchronous machine

% Author(s): Yue Zhu
% Reference Book: Power System Dynamic and Stability, P54

%% Class

classdef SynchronousMachineFull_SM < SimplusGT.Class.ModelAdvance
    
    properties(Access = protected)
        psi_f;
        ws;
        L;

        avr_i;%indicator
        gov_i;
        pss_i;
    end
    
    methods
        % constructor
        function obj = SynchronousMachineFull_SM(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        %% signal list
        
        function [State,Input,Output] = SignalList(obj)
            type=obj.ApparatusType-60; % decoder, type= PSS_AVR_GOV, binary
            if mod(type,2)==1
                obj.gov_i=1;
            else 
                obj.gov_i=0;
            end
            type=floor(type/2);
            if mod(type,2)==1
                obj.avr_i=1;
            else 
                obj.avr_i=0;
            end
            type=floor(type/2);
            if mod(type,2)==1
                obj.pss_i=1;
            else 
                obj.pss_i=0;
            end
            State = {'id', 'iq', 'w', 'theta','Ed1','Eq1','Psi1d','Psi2q'};
            if obj.avr_i==1
                State=[State,{'Efd','Vr','Rf','Va','Vx1','Vx2'}];
            end
            if obj.pss_i==1
                State=[State,{'Vss1','Vss2','Vss3','Vss'}];
            end
            if obj.gov_i==1
                State=[State,{'G2', 'G3', 'T_m'}];
            end
                Input	 = {'v_d','v_q','T_m','Efd'};
            if obj.avr_i==1 && obj.pss_i==0
                Input	 = {'v_d','v_q','T_m','Efd','Vss'};
            end
                Output = {'i_d','i_q','w','theta'};
        end
        %% Equilibrium point
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);            
            % Get parameters
            % Synchronous machine
            X=obj.Para(1);
            R=obj.Para(2);
            Xd=obj.Para(3); %synchronous reactance in d axis
            Xd1=obj.Para(4); %transient reactance
            Xd2=obj.Para(5); %subtransient reactance
            Td1=obj.Para(6); %d-axis open circuit transient time constant
            Td2=obj.Para(7); %d-axis open circuit sub-transient time constant
            Xq=obj.Para(8);
            Xq1=obj.Para(9);
            Xq2=obj.Para(10);
            Tq1=obj.Para(11);
            Tq2=obj.Para(12);
            H=obj.Para(13);
            D=obj.Para(14);
            if obj.gov_i==1
                Rgov = obj.Para(42);
                T1gov = obj.Para(43);
                T2gov = obj.Para(44);
                T3gov = obj.Para(45);
                Dtgov = obj.Para(46);
            end
            if obj.avr_i==1
                TR=obj.Para(15);
                KA=obj.Para(16);
                TA=obj.Para(17);
                VRmax=obj.Para(18);
                VRmin=obj.Para(19);
                KE=obj.Para(20);
                TE=obj.Para(21);
                E1=obj.Para(22);
                SE1=obj.Para(23);
                E2=obj.Para(24);
                SE2=obj.Para(25);
                KF=obj.Para(26);
                TF=obj.Para(27);
                KP=obj.Para(28);
                KI=obj.Para(29);
                KD=obj.Para(30);
                TD=obj.Para(31);
                Bex=log(SE1/SE2)/(E1-E2);
                Aex=SE1*exp(-Bex*E1);
            end
            if obj.pss_i==1

            end
            % Calculate Equilibriums
            i_DQ = (conj(P+1j*Q)/V);
            Eint = V - i_DQ*(R+1j*Xq); %internal voltage
            sigma = angle(Eint);           
            i_dq = i_DQ * exp(1j*(-sigma+pi/2));
            i_d = real(i_dq);
            i_q = imag(i_dq);
            v_dq = V * exp(1j*(-sigma+pi/2));
            v_d = real(v_dq);
            v_q = imag(v_dq);           
            Efd = abs(Eint)+(Xq-Xd)*i_d;
            Eq1 = Efd+(Xd-Xd1)*i_d;           
            Ed1 = -(Xq-Xq1)*i_q;
            Psi1d=Eq1+(Xd1-X)*i_d;
            Psi2q=-Ed1+(Xq1-X)*i_q;

            obj.ws = w; % record synchronous rotor speed.
            obj.L = X/w; % record inductor value

            Psi_q = R*i_d-v_d;
            Psi_d = -R*i_q+v_q;            
            T_e = Psi_d*i_q-Psi_q*i_d;%Psi_q*i_d-Psi_d*i_q;
            T_m = T_e;
            if obj.gov_i==1 
                G3=T_m;
                G2=G3;
                G1=G2;
                Tref=G1*Rgov;
            end
            Vss1=0;
            Vss2=0;
            Vss3=0;
            Vss=0;
            if obj.avr_i==1 
                Vt=abs(v_d+1j*v_q);
                Vref=Vt;
                Vr = Vt;
                Va=KE*Efd+Efd*Aex*exp(Bex*Efd);  % not 0
                if Va>=Vt*VRmax
                    Va = Vt*VRmax;
                elseif Va<= Vt*VRmin
                    Va = Vt*VRmin;
                end  
                Rf=-KF/TF*Efd;  % not 0
                Vf=Rf+KF/TF*Efd; %=0
                Vx=Vss+Vref-Vr-Vf; %=0
                Vx2=-KD/TD*Vx; %=0
                Vpid=Va/KA; % not 0
                if Vpid>=VRmax
                    Vpid = VRmax;
                elseif Vpid<= VRmin
                    Vpid = VRmin;
                end            
                Vx1=Vpid-Vx2-Vx*(KP+KD/TD); % not 0
                % obj.Vpid0=Vpid; % save the steady value for later saturation.
                % obj.Va0=Va;
            end
            xi = xi+sigma-pi/2;
            theta = xi;
           
            % Get equilibrium
            x_e = [i_d; i_q; w; theta; Ed1; Eq1; Psi1d; Psi2q];
            u_e = [v_d; v_q; T_m; Efd];
            % Get equilibrium
            
            if obj.avr_i==1
                x_e=[x_e;[Efd; Vr; Rf; Va; Vx1; Vx]];
                u_e(4)=Vref;
                if obj.pss_i==0
                    u_e(5)=Vss;
                end
            end
            if obj.pss_i==1
                x_e=[x_e; Vss1; Vss2; Vss3; Vss];
            end
            if obj.gov_i==1
                x_e=[x_e;[G2; G3; T_m]];
                u_e(3)=Tref;
            end
            %u_e = [v_d; v_q; Tref; Vref];
            % Get equilibrium
            xi  = [xi];    %connection angle to the whole grid.
        end
        
        %% State-space
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
                        % Get parameters
            L=obj.L;
            ws=obj.ws;
            X=obj.Para(1);
            R=obj.Para(2);
            Xd=obj.Para(3); %synchronous reactance in d axis
            Xd1=obj.Para(4); %transient reactance
            Xd2=obj.Para(5); %subtransient reactance
            Td1=obj.Para(6); %d-axis open circuit transient time constant
            Td2=obj.Para(7); %d-axis open circuit sub-transient time constant
            Xq=obj.Para(8);
            Xq1=obj.Para(9);
            Xq2=obj.Para(10);
            Tq1=obj.Para(11);
            Tq2=obj.Para(12);
            H=obj.Para(13);
            D=obj.Para(14);
            % if obj.Timer <= 5 && obj.Timer>0
            %     D = H*3;
            % else D=obj.Para(14);
            % end
            if obj.avr_i==1
                TR=obj.Para(15);
                KA=obj.Para(16);
                TA=obj.Para(17);
                VRmax=obj.Para(18);
                VRmin=obj.Para(19);
                KE=obj.Para(20);
                TE=obj.Para(21);
                E1=obj.Para(22);
                SE1=obj.Para(23);
                E2=obj.Para(24);
                SE2=obj.Para(25);
                KF=obj.Para(26);
                TF=obj.Para(27);
                KP=obj.Para(28);
                KI=obj.Para(29);
                KD=obj.Para(30);
                TD=obj.Para(31);
                Bex=log(SE1/SE2)/(E1-E2);
                Aex=SE1*exp(-Bex*E1);
            end
            if obj.pss_i==1
                KPSS=obj.Para(32);
                TW=obj.Para(33);
                T11=obj.Para(34);
                T12=obj.Para(35);
                T21=obj.Para(36);
                T22=obj.Para(37);
                T31=obj.Para(38);
                T32=obj.Para(39);
                VSSmax=obj.Para(40);
                VSSmin=obj.Para(41);
            end
            if obj.gov_i==1
                Rgov = obj.Para(42);
                T1gov = obj.Para(43);
                T2gov = obj.Para(44);
                T3gov = obj.Para(45);
                Dtgov = obj.Para(46);
            end
            % Get states
            i_d = x(1);
            i_q = x(2);
            w = x(3);
            theta = x(4);
            Ed1 = x(5);
            Eq1 = x(6);
            Psi1d = x(7);
            Psi2q = x(8);
            if obj.avr_i==1
                Efd = x(9);
                Vr = x(10);
                Rf = x(11);
                Va = x(12);
                Vx1= x(13);
                Vx2=x(14);
                if obj.pss_i==0
                    Vss = u(5);
                end
            end
            if obj.pss_i==1
                if obj.avr_i==1
                    Vss1=x(15);
                    Vss2=x(16);
                    Vss3=x(17);
                    Vss=x(18);
                else
                    Vss1=x(9);
                    Vss2=x(10);
                    Vss3=x(11);
                    Vss=x(12);
                end
                if Vss>=VSSmax
                    Vss = VSSmax;
                elseif Vss<=VSSmin
                    Vss = VSSmin;
                end
            end
            if obj.gov_i==1
                if obj.avr_i==1 && obj.pss_i==1
                    G2=x(19);
                    G3=x(20);
                    T_m=x(21);
                elseif obj.avr_i==0 && obj.pss_i==1
                    G2=x(13);
                    G3=x(14);
                    T_m=x(15);
                elseif obj.avr_i==1 && obj.pss_i==0
                    G2=x(15);
                    G3=x(16);
                    T_m=x(17);
                else
                    G2=x(9);
                    G3=x(10);
                    T_m=x(11);
                end
            end

            % Get input signals
            v_d  = u(1);
            v_q  = u(2);
            if obj.gov_i==1
                Tref  = u(3);
            else
                T_m  = u(3);
            end
            if obj.avr_i==1
                Vref  = u(4);
            else
                Efd  = u(4);
            end
            if CallFlag == 1        
            % ### Call state equation: dx/dt = f(x,u)
                % Auxiliary equations
                Psi_d = Xd2*i_d+(Xd2-X)/(Xd1-X)*Eq1 + (Xd1-Xd2)/(Xd1-X)*Psi1d; %book 
                Psi_q = Xq2*i_q-(Xq2-X)/(Xq1-X)*Ed1 + (Xq1-Xq2)/(Xq1-X)*Psi2q;
                T_e = Psi_d*i_q-Psi_q*i_d;%T_e = (Psi_q*i_d-Psi_d*i_q);
                ws=obj.ws;
                % SM
                dtheta = w;
                dw     = ws*(T_e - T_m - (D/ws)*(w-ws))/(2*H); % original
                dEd1 = (-Ed1+(Xq-Xq1)*(-i_q+(Xq1-Xq2)/(Xq1-X)^2*(-Psi2q+(Xq1-X)*i_q-Ed1)))/Tq1;
                dEq1 = (Efd-Eq1+(Xd-Xd1)*(i_d+(Xd1-Xd2)/(Xd1-X)^2*(Psi1d-(Xd1-X)*i_d-Eq1)))/Td1;                 
                dPsi1d = (-Psi1d+Eq1+(Xd1-X)*i_d)/Td2;
                dPsi2q = (-Psi2q-Ed1+(Xq1-X)*i_q)/Tq2;
                di_d = (ws*(v_d-R*i_d+w/ws*Psi_q)-(Xd2-X)/(Xd1-X)*dEq1-(Xd1-Xd2)/(Xd1-X)*dPsi1d)/Xd2;
                di_q = (ws*(v_q-R*i_q-w/ws*Psi_d)+(Xq2-X)/(Xq1-X)*dEd1-(Xq1-Xq2)/(Xq1-X)*dPsi2q)/Xq2;
                % PSS   %whash out time is negelected.
                if obj.pss_i==1
                    Sm=(w-ws)/ws;
                    dVss1=(-Vss1+Sm*KPSS+T11*KPSS*dw/ws)/T12;
                    dVss2=(-Vss2+Vss1+T21*dVss1)/T22;
                    dVss3=(-Vss3+Vss2+T31*dVss2)/T32;
                    dVss =(-Vss+TW*dVss3)/TW;
                end
                % Exciter     
                if obj.avr_i==1
                    Vf=Rf+KF/TF*Efd;
                    Vx=Vss+Vref-Vr-Vf;
                    Vpid=Vx1+Vx2+(KP+KD/TD)*Vx;
                    Vt=abs(v_d+1j*v_q);
                    if Va>=Vt*VRmax
                        Va = Vt*VRmax;
                    elseif Va<=Vt*VRmin
                        Va = Vt*VRmin;
                    end
                    if Vpid>=VRmax
                        Vpid = VRmax;
                    elseif Vpid<=VRmin
                        Vpid = VRmin;
                    end               
                    dEfd = (Va-(KE*Efd+Efd*Aex*exp(Bex*Efd)))/TE;
                    dVr  = (Vt-Vr)/TR;
                    dRf  = (-Rf-KF/TF*Efd)/TF;
                    dVa  = (-Va+KA*Vpid)/TA;
                    dVx1 = KI*Vx;
                    dVx2 = (-Vx2-KD/TD*Vx)/TD;
                end
                %governor
                if obj.gov_i==1
                    Delta_w = -(w - ws);
                    G1 = (Tref-Delta_w)/Rgov;
                    dG2 = (G1-G2)/T1gov;
                    dG3=(G2+T2gov*dG2-G3)/T3gov;
                    dT_m=(dG3-(-dw*Dtgov));
                end
                % end of state calculation
                f_xu = [di_d; di_q; dw; dtheta; dEd1; dEq1; dPsi1d; dPsi2q];
                if obj.avr_i==1
                    f_xu = [f_xu;dEfd; dVr; dRf; dVa; dVx1; dVx2];
                end
                if obj.pss_i==1
                    f_xu = [f_xu; dVss1; dVss2; dVss3; dVss];
                end
                if obj.gov_i==1
                    f_xu=[f_xu;[dG2; dG3; dT_m]];
                end
                Output = f_xu;
            elseif CallFlag == 2    
            % ### Call output equation: y = g(x,u)
                %i_ex = 0;
                g_xu = [i_d; i_q; w; theta];
                Output = g_xu;
            end                       
        end
        
    end

end     % End class definition
