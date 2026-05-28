classdef TwoPort < handle
    properties 
        ABCD = {};
        abcd = [];
        n = 0;
        f = [];
    end

    methods (Access = 'public')
        function self = TwoPort(f)
            %Add frequency in Hz
            self.f=f(:);
        end

        function Y = ABCD_to_Y(self,abcd)
            Y=zeros(size(abcd));
            dabcd=abcd(:,1).*abcd(:,4)-abcd(:,2).*abcd(:,3);
            Y(:,1)=abcd(:,4)./abcd(:,2);
            Y(:,2)=-dabcd./abcd(:,2);
            Y(:,3)=-1./abcd(:,2);
            Y(:,4)=abcd(:,1)./abcd(:,2);
        end

        function abcd = Y_to_ABCD(self,Y)
            abcd=zeros(size(Y));
            abcd(:,1)=-1.*Y(:,4)./Y(:,3);
            abcd(:,2)=-1./Y(:,3);
            abcd(:,3)=-1.*(Y(:,1).*Y(:,4)-Y(:,2).*Y(:,3))./Y(:,3);
            abcd(:,4)=-Y(:,1)./Y(:,3);
        end

        function abcd=Z_to_ABCD(self,Z)
            abcd=zeros(size(Z));
            abcd(:,1)=Z(:,1)./Z(:,3);
            abcd(:,2)=(Z(:,1).*Z(:,4) - Z(:,2).*Z(:,3))./Z(:,3);
            abcd(:,3)=1./Z(:,3);
            abcd(:,4)=Z(:,4)./Z(:,3);
        end

        function Z = ABCD_to_Z(self,abcd)
            Z=zeros(size(abcd));
            dabcd=abcd(:,1).*abcd(:,4)-abcd(:,2).*abcd(:,3);
            Z(:,1)=abcd(:,1)./abcd(:,3);
            Z(:,2)=-dabcd./abcd(:,3);
            Z(:,3)=1./abcd(:,3);
            Z(:,4)=abcd(:,4)./abcd(:,3);
        end

        function Z = Zsc(self,abcd)
            Z=abcd(:,2)./abcd(:,4);
        end

        function Z = Zoc(self,abcd)
            Z=abcd(:,1)./abcd(:,3);            
        end

        function abcd = cascade(self)
            abcd=self.ABCD{1};

            if self.n==1
                self.abcd=abcd;
                return
            end

            for i=2:self.n
                next=self.ABCD{i};
                abcd_f(:,1)=abcd(:,1).*next(:,1)+abcd(:,2).*next(:,3);
                abcd_f(:,2)=abcd(:,1).*next(:,2)+abcd(:,2).*next(:,4);
                abcd_f(:,3)=abcd(:,3).*next(:,1)+abcd(:,4).*next(:,3);
                abcd_f(:,4)=abcd(:,3).*next(:,2)+abcd(:,4).*next(:,4);
    
                abcd=abcd_f;
            end
            self.abcd=abcd;
        end

        function abcd = parallel_C(self,value)
            abcd=repmat([1 0 -99 1],length(self.f),1);
            abcd(:,3)=1i*self.f*2*pi*value;
        end

        function abcd = parallel_L(self,value)
            abcd=repmat([1 0 -99 1],length(self.f),1);
            abcd(:,3)=1./(1i*self.f*2*pi*value);
        end

        function abcd = parallel_r(self,value)
            abcd=repmat([1 0 1./value 1],length(self.f),1);
        end

        function abcd = parallel_series_lc(self,L,C)
            %Capacitance cannot be 0
            abcd=repmat([1 0 -99 1],length(self.f),1);
            XL=(1i*self.f*2*pi*L);
            XC=1./(1i*self.f*2*pi*C);
            value = XL+XC;
            abcd(:,3)=1./value;
        end             

        function abcd = parallel_series_rlc(self,R,L,C)
            %Capacitance cannot be 0
            abcd=repmat([1 0 -99 1],length(self.f),1);
            XL=(1i*self.f*2*pi*L);
            XC=1./(1i*self.f*2*pi*C);
            value = R+XL+XC;
            abcd(:,3)=1./value;
        end             

        function abcd = series_C(self,value)
            abcd=repmat([1 -99 0 1],length(self.f),1);
            abcd(:,2)=1./(1i*self.f*2*pi*value);
        end

        function abcd = series_L(self,value)
            abcd=repmat([1 -99 0 1],length(self.f),1);
            abcd(:,2)=(1i*self.f*2*pi*value);
        end

        function abcd = series_r(self,value)
            abcd=repmat([1 value 0 1],length(self.f),1);
        end

        function abcd = series_parallel_rlc(self,R,L,C)
            abcd=repmat([1 -99 0 1],length(self.f),1);
            %Values cannot be zero - not check
            XL=(1i*self.f*2*pi*L);
            XC=1./(1i*self.f*2*pi*C);
            value = 1./ ( 1./R + 1./XL + 1./XC);
            abcd(:,2)=value;
        end

        function abcd = series_parallel_rl(self,R,L)
            %Values cannot be zero - not check
            abcd=repmat([1 -99 0 1],length(self.f),1);
            XL=(1i*self.f*2*pi*L);
            value = 1./ ( 1./R + 1./XL);
            abcd(:,2)=value;
        end

        function abcd = series_parallel_rc(self,R,C)
            abcd=repmat([1 -99 0 1],length(self.f),1);
            %Values cannot be zero - not check
            XC=1./(1i*self.f*2*pi*C);
            value = 1./ ( 1./R + 1./XC);
            abcd(:,2)=value;
        end

        function abcd = series_parallel_lc(self,L,C)
            abcd=repmat([1 -99 0 1],length(self.f),1);
            %Values cannot be zero - not check
            XL=(1i*self.f*2*pi*L);
            XC=1./(1i*self.f*2*pi*C);
            value = 1./ ( 1./XL + 1./XC);
            abcd(:,2)=value;
        end        

        function abcd = transformer(self,n)
            abcd=repmat([n 0 0 1/n],length(self.f),1);
        end

        function self = addComponent(self,abcd)
            self.n=self.n+1;
            self.ABCD{self.n}=abcd;
        end

        function [RL_dB, ML_dB, Pnet_dB, PL_dB] = networkLoss(self, Zs, abcd, ZL)
            Vin=1;
            A=abcd(:,1); B=abcd(:,2); C=abcd(:,3); D=abcd(:,4);
            if(~iscolumn(ZL)) ZL=ZL.'; end

            Zin = (A.*ZL+B)./(C.*ZL+D); %Input impedance            
            %Assume we have a voltage at the input, not at source
            %Vin = Vs.*Zin./(Zin+Zs);    %Voltage at input to network
            VL  = Vin./(A+B./ZL);       %Voltage at Load
            IL  = VL./ZL;               %Current at Load

            PL=0.5*real(VL.*conj(IL));  %Power dissipated in Load
            Iin=Vin./Zin;                 %Current delivered by source
            Pin=0.5*real(Vin.*conj(Iin));  %Power delivered by source

            S11 =(Zin-Zs)./(Zin+Zs);
            RL_dB=20*log10(abs(S11));
            ML = 1-abs(S11).^2;
            ML_dB=10*log10(ML);

            Pnet=Pin-PL;                   %Power in network
            Ps=Pin./ML;                    %Power Delivered by source
            Pr=Ps-Pin;                     %Power Reflected
            Pnet_dB=10*log10(Pnet./Ps);
            PL_dB=10*log10(PL./Ps);
        end

    end
end