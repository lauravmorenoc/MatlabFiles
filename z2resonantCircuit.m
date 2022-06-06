function [R, L, C] = z2resonantCircuit(varargin)
    %{
    returns R, L and C values and (op) prints circuit description
    [R, L, C]=z2resonantCircuit(z, f), where z is the complex value of impedance and f the frequency range in hertz
    [R, L, C]=z2resonantCircuit(z, f, 'no print') for keeping the function from printing in console
    %} 
    
    [aux1, aux2, aux3]=deal(0,0,0);
    z=varargin{1};
    freq=varargin{2};
   	x=imag(z);
    type='nresonant';
    fc=0;
    
    if length(z)~=length(freq)
        fprintf('Z and frequency vectors must have the same length')        
    else
   		found=find(x==0);
        
   		if ~isempty(found)
            fc=freq(found(1));
            %zc=z(found(1));
            pos=found(1);
        else
            
       		for k=2:length(z) 
                if fc==0
                    if ((x(k)>0)&&(x(k-1)<0))||((x(k)<0)&&(x(k-1)>0))
                        fc=freq(k)
                        zc=z(k);
                        pos=k;
                        type='yresonant';
                    end
                end
            end
        end
        if type=='yresonant'
            R=abs(real(zc));
            if fc<freq(ceil(length(freq)/2))
                f2=freq(pos+ceil(length(freq)/2));
                x2=x(pos+ceil(length(freq)/2));
            else
                f2=freq(pos-ceil(length(freq)/5));
                x2=x(pos-ceil(length(freq)/5));
            end
           
 
     		syms L C
        	eq1=fc==1/(2*pi*sqrt(L*C));
            
            if (nargin>2)&&(length(varargin{3})==1)&&(varargin{3}=='X')
                eq2=x2==(2*pi*f2*L)*(-1/(2*pi*f2*C))/(2*pi*f2*L-(1/(2*pi*f2*C)));
                circ='X';
            elseif (nargin>2)&&(length(varargin{3})==1)&&(varargin{3}=='R') % Resistance model
                eq1=(2*pi*fc)^2*L*C==1/(1+(1/((2*pi*fc*C*R)^2))); % Resonance changes
                eq2=x2==(-1/(2*pi*f2*C))*R^2/(R^2+(1/(2*pi*f2*C))^2)+2*pi*f2*L;
                circ='R';
            elseif (nargin>2)&&(length(varargin{3})==1)&&(varargin{3}=='A')
                eq1=(2*pi*fc)^2*L*C==((2*pi*fc*L)^2+R^2)/(R^2); % Resonance changes
                eq2=x2==(-1/(2*pi*f2*C))+2*pi*f2*L*R^2/(R^2+(2*pi*f2*L)^2);
                %eq=x2==(-1/(2*pi*f2*C))*R^2/(R^2+(1/(2*pi*f2*C))^2);
                circ='A';
            elseif (nargin>2)&&(length(varargin{3})==1)&&(varargin{3}=='P')
                XC=-1/(2*pi*f2*C);
                XL=2*pi*f2*L;
                Z1=XL*XC/(XL+XC);
                Z2=Z1*R/(Z1+R);
                eq2=x2==imag(Z2);
                circ='P';
            else
                eq2=x2==2*pi*f2*L-(1/(2*pi*f2*C));
                circ='S';
            end
            
            
            assume(L>0);
            assume(C>0);
            S=solve([eq1,eq2], [L,C],'ReturnConditions',true);
            %S=solve(eq, C,'ReturnConditions',true);
            
            L=subs(S.L);
            %L=0;
            C=subs(S.C);
            
           
        else
            fprintf('No resonant')
            % The circuit doesn't have resonance in the freq range
            f1=freq(ceil(length(freq)/3));
            x1=x(ceil(length(freq)/3));
            f2=freq(ceil(2*length(freq)/3));
            x2=x(ceil(2*length(freq)/3));
            
            syms L C
            XL=2*pi*freq*L;
            XC=1/(2*pi*freq*C);
            eq1=x1==2*pi*f1*L-(1/(2*pi*f1*C));
            eq2=x2==2*pi*f2*L-(1/(2*pi*f2*C));
            if all(x>0)
                assume(XL>XC);
            else
                assume(XL<XC);
            end
            S=solve([eq1, eq2],[L, C],'ReturnConditions',true); 
            L=double(S.L);
            C=double(S.C);
            R=(real(z(ceil(length(freq)/3)))+real(z(ceil(2*length(freq)/3))))/2;
            
        end
        
        [aux1, aux2, aux3]=deal(R, L, C);
        
        if (nargin==3)&&(length(varargin{3})==8)&&varargin{3}~='no print'
            fprintf('The third parameter must be "no print" to keep the function from printing in console')
        else
            if circ=='S'
                fprintf('\nSeries resonant circuit found:\n')
            elseif circ=='X'
                fprintf('\nResistance in series with parallel capacitor and inductor circuit found:\n')
            elseif circ=='R'
                fprintf('\nInductor in series with parallel capacitor and resistor circuit found:\n')
            elseif circ=='P'
                fprintf('\nRLC parallel circuit found:\n')
            elseif circ=='A'
                fprintf('\nCapacitor in series with parallel inductor and resistor circuit found:\n')
            end
            fprintf('R=%g, L=%g, C=%g\n', R, L, C)
        end
        %}
    end
        
    [R, L, C]=deal(aux1, aux2, aux3);
    
end