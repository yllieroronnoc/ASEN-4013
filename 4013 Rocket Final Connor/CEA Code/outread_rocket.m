
function [of, T, Cp, gamma] = outread(dummy)

    fid1=fopen('Detn.out','r');

while 1
    
    % read in single line from file
    s1=fgetl(fid1);
    
    % if EOF encountered exit loop
    
    if ~ischar(s1), break, end

    % if the current line is not blank perform processing operations
    
    if (length(deblank(s1)) ~= 0)
   
    % look for line that lists the O/F ratio this defines the beginning of a zone
    % for either equilibrium or frozen flow conditions

%     if (strncmpi(s1,'O/F=',4)==1)
%         keyboard
% %        if (col == 2) 
% %            counter = counter - 1;
% %        end
% % 
% %        if (col > 2)
% %            col = 1;
% %        end
%        % construct the O/F data array
%        of=str2num(s1(6:16));
%        
%        
%     end
    
    t1=findstr(s1,'T, K  ');
    t2=findstr(s1,'Cp, KJ/(KG)(K)');
    t3=findstr(s1,'GAMMAs');
    t4=findstr(s1,'O/F =');
%     t4=findstr(s1,'I, LBF-SEC/LBM');
%     t5=findstr(s1,'RHO*I, LBF-SEC');
% 
%     a1=findstr(s1,'THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION DURING EXPANSION');
%     a2=findstr(s1,'THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM COMPOSITION DURING EXPANSION');

%     if (~isempty(a1))
%         disp(s1);
%     end
%     if (~isempty(a2))
%         disp(s1)
%     end
    
    % extract CSTAR @ nozzle exit from string
    if (~isempty(t1))
        dat = sscanf(s1,'%s %s %f %f');
        T=dat(length(dat)-1);
    end
    
    % extract RHO*IVAC @ nozzle exit from string
    if (~isempty(t2))
        dat = sscanf(s1,'%s %s %f %f');
        Cp=dat(length(dat)-1);   
    end
    
    % extract IVAC @ nozzle exit from string
    if (~isempty(t3))
        dat = sscanf(s1,'%s %f %f');
        gamma=dat(length(dat)-1);    
    end
    
    % extract Isp @ nozzle exit from string
    if (~isempty(t4))
        dat = sscanf(s1,'%s %s %f');
        of=dat(length(dat));
    end
%     
%     % extract rho*Isp @ nozzle exit from string
%     if (~isempty(t5))
%         disp(s1)
%         dat = sscanf(s1,'%s %s %f %f %f');
%         rhoI(counter,col)=dat(length(dat));
%     end
    end
end


fclose(fid1);

%%% convert si unit to english
%%% Cp will remain KJ/(Kg-K)
T = 9/5*T;  %  K to R
Cp = Cp*1.002*778.169; % [cal/g-K] to [ft-lb/lb-R]

