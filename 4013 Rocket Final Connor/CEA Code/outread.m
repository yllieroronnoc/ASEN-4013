%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This file reads output from single Detonation solution from the CEC600
%
% Yu Matsutomi, April 29, 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [of, mol1, gam1, son1, P2, T2, rho2, Cp2, gam2, son2, Mdet, Vdet, mol2] = outread(dummy)
% open CEC600 output file
fid1=fopen('Detn.out','r');
Cpcount = 1;
while 1
    
    % read in single line from file
    s1=fgetl(fid1);
    
    % if EOF encountered exit loop
    
        if ~ischar(s1), break, end

    % if the current line is not blank perform processing operations
    
    if (length(deblank(s1)) ~= 0)
   
    % look for line that lists the O/F ratio this defines the beginning of a zone
    % for either equilibrium or frozen flow conditions
    
        if (strncmpi(s1,' O/F=',5)==1)
        of=str2num(s1(6:16));
        end
    
    t1=findstr(s1,'M1, (1/n)');
    t2=findstr(s1,'GAMMA1');
    t3=findstr(s1,'SON VEL1,M/SEC');
    t4=findstr(s1,'P, ATM');
    t5=findstr(s1,'T, K');
    t6=findstr(s1,'RHO, G/CC');
    t7=findstr(s1,'Cp, CAL/(G)(K)');
    t8=findstr(s1,'GAMMAs');
    t9=findstr(s1,'SON VEL,M/SEC');
    t10=findstr(s1,'DET MACH NUMBER');
    t11=findstr(s1,'DET VEL,M/SEC');
    t12=findstr(s1,'M, (1/n)');

      
        % extract molecular weight of unburned gas molecular weight
        if (~isempty(t1))
        dat = sscanf(s1,'%s %s %f');
        mol1=dat(length(dat));  %symbol is alphabet mol and number 1.
        end
    
        % extract gamma of unburned gas
        if (~isempty(t2))
        dat = sscanf(s1,'%s %f');
        gam1=dat(length(dat));   
        end
    
        % extract sonic velocity of unburned gas
        if (~isempty(t3))
        dat = sscanf(s1,'%s %s %f');
        son1=dat(length(dat));    
        end
    
        % extract pressure of burned gas
        if (~isempty(t4))
        dat = sscanf(s1,'%s %s %f');
        P2=dat(length(dat));
        end
    
        % extract temperature of burned gas
        if (~isempty(t5))
        dat = sscanf(s1,'%s %s %f');
        T2=dat(length(dat));
        end
        
        % extract density of burned gas
        % this extract is very unstable.  the code assume the output is in form of 4.00-3, for 4.00e-3
        if (~isempty(t6))
        dat = sscanf(s1,'%s %s %f %f');
        rho2=dat(length(dat)-1)*10^dat(length(dat));  
        end
    
        % extract Cp of unburned gas
        if (~isempty(t7) & Cpcount)
        dat = sscanf(s1,'%s %s %f');
        Cp2=dat(length(dat)); 
        Cpcount = 0;
        end
    
        % extract gamma of burned gas
        if (~isempty(t8))
        dat = sscanf(s1,'%s %f');
        gam2=dat(length(dat));    
        end
    
        % extract sonic velocity of burned gas
        if (~isempty(t9))
        dat = sscanf(s1,'%s %s %f');
        son2=dat(length(dat));
        end
    
        % extract detonation mach number
        if (~isempty(t10))
        dat = sscanf(s1,'%s %s %s %f');
        Mdet=dat(length(dat));
        end
        
        % extract detonation velocity
        if (~isempty(t11))
        dat = sscanf(s1,'%s %s %f');
        Vdet=dat(length(dat));
        end
        
        % extract burned gas molecular weight
        if (~isempty(t12))
        dat = sscanf(s1,'%s %s %f');
        mol2 = dat(length(dat));
        end
    end
end

fclose(fid1);

%%% convert si unit to english
%%% Cp will remain cal/(g-K)
son1 = 3.28084*son1;
P2 = 14.695949*P2;
T2 = 9/5*T2;
rho2 = 62.427962*rho2;
son2 = 3.28084*son2;
Vdet = 3.28084*Vdet;
