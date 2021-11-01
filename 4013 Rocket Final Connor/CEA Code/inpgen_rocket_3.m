%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This file makes CEA rocket case input file for the CEA600
%  Created by:  Yu Matsutomi, April 29, 2003
%  Modified by: Kevin Miller, August 27, 2003
%  Modified by: Steven Shark, September 17,2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function inpgen_rocket_3a(ox1,ox2,ox1wt,ox2wt,ox1T,ox2T,fu1,fu2,fu1wt,fu2wt,fu1T,fu2T,Pc,OF,PR,subar,supar,CR,flow,out)
function inpgen_rocket_3(ox1,ox2,ox1wt,ox2wt,ox1T,ox2T,ox1chem,ox2chem,ox1H,ox2H,fu1,fu2,fu1wt,fu2wt,fu1T,fu2T,fu1chem,fu2chem,fu1H,fu2H,Pc,OF,Phi,PR,subar,supar,CR,flow,out);
% Inputs:
% ox1     = primary oxidizer (string)
% ox2     = secondary oxidizer (string)
% ox1wt   = wt. fraction of primary oxidizer (%mass)
% ox2wt   = wt. fraction of secondary oxidizer (%mass)
% ox1T    = initial (injection) temperature of primary oxidizer (degR)
% ox2T    = initial (injection) temperature of secondary oxidizer (degR)
% ox1chem = primary oxidizer chemical formula if required (captilize chemical symbols)
% ox2chem = secondary oxidizer chemical formula if required (captilize chemical symbols)
% ox1H    = primary oxid enthalpy of formation [cal/mol]
% ox2H    = secondary oxid enthalpy of formation [cal/mol]
% fu1     = primary fuel (string)
% fu2     = secondary oxidizer (string)
% fu1wt   = wt. fraction of primary fuel (%mass)
% fu2wt   = wt. fraction of secondary oxidizer (%mass)
% fu1T    = initial (injection) temperature of primary fuel (degR)
% fu2T    = initial (injection) temperature of secondary fuel (degR)
% fu1che  = primary fuel chemical formula if requied (captilize chemical symbols)
% fu2chem = secondary fuel chemical formula if requied (captilize chemical symbols)
% fu1H    = primary fuel enthalpy of formation [cal/mol]
% fu2H    = secondary fuel enthalpy of formation [cal/mol]
% Pc      = chamber pressure (psia)
% OF      = mixture ratio (total mass oxidizer/total mass fuel)
% PR      = pressure ratios (1)
% supar   = supersonic area ratios (1)
% subar   = subsonic area ratios (1)
% CR      = contraction ratio (1)
% flow    = type of flow assumed (frozen (fz) or equilibrium (eq))
% out     = array of output quantities

% check plot file matrix size
maxrows = 100;
numrows = (length(Pc)*(2+length(PR)+length(subar)+length(supar)));
if (numrows > maxrows)
    warn = 'Plot file will be too large, reduce number of points to calculate (areas, pressures)';
    disp(warn)
end

% Open input file in order to write to it
FID = fopen('Detn.inp','w');   %This file name is hardcoded into the CEA executable Im running

% Problem namelist:
fprintf(FID,'problem  rocket %s\n ',flow);     %problem type and flow typr
if (CR>0)
    fprintf(FID,'fac ');                     %finite chamber
end            

if(isempty(OF)==0)
    fprintf(FID,'o/f=');                        %mixture ratio
    for i = 1:length(OF)
        fprintf(FID,'%9.4f,',OF(i));
    end
    fprintf(FID,'\n');
end

if(isempty(Phi)==0)
    fprintf(FID,'phi,eq.ratio=');                  %equivalence ratio
    for i = 1:length(Phi)
        fprintf(FID,'%9.4f,',Phi(i));
    end
    fprintf(FID,'\n');
end

if (CR>0)
    fprintf(FID,'ac/at=%9.4f ',CR);             %finite chamber
end

fprintf(FID,'\n');
fprintf(FID,'p,psia=');                      %chamber pressures
for i = 1:length(Pc)
   fprintf(FID,'%9.2f  ',Pc(i));
end
fprintf(FID,'\n');

if(isempty(PR)==0)
    fprintf(FID,'pi/p=');                        %pressure ratios
    for i = 1:length(PR)
        fprintf(FID,'%9.4f  ',PR(i));
    end
    fprintf(FID,'\n');
end

if(isempty(subar)==0)
    fprintf(FID,'subar=');                       %subsonic area ratios
    for i = 1:length(subar)
        fprintf(FID,'%9.4f  ',subar(i));
    end
    fprintf(FID,'\n');
end

if(isempty(supar)==0)
    fprintf(FID,'supar=');                       %supersonic area ratios
    for i = 1:length(supar)
        fprintf(FID,'%9.4f  ',supar(i));
    end
    fprintf(FID,'\n');
end

% Reactants namelist

fprintf(FID,'reactants \n');

if ox1wt~=0
    fprintf(FID,'reac  oxid %s  wt%%=%3f   t(r)=%7.2f    h(cal/mol)=%7.2f\n           %s\n',ox1,ox1wt,ox1T,ox1H,ox1chem);
end
if ox2wt~=0
    fprintf(FID,'      oxid %s  wt%%=%3f   t(r)=%7.2f    h(cal/mol)=%7.2f\n           %s\n',ox2,ox2wt,ox2T,ox2H,ox2chem);
end
if fu1wt~=0
    fprintf(FID,'      fuel %s  wt%%=%3f   t(r)=%7.2f    h(cal/mol)=%7.2f\n           %s\n',fu1,fu1wt,fu1T,fu1H,fu1chem);
end
if fu2wt~=0
    fprintf(FID,'      fuel %s  wt%%=%3f   t(r)=%7.2f    h(cal/mol)=%7.2f\n           %s\n',fu2,fu2wt,fu2T,fu2H,fu2chem);
end

fprintf(FID,'\n');

% Output options:
% joules instead of calories (calories)
% mole fractions of mass fractions (massf or leave blank!)
% return transport properties
% plot file includes: area ratios, temps, pressures, etc.
%fprintf(FID,'output calories transport massf plot '); 
fprintf(FID,'output transport massf plot '); 
fprintf(FID,'%s \n',out);

% End input file
fprintf(FID,'end');
fclose(FID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%   rocket problem sample format  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  problem  rocket  equilibrium  o/f=5.55157  
% case=8  p,bar=53.3172 subar=1.58,pi/p=10,100,1000,supar=25,50,75
%       reactants   
% fuel = H2(L)  wt% 100.   t(k) 20.27 
% oxid = O2(L)  wt% 100.   t(k) 90.17 
%       output  siunits
%       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%