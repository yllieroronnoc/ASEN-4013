%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This file makes CEA rocket case input file for the CEA600
%  Created by:  Yu Matsutomi, April 29, 2003
%  Modified by: Kevin Miller, August 27, 2003
%  Modified " " April 20, 2004 to allow for "exploded" reactant input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cea_input_gen(ox1,ox2,ox1wt,ox2wt,ox1T,ox2T,ox1st,ox2st,fu1,fu2,fu1wt,fu2wt,fu1T,fu2T,fu1st,fu2st,...
                       ox1hf,ox2hf,fu1hf,fu2hf,Pc,OF,PR,subar,supar,CR,flow,out)
% Inputs:
% ox1   = primary oxidizer (string)
% ox2   = secondary oxidizer (string)
% ox1wt = wt. fraction of primary oxidizer (%mass)
% ox2wt = wt. fraction of secondary oxidizer (%mass)
% ox1T  = initial (injection) temperature of primary oxidizer (degR)
% ox2T  = initial (injection) temperature of secondary oxidizer (degR)
% fu1   = primary fuel (string)
% fu2   = secondary oxidizer (string)
% fu1wt = wt. fraction of primary fuel (%mass)
% fu2wt = wt. fraction of secondary oxidizer (%mass)
% fu1T  = initial (injection) temperature of primary fuel (degR)
% fu2T  = initial (injection) temperature of secondary fuel (degR)
% Pc    = chamber pressure (psia)
% OF    = mixture ratio (total mass oxidizer/total mass fuel)
% PR    = pressure ratios (1)
% supar = supersonic area ratios (1)
% subar = subsonic area ratios (1)
% CR    = contraction ratio (1)
% flow  = type of flow assumed (frozen (fz) or equilibrium (eq))
% out   = array of output quantities

% check plot file matrix size
maxrows = 100;
numrows = (length(Pc)*(2+length(PR)+length(subar)+length(supar)));
if (numrows > maxrows)
    warn = 'Plot file will be too large, reduce number of points to calculate (areas, pressures)';
    disp(warn)
end

% Open input file in order to write to it
FID = fopen('Detn.inp','w');   %This file name is hardcoded into the CEA executable I am running

fprintf(FID,'\n');
fprintf(FID,'input thermo');              %mixture ratio
fprintf(FID,'\n');

% Problem namelist:
fprintf(FID,'problem  exploded %s ',flow);     %problem type and flow typr
if (CR>0)
    fprintf(FID,'fac ');                     %finite chamber
end

fprintf(FID,'\n');
fprintf(FID,'o/f=%9.4f ',OF);              %mixture ratio
fprintf(FID,'\n');

if (CR>0)
    fprintf(FID,'ac/at=%9.4f ',CR);          %finite chamber
end

fprintf(FID,'p,psia=');                      %chamber pressures
for i = 1:length(Pc)
   fprintf(FID,'%9.2f  ',Pc(i));
end
fprintf(FID,'\n');

if(isempty(PR)==0)
    fprintf(FID,'pi/p=');                    %pressure ratios
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
fprintf(FID,'reac  name %s%s wt%%=%3f  t(r)=%7.2f  h,cal=%10.1f\n',ox1,ox1st,ox1wt,ox1T,ox1hf);
if ox2wt~=0
fprintf(FID,'      name %s%s wt%%=%3f  t(r)=%7.2f  h,cal=%10.1f\n',ox2,ox2st,ox2wt,ox2T,ox2hf);
end

fprintf(FID,'      name %s%s wt%%=%3f  t(r)=%7.2f  h,cal=%10.1f\n',fu1,fu1st,fu1wt,fu1T,fu1hf);
if fu2wt~=0
%fprintf(FID,'      name %s %s wt%%=%3f  t(r)=%7.2f  h,cal=%10.1f\n',fu2,fu2st,fu2wt,fu2T,fu2hf);
fprintf(FID,'      fuel MAT(S) Mn 1.0 C 4.0 H 14.0 O 8.0 wt%%=%3f  t(r)=%7.2f  h,cal=%10.1f\n',fu2wt,fu2T,fu2hf);
end
fprintf(FID,'\n');


% Output options:
% joules instead of calories (calories)
% mole fractions of mass fractions (massf)
% return transport properties
% plot file includes: area ratios, temps, pressures, etc.
fprintf(FID,'output calories transport molef plot '); 
fprintf(FID,'%s \n',out);

% fprintf(FID,'\n');
% fprintf(FID,'insert  MnO(L)');              %mixture ratio
% fprintf(FID,'\n\n');


% End input file
fprintf(FID,'end');
fclose(FID);
