% CEA input class for easier use of NASA Chemical Equilibrium Analysis
% 
% Drew Sherman
% Purdue University
% Made 4/12/2017
% 
% Mitch Woolever 
% University of Colorado Boulder
% Modified 4/20/2017
%
% Defines a CEA object with all the properties required to run CEA
%   (SEE BELOW)
%
% METHODS:
%   runCEA(obj): used by a CEA object to run CEA for the values (usually
%       defined in a parent code using '.' notation
%   
%   readCEAout(obj,location): reads the Detn.out file of CEA for the values the
%       user requested using obj.out command. The function returns a structure
%       with all relevant data from the CEA run


classdef CEAinput
   properties
      Filename = ''; %specify name for output file(.txt) 
        
        % specify oxidizer & fuel conditions inputs
        ox1   = '';             %% primary oxidizer
        ox2   = '';             %% secondary oxidizer
        fu1   = '';             %% primary fuel: RP-1
        fu2   = '';             %% secondary fuel
        ox1chem   = '';         %% optional primary oxidizer chemical formula if required (captilize chemical symbols)
        ox2chem   = '';         %% optional secondary oxidizer chemical formula if required (captilize chemical symbols)
        fu1chem   = '';         %% optional primary fuel chemical formula if requied (captilize chemical symbols)
        fu2chem   = '';         %% optional secondary fuel chemical formula if requied (captilize chemical symbols)

        ox1wt = 0;            %% wt fraction of primary oxid in total oxid [1]
        ox2wt = 0               %% wt fraction of secondary oxid in total oxid [1]
        fu1wt = 0;            %% wt fraction of primary fuel in total fuel [1]
        fu2wt = 0;              %% wt fraction of secondary fuel in total fuel [1]
        ox1T  = 0;              %% optional input of primary oxid temperature which enthalpy is evaluated [degR]
        ox2T  = 0;              %% optional input of secondary oxid temperature which enthalpy is evaluated [degR]
        fu1T  = 0;              %% optional input of primary fuel temperature which enthalpy is evaluated [degR]
        fu2T  = 0;              %% optional input of secondary fuel temperature which enthalpy is evaluated [degR]
        ox1H  = 0;              %% optional input of primary oxid enthalpy of formation [cal/mol]
        ox2H  = 0;              %% optional input of secondary oxid enthalpy of formation [cal/mol]          
        fu1H  = 0;              %% optional input of primary fuel enthalpy of formation [cal/kg]
        fu2H  = 0;              %% optional input of secondary fuel enthalpy of formation [cal/mol]
        
       
        Pc    = [];             %% Chamber pressure [psia]
        OF    = [];             %% mixture ratio [wt oxid/wt fuel]
        Phi   = [];             %% equivalence ratio
        
        Pe    = []              %% optional input for exit pressure
        PR    = [];             %% pressure ratios [Pc/Pe]
        subar = [];             %% subsonic area ratios [A/At]
        supar = [];             %% supersonic area ratios [A/At]



        CR    = 0;              %% chamber contraction ratio [Ac/At]
        flow  = 'eq';           %% flow type [eq or fz]
        out   = 'isp ivac cp p t gam m'; %%'aeat p t';  %maximum eight output parameters in one call to CEA
        %outlab= ['Area ratio' 'Temperature' 'Pressure' 'xH2O' 'xCO2' 'xCO'];
        
   end
   %%
   methods
       function runCEA(obj)
         %finds and runs CEA given the properties of the CEA object defined
         
         
         cd './CEA Code'                %find the CEA code directory (should be one level down from this class def)
         
         %call input generation function
         inpgen_rocket_3(obj.ox1,obj.ox2,obj.ox1wt,obj.ox2wt,obj.ox1T,...
             obj.ox2T,obj.ox1chem,obj.ox2chem,obj.ox1H,obj.ox2H,obj.fu1,...
             obj.fu2,obj.fu1wt,obj.fu2wt,obj.fu1T,obj.fu2T,obj.fu1chem,...
             obj.fu2chem,obj.fu1H,obj.fu2H,obj.Pc,obj.OF,obj.Phi,obj.PR,...
             obj.subar,obj.supar,obj.CR,obj.flow,obj.out);
         
         % execute CEA600
         system('CEA600.exe');
         
         cd '..'                        %return to the parent directory
       end
     
      function CEA_out = getCEAresults(obj, dataLoc, units)
         %location = 1: chamber, location = 2: throat, location = 3: exit
         %units: 'si' or 'en' for SI or English units returned
         %
         % this function returns a structure with all the relevant output
         %  data from CEA
         
         % INPUTS:
        %   dataLoc: 1,2,3: chamber, throat, exit
        %   units: 'si' or 'en' ONLY!
        %
        % OUTPUTS:
        %   CEA_out: structure with Gamma, P0, T0, rho, Cp, mu, k, Pr, a
        %

        %% Setup output

        % output file
        directory = '.\CEA Code';
        FID = fopen(fullfile(directory, 'Detn.out'));
        CEA_output = fscanf(FID,'%c');
        fclose(FID);

        %---------- Extract Output from File -------------%

        %% Stagnation Properties

            % find Mach
            Mach_index = strfind(CEA_output, 'MACH NUMBER') + 15; 
            end_ofline_index = strfind(CEA_output, 'TRANSPORT PROPERTIES (GASES ONLY)') - 1;           %finds endline
            line_of_Mach = CEA_output(Mach_index:end_ofline_index);              %makes array with only the mach values
            Mach_array = sscanf(line_of_Mach, '%*c %f');                           %gets the mach values
            
            
            % find P0
            P0_index = strfind(CEA_output, 'P, BAR') + 15; 
            end_ofline_index = strfind(CEA_output, 'T, K') - 1;           %finds endline
            line_of_P = CEA_output(P0_index:end_ofline_index);              %makes array with only the GAMMA values
            P_array = sscanf(line_of_P, '%*c %f');                           %gets the gamma values for eq flow
            P0 = P_array(1);


            % find T0
            T0_index = strfind(CEA_output, 'T, K') + 15; 
            end_ofline_index = strfind(CEA_output, 'RHO,') - 1;           %finds endline
            line_of_T = CEA_output(T0_index:end_ofline_index);              %makes array with only the GAMMA values
            T_array = sscanf(line_of_T, '%*c %f');                           %gets the gamma values for eq flow
            T0 = T_array(1);

            % find C*
            C_index = strfind(CEA_output, 'CSTAR') + 15;
            C_index = C_index(1);
            end_ofline_index = strfind(CEA_output, 'CF') - 1;           %finds endline
            line_of_C = CEA_output(C_index:end_ofline_index);              %makes array with only the GAMMA values
            C_array = sscanf(line_of_C, '%*c %f');                           %gets the gamma values for eq flow
            C_star = C_array(1);


            % find specific impulse in VACUUM
            Isp_index = strfind(CEA_output, 'Ivac') + 15;
            end_ofline_index = strfind(CEA_output, 'Isp') - 1;           %finds endline
            line_of_Isp = CEA_output(Isp_index:end_ofline_index);              %makes array with only the GAMMA values
            Isp_array = sscanf(line_of_Isp, '%*c %f');                           %gets the gamma values for eq flow



            % find specific impulse
            Isp_index = strfind(CEA_output, 'Isp') + 15;
            end_ofline_index = strfind(CEA_output, 'MASS FRACTIONS') - 3;           %finds endline
            line_of_Isp = CEA_output(Isp_index:end_ofline_index);              %makes array with only the impulse values
            Isp_atm_array = sscanf(line_of_Isp, '%*c %f');                           %gets the impulse values for eq flow



            % find Cf
            Cf_index = strfind(CEA_output, 'CF') + 15; 
            end_ofline_index = strfind(CEA_output, 'Ivac') - 1;           %finds endline
            line_of_Cf = CEA_output(Cf_index:end_ofline_index);              %makes array with only the GAMMA values
            Cf_array = sscanf(line_of_Cf, '%*c %f');                           %gets the gamma values for eq flow


            % find Ae/At
            eps_index = strfind(CEA_output, 'Ae/At') + 15; 
            end_ofline_index = strfind(CEA_output, 'CSTAR') - 1;           %finds endline
            line_of_eps = CEA_output(eps_index:end_ofline_index);              %makes array with only the GAMMA values
            eps_array = sscanf(line_of_eps, '%*c %f');                           %gets the gamma values for eq flow


            %If the user requests chamber conditions, all the performance values
            %   are returned at the throat
            if dataLoc == 1 
                Isp = Isp_atm_array(1);             
                Isp_vac = Isp_array(1);
                Cf = Cf_array(1);
                eps = 1;
            else
                Isp = Isp_atm_array(2);
                Isp_vac = Isp_array(2);
                Cf = Cf_array(2);
                if length(eps_array) <2     %in case user asks for exit and the exit doesnt exist
                    eps = eps_array(1);
                else 
                    eps = eps_array(2);
                end
            end


        %% Find the non_stagnation values
            % find GAMMA
            gamma_index = strfind(CEA_output, 'GAMMAs') + 15; 
            end_ofline_index = strfind(CEA_output, 'SON VEL,M/SEC') - 1;           %finds endline
            line_of_gamma = CEA_output(gamma_index:end_ofline_index);              %makes array with only the GAMMA values
            gamma_array = sscanf(line_of_gamma, '%*c %f');                           %gets the gamma values for eq flow



            % find rho
            rho_index = strfind(CEA_output, 'RHO, KG/CU M') + 15; 
            end_ofline_index = strfind(CEA_output, 'H, KJ/KG') - 1;           %finds endline
            line_of_rho = CEA_output(rho_index:end_ofline_index);              %makes array with only the GAMMA values
            rho_array = sscanf(line_of_rho, '%f');                           %gets the gamma values for eq flow
            
            if length(gamma_array) > 1 
                rho_array(1) = rho_array(1) * 10^(rho_array(2));
                rho_array(2) = rho_array(3) * 10^(rho_array(4));
                
                if length(gamma_array) > 2
                    rho_array(3) = rho_array(5) * 10^(rho_array(6));
                end
            end


            % find Cp
            cp_index = strfind(CEA_output, 'Cp,') + 15;
            cp_index = cp_index(1);
            end_ofline_index = strfind(CEA_output, 'GAMMAs') - 1;           %finds endline
            line_of_cp = CEA_output(cp_index:end_ofline_index);              %makes array with only the GAMMA values
            cp_array = sscanf(line_of_cp, '%*c %f');                           %gets the gamma values for eq flow



            % find Viscosity
            mu_index = strfind(CEA_output, 'VISC,') + 15;
            mu_index = mu_index(1);
            end_ofline_index = strfind(CEA_output, ' WITH EQUILIBRIUM') - 1;           %finds endline
            line_of_mu = CEA_output(mu_index:end_ofline_index);              %makes array with only the GAMMA values
            mu_array = sscanf(line_of_mu, '%*c %f');                           %gets the gamma values for eq flow



            % find conductivity
            k_index = strfind(CEA_output, 'CONDUCTIVITY') + 15;
            k_index = k_index(2);
            end_ofline_index = strfind(CEA_output, 'PRANDTL') - 1;           %finds endline
            end_ofline_index = end_ofline_index(1);
            line_of_k = CEA_output(k_index:end_ofline_index);              %makes array with only the GAMMA values
            k_array = sscanf(line_of_k, '%*c %f');                           %gets the gamma values for eq flow



            % find Prandtl Number
            Pr_index = strfind(CEA_output, 'PRANDTL') + 15;
            Pr_index = Pr_index(1);
            end_ofline_index = strfind(CEA_output, ' WITH FROZEN') - 1;           %finds endline
            line_of_Pr = CEA_output(Pr_index:end_ofline_index);              %makes array with only the GAMMA values
            Pr_array = sscanf(line_of_Pr, '%*c %f');                           %gets the gamma values for eq flow



            % find sonic velocity
            a_index = strfind(CEA_output, 'SON VEL') + 15;
            end_ofline_index = strfind(CEA_output, 'MACH NUMBER') - 1;          %finds endline
            line_of_a = CEA_output(a_index:end_ofline_index);              %makes array with only the GAMMA values
            a_array = sscanf(line_of_a, '%*c %f');                           %gets the gamma values for eq flow



            %find P
            P_index = strfind(CEA_output, 'P, BAR') + 15; 
            end_ofline_index = strfind(CEA_output, 'T, K') - 1;           %finds endline
            line_of_P = CEA_output(P_index:end_ofline_index);              %makes array with only the GAMMA values
            P_array = sscanf(line_of_P, '%*c %f');                           %gets the gamma values for eq flow



            % find T
            T_index = strfind(CEA_output, 'T, K') + 15; 
            end_ofline_index = strfind(CEA_output, 'RHO,') - 1;           %finds endline
            line_of_T = CEA_output(T_index:end_ofline_index);              %makes array with only the GAMMA values
            T_array = sscanf(line_of_T, '%*c %f');                           %gets the gamma values for eq flow


            %
            % still need to grab mole/mass fractions
            if isempty(strfind(CEA_output, 'MASS FRACTIONS')) ~= 1
    
               massf_index = strfind(CEA_output, 'MASS FRACTION') + 19;
               end_ofline_index = strfind(CEA_output, '* THERMODYNAMIC') - 6;
               line_of_massf = CEA_output(massf_index:end_ofline_index);
%                
               
               split = strsplit(line_of_massf, '\n');       %get number of compounds by splitting at endline
               expression = '\w*';                          %arugment for regexp (find word or numeric)    
               value_array = zeros(length(split));
               names = zeros(1,length(split));
            
               massf = struct();
               for i = 1:length(split)
                 matchStr = regexp(split(i), expression, 'match');  %find all the useful data in each line of mass/mole fractions (string and numbers)
                 compounds = matchStr{1,1};         %find the compound names
                 comp_array(i) = compounds(1,1);    %save the compound names
                 
                 
              
                 for j  = 1:length(compounds(1,2:end))
                    value_array(i,j) = str2double(compounds(1,j+1))/1e5;    %grab the mass/mole fraction values, regexp destroys decimal so this is how to get it back (CEA displays in .xxxxx)
                 end
                 
                
                 massf.(comp_array{i}) = value_array(i, dataLoc);        %dynamically name struct fields as the compound
               end

               molef = [];
            else
               molef_index = strfind(CEA_output, 'MOLE FRACTION') + 19;
               end_ofline_index = strfind(CEA_output, ' * THERMODYNAMIC') - 6;
               line_of_molef = CEA_output(molef_index:end_ofline_index);
                
               split = strsplit(line_of_molef, '\n');       %get number of compounds by splitting at endline
               expression = '\w*';                          %arugment for regexp (find word or numeric)    
               value_array = zeros(length(split));
               
            
               molef = struct();
               for i = 1:length(split)
                 matchStr = regexp(split(i), expression, 'match');  %find all the useful data in each line of mass/mole fractions (string and numbers)
                 compounds = matchStr{1,1};         %find the compound names
                 comp_array(i) = compounds(1,1);    %save the compound names
  
              
                 for j  = 1:length(compounds(1,2:end))
                    value_array(i,j) = str2double(compounds(1,j+1))/1e5;    %grab the mass/mole fraction values, regexp destroys decimal so this is how to get it back (CEA displays in .xxxxx)
                 end
                 
                
                 molef.(comp_array{i}) = value_array(i, dataLoc);        %dynamically name struct fields as the compound
               end
               massf= [];
            end



        %% Convert Units
            if strcmp(units,'en') == 1

                %convert Cp, mu, k, a, rho
                cp_array = cp_array * 1000 * 2.388e-4;              %J/Kg K to BTU/lbm*F
                mu_array = mu_array * 0.000067197 / 12;              %mPoise to lb/in*s
                k_array = k_array * 100/1000 * 0.001927  / 144;     %mW/cmK to BTU/in*s*F
                a_array = a_array * 3.28084 * 12;                     %m/s to in/s
                rho_array = rho_array * 3.61273e-5;                 %kg/m3 to lbm/in^3
                T_array = T_array * 1.8;                            %K to R
                P_array = P_array * 14.5038;                        %bar to psi


                %convert units
                P0 = P0 * 14.5038;                      %bar to psi
                T0 = T0 * 1.8;                          %K to R
                C_star = C_star * 3.28084;              %m/s to ft/s
                Isp_vac = Isp_vac * 3.28084;            %m/s to ft/s
                Isp = Isp * 3.28084;                    %m/s to ft/s

            end

        %% Output struct

            CEA_out = struct('Mach', Mach_array(dataLoc),'P0', P0, 'T0', T0, ...
                'Cstar', C_star, 'Ivac', Isp_vac, ... 
                'Isp', Isp, 'CF', Cf, 'GAMMA', gamma_array(dataLoc), 'Cp',...
                cp_array(dataLoc), 'mu', mu_array(dataLoc), 'k', ...
                k_array(dataLoc), 'Pr', Pr_array(dataLoc),'a', a_array(dataLoc),...
                'rho', rho_array(dataLoc), 'T', T_array(dataLoc), 'P', ...
                P_array(dataLoc), 'eps', eps, 'massf', massf, ...
                'molef', molef);
                
                CEA_out.compounds = comp_array; %note this is a cell array



       end
   end
end