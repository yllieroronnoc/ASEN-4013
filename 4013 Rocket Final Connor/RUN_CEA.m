function[SM_results1,SM_results3] = RUN_CEA(Pc_en, AR_sup)
    SM_inputs = CEAinput();           %cea input class object for easy input and definition
    
    % set conditions for the CEA run of H202 
    SM_inputs.ox1      = 'NH4CLO4(I)';                 %primary oxidizer
    SM_inputs.ox1T     = 536;                    %primary ox temp (R)

    SM_inputs.ox2      = 'AL(cr)';                 %primary oxidizer
    SM_inputs.ox2T     = 536;                    %primary ox temp (R)

    SM_inputs.fu1      = 'C4H6,butadiene';                 %primary fuel
    SM_inputs.fu1T     = 536;                    %primary fuel temp (R)
    
    SM_inputs.ox1wt    = 72;                      %primary oxidizer weight (by mass) in total 
    SM_inputs.ox2wt    = 10;                        % secondary oxidizer weight (by mass) in total
    SM_inputs.fu1wt    = 18;                        % primary fuel weight (by mass) in total
    
    SM_inputs.Pc       = Pc_en;                      %chamber pressure (psi)
    
    SM_inputs.supar    = AR_sup;                    % nozzle expansion ratio
    
    %% Run CEA

    SM_inputs.runCEA();                             %execute CEA for the above conditions
    SM_results1 = SM_inputs.getCEAresults(1,'si');
    SM_results3 = SM_inputs.getCEAresults(3,'si');           %1 for chamber conditions (2 for throat, 3 for nozzle exit: requires arg for area ratio or exit pressure), 'en' for english units
    