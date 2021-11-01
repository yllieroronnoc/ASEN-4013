function [Ab] = burn_geometry(ri,r0,h,rb)
    
    if rb >= r0 % motor is burnt out
        Ab = 0; % [m^2] 
    else % there is grain remaining
        %% BURN AREA
        Ab = 2*pi()*(r0^2-(ri+rb)^2) + 2*pi()*(h-2*rb)*(ri+rb); % [m^2] total burn area

        %% BURN VOLUME
        % Vb = ; % [m^3]
    end