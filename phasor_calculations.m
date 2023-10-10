function [phasor_coords, percentages] = phasor_calculations(pure_spectra, exp_spectrum, lambda_min, lambda_max)
    %phasor_calculations() calculates and returns the phasor coordinates for two or three pure component 
        %spectra and one experimental spectrum. The percent of the experimental
        %spectrum composed by each pure spectrum is also calculated and returned.
    
    %Input Variables:
        %pure_spectra: nx4 (2 pure) or nx6 (3 pure) matrix where odd 
            %columns are the wavelengths and even columns are the intensity
        %exp_spctrum: mx2 matrix where column 1 is the wavelength and column 2 is the
            %intensity
        %lambda_min: minimum wavelength for calculations
        %lambda_max: maximum wavelength for calculations
            % All phasors are calculated using the same lambda_min and lambda_max
    
    %Returned Variables:
        %phasor_coords: 4x2 (3 pure) or 3x2 (2 pure) matrix containing
            %the phasor coordinates for all pure and experimental spectra.
            %Column 1 contains G values and column 2 contains S values. The
            %phasors are returned in the order [pure_1; pure_2; pure_3(optional); exp] 
            %where the order of the pure elements is the same
            %as input in pure_spectra
        %percentages: 3x1 (3 pure) or 2x1 (2 pure) matrix representing the
            %percent composition of each pure element in the experimental
            %spectrum as determined by phasor analysis in the same as order as
            %input in pure_spectra

    size_ = size(pure_spectra);
    phasor_coords = 0;
    percentages = 0;
    
    %3 pure components
    if size_(2) == 6
        G_1 = 0; S_1 = 0; I_1 = 0;
        G_2 = 0; S_2 = 0; I_2 = 0;
        G_3 = 0; S_3 = 0; I_3 = 0;
        G_exp = 0; S_exp = 0; I_exp = 0;
        %summation of spectra for G and S transforms
        for n = 1:size_(1)
            if pure_spectra(n,1) >= lambda_min && pure_spectra(n,1) <= lambda_max
                G_1 = G_1+ pure_spectra(n,2)*cos(2*pi*(pure_spectra(n,1) - lambda_min)/(lambda_max - lambda_min)); 
                S_1 = S_1+ pure_spectra(n,2)*sin(2*pi*(pure_spectra(n,1) - lambda_min)/(lambda_max - lambda_min));
                I_1 = I_1+ pure_spectra(n,2);
            end
            if pure_spectra(n,3) >= lambda_min && pure_spectra(n,3) <= lambda_max
                G_2 = G_2+ pure_spectra(n,4)*cos(2*pi*(pure_spectra(n,3) - lambda_min)/(lambda_max - lambda_min)); 
                S_2 = S_2+ pure_spectra(n,4)*sin(2*pi*(pure_spectra(n,3) - lambda_min)/(lambda_max - lambda_min));
                I_2 = I_2+ pure_spectra(n,4);
            end
            if pure_spectra(n,5) >= lambda_min && pure_spectra(n,5) <= lambda_max
                G_3 = G_3+ pure_spectra(n,6)*cos(2*pi*(pure_spectra(n,5) - lambda_min)/(lambda_max - lambda_min)); 
                S_3 = S_3+ pure_spectra(n,6)*sin(2*pi*(pure_spectra(n,5) - lambda_min)/(lambda_max - lambda_min));
                I_3 = I_3+ pure_spectra(n,6);
            end
        end
        for n = 1:length(exp_spectrum)
            if exp_spectrum(n,1) >= lambda_min && exp_spectrum(n,1) <= lambda_max
                G_exp = G_exp + exp_spectrum(n,2)*cos(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min)); 
                S_exp = S_exp + exp_spectrum(n,2)*sin(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min));
                I_exp = I_exp + exp_spectrum(n,2);
            end
        end
        
        %normalization of G and S coords to summed intensity
        G_1 = G_1/I_1;  S_1 = S_1/I_1;
        G_2 = G_2/I_2;  S_2 = S_2/I_2;
        G_3 = G_3/I_3;  S_3 = S_3/I_3;
        G_exp = G_exp/I_exp;  S_exp = S_exp/I_exp;
        
        
        %solving the linear phasor equation
        A = [G_1, G_2, G_3;
            S_1, S_2, S_3;
            1, 1, 1];
        B = [G_exp; S_exp; 1];
        percentages = lsqnonneg(A,B);
        phasor_coords = [G_1, S_1;
                         G_2, S_2;
                         G_3, S_3;
                         G_exp, S_exp];
    else
        %2 pure components
        if size_(2) == 4
            G_1 = 0; S_1 = 0; I_1 = 0;
            G_2 = 0; S_2 = 0; I_2 = 0;
            G_exp = 0; S_exp = 0; I_exp = 0;
            %summation of spectra for G and S transforms
            for n = 1:size_(1)
                if pure_spectra(n,1) >= lambda_min && pure_spectra(n,1) <= lambda_max
                    G_1 = G_1+ pure_spectra(n,2)*cos(2*pi*(pure_spectra(n,1) - lambda_min)/(lambda_max - lambda_min)); 
                    S_1 = S_1+ pure_spectra(n,2)*sin(2*pi*(pure_spectra(n,1) - lambda_min)/(lambda_max - lambda_min));
                    I_1 = I_1+ pure_spectra(n,2);
                end
                if pure_spectra(n,3) >= lambda_min && pure_spectra(n,3) <= lambda_max
                    G_2 = G_2+ pure_spectra(n,4)*cos(2*pi*(pure_spectra(n,3) - lambda_min)/(lambda_max - lambda_min)); 
                    S_2 = S_2+ pure_spectra(n,4)*sin(2*pi*(pure_spectra(n,3) - lambda_min)/(lambda_max - lambda_min));
                    I_2 = I_2+ pure_spectra(n,4);
                end
            end
            for n = 1:length(exp_spectrum)
                if exp_spectrum(n,1) >= lambda_min && exp_spectrum(n,1) <= lambda_max
                    G_exp = G_exp + exp_spectrum(n,2)*cos(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min)); 
                    S_exp = S_exp + exp_spectrum(n,2)*sin(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min));
                    I_exp = I_exp + exp_spectrum(n,2);
                end
            end

            %normalization of G and S coords to summed intensity
            G_1 = G_1/I_1;  S_1 = S_1/I_1;
            G_2 = G_2/I_2;  S_2 = S_2/I_2;
            G_exp = G_exp/I_exp;  S_exp = S_exp/I_exp;

            %solving the linear phasor equation
            A = [G_1, G_2;
                S_1, S_2];
            B = [G_exp; S_exp];
            percentages = lsqnonneg(A,B); 
            phasor_coords = [G_1, S_1;
                             G_2, S_2;
                             G_exp, S_exp];
        else
            if size_(1) == 2
                G_exp = 0; S_exp = 0; I_exp = 0;
                for n = 1:length(exp_spectrum)
                    if exp_spectrum(n,1) >= lambda_min && exp_spectrum(n,1) <= lambda_max
                        G_exp = G_exp + exp_spectrum(n,2)*cos(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min)); 
                        S_exp = S_exp + exp_spectrum(n,2)*sin(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min));
                        I_exp = I_exp + exp_spectrum(n,2);
                    end
                end
                G_exp = G_exp/I_exp;  S_exp = S_exp/I_exp;
                
                %solving the linear phasor equation
                A = [pure_spectra(1,1), pure_spectra(2,1);
                    pure_spectra(1,2), pure_spectra(2,2)];
                B = [G_exp; S_exp];
                percentages = lsqnonneg(A,B); 
                phasor_coords = [pure_spectra(1,1), pure_spectra(1,2);
                                 pure_spectra(2,1), pure_spectra(2,2);
                                 G_exp, S_exp];
            else
                if size_(1) == 3
                    G_exp = 0; S_exp = 0; I_exp = 0;
                    for n = 1:length(exp_spectrum)
                        if exp_spectrum(n,1) >= lambda_min && exp_spectrum(n,1) <= lambda_max
                            G_exp = G_exp + exp_spectrum(n,2)*cos(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min)); 
                            S_exp = S_exp + exp_spectrum(n,2)*sin(2*pi*(exp_spectrum(n,1) - lambda_min)/(lambda_max - lambda_min));
                            I_exp = I_exp + exp_spectrum(n,2);
                        end
                    end
                    G_exp = G_exp/I_exp;  S_exp = S_exp/I_exp;
                    
                    %solving the linear phasor equation
                    A = [pure_spectra(1,1), pure_spectra(2,1), pure_spectra(3,1)
                        pure_spectra(1,2), pure_spectra(2,2), pure_spectra(3,2)];
                    B = [G_exp; S_exp];
                    percentages = lsqnonneg(A,B); 
                    phasor_coords = [pure_spectra(1,1), pure_spectra(1,2);
                                     pure_spectra(2,1), pure_spectra(2,2);
                                     pure_spectra(3,1), pure_spectra(3,2);
                                     G_exp, S_exp];
                end
            end
    end
end

