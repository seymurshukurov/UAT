classdef UAT_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        CorrugatedAntennaTab            matlab.ui.container.Tab
        VariablepwslotMCDeltaEditField  matlab.ui.control.NumericEditField
        VariablepwslotMCDeltaEditFieldLabel  matlab.ui.control.Label
        RingLoadedMCDeltaEditField      matlab.ui.control.NumericEditField
        RingLoadedMCDeltaEditFieldLabel  matlab.ui.control.Label
        NumberofModeConvertersEditField  matlab.ui.control.NumericEditField
        NumberofModeConvertersEditFieldLabel  matlab.ui.control.Label
        NumberofCorrugationsEditField   matlab.ui.control.NumericEditField
        NumberofCorrugationsEditFieldLabel  matlab.ui.control.Label
        InputWaveguideLengthmmEditField  matlab.ui.control.NumericEditField
        InputWaveguideLengthmmEditFieldLabel  matlab.ui.control.Label
        SigmaEditField                  matlab.ui.control.NumericEditField
        SigmaEditFieldLabel             matlab.ui.control.Label
        PitchtoWidthpWRatioEditField    matlab.ui.control.NumericEditField
        PitchtoWidthpWRatioEditFieldLabel  matlab.ui.control.Label
        PitchEditField                  matlab.ui.control.NumericEditField
        PitchEditFieldLabel             matlab.ui.control.Label
        OutputRadiusCoefficientEditField  matlab.ui.control.NumericEditField
        OutputRadiusCoefficientEditFieldLabel  matlab.ui.control.Label
        MaximumFrequencyGHzEditField    matlab.ui.control.NumericEditField
        MaximumFrequencyGHzEditFieldLabel  matlab.ui.control.Label
        MinimumFrequencyGHzEditField    matlab.ui.control.NumericEditField
        MinimumFrequencyGHzEditFieldLabel  matlab.ui.control.Label
        ModeConverterTypeDropDown       matlab.ui.control.DropDown
        ModeConverterTypeDropDownLabel  matlab.ui.control.Label
        HornProfileDropDown             matlab.ui.control.DropDown
        HornProfileDropDownLabel        matlab.ui.control.Label
        ConstructButton                 matlab.ui.control.Button
        PreviewButton                   matlab.ui.control.Button
        UIAxes                          matlab.ui.control.UIAxes
        FSSTab                          matlab.ui.container.Tab
        ElementSpacingmmEditField       matlab.ui.control.NumericEditField
        ElementSpacingmmEditFieldLabel  matlab.ui.control.Label
        InnerRadiusLengthmmEditField    matlab.ui.control.NumericEditField
        InnerRadiusLengthmmEditFieldLabel  matlab.ui.control.Label
        OuterRadiusLengthmmEditField    matlab.ui.control.NumericEditField
        OuterRadiusLengthmmEditFieldLabel  matlab.ui.control.Label
        ConstructFSSButton              matlab.ui.control.Button
        PreviewFSSButton                matlab.ui.control.Button
        LatticeStructureDropDown        matlab.ui.control.DropDown
        LatticeStructureDropDownLabel   matlab.ui.control.Label
        UnitCellTypeDropDown            matlab.ui.control.DropDown
        UnitCellTypeDropDownLabel       matlab.ui.control.Label
        UIAxes_FSS                      matlab.ui.control.UIAxes
        ReflectorAntennaTab             matlab.ui.container.Tab
        SeptumPolarizerTab              matlab.ui.container.Tab
        WaveguideFilterTab              matlab.ui.container.Tab
        UltimateRFDesignerToolLabel     matlab.ui.control.Label
    end


    properties (Access = public)
        %         fmin;       %minimum frequency
        %         fmax;       %maximum frequency
        %         ao;         %Output Radius Coefficient
        %         p;          %pitch
        %         delta;      %pitch to width ratio
        %         sigma;      %sigma
        %         wgl;        %waveguide length
        %         N;          %number of corrugations
        %         n_mc;       %Number of mode converters
        %         delta2;     %delta2 is only used for case 2, ring loaded slot mode converter
        %         delta_min;  % delta_min is only used for case 3,variable pitch to width slot mode converter
        %         horn_pro;    %horn antenna profile
        aaa;
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: PreviewButton
        function PreviewButtonPushed(app, event)


            fmin = app.MinimumFrequencyGHzEditField.Value; % in ghz
            fmax = app.MaximumFrequencyGHzEditField.Value; % in ghz

            scl=fmax/fmin;

            if scl < 1.4
                fc=sqrt(fmin*fmax);
                fo=1.02*fc;
            else
                fc=1.2*fmin;
                fo=1.1*fc;
            end

            lam_c=physconst('LightSpeed')/(fc*10^9);
            lam_c_mm = lam_c * 1000;
            lam_o=physconst('LightSpeed')/(fo*10^9);
            lam_o_mm = lam_o * 1000;

            ai=3*lam_c_mm/(2*pi); %in mm
            %from fig 3
            ao=app.OutputRadiusCoefficientEditField.Value*lam_c_mm; %in mm
            %pitch between 5 (narrowband) - 10 (broadband)
            p=lam_c_mm/app.PitchEditField.Value;
            % pitch to width ratio (affect XPD)
            delta = app.PitchtoWidthpWRatioEditField.Value;
            sigma = app.SigmaEditField.Value; % sigma 0.4 - 0.5
            kc = (2*pi)/lam_c_mm;            % Wave number at center frequency
            ko = (2*pi)/lam_o_mm;            % Wave number at output frequency
            wgl = app.InputWaveguideLengthmmEditField.Value; % input waveguide length

            %Number of total corrugations
            N = app.NumberofCorrugationsEditField.Value;
            len = N*p;
            %Number of mode converters
            n_mc = app.NumberofModeConvertersEditField.Value;
            z=0:p:len;

            % delta2 is only used for case 2, ring loaded slot mode converter
            delta2 = app.RingLoadedMCDeltaEditField.Value;

            % delta_min is only used for case 3,variable pitch to width slot mode converter
            delta_min = app.VariablepwslotMCDeltaEditField.Value; % Should be greater than 0.125 and less than delta

            %1=LINEAR, 2=SINUSOID, 3=ASYMMETRIC SINE-SQUARED, 4=TANGENTIAL,
            % 5=x.rho, 6=EXPONENTIAL, 7=HYPERBOLIC, 8=POLYNOMIAL

            popupmenu1value = app.HornProfileDropDown.Value;

            switch popupmenu1value
                case 'Linear'
                    horn_profile = 1;
                case 'Sinusoidal'
                    horn_profile = 2;
                case 'Asymmetric Sine-Squared'
                    horn_profile = 3;
                case 'Tangential'
                    horn_profile = 4;
                case 'x.rho'
                    horn_profile = 5;
                case 'Exponential'
                    horn_profile = 6;
                case 'Hyperbolic'
                    horn_profile = 7;
                case 'Poynomial'
                    horn_profile = 8;
            end

            %1=VARIABLE SLOT DEPTH MC,
            % 2=RING LOADED SLOTS MC, 3=VARIABLE PITCH TO WIDTH SLOT MC
            popupmenu2value = app.ModeConverterTypeDropDown.Value;
            switch popupmenu2value
                case 'Variable Slot Depth'
                    mode_converter_type = 1;
                case 'Ring Loaded Slots'
                    mode_converter_type = 2;
                case 'Variable p/w Slot'
                    mode_converter_type = 3;
            end

            switch(horn_profile)
                case 1
                    %%% Linear profile %%%
                    a = ai+(ao-ai)*z/len;

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Linear Horn Profile', 'FontSize', 16 );

                case 2
                    %%% Sinusoid profile %%%
                    A = 1;     % Amplitude factor 'A' should be between 0 and 1
                    rho = 2;     % rho should be between 0.5 and 5, default is 2
                    a = ai+(ao-ai)*((1-A)*(z/len)+A*power(sin((pi*z)/(2*len)),rho));

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Sinusoidal Horn Profile', 'FontSize', 16 );

                case 3
                    %%% Asymmetric Sine-Squared profile %%%
                    L1 = len/3;    % Choose a value for L1, must be less than the horn length
                    L2 = len-L1;   % L2 is the length between L1 and the end of the horn
                    gamma = L2/L1;
                    idx = find(z <= L1);     % Find the index of z corresponding to L1
                    zelements = size(z,2);   % Total number of points in z axis of the horn

                    za = z(1: max(idx));
                    aa = ai+((2*(ao-ai))/(1+gamma))*sin((pi*za)/(4*L1)).^2;

                    zb = z(max(idx)+1 : zelements);
                    ab = ai+((2*(ao-ai))/(1+gamma))*(gamma*sin(((pi*(zb+L2-L1))/(4*L2))).^2+((1-gamma)/2));

                    a = [aa,ab];
                    z = [za,zb];

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Asymmetric Sine Squared Horn Profile', 'FontSize', 16 );

                case 4
                    %%% Tangential profile %%%
                    A = 1;
                    rho = 2;
                    a = ai+(ao-ai)*((1-A)*(z/len)+A*power(tan((pi*z)/(4*len)),rho));

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Tangential Horn Profile', 'FontSize', 16 );

                case 5
                    %%% x.rho profile %%%
                    A = 1;
                    rho = 2;
                    a = ai+(ao-ai)*((1-A)*(z/len)+A*power(z/len,rho));

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'xp Horn Profile', 'FontSize', 16 );

                case 6
                    %%% Exponential profile %%%
                    a=ai*exp(log(ao/ai)*(z/len));

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Exponential Horn Profile', 'FontSize', 16 );

                case 7
                    %%% Hyperbolic profile %%%
                    a = sqrt(ai^2 + (power(z,2) * (ao^2-ai^2) / len^2));

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Hyperbolic Horn Profile', 'FontSize', 16 );

                case 8
                    %%% POLYNOMIAL Profile %%%
                    rho = 3;
                    a=ai+(rho+1)*(ao-ai)*(1-((rho*z)/((rho+1)*len))).*power(z/len,rho);

                    plot(app.UIAxes, z, a);
                    set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                    xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
                    ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
                    title(app.UIAxes, 'Polynomial Horn Profile', 'FontSize', 16 );
            end

            app.aaa=a;
            switch(mode_converter_type)

                case 1                            % case 1=VARIABLE SLOT DEPTH MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                     % Index range for mode converter
                    idx = 1:n_mc;
                    djmc = (sigma-((idx-1)./n_mc).*(sigma-(0.25.*exp(1./(2.114.*(kc*ajmc).^1.134)))))*lam_c_mm;
                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as len and rad vector
                    n = 0;
                    lent(1) = 0;
                    lent(2) = 0;
                    for i = 1:N
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+delta*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
                        lent(i+n+5) = lent(i+n+4);
                        n = n+3;
                    end

                    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
                    lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length

                case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                 % Index range for mode converter
                    idx = 1:n_mc;
                    djmc = (lam_c_mm/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
                    % Width of bjth slot for mode converter
                    bj = (0.1+(idx-1).*((delta2-0.1)./n_mc)).*p;
                    % Height of hjth slot for mode converter
                    hj = (2/3).*djmc;

                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as lent,rad vector
                    n = 5;
                    lent = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
                    rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
                    for i = 2:n_mc
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i)-hj(i);
                        rad(i+n+2) = a(i)+d(i)-hj(i);
                        rad(i+n+3) = a(i)+d(i);
                        rad(i+n+4) = a(i)+d(i);
                        rad(i+n+5) = a(i+1);
                        lent(i+n) = i*p-p;
                        lent(i+n+1) = i*p-p;
                        lent(i+n+2) = lent(i+n+1)-(delta*p)+bj(i);
                        lent(i+n+3) = lent(i+n+1)-(delta*p)+bj(i);
                        lent(i+n+4) = lent(i+n+3)+(delta*p)+bj(i);
                        lent(i+n+5) = lent(i+n+3)+(delta*p)+bj(i);
                        n = n+5;
                    end
                    % Add extra coordinate points before remaining corrugations
                    lent(n_mc*(n_mc+1)+1) = lent(n_mc*(n_mc+1))+(1-delta)*p;
                    rad(n_mc*(n_mc+1)+1) = a(n_mc+1);

                    n = n+n_mc+1;
                    for i = n_mc+1:N
                        rad(n) = a(i);
                        rad(n+1) = a(i)+d(i);
                        rad(n+2) = a(i)+d(i);
                        rad(n+3) = a(i+1);
                        rad(n+4) = a(i+1);
                        lent(n+1) = lent(n);
                        lent(n+2) = lent(n+1)+delta*p;
                        lent(n+3) = lent(n+2);
                        lent(n+4) = lent(n+3)+(1-delta)*p;
                        n = n+4;
                    end

                    z_number = (n_mc*2)+(N*4)+1; % Number of coordinate points for corrugated length of horn


                case 3                 % case 3=VARIABLE PITCH TO WIDTH SLOT MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                     % First indexes for mode converter
                    idx = 1:n_mc;
                    djmc = (sigma*(lam_c_mm/1.15)+((idx-1)./(n_mc-1)).*(lam_c_mm/4-(sigma*lam_c_mm/1.15))).*exp(1./(2.114.*(kc*ajmc).^1.134));
                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as lent,rad vector
                    n = 0;
                    lent(1) = 0;
                    lent(2) = 0;
                    for i = 1:n_mc
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min)))*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min))))*p;
                        lent(i+n+5) = lent(i+n+4);
                        n = n+3;
                    end

                    for i = n_mc+1:N
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+delta*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
                        lent(i+n+5) = lent(i+n+4);
                        n = n+3;
                    end

                    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
                    lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length
            end

            % Add the rest of the geometry to create a closed path
            % a_offset is the inner horn profile shifted up to give the horn a thickness
            a_offset = a+(lam_c_mm/2+2);

            % Add vertical surface at horn aperture
            lent = [lent, lent(z_number)];
            rad = [rad, a_offset(N)];
            radmsh=rad;                 % radmesh to fix mesh lines to corrugations

            % Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
            outer_surface = fliplr(a_offset);
            z_flip = fliplr(z);
            extent = lent(end);  % Fudge to make horn aperture planar for ring loaded slot MC
            z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
            % Add outer profile and circular waveguide to horn
            lent = [lent, z_flip, -wgl, -wgl, 0];
            rad = [rad, outer_surface, ai+(lam_c_mm/2+2), ai, ai];

            plot(app.UIAxes, lent,rad);
            xlim(app.UIAxes, [-wgl-10 N*p+10])
            ylim(app.UIAxes, [0 Inf])
            set(app.UIAxes, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
            xlabel(app.UIAxes, 'Dimension in z Direction (mm)', 'FontSize', 14 );
            ylabel(app.UIAxes, 'Dimension in y Direction (mm)', 'FontSize', 14 );
            title(app.UIAxes, 'Complete Corrugated Horn Profile', 'FontSize', 16 );

        end

        % Button pushed function: ConstructButton
        function ConstructButtonPushed(app, event)
            CST = CST_MicrowaveStudio(cd,'dipole');
            a=app.aaa;

            fmin = app.MinimumFrequencyGHzEditField.Value; % in ghz
            fmax = app.MaximumFrequencyGHzEditField.Value; % in ghz
            CST.setFreq(fmin,fmax);
            CST.setBoundaryCondition('xmin','open add space','xmax','open add space','ymin','open add space',...
                'ymax','open add space','zmin','open add space','zmax','open add space')

            scl=fmax/fmin;

            if scl < 1.4
                fc=sqrt(fmin*fmax);
                fo=1.02*fc;
            else
                fc=1.2*fmin;
                fo=1.1*fc;
            end

            lam_c=physconst('LightSpeed')/(fc*10^9);
            lam_c_mm = lam_c * 1000;
            lam_o=physconst('LightSpeed')/(fo*10^9);
            lam_o_mm = lam_o * 1000;

            ai=3*lam_c_mm/(2*pi); %in mm
            %from fig 3
            ao=app.OutputRadiusCoefficientEditField.Value*lam_c_mm; %in mm
            %pitch between 5 (narrowband) - 10 (broadband)
            p=lam_c_mm/app.PitchEditField.Value;
            % pitch to width ratio (affect XPD)
            delta = app.PitchtoWidthpWRatioEditField.Value;
            sigma = app.SigmaEditField.Value; % sigma 0.4 - 0.5
            kc = (2*pi)/lam_c_mm;            % Wave number at center frequency
            ko = (2*pi)/lam_o_mm;            % Wave number at output frequency
            wgl = app.InputWaveguideLengthmmEditField.Value; % input waveguide length

            %Number of total corrugations
            N = app.NumberofCorrugationsEditField.Value;
            len = N*p;
            %Number of mode converters
            n_mc = app.NumberofModeConvertersEditField.Value;
            z=0:p:len;

            % delta2 is only used for case 2, ring loaded slot mode converter
            delta2 = app.RingLoadedMCDeltaEditField.Value;

            % delta_min is only used for case 3,variable pitch to width slot mode converter
            delta_min = app.VariablepwslotMCDeltaEditField.Value; % Should be greater than 0.125 and less than delta

            %1=LINEAR, 2=SINUSOID, 3=ASYMMETRIC SINE-SQUARED, 4=TANGENTIAL,
            % 5=x.rho, 6=EXPONENTIAL, 7=HYPERBOLIC, 8=POLYNOMIAL

            popupmenu1value = app.HornProfileDropDown.Value;

            switch popupmenu1value
                case 'Linear'
                    horn_profile = 1;
                case 'Sinusoidal'
                    horn_profile = 2;
                case 'Asymmetric Sine-Squared'
                    horn_profile = 3;
                case 'Tangential'
                    horn_profile = 4;
                case 'x.rho'
                    horn_profile = 5;
                case 'Exponential'
                    horn_profile = 6;
                case 'Hyperbolic'
                    horn_profile = 7;
                case 'Poynomial'
                    horn_profile = 8;
            end

            %1=VARIABLE SLOT DEPTH MC,
            % 2=RING LOADED SLOTS MC, 3=VARIABLE PITCH TO WIDTH SLOT MC
            popupmenu2value = app.ModeConverterTypeDropDown.Value;
            switch popupmenu2value
                case 'Variable Slot Depth'
                    mode_converter_type = 1;
                case 'Ring Loaded Slots'
                    mode_converter_type = 2;
                case 'Variable p/w Slot'
                    mode_converter_type = 3;
            end

            switch(horn_profile)
                case 1
                    %%% Linear profile %%%
                    a = ai+(ao-ai)*z/len;
                case 2
                    %%% Sinusoid profile %%%
                    A = 1;     % Amplitude factor 'A' should be between 0 and 1
                    rho = 2;     % rho should be between 0.5 and 5, default is 2
                    a = ai+(ao-ai)*((1-A)*(z/len)+A*power(sin((pi*z)/(2*len)),rho));
                case 3
                    %%% Asymmetric Sine-Squared profile %%%
                    L1 = len/3;    % Choose a value for L1, must be less than the horn length
                    L2 = len-L1;   % L2 is the length between L1 and the end of the horn
                    gamma = L2/L1;
                    idx = find(z <= L1);     % Find the index of z corresponding to L1
                    zelements = size(z,2);   % Total number of points in z axis of the horn

                    za = z(1: max(idx));
                    aa = ai+((2*(ao-ai))/(1+gamma))*sin((pi*za)/(4*L1)).^2;

                    zb = z(max(idx)+1 : zelements);
                    ab = ai+((2*(ao-ai))/(1+gamma))*(gamma*sin(((pi*(zb+L2-L1))/(4*L2))).^2+((1-gamma)/2));

                    a = [aa,ab];
                    z = [za,zb];
                case 4
                    %%% Tangential profile %%%
                    A = 1;
                    rho = 2;
                    a = ai+(ao-ai)*((1-A)*(z/len)+A*power(tan((pi*z)/(4*len)),rho));
                case 5
                    %%% x.rho profile %%%
                    A = 1;
                    rho = 2;
                    a = ai+(ao-ai)*((1-A)*(z/len)+A*power(z/len,rho));
                case 6
                    %%% Exponential profile %%%
                    a=ai*exp(log(ao/ai)*(z/len));
                case 7
                    %%% Hyperbolic profile %%%
                    a = sqrt(ai^2 + (power(z,2) * (ao^2-ai^2) / len^2));
                case 8
                    %%% POLYNOMIAL Profile %%%
                    rho = 3;
                    a=ai+(rho+1)*(ao-ai)*(1-((rho*z)/((rho+1)*len))).*power(z/len,rho);
            end

            switch(mode_converter_type)

                case 1                            % case 1=VARIABLE SLOT DEPTH MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                     % Index range for mode converter
                    idx = 1:n_mc;
                    djmc = (sigma-((idx-1)./n_mc).*(sigma-(0.25.*exp(1./(2.114.*(kc*ajmc).^1.134)))))*lam_c_mm;
                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as len and rad vector
                    n = 0;
                    lent(1) = 0;
                    lent(2) = 0;
                    for i = 1:N
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+delta*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
                        lent(i+n+5) = lent(i+n+4);
                        n = n+3;
                    end

                    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
                    lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length

                case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                 % Index range for mode converter
                    idx = 1:n_mc;
                    djmc = (lam_c_mm/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
                    % Width of bjth slot for mode converter
                    bj = (0.1+(idx-1).*((delta2-0.1)./n_mc)).*p;
                    % Height of hjth slot for mode converter
                    hj = (2/3).*djmc;

                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as lent,rad vector
                    n = 5;
                    lent = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
                    rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
                    for i = 2:n_mc
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i)-hj(i);
                        rad(i+n+2) = a(i)+d(i)-hj(i);
                        rad(i+n+3) = a(i)+d(i);
                        rad(i+n+4) = a(i)+d(i);
                        rad(i+n+5) = a(i+1);
                        lent(i+n) = i*p-p;
                        lent(i+n+1) = i*p-p;
                        lent(i+n+2) = lent(i+n+1)-(delta*p)+bj(i);
                        lent(i+n+3) = lent(i+n+1)-(delta*p)+bj(i);
                        lent(i+n+4) = lent(i+n+3)+(delta*p)+bj(i);
                        lent(i+n+5) = lent(i+n+3)+(delta*p)+bj(i);
                        n = n+5;
                    end
                    % Add extra coordinate points before remaining corrugations
                    lent(n_mc*(n_mc+1)+1) = lent(n_mc*(n_mc+1))+(1-delta)*p;
                    rad(n_mc*(n_mc+1)+1) = a(n_mc+1);

                    n = n+n_mc+1;
                    for i = n_mc+1:N
                        rad(n) = a(i);
                        rad(n+1) = a(i)+d(i);
                        rad(n+2) = a(i)+d(i);
                        rad(n+3) = a(i+1);
                        rad(n+4) = a(i+1);
                        lent(n+1) = lent(n);
                        lent(n+2) = lent(n+1)+delta*p;
                        lent(n+3) = lent(n+2);
                        lent(n+4) = lent(n+3)+(1-delta)*p;
                        n = n+4;
                    end

                    z_number = (n_mc*2)+(N*4)+1; % Number of coordinate points for corrugated length of horn

                case 3                 % case 3=VARIABLE PITCH TO WIDTH SLOT MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                     % First indexes for mode converter
                    idx = 1:n_mc;
                    djmc = (sigma*(lam_c_mm/1.15)+((idx-1)./(n_mc-1)).*(lam_c_mm/4-(sigma*lam_c_mm/1.15))).*exp(1./(2.114.*(kc*ajmc).^1.134));
                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as lent,rad vector
                    n = 0;
                    lent(1) = 0;
                    lent(2) = 0;
                    for i = 1:n_mc
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min)))*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min))))*p;
                        lent(i+n+5) = lent(i+n+4);
                        n = n+3;
                    end

                    for i = n_mc+1:N
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+delta*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
                        lent(i+n+5) = lent(i+n+4);
                        n = n+3;
                    end

                    z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
                    lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length

            end

            % Add the rest of the geometry to create a closed path
            % a_offset is the inner horn profile shifted up to give the horn a thickness
            a_offset = a+(lam_c_mm/2+2);
            %figure;                    % Uncomment these three lines for debugging
            %plot(z, a_offset);
            %axis equal;

            % Add vertical surface at horn aperture
            lent = [lent, lent(z_number)];
            rad = [rad, a_offset(N)];
            radmsh=rad;                 % radmesh to fix mesh lines to corrugations
            %figure;                    % Uncomment these three lines for debugging
            %plot(lent, rad);
            %axis equal;

            % Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
            outer_surface = fliplr(a_offset);
            z_flip = fliplr(z);
            extent = lent(end);  % Fudge to make horn aperture planar for ring loaded slot MC
            z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
            % Add outer profile and circular waveguide to horn
            lent = [lent, z_flip, -wgl, -wgl, 0];
            rad = [rad, outer_surface, ai+(lam_c_mm/2+2), ai, ai];

            material = 'Copper (annealed)';
            Name = 'Waveguide';
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = ai;
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [-wgl 0];
            CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

            switch(mode_converter_type)

                case 1

                    n = 0;
                    lent(1) = 0;
                    lent(2) = 0;
                    for i=1:N
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+delta*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
                        lent(i+n+5) = lent(i+n+4);

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'b'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+1);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n) lent(i+n+2)];
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'a'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+3);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n+3) lent(i+n+5)];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        n = n+3;
                    end

                case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
                    % Mode Converter depths for element j
                    ajmc = a(1:n_mc);                 % Index range for mode converter
                    idx = 1:n_mc;
                    djmc = (lam_c_mm/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
                    % Width of bjth slot for mode converter
                    bj = (0.1+(idx-1).*((delta2-0.1)./n_mc)).*p;
                    % Height of hjth slot for mode converter
                    hj = (2/3).*djmc;

                    % Depth of remaining corrugations
                    aj = a(n_mc+1:end);
                    idx = n_mc+1:N+1;
                    dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
                    d = [djmc, dj];       % Combining the mode converter and horn depth values

                    % Generate z,y coordinates as lent,rad vector
                    n = 5;
                    lent = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
                    rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
                    i=1;
                    material = 'Copper (annealed)';
                    Name = ['Coor' num2str(1) 'c'];
                    OuterRadius = rad(i+1);
                    InnerRadius = rad(i);
                    Xcenter = 0;
                    Ycenter = 0;
                    Zrange = [lent(i+1) lent(i+2)];
                    % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                    CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)


                    material = 'Copper (annealed)';
                    Name = ['Coor' num2str(1) 'd'];
                    OuterRadius = ao+(lam_c_mm/2+2);
                    InnerRadius = rad(i+3);
                    Xcenter = 0;
                    Ycenter = 0;
                    Zrange = [lent(i+3) lent(i+4)];
                    % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                    CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)


                    material = 'Copper (annealed)';
                    Name = ['Coor' num2str(1) 'e'];
                    OuterRadius = ao+(lam_c_mm/2+2);
                    InnerRadius = rad(i+5);
                    Xcenter = 0;
                    Ycenter = 0;
                    Zrange = [lent(i+4) lent(i+4)+p-(delta*p)];
                    % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                    CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)


                    material = 'Copper (annealed)';
                    Name = ['Coor' num2str(1) 'f'];
                    OuterRadius = rad(i+3);
                    InnerRadius = rad(i+1);
                    Xcenter = 0;
                    Ycenter = 0;
                    Zrange = [lent(i+1) lent(i+2)];
                    % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                    CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)


                    component1 = 'component1:Waveguide';
                    component2 = 'component1:Coor1f';
                    % CstSubtract(mws,component1,component2)
                    CST.subtractObject(component1,component2)

                    for i = 2:n_mc
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i)-hj(i);
                        rad(i+n+2) = a(i)+d(i)-hj(i);
                        rad(i+n+3) = a(i)+d(i);
                        rad(i+n+4) = a(i)+d(i);
                        rad(i+n+5) = a(i+1);
                        lent(i+n) = i*p-p;
                        lent(i+n+1) = i*p-p;
                        lent(i+n+2) = lent(i+n+1)-(delta*p)+bj(i);
                        lent(i+n+3) = lent(i+n+1)-(delta*p)+bj(i);
                        lent(i+n+4) = lent(i+n+3)+(delta*p)+bj(i);
                        lent(i+n+5) = lent(i+n+3)+(delta*p)+bj(i);

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'c'];
                        OuterRadius = rad(i+n+1);
                        InnerRadius = rad(i+n);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n+1) lent(i+n+2)];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)


                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'd'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+3);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n+3) lent(i+n+4)];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'e'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+5);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n+4) lent(i+n+4)+(1-delta)*p];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        n = n+5;
                    end
                    % Add extra coordinate points before remaining corrugations
                    lent(n_mc*(n_mc+1)+1) = lent(n_mc*(n_mc+1))+(1-delta)*p;
                    rad(n_mc*(n_mc+1)+1) = a(n_mc+1);

                    n = n+n_mc+1;
                    for i = n_mc+1:N
                        rad(n) = a(i);
                        rad(n+1) = a(i)+d(i);
                        rad(n+2) = a(i)+d(i);
                        rad(n+3) = a(i+1);
                        rad(n+4) = a(i+1);
                        lent(n+1) = lent(n);
                        lent(n+2) = lent(n+1)+delta*p;
                        lent(n+3) = lent(n+2);
                        lent(n+4) = lent(n+3)+(1-delta)*p;

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'b'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(n+1);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(n) lent(n+2)];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'a'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(n+3);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(n+2) lent(n+4)];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)
                        n = n+4;
                    end


                case 3
                    n = 0;
                    lent(1) = 0;
                    lent(2) = 0;
                    for i = 1:n_mc
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min)))*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min))))*p;
                        lent(i+n+5) = lent(i+n+4);

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'b'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+1);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n) lent(i+n+2)];
                        % Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'a'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+3);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n+3) lent(i+n+5)];
                        %                         Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)
                        n = n+3;
                    end

                    for i = n_mc+1:N
                        rad(i+n) = a(i);
                        rad(i+n+1) = a(i)+d(i);
                        rad(i+n+2) = a(i)+d(i);
                        rad(i+n+3) = a(i+1);
                        rad(i+n+4) = a(i+1);
                        lent(i+n+2) = lent(i+n)+delta*p;
                        lent(i+n+3) = lent(i+n+2);
                        lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
                        lent(i+n+5) = lent(i+n+4);

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'b'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+1);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n) lent(i+n+2)];
                        %                         Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)

                        material = 'Copper (annealed)';
                        Name = ['Coor' num2str(i) 'a'];
                        OuterRadius = ao+(lam_c_mm/2+2);
                        InnerRadius = rad(i+n+3);
                        Xcenter = 0;
                        Ycenter = 0;
                        Zrange = [lent(i+n+3) lent(i+n+5)];
                        %                         Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
                        CST.addCylinder(OuterRadius,InnerRadius,'z',Xcenter,Ycenter,Zrange,Name,'component1',material)
                        n = n+3;
                    end
            end
        end

        % Button pushed function: PreviewFSSButton
        function PreviewFSSButtonPushed(app, event)
            celltype = app.UnitCellTypeDropDown.Value;
            lattice = app.LatticeStructureDropDown.Value;

            switch celltype
                case 'Circular'
                    switch lattice
                        case 'Triangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;

                            th = 0:pi/50:2*pi;
                            xunit_out = r_out * cos(th);
                            yunit_out = r_out * sin(th);

                            plot(app.UIAxes_FSS, xunit_out, yunit_out,'r');
                            hold(app.UIAxes_FSS,'on')
                            set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                            fill(app.UIAxes_FSS, [0* cos(th) flip(xunit_out)],[0*sin(th) flip(yunit_out)],'b')

                            for ii=0:60:300
                                plot(app.UIAxes_FSS, xunit_out+dist*sind(ii), yunit_out+dist*cosd(ii),'r');
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                fill(app.UIAxes_FSS, [0* sin(th)+dist*sind(ii) flip(xunit_out+dist*sind(ii))],...
                                    [0*cos(th)+dist*cosd(ii) flip(yunit_out+dist*cosd(ii))],'b')
                            end
                        case 'Rectangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;

                            th = 0:pi/50:2*pi;
                            xunit_out = r_out * cos(th);
                            yunit_out = r_out * sin(th);

                            for ii=45:90:315
                                plot(app.UIAxes_FSS, xunit_out+dist*sqrt(2)*sind(ii)/2, yunit_out+dist*sqrt(2)*cosd(ii)/2,'b');
                                hold(app.UIAxes_FSS,'on')
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                fill(app.UIAxes_FSS, [0*cos(th)+dist*sqrt(2)*sind(ii)/2 flip(xunit_out+dist*sqrt(2)*sind(ii)/2)] ...
                                    ,[0*sin(th)+dist*sqrt(2)*cosd(ii)/2 flip(yunit_out+dist*sqrt(2)*cosd(ii)/2)],'b')
                            end
                    end
                case 'Circular Loop'
                    switch lattice
                        case 'Triangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            r_in=app.InnerRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;

                            th = 0:pi/50:2*pi;
                            xunit_out = r_out * cos(th);
                            yunit_out = r_out * sin(th);
                            plot(app.UIAxes_FSS, xunit_out, yunit_out);
                            hold(app.UIAxes_FSS,'on')
                            set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])

                            xunit_in = r_in * cos(th);
                            yunit_in = r_in * sin(th);
                            plot(app.UIAxes_FSS, xunit_in, yunit_in);
                            set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])

                            fill(app.UIAxes_FSS, [xunit_in flip(xunit_out)],[yunit_in flip(yunit_out)],'b')
                            for ii=0:60:300
                                plot(app.UIAxes_FSS, xunit_out+dist*sind(ii), yunit_out+dist*cosd(ii));
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                plot(app.UIAxes_FSS, xunit_in+dist*sind(ii), yunit_in+dist*cosd(ii));
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                fill(app.UIAxes_FSS, [xunit_in+dist*sind(ii) flip(xunit_out+dist*sind(ii))], ...
                                    [yunit_in+dist*cosd(ii) flip(yunit_out+dist*cosd(ii))],'b')
                            end
                        case 'Rectangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            r_in=app.InnerRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;

                            th = 0:pi/50:2*pi;
                            xunit_out = r_out * cos(th);
                            yunit_out = r_out * sin(th);
                            xunit_in = r_in * cos(th);
                            yunit_in = r_in * sin(th);

                            for ii=45:90:315
                                plot(app.UIAxes_FSS, xunit_out+dist*sqrt(2)*sind(ii)/2, yunit_out+dist*sqrt(2)*cosd(ii)/2,'b');
                                hold(app.UIAxes_FSS,'on')
                                plot(app.UIAxes_FSS, xunit_in+dist*sqrt(2)*sind(ii)/2, yunit_in+dist*sqrt(2)*cosd(ii)/2,'b');
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                fill(app.UIAxes_FSS, [xunit_in+dist*sqrt(2)*sind(ii)/2 flip(xunit_out+dist*sqrt(2)*sind(ii)/2)] ...
                                    ,[yunit_in+dist*sqrt(2)*cosd(ii)/2 flip(yunit_out+dist*sqrt(2)*cosd(ii)/2)],'b')
                            end
                    end
                case 'Rectangular'
                    switch lattice
                        case 'Triangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;
                            plot(app.UIAxes_FSS,0,0) % to clear the precious plot, hold off is not working

                            rectangle(app.UIAxes_FSS,'Position',[-r_out/2 -r_out/2 r_out r_out],'FaceColor',[0 0 1]);
                            hold(app.UIAxes_FSS,'on')
                            set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])

                            for ii=0:60:300
                                rectangle(app.UIAxes_FSS,'Position',[-r_out/2+dist*sind(ii) -r_out/2+dist*cosd(ii) r_out r_out],'FaceColor',[0 0 1]);
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                            end
                        case 'Rectangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;
                            plot(app.UIAxes_FSS,0,0) % to clear the precious plot, hold off is not working

                            for ii=45:90:315
                                rectangle(app.UIAxes_FSS,'Position',[-r_out/2+dist*sqrt(2)*sind(ii)/2 -r_out/2+dist*sqrt(2)*cosd(ii)/2 r_out r_out],'FaceColor',[0 0 1]);
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                hold(app.UIAxes_FSS,'on')
                            end
                    end
                case 'Ractangular Loop'
                    switch lattice
                        case 'Triangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            r_in=app.InnerRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;
                            plot(app.UIAxes_FSS,0,0) % to clear the precious plot, hold off is not working

                            rectangle(app.UIAxes_FSS,'Position',[-r_out/2 -r_out/2 r_out r_out],'FaceColor',[0 0 1]);
                            hold(app.UIAxes_FSS,'on')
                            rectangle(app.UIAxes_FSS,'Position',[-r_in/2 -r_in/2 r_in r_in],'FaceColor',[1 1 1]);
                            set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])

                            for ii=0:60:300
                                rectangle(app.UIAxes_FSS,'Position',[-r_out/2+dist*sind(ii) -r_out/2+dist*cosd(ii) r_out r_out],'FaceColor',[0 0 1]);
                                rectangle(app.UIAxes_FSS,'Position',[-r_in/2+dist*sind(ii) -r_in/2+dist*cosd(ii) r_in r_in],'FaceColor',[1 1 1]);
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                            end
                        case 'Rectangular'
                            hold(app.UIAxes_FSS,'off')
                            r_out=app.OuterRadiusLengthmmEditField.Value;
                            dist=app.ElementSpacingmmEditField.Value;
                            r_in=app.InnerRadiusLengthmmEditField.Value;
                            plot(app.UIAxes_FSS,0,0) % to clear the precious plot, hold off is not working

                            for ii=45:90:315
                                rectangle(app.UIAxes_FSS,'Position',[-r_out/2+dist*sqrt(2)*sind(ii)/2 -r_out/2+dist*sqrt(2)*cosd(ii)/2 r_out r_out],'FaceColor',[0 0 1]);
                                rectangle(app.UIAxes_FSS,'Position',[-r_in/2+dist*sqrt(2)*sind(ii)/2 -r_in/2+dist*sqrt(2)*cosd(ii)/2 r_in r_in],'FaceColor',[1 1 1]);
                                set(app.UIAxes_FSS, "linewidth",2, "fontsize", 14 ,'DataAspectRatio',[1 1 1])
                                hold(app.UIAxes_FSS,'on')
                            end
                    end
                case 'Triangular'
                    horn_profile = 5;
                case 'Triangular Loop'
                    horn_profile = 6;
                case 'Hexagonal'
                    horn_profile = 7;
                case 'Hexagonal Loop'
                    horn_profile = 8;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1037 598];
            app.UIFigure.Name = 'MATLAB App';

            % Create UltimateRFDesignerToolLabel
            app.UltimateRFDesignerToolLabel = uilabel(app.UIFigure);
            app.UltimateRFDesignerToolLabel.HorizontalAlignment = 'center';
            app.UltimateRFDesignerToolLabel.FontSize = 20;
            app.UltimateRFDesignerToolLabel.FontWeight = 'bold';
            app.UltimateRFDesignerToolLabel.FontColor = [0.149 0.149 0.149];
            app.UltimateRFDesignerToolLabel.Position = [367 574 311 25];
            app.UltimateRFDesignerToolLabel.Text = 'Ultimate RF Designer Tool';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [16 17 1012 524];

            % Create CorrugatedAntennaTab
            app.CorrugatedAntennaTab = uitab(app.TabGroup);
            app.CorrugatedAntennaTab.Title = 'Corrugated Antenna';
            app.CorrugatedAntennaTab.BackgroundColor = [0.9412 0.9412 0.9412];

            % Create UIAxes
            app.UIAxes = uiaxes(app.CorrugatedAntennaTab);
            title(app.UIAxes, 'Complete Corrugated Horn Profile')
            xlabel(app.UIAxes, 'Dimension in z Direction (mm)')
            ylabel(app.UIAxes, 'Dimension in y Direction (mm)')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.FontWeight = 'bold';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.FontSize = 12;
            app.UIAxes.Position = [284 64 713 427];

            % Create PreviewButton
            app.PreviewButton = uibutton(app.CorrugatedAntennaTab, 'push');
            app.PreviewButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewButtonPushed, true);
            app.PreviewButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.PreviewButton.Position = [583 23 100 22];
            app.PreviewButton.Text = 'Preview';

            % Create ConstructButton
            app.ConstructButton = uibutton(app.CorrugatedAntennaTab, 'push');
            app.ConstructButton.ButtonPushedFcn = createCallbackFcn(app, @ConstructButtonPushed, true);
            app.ConstructButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.ConstructButton.Position = [746 23 100 22];
            app.ConstructButton.Text = 'Construct';

            % Create HornProfileDropDownLabel
            app.HornProfileDropDownLabel = uilabel(app.CorrugatedAntennaTab);
            app.HornProfileDropDownLabel.HorizontalAlignment = 'right';
            app.HornProfileDropDownLabel.Position = [47 64 69 22];
            app.HornProfileDropDownLabel.Text = 'Horn Profile';

            % Create HornProfileDropDown
            app.HornProfileDropDown = uidropdown(app.CorrugatedAntennaTab);
            app.HornProfileDropDown.Items = {'Linear', 'Sinusoidal', 'Asymmetric Sine-Squared', 'Tangential', 'x.rho', 'Exponential', 'Hyperbolic', 'Poynomial'};
            app.HornProfileDropDown.Position = [131 64 100 22];
            app.HornProfileDropDown.Value = 'Linear';

            % Create ModeConverterTypeDropDownLabel
            app.ModeConverterTypeDropDownLabel = uilabel(app.CorrugatedAntennaTab);
            app.ModeConverterTypeDropDownLabel.HorizontalAlignment = 'right';
            app.ModeConverterTypeDropDownLabel.Position = [22 31 121 22];
            app.ModeConverterTypeDropDownLabel.Text = 'Mode Converter Type';

            % Create ModeConverterTypeDropDown
            app.ModeConverterTypeDropDown = uidropdown(app.CorrugatedAntennaTab);
            app.ModeConverterTypeDropDown.Items = {'Variable Slot Depth', 'Ring Loaded Slots', 'Variable p/w Slot'};
            app.ModeConverterTypeDropDown.Position = [158 31 139 22];
            app.ModeConverterTypeDropDown.Value = 'Variable Slot Depth';

            % Create MinimumFrequencyGHzEditFieldLabel
            app.MinimumFrequencyGHzEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.MinimumFrequencyGHzEditFieldLabel.HorizontalAlignment = 'right';
            app.MinimumFrequencyGHzEditFieldLabel.Position = [35 449 150 22];
            app.MinimumFrequencyGHzEditFieldLabel.Text = 'Minimum Frequency (GHz)';

            % Create MinimumFrequencyGHzEditField
            app.MinimumFrequencyGHzEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.MinimumFrequencyGHzEditField.Limits = [0 Inf];
            app.MinimumFrequencyGHzEditField.Position = [218 449 44 22];
            app.MinimumFrequencyGHzEditField.Value = 10.7;

            % Create MaximumFrequencyGHzEditFieldLabel
            app.MaximumFrequencyGHzEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.MaximumFrequencyGHzEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumFrequencyGHzEditFieldLabel.Position = [35 418 153 22];
            app.MaximumFrequencyGHzEditFieldLabel.Text = 'Maximum Frequency (GHz)';

            % Create MaximumFrequencyGHzEditField
            app.MaximumFrequencyGHzEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.MaximumFrequencyGHzEditField.Limits = [0 Inf];
            app.MaximumFrequencyGHzEditField.Position = [218 418 44 22];
            app.MaximumFrequencyGHzEditField.Value = 14.5;

            % Create OutputRadiusCoefficientEditFieldLabel
            app.OutputRadiusCoefficientEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.OutputRadiusCoefficientEditFieldLabel.HorizontalAlignment = 'right';
            app.OutputRadiusCoefficientEditFieldLabel.Position = [35 387 142 22];
            app.OutputRadiusCoefficientEditFieldLabel.Text = 'Output Radius Coefficient';

            % Create OutputRadiusCoefficientEditField
            app.OutputRadiusCoefficientEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.OutputRadiusCoefficientEditField.Limits = [0 Inf];
            app.OutputRadiusCoefficientEditField.Position = [218 387 44 22];
            app.OutputRadiusCoefficientEditField.Value = 1.95;

            % Create PitchEditFieldLabel
            app.PitchEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.PitchEditFieldLabel.HorizontalAlignment = 'right';
            app.PitchEditFieldLabel.Position = [35 356 32 22];
            app.PitchEditFieldLabel.Text = 'Pitch';

            % Create PitchEditField
            app.PitchEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.PitchEditField.Limits = [5 10];
            app.PitchEditField.Position = [218 356 44 22];
            app.PitchEditField.Value = 8;

            % Create PitchtoWidthpWRatioEditFieldLabel
            app.PitchtoWidthpWRatioEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.PitchtoWidthpWRatioEditFieldLabel.HorizontalAlignment = 'right';
            app.PitchtoWidthpWRatioEditFieldLabel.Position = [35 325 144 22];
            app.PitchtoWidthpWRatioEditFieldLabel.Text = 'Pitch to Width (p/W) Ratio';

            % Create PitchtoWidthpWRatioEditField
            app.PitchtoWidthpWRatioEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.PitchtoWidthpWRatioEditField.Limits = [0 1];
            app.PitchtoWidthpWRatioEditField.Position = [218 325 44 22];
            app.PitchtoWidthpWRatioEditField.Value = 0.8;

            % Create SigmaEditFieldLabel
            app.SigmaEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.SigmaEditFieldLabel.HorizontalAlignment = 'right';
            app.SigmaEditFieldLabel.Position = [35 295 40 22];
            app.SigmaEditFieldLabel.Text = 'Sigma';

            % Create SigmaEditField
            app.SigmaEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.SigmaEditField.Limits = [0.4 0.5];
            app.SigmaEditField.Position = [218 295 44 22];
            app.SigmaEditField.Value = 0.42;

            % Create InputWaveguideLengthmmEditFieldLabel
            app.InputWaveguideLengthmmEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.InputWaveguideLengthmmEditFieldLabel.HorizontalAlignment = 'right';
            app.InputWaveguideLengthmmEditFieldLabel.Position = [35 265 166 22];
            app.InputWaveguideLengthmmEditFieldLabel.Text = 'Input Waveguide Length (mm)';

            % Create InputWaveguideLengthmmEditField
            app.InputWaveguideLengthmmEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.InputWaveguideLengthmmEditField.Limits = [0 Inf];
            app.InputWaveguideLengthmmEditField.Position = [218 265 44 22];
            app.InputWaveguideLengthmmEditField.Value = 30;

            % Create NumberofCorrugationsEditFieldLabel
            app.NumberofCorrugationsEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.NumberofCorrugationsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofCorrugationsEditFieldLabel.Position = [35 235 134 22];
            app.NumberofCorrugationsEditFieldLabel.Text = 'Number of Corrugations';

            % Create NumberofCorrugationsEditField
            app.NumberofCorrugationsEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.NumberofCorrugationsEditField.Limits = [0 Inf];
            app.NumberofCorrugationsEditField.Position = [218 235 44 22];
            app.NumberofCorrugationsEditField.Value = 60;

            % Create NumberofModeConvertersEditFieldLabel
            app.NumberofModeConvertersEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.NumberofModeConvertersEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofModeConvertersEditFieldLabel.Position = [35 205 157 22];
            app.NumberofModeConvertersEditFieldLabel.Text = 'Number of Mode Converters';

            % Create NumberofModeConvertersEditField
            app.NumberofModeConvertersEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.NumberofModeConvertersEditField.Limits = [0 Inf];
            app.NumberofModeConvertersEditField.Position = [218 205 44 22];
            app.NumberofModeConvertersEditField.Value = 5;

            % Create RingLoadedMCDeltaEditFieldLabel
            app.RingLoadedMCDeltaEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.RingLoadedMCDeltaEditFieldLabel.HorizontalAlignment = 'right';
            app.RingLoadedMCDeltaEditFieldLabel.Position = [35 175 127 22];
            app.RingLoadedMCDeltaEditFieldLabel.Text = 'Ring Loaded MC Delta';

            % Create RingLoadedMCDeltaEditField
            app.RingLoadedMCDeltaEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.RingLoadedMCDeltaEditField.Limits = [0 1];
            app.RingLoadedMCDeltaEditField.Position = [218 175 44 22];
            app.RingLoadedMCDeltaEditField.Value = 0.1;

            % Create VariablepwslotMCDeltaEditFieldLabel
            app.VariablepwslotMCDeltaEditFieldLabel = uilabel(app.CorrugatedAntennaTab);
            app.VariablepwslotMCDeltaEditFieldLabel.HorizontalAlignment = 'right';
            app.VariablepwslotMCDeltaEditFieldLabel.Position = [35 145 146 22];
            app.VariablepwslotMCDeltaEditFieldLabel.Text = 'Variable p/w slot MC Delta';

            % Create VariablepwslotMCDeltaEditField
            app.VariablepwslotMCDeltaEditField = uieditfield(app.CorrugatedAntennaTab, 'numeric');
            app.VariablepwslotMCDeltaEditField.Limits = [0 1];
            app.VariablepwslotMCDeltaEditField.Position = [218 145 44 22];
            app.VariablepwslotMCDeltaEditField.Value = 0.125;

            % Create FSSTab
            app.FSSTab = uitab(app.TabGroup);
            app.FSSTab.Title = 'FSS';

            % Create UIAxes_FSS
            app.UIAxes_FSS = uiaxes(app.FSSTab);
            title(app.UIAxes_FSS, 'Unit Cell Configuration')
            xlabel(app.UIAxes_FSS, 'Dimension in z Direction (mm)')
            ylabel(app.UIAxes_FSS, 'Dimension in y Direction (mm)')
            zlabel(app.UIAxes_FSS, 'Z')
            app.UIAxes_FSS.FontWeight = 'bold';
            app.UIAxes_FSS.XGrid = 'on';
            app.UIAxes_FSS.YGrid = 'on';
            app.UIAxes_FSS.FontSize = 12;
            app.UIAxes_FSS.Position = [439 64 558 427];

            % Create UnitCellTypeDropDownLabel
            app.UnitCellTypeDropDownLabel = uilabel(app.FSSTab);
            app.UnitCellTypeDropDownLabel.Position = [23 449 78 22];
            app.UnitCellTypeDropDownLabel.Text = 'Unit Cell Type';

            % Create UnitCellTypeDropDown
            app.UnitCellTypeDropDown = uidropdown(app.FSSTab);
            app.UnitCellTypeDropDown.Items = {'Circular', 'Circular Loop', 'Rectangular', 'Ractangular Loop', 'Hexagonal', 'Hexagonal Loop'};
            app.UnitCellTypeDropDown.Position = [132 449 139 22];
            app.UnitCellTypeDropDown.Value = 'Circular';

            % Create LatticeStructureDropDownLabel
            app.LatticeStructureDropDownLabel = uilabel(app.FSSTab);
            app.LatticeStructureDropDownLabel.Position = [23 418 94 22];
            app.LatticeStructureDropDownLabel.Text = 'Lattice Structure';

            % Create LatticeStructureDropDown
            app.LatticeStructureDropDown = uidropdown(app.FSSTab);
            app.LatticeStructureDropDown.Items = {'Triangular', 'Rectangular'};
            app.LatticeStructureDropDown.Position = [132 418 139 22];
            app.LatticeStructureDropDown.Value = 'Triangular';

            % Create PreviewFSSButton
            app.PreviewFSSButton = uibutton(app.FSSTab, 'push');
            app.PreviewFSSButton.ButtonPushedFcn = createCallbackFcn(app, @PreviewFSSButtonPushed, true);
            app.PreviewFSSButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.PreviewFSSButton.Position = [612 23 100 22];
            app.PreviewFSSButton.Text = 'Preview';

            % Create ConstructFSSButton
            app.ConstructFSSButton = uibutton(app.FSSTab, 'push');
            app.ConstructFSSButton.BackgroundColor = [0.9412 0.9412 0.9412];
            app.ConstructFSSButton.Position = [775 23 100 22];
            app.ConstructFSSButton.Text = 'Construct';

            % Create OuterRadiusLengthmmEditFieldLabel
            app.OuterRadiusLengthmmEditFieldLabel = uilabel(app.FSSTab);
            app.OuterRadiusLengthmmEditFieldLabel.Position = [23 378 148 22];
            app.OuterRadiusLengthmmEditFieldLabel.Text = 'Outer Radius/Length (mm)';

            % Create OuterRadiusLengthmmEditField
            app.OuterRadiusLengthmmEditField = uieditfield(app.FSSTab, 'numeric');
            app.OuterRadiusLengthmmEditField.Limits = [0 Inf];
            app.OuterRadiusLengthmmEditField.Position = [205 378 44 22];
            app.OuterRadiusLengthmmEditField.Value = 5;

            % Create InnerRadiusLengthmmEditFieldLabel
            app.InnerRadiusLengthmmEditFieldLabel = uilabel(app.FSSTab);
            app.InnerRadiusLengthmmEditFieldLabel.Position = [23 347 145 22];
            app.InnerRadiusLengthmmEditFieldLabel.Text = 'Inner Radius/Length (mm)';

            % Create InnerRadiusLengthmmEditField
            app.InnerRadiusLengthmmEditField = uieditfield(app.FSSTab, 'numeric');
            app.InnerRadiusLengthmmEditField.Limits = [0 Inf];
            app.InnerRadiusLengthmmEditField.Position = [205 347 44 22];
            app.InnerRadiusLengthmmEditField.Value = 4.5;

            % Create ElementSpacingmmEditFieldLabel
            app.ElementSpacingmmEditFieldLabel = uilabel(app.FSSTab);
            app.ElementSpacingmmEditFieldLabel.Position = [23 316 128 22];
            app.ElementSpacingmmEditFieldLabel.Text = 'Element Spacing (mm)';

            % Create ElementSpacingmmEditField
            app.ElementSpacingmmEditField = uieditfield(app.FSSTab, 'numeric');
            app.ElementSpacingmmEditField.Limits = [0 Inf];
            app.ElementSpacingmmEditField.Position = [205 316 44 22];
            app.ElementSpacingmmEditField.Value = 15;

            % Create ReflectorAntennaTab
            app.ReflectorAntennaTab = uitab(app.TabGroup);
            app.ReflectorAntennaTab.Title = 'Reflector Antenna';

            % Create SeptumPolarizerTab
            app.SeptumPolarizerTab = uitab(app.TabGroup);
            app.SeptumPolarizerTab.Title = 'Septum Polarizer';

            % Create WaveguideFilterTab
            app.WaveguideFilterTab = uitab(app.TabGroup);
            app.WaveguideFilterTab.Title = 'Waveguide Filter';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = UAT_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end