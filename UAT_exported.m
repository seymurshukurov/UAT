classdef UAT_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        ReflectorAntennaDesignButton   matlab.ui.control.Button
        CorrugatedAntennaDesignButton  matlab.ui.control.Button
        UltimateAntennaToolLabel       matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: CorrugatedAntennaDesignButton
        function CorrugatedAntennaDesignButtonPushed(app, event)
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 632 244];
            app.UIFigure.Name = 'MATLAB App';

            % Create UltimateAntennaToolLabel
            app.UltimateAntennaToolLabel = uilabel(app.UIFigure);
            app.UltimateAntennaToolLabel.HorizontalAlignment = 'center';
            app.UltimateAntennaToolLabel.FontSize = 20;
            app.UltimateAntennaToolLabel.FontWeight = 'bold';
            app.UltimateAntennaToolLabel.Position = [213 220 216 25];
            app.UltimateAntennaToolLabel.Text = 'Ultimate Antenna Tool';

            % Create CorrugatedAntennaDesignButton
            app.CorrugatedAntennaDesignButton = uibutton(app.UIFigure, 'push');
            app.CorrugatedAntennaDesignButton.ButtonPushedFcn = createCallbackFcn(app, @CorrugatedAntennaDesignButtonPushed, true);
            app.CorrugatedAntennaDesignButton.FontWeight = 'bold';
            app.CorrugatedAntennaDesignButton.Position = [29 135 175 40];
            app.CorrugatedAntennaDesignButton.Text = 'Corrugated Antenna Design';

            % Create ReflectorAntennaDesignButton
            app.ReflectorAntennaDesignButton = uibutton(app.UIFigure, 'push');
            app.ReflectorAntennaDesignButton.FontWeight = 'bold';
            app.ReflectorAntennaDesignButton.Position = [236 135 162 40];
            app.ReflectorAntennaDesignButton.Text = 'Reflector Antenna Design';

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