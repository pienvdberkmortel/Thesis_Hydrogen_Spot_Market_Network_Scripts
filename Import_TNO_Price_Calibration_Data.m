%% This function is used to calibrate hydrogen prices. Data is used form the TNO I-ELGAS Model


function Calibration_prices_2035 = Import_TNO_Price_Calibration_Data()
    try
        % Set up the Import Options and import the data
        opts = spreadsheetImportOptions("NumVariables", 19);
        opts.Sheet = 'Prijzen';
        opts.DataRange = 'A1:S8766';
        opts.SelectedVariableNames = {'Var16', 'Var17', 'Var18', 'Var19'};

        % Import the data
        Calibration_prices_2035 = readtable("C:\Users\pienv\OneDrive\Documenten\Master S&C\S&C Thesis Economic engineering\MATLAB\Matlab and Simulink files\Thesis_Repository\TNO_Price_Calibration_Data.xlsx", opts, 'UseExcel', false);

        % Delete the first three rows
        Calibration_prices_2035(1:3, :) = [];

        % Delete the last three rows
        Calibration_prices_2035(end-2:end, :) = [];

        % Convert selected variables to double
        Calibration_prices_2035.Var16 = str2double(Calibration_prices_2035.Var16);
        Calibration_prices_2035.Var17 = str2double(Calibration_prices_2035.Var17);
        Calibration_prices_2035.Var18 = str2double(Calibration_prices_2035.Var18);
        Calibration_prices_2035.Var19 = str2double(Calibration_prices_2035.Var19);

        % Suppress table display
        disp('');

    catch ME
        % Display the error message
        disp(getReport(ME));
    end
end
