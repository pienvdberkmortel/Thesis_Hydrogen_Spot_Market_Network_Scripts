%% This function loads wind data from an excel file containing wind measurements from the KNMI from 2019. 

function Wind_data_2019 = Import_KNMI_Wind_Data_2019()
    try
        % Set up the Import Options and import the data
        opts = spreadsheetImportOptions("NumVariables", 5);

        % Specify sheet and range
        opts.Sheet = "Blad1";
        opts.DataRange = "A2:E8762";

        % Specify column names and types
        opts.VariableNames = ["Var1", "Var2", "FH", "FF", "windspeed"];
        opts.SelectedVariableNames = ["FH", "FF", "windspeed"];
        opts.VariableTypes = ["char", "char", "double", "double", "double"];

        % Specify file level properties
        opts.ImportErrorRule = "omitrow";
        opts.MissingRule = "omitrow";

        % Specify variable properties
        opts = setvaropts(opts, ["Var1", "Var2"], "WhitespaceRule", "preserve");
        opts = setvaropts(opts, ["Var1", "Var2"], "EmptyFieldRule", "auto");
        opts = setvaropts(opts, ["FH", "FF", "windspeed"], "TreatAsMissing", '');

        % Import the data
        % Make sure the path is correct
        Wind_data_2019 = readtable("C:\Users\pienv\OneDrive\Documenten\Master S&C\S&C Thesis Economic engineering\MATLAB\Matlab and Simulink files\Thesis_Repository\KNMI_Wind_Data_2019.xlsx", opts, "UseExcel", false); 

        % Suppress table display
        disp('');

    catch ME
        % Display the error message
        disp(getReport(ME));
    end
end