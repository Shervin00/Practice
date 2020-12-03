clc
% Initialize variables
TimeClose1 = zeros(40,1);
sr=50;
TrialCloseTimeP1 = zeros(137,40);
TrialCloseTimeP2 = zeros(135,40);
SquareCloseTime1 = zeros(137,40);
QuestionCloseTimeP2 = zeros(135,40);
AnswerCloseTimeP1 = zeros(137,40);
AnswerCloseTimeP2 = zeros(135,40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this section read in data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define file path 
% EyeFolder = 'M:\eyetracker\txt_files';
%EyeFolder = 'P:\3013070.01\Physiology_Analysis\Iwan_Analysis\Scripts_Iwan\eyetracker\txt_files';
EyeFolder = 'M:\Documents\MATLAB\Shervin_Scripts\txt_filesSherv';
if exist(EyeFolder, 'dir') ~= 7
  Message = sprintf('Error: The following folder does not exist:\n%s', EyeFolder);
  uiwait(warndlg(Message));
  return;
end
filePattern = fullfile(EyeFolder, '*.txt');
Samplefiles = dir(filePattern);

%% Start loop for reading all the data
for k = 1:length(Samplefiles); EyeFolder;
    
%% Read in all data
%% Read in Eyetracking data
 baseFileName = Samplefiles(k).name;
  Eyetrackingfiles = fullfile(EyeFolder, baseFileName);
  fprintf('Now reading %s\n', Eyetrackingfiles);
  

%% Read in first part of the study output
%StudyFolder1 = 'P:\3013070.01\Analysis\Scripts_Iwan\LC';
StudyFolder1 = 'M:\Documents\MATLAB\Shervin_Scripts\LC_logfiles';
if exist(StudyFolder1, 'dir') ~= 7
  Message = sprintf('Error: The following folder does not exist:\n%s', StudyFolder1);
  uiwait(warndlg(Message));
  return;
end

filePattern1 = fullfile(StudyFolder1, '*.log');
StudytxtFiles1   = dir(filePattern1);

  baseFileName1 = StudytxtFiles1(k).name;
  Studyfiles1 = fullfile(StudyFolder1, baseFileName1);
  fprintf('Now reading %s\n', Studyfiles1);
  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section opens the files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Open Eyetracking data
startRow = 47; %39
endRow = inf;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
delimiter = '\t';

 fileID(k)=fopen(strcat(EyeFolder,filesep,Samplefiles(k).name));
  EyetrackingData =textscan(fileID(k), formatSpec,'delimiter', delimiter, 'ReturnOnError', false);
  for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(EyetrackingData)
        EyetrackingData{col} = [EyetrackingData{col};dataArrayBlock{col}];
    end
  end
  
% Allocate imported array to column variable names
VarName1 = EyetrackingData{:, 1};
VarName2 = EyetrackingData{:, 2};
VarName3 = EyetrackingData{:, 3};
VarName4 = EyetrackingData{:, 4};
VarName5 = EyetrackingData{:, 5};
VarName6 = EyetrackingData{:, 6};
VarName7 = EyetrackingData{:, 7};
VarName8 = EyetrackingData{:, 8};
VarName9 = EyetrackingData{:, 9};
VarName10 = EyetrackingData{:, 10};
VarName11 = EyetrackingData{:, 11};   
VarName12 = EyetrackingData{:, 12};
VarName13 = EyetrackingData{:, 13};
VarName14 = EyetrackingData{:, 18};

%change variable names
time = VarName1;
smp = VarName2;
third = VarName3;
raw_x = VarName4;
raw_y = VarName5;
dia_x = VarName6;
dia_y = VarName7;
cr1_x = VarName8;
cr1_y = VarName9;
por_x = VarName10;
por_y = VarName11;
timing = VarName12;
pupil = VarName13;
trigger = VarName14; 

% clear variables 
clearvars filename delimiter startRow formatSpec fileID EyetrackingData ans;

% Index Scanning triggers
TriggerIndex=find(strcmp(trigger,'256		'));

%% Open Studyoutput1
% Initialize variables
delimiter = '\t';
startRow = 7;
formatSpec = '%s%s%s%s%s%*s%[^\n\r]';

fileID(k)=fopen(strcat(StudyFolder1,filesep,StudytxtFiles1(k).name));
StudyData1 = textscan(fileID(k), formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Allocate imported array to column variable names
subject = StudyData1{1, 1}(1,1);
pulse = StudyData1{:, 2};
event = StudyData1{:, 3};
code = StudyData1{:, 4};
timeLC = StudyData1{:, 5};
uncertainty = StudyData1{:, 6};

% clear variables 
clearvars filename delimiter startRow formatSpec fileID StudyData1 ans;
 
%% Get onset timing of red X, red O, green X and green O squares

RedXTime = zeros(16,1);
GreenXTime = zeros(16,1);
RedOTime = zeros(16,1);
GreenOTime= zeros(16,1);
color = 0;
RedXAmount = 1;
GreenXAmount = 1;
RedOAmount = 1;
GreenOAmount = 1;
timeLCmatrix = cellfun(@str2double,timeLC);
for m=1:length(code)
    if color == 0 || color == 1
        if contains(code{m},'redx')
            code{m} = 'fredx';
            RedXTime(RedXAmount,1) = timeLCmatrix(m,1) - timeLCmatrix(4,1);
            RedXAmount = RedXAmount + 1;
            color = 2;
        end
    end
    if color == 0 || color == 2
        if contains(code{m},'greenx')
            code{m} = 'fgreenx';
            GreenXTime(GreenXAmount,1) = timeLCmatrix(m,1) - timeLCmatrix(4,1);
            GreenXAmount = GreenXAmount + 1;
            color = 1;
        end
    end
     if color == 0 || color == 1
        if contains(code{m},'redo')
            code{m} = 'fredo';
            RedOTime(RedOAmount,1) = timeLCmatrix(m,1) - timeLCmatrix(4,1);
            RedOAmount = RedOAmount + 1;
            color = 2;
        end
     end
                if color == 0 || color == 1
        if contains(code{m},'greeno')
            code{m} = 'greeno';
            GreenOTime(GreenOAmount,1) = timeLCmatrix(m,1) - timeLCmatrix(4,1);
            GreenOAmount = GreenOAmount + 1;
            color = 2;
        end
                end
                
end
EndTime(1,1) = timeLCmatrix(m) - timeLCmatrix(4);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section uses the triggers to find the right task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear TriggerBlock

test = 1;
TriggerBlock(test,1) = TriggerIndex(1);

for jj=1:length(TriggerIndex)-1
    if TriggerIndex(jj+1)-TriggerIndex(jj)> 1000
        TriggerBlock(test,2) = TriggerIndex(jj);
        TriggerBlock(test,3) = TriggerBlock(test,2) - TriggerBlock(test,1);        
        test = test + 1;
        TriggerBlock(test,1) = TriggerIndex(jj+1);
    end
end
TriggerBlock(test,2) = TriggerIndex(jj);
TriggerBlock(test,3) = TriggerBlock(test,1) - TriggerBlock(test-1,1);

temp = TriggerBlock(1:end-1,:);
temp = sortrows(temp, 3);
temp = temp(end-1,:);

TriggerTask{k} = TriggerBlock;
if temp(1,3) < 53000
    if temp(1,3) > 45000
        TriggerLC{k} = temp;
    else
    TriggerLC{k} = 0;
    end
else
    TriggerLC{k} = 0;
end
    


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Eye Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eye = dia_y;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section calculates the eye close time in the Question part of each
% trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Study Datapart 1

% Find all lines with zeros in them
ZeroLines=find(strcmp(Eye,'0.000'));

for m=1:length(RedXTime)  
    % Define start of the squares of each trial in the study data of presentation
    StartRedXSquare = RedXTime(m)/10000;
    if RedXTime(1) < GreenXTime(1)
        EndRedXSquare = GreenXTime(m)/10000;
    else
        if m == 16
            EndRedXSquare = EndTime(1)/10000;
        else
            EndRedXSquare = GreenXTime(m+1)/10000;
        end
    end  
    %eyetracker
    StartRedXSquareEye = StartRedXSquare * sr + temp(1,1);
    EndRedSquareEye = EndRedXSquare * sr + temp(1,1);
    %zerolines
    ZeroRedXSquareStart = find(ZeroLines > StartRedXSquareEye);
    ZeroRedXSquareEnd = find (ZeroLines < EndRedSquareEye);
    ZeroSquares = intersect(ZeroRedXSquareStart,ZeroRedXSquareEnd);
    EyeCloseTimeSquares = length(ZeroSquares)/sr;
    %save
    SquareCloseTimeRedX(k,m) = EyeCloseTimeSquares;

    %presentation
    StartGreenXSquare = GreenXTime(m)/10000;
    if GreenXTime(1) < RedXTime(1)
        EndGreenXSquare = RedXTime(m)/10000;
    else
        if m ==16
            EndGreenXSquare = EndTime(1)/10000;
        else
            EndGreenXSquare = RedXTime(m+1)/10000;
        end
    end
    %eyetracker
    StartGreenXSquareEye = StartGreenXSquare * sr + temp(1,1);
    EndGreenXSquareEye = EndGreenXSquare * sr + temp(1,1);
    %zerolines
    ZeroGreenXSquareStart = find(ZeroLines > StartGreenXSquareEye);
    ZeroGreenXSquareEnd = find (ZeroLines < EndGreenXSquareEye);
    ZeroSquares = intersect(ZeroGreenXSquareStart,ZeroGreenXSquareEnd);
    EyeCloseTimeSquares = length(ZeroSquares)/50;
    %save
    SquareCloseTimeGreenX(k,m) = EyeCloseTimeSquares;

    % Define start of the squares of each trial in the study data of presentation
    StartRedOSquare = RedOTime(m)/10000;
    if RedOTime(1) < GreenOTime(1)
        EndRedOSquare = GreenOTime(m)/10000;
    else
        if m == 16
            EndRedOSquare = EndTime(1)/10000;
        else
            EndRedOSquare = GreenOTime(m+1)/10000;
        end
    end  
%eyetracker
    StartRedOSquareEye = StartRedOSquare * sr + temp(1,1);
    EndRedOSquareEye = EndRedOSquare * sr + temp(1,1);
    %zerolines
    ZeroRedOSquareStart = find(ZeroLines > StartRedOSquareEye);
    ZeroRedOSquareEnd = find (ZeroLines < EndRedOSquareEye);
    ZeroSquares = intersect(ZeroRedOSquareStart,ZeroRedOSquareEnd);
    EyeCloseTimeSquares = length(ZeroSquares)/sr;
    %save
    SquareCloseTimeRedO(k,m) = EyeCloseTimeSquares;

    StartGreenOSquare = GreenOTime(m)/10000;
    if GreenOTime(1) < RedOTime(1)
        EndGreenOSquare = RedOTime(m)/10000;
    else
        if m ==16
            EndGreenOSquare = EndTime(1)/10000;
        else
            EndGreenOSquare = RedOTime(m+1)/10000;
        end
    end
    %eyetracker
    StartGreenOSquareEye = StartGreenOSquare * sr + temp(1,1);
    EndGreenOSquareEye = EndGreenOSquare * sr + temp(1,1);
    %zerolines
    ZeroGreenOSquareStart = find(ZeroLines > StartGreenXSquareEye);
    ZeroGreenOSquareEnd = find (ZeroLines < EndGreenXSquareEye);
    ZeroSquares = intersect(ZeroGreenOSquareStart,ZeroGreenOSquareEnd);
    EyeCloseTimeSquares = length(ZeroSquares)/50;
    %save
    SquareCloseTimeGreenO(k,m) = EyeCloseTimeSquares;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section calculates the mean pupil dilation during red X squares, red O or
% green X or green O squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eye_double = cellfun(@str2double,Eye);

for m=1:length(RedXTime)  
    % Define start of the squares of each trial in the study data of presentation
    StartRedXSquare = RedXTime(m)/10000;
    if RedXTime(1) < GreenXTime(1)
        EndRedXSquare = GreenXTime(m)/10000;
    else
        if m == 16
            EndRedXSquare = EndTime(1)/10000;
        else
            EndRedXSquare = GreenTime(m+1)/10000;
        end
    end  
    %eyetracker
    StartRedXSquareEye = round(StartRedXSquare * sr + temp(1,1));
    EndRedSquareEye = round(EndRedXSquare * sr + temp(1,1));
    %zerolines
    Dilation = 0;
    ZeroDilation = 0;
    kk = 1;
    for mm=StartRedXSquareEye:EndRedSquareEye
        if Eye_double(mm,1) == 0
            ZeroDilation = ZeroDilation + 1;
        else
            Dilation(kk,1) = Eye_double(mm,1);
            kk = kk + 1;
        end
    end
    MeanDilation = mean(Dilation);
    Std = std(Dilation);
    kk = 1;
    error = 0;
    for mm=1:length(Dilation(:,1))
        if Dilation(mm,1) > MeanDilation-(4*Std)
            if Dilation(mm,1) < MeanDilation+(4*Std)
                RealDilation(kk,1) = Dilation(mm,1);
                kk = kk + 1;
            else 
                error = error + 1;
            end
        else
            error = error + 1;
        end
    end
    MeanDilationRedX(k,m) = mean(RealDilation);
    ErrorRedX(k,m) = error;

    %presentation
    StartGreenXSquare = GreenXTime(m)/10000;
    if GreenXTime(1) < RedXTime(1)
        EndGreenXSquare = RedXTime(m)/10000;
    else
        if m ==16
            EndGreenXSquare = EndTime(1)/10000;
        else
            EndGreenXSquare = RedXTime(m+1)/10000;
        end
    end
    %eyetracker
    StartGreenXSquareEye = round(StartGreenXSquare * sr + temp(1,1));
    EndGreenXSquareEye = round(EndGreenXSquare * sr + temp(1,1));
     %zerolines
    Dilation = 0;
    ZeroDilation = 0;
    kk = 1;
    for mm=StartGreenXSquareEye:EndGreenXSquareEye
        if Eye_double(mm,1) == 0
            ZeroDilation = ZeroDilation + 1;
        else
            Dilation(kk,1) = Eye_double(mm,1);
            kk = kk + 1;
        end
    end
    MeanDilation = mean(Dilation);
    Std = std(Dilation);
    kk = 1;
    error = 0;
    for mm=1:length(Dilation(:,1))
        if Dilation(mm,1) > MeanDilation-(4*Std)
            if Dilation(mm,1) < MeanDilation+(4*Std)
                RealDilation(kk,1) = Dilation(mm,1);
                kk = kk + 1;
            else 
                error = error + 1;
            end
        else
            error = error + 1;
        end
    end
    MeanDilationGreenX(k,m) = mean(RealDilation);
    ErrorGreenX(k,m) = error;
    
       % Define start of the squares of each trial in the study data of presentation
    StartRedOSquare = RedOTime(m)/10000;
    if RedOTime(1) < GreenOTime(1)
        EndRedOSquare = GreenOTime(m)/10000;
    else
        if m == 16
            EndRedOSquare = EndTime(1)/10000;
        else
            EndRedOSquare = GreenTime(m+1)/10000;
        end
    end  
    %eyetracker
    StartRedXSquareEye = round(StartRedXSquare * sr + temp(1,1));
    EndRedSquareEye = round(EndRedOSquare * sr + temp(1,1));
    %zerolines
    Dilation = 0;
    ZeroDilation = 0;
    kk = 1;
    for mm=StartRedOSquareEye:EndRedSquareEye
        if Eye_double(mm,1) == 0
            ZeroDilation = ZeroDilation + 1;
        else
            Dilation(kk,1) = Eye_double(mm,1);
            kk = kk + 1;
        end
    end
    MeanDilation = mean(Dilation);
    Std = std(Dilation);
    kk = 1;
    error = 0;
    for mm=1:length(Dilation(:,1))
        if Dilation(mm,1) > MeanDilation-(4*Std)
            if Dilation(mm,1) < MeanDilation+(4*Std)
                RealDilation(kk,1) = Dilation(mm,1);
                kk = kk + 1;
            else 
                error = error + 1;
            end
        else
            error = error + 1;
        end
    end
    MeanDilationRedO(k,m) = mean(RealDilation);
    ErrorRedO(k,m) = error;
    
    %presentation
    StartGreenOSquare = GreenOTime(m)/10000;
    if GreenOTime(1) < RedOTime(1)
        EndGreenOSquare = RedOTime(m)/10000;
    else
        if m ==16
            EndGreenOSquare = EndTime(1)/10000;
        else
            EndGreenOSquare = RedOTime(m+1)/10000;
        end
    end
    %eyetracker
    StartGreenOSquareEye = round(StartGreenOSquare * sr + temp(1,1));
    EndGreenOSquareEye = round(EndGreenOSquare * sr + temp(1,1));
     %zerolines
    Dilation = 0;
    ZeroDilation = 0;
    kk = 1;
    for mm=StartGreenOSquareEye:EndGreenOSquareEye
        if Eye_double(mm,1) == 0
            ZeroDilation = ZeroDilation + 1;
        else
            Dilation(kk,1) = Eye_double(mm,1);
            kk = kk + 1;
        end
    end
    MeanDilation = mean(Dilation);
    Std = std(Dilation);
    kk = 1;
    error = 0;
    for mm=1:length(Dilation(:,1))
        if Dilation(mm,1) > MeanDilation-(4*Std)
            if Dilation(mm,1) < MeanDilation+(4*Std)
                RealDilation(kk,1) = Dilation(mm,1);
                kk = kk + 1;
            else 
                error = error + 1;
            end
        else
            error = error + 1;
        end
    end
    MeanDilationGreenO(k,m) = mean(RealDilation);
    ErrorGreenO(k,m) = error;
end


%%
% Clear variables
clearvars filename delimiter startRow formatSpec fileID EyetrackingData ans;
end
%% Show the times that the eyes are closed in each participant and export to an excel file
SquareCloseTimeRedX;
xlswrite('SquareCloseTimeRedXLArea.xls',SquareCloseTimeRedX);

SquareCloseTimeGreenX;
xlswrite('SquareCloseTimeGreenXLArea.xls',SquareCloseTimeGreenX);

SquareCloseTimeRedO;
xlswrite('SquareCloseTimeRedXLArea.xls',SquareCloseTimeRedX);

SquareCloseTimeGreenO;
xlswrite('SquareCloseTimeGreenXLArea.xls',SquareCloseTimeGreenX);


%Show the mean dilation for each participant for each condition for each
%block
MeanDilationRedX;
xlswrite('MeanDilationRedOLArea.xls',MeanDilationRedO);

MeanDilationGreenX;
xlswrite('MeanDilationGreenOLArea.xls', MeanDilationGreenO);

MeanDilationRedO;
xlswrite('MeanDilationRedOLArea.xls',MeanDilationRedO);

MeanDilationGreenO;
xlswrite('MeanDilationGreenOLArea.xls', MeanDilationGreenO);

% MeanDilation = [MeanDilationRed;MeanDilationGreen];
% xlswrite('MeanDilation.xls', MeanDilation);