# SAP2000-modeling-automization-using-MATLAB
## Introduction
Modeling buildings and extracting data using FEA programs is a time-consuming process. Additionally, significant effort is required for data storage, analysis, and parsing. However, in situations where parametric studies are permitted, various buildings can be automatically modeled through API, and data output and processing can be easily done using API as well. SAP2000 is considered the optimal FEA program due to its diversity in FEA and API support for all functions. The API code was developed using MATLAB, which is the most suitable interpreter for research environments.

>  SAP2000 provides APIs for various environments including MATLAB, C, and Python. However, due to incomplete documentation in the technical manual, some trial and error is necessary. ETABS, another product from the same company, offers only partial API functionality, resulting in many features of the code developed for SAP2000 being unusable when adapted for ETABS.

## Usage
### run_sap.m
First run `SAP2000` using `MATLAB` API

```matlab
% Set flags for program instance and path
AttachToInstance = false;
SpecifyPath = false;

% Define paths (replace with your actual paths)
ProgramPath = 'path/to/SAP2000.exe';
APIDLLPath = 'path/to/SAP2000v1.dll';
ModelDirectory = 'path/to/model/directory';

% Create model directory if it doesn't exist
if ~exist(ModelDirectory, 'dir')
    mkdir(ModelDirectory);
end

ModelName = "model_name.sdb";
ModelPath = fullfile(ModelDirectory, ModelName);

% Add SAP2000 API assembly
a = NET.addAssembly(APIDLLPath);
helper = SAP2000v1.Helper;
helper = NET.explicitCast(helper,'SAP2000v1.cHelper');

% Create or attach to SAP2000 instance
if AttachToInstance
    SapObject = helper.GetObject('CSI.SAP2000.API.SapObject');
else
    if SpecifyPath
        SapObject = helper.CreateObject(ProgramPath);
    else
        SapObject = helper.CreateObjectProgID('CSI.SAP2000.API.SapObject');
    end
    SapObject.ApplicationStart;
end

SapObject = NET.explicitCast(SapObject,'SAP2000v1.cOAPI');
SapModel = NET.explicitCast(SapObject.SapModel,'SAP2000v1.cSapModel');

run('main_body.m')
```
### main_body.m

Make the blank space(values are unimportant) and set the unit
```matlab
ret = SapModel.SetModelIsLocked(false());
ret = SapModel.InitializeNewModel(SAP2000v1.eUnits.kN_m_C);
File = NET.explicitCast(SapModel.File,'SAP2000v1.cFile');
ret = File.NewBlank;
Options = NET.explicitCast(SapModel.Options,'SAP2000v1.cOptions');
ret = Options.SetDimensions(12, 2, 1, 4,1,1,1,1,1,1,120);
Analyze = NET.explicitCast(SapModel.Analyze,'SAP2000v1.cAnalyze');
ret=Analyze.SetSolverOption_3(2, 2, 0, -1, NumThreads, "DEAD");
ret=Analyze.SetActiveDOF(DOF_option);
```


```matlab

-----------------------------------------------
% set the parameter or get the plan informations
-----------------------------------------------
% modeling
if model_mode=="MDOF"
    run('MDOF_modeling.m')
else
    run('modeling.m')
end
% get the design wind speed
Kzr=1.87566687625052*exp(-6.4856276*alpha)*H_bld^alpha;
V_H_bld=V_0*Kzr;


% run the model to get the modal informations (period, mode shape)
run('MDOF_pre_sap_running.m');
run('find_beta_MDOF.m')
run('get_period.m')

% run time history analysis
if model_mode=="MDOF"
    run('MDOF_sap_running.m');
else
    run('sap_running.m')
end

% correct wind tunnel test environment based on KDS code
run('MDOF_correction_run.m');

% get the response data
```
### modeling.m
```matlab
% Initialize model and set basic properties
ret = SapModel.SetModelIsLocked(false); % unlock model file after analysis
ret = SapModel.InitializeNewModel(SAP2000v1.eUnits.kN_m_C);
File = NET.explicitCast(SapModel.File, 'SAP2000v1.cFile');
ret = File.NewBlank; % reset the page
ret = NET.explicitCast(SapModel.Options, 'SAP2000v1.cOptions').SetDimensions(12, 2, 1, 4, 1, 1, 1, 1, 1, 1, 120);
ret = NET.explicitCast(SapModel.Analyze, 'SAP2000v1.cAnalyze').SetSolverOption_3(2, 0, 0, 0, 0, "DEAD");

% Define grid
N_X_line = 13;
N_Y_line = 13;
X_spacing = B_bld / (N_Y_line - 1);
Y_spacing = D_bld / (N_X_line - 1);

% Set up structural properties
setupStructuralProperties(SapModel);

% Create story constraints
for i = 1:N_story
    slName = "Slab" + num2str(i);
    ret = NET.explicitCast(SapModel.ConstraintDef, 'SAP2000v1.cConstraint').SetDiaphragm(slName, SAP2000v1.eConstraintAxis.Z);
end

% Create structural elements
createStructuralElements(SapModel, N_story, typ_story_H, SlabP, CB_str, CB_end, MCW_str, MCW_end, SCW_str, SCW_end, Col_str, Col_end, MB_str, MB_end, OB_str, OB_end);

% Set up groups and assignments
setupGroupsAndAssignments(SapModel, N_story, D_bld, B_bld, typ_story_H);

% Set up diaphragms
setupDiaphragms(SapModel,  typ_story_H);

% Save the model
ret = File.Save(ModelPath);

function setupStructuralProperties(SapModel)
    PropFrame = NET.explicitCast(SapModel.PropFrame, 'SAP2000v1.cPropFrame');
    PropArea = NET.explicitCast(SapModel.PropArea, 'SAP2000v1.cPropArea');
    
    % Set up frame properties
    frameProps = struct('H_COL', [1, 1], 'L_COL', [1, 1], 'H_MB', [1, 1], 'L_MB', [1, 1], 'H_CB', [1, 1], 'L_CB', [1, 1], 'H_OB', [1, 1], 'L_OB', [1, 1]);
    fields = fieldnames(frameProps);
    for i = 1:length(fields)
        ret = PropFrame.SetRectangle(fields{i}, "4000Psi", frameProps.(fields{i})(1), frameProps.(fields{i})(2), 255);
    end
    
    % Set up area properties
    areaProps = struct('SLAB', slab_t, 'LMW', core_main_t, 'LSW', core_sub_t, 'HMW', core_main_t, 'HSW', core_sub_t);
    fields = fieldnames(areaProps);
    for i = 1:length(fields)
        ret = PropArea.SetShell(fields{i}, 1, "4000psi", 0, areaProps.(fields{i}), areaProps.(fields{i}));
    end
    
    % Set modifiers (example for one property, repeat for others as needed)
    MyValue = NET.createArray('System.Double', 8);
    MyValue(1:4) = 1; MyValue(5:6) = 0.7; MyValue(7:8) = 1;
    ret = PropFrame.SetModifiers("H_COL", MyValue);
end

function createStructuralElements(SapModel, N_story, typ_story_H, SlabP, CB_str, CB_end, MCW_str, MCW_end, SCW_str, SCW_end, Col_str, Col_end, MB_str, MB_end, OB_str, OB_end)
    AreaObj = NET.explicitCast(SapModel.AreaObj, 'SAP2000v1.cAreaObj');
    FrameObj = NET.explicitCast(SapModel.FrameObj, 'SAP2000v1.cFrameObj');
    
    for i = 1:N_story
        % Create slab
        slabCoords = [SlabP(:,1), SlabP(:,2), repmat(i*typ_story_H, size(SlabP, 1), 1)];
        ret = AreaObj.AddByCoord(4, slabCoords(:,1), slabCoords(:,2), slabCoords(:,3), "1", "SLAB", "Slab" + num2str(i));
        
        % Create other elements (walls, columns, beams)
        createWalls(AreaObj, i, typ_story_H, MCW_str, MCW_end, SCW_str, SCW_end);
        createFrameElements(FrameObj, i, typ_story_H, CB_str, CB_end, Col_str, Col_end, MB_str, MB_end, OB_str, OB_end);
    end
end

function createWalls(AreaObj, i, typ_story_H, MCW_str, MCW_end, SCW_str, SCW_end)
    wallTypes = {'LMW', 'HMW'};
    for t = 1:2
        wallType = wallTypes{t};
        if (i/N_story < 1/6 && t == 1) || (i/N_story >= 1/6 && t == 2)
            for k = 1:size(MCW_str, 1)
                coords = [MCW_str(k,:), (i-1)*typ_story_H;
                          MCW_str(k,:), i*typ_story_H;
                          MCW_end(k,:), i*typ_story_H;
                          MCW_end(k,:), (i-1)*typ_story_H];
                ret = AreaObj.AddByCoord(4, coords(:,1), coords(:,2), coords(:,3), "1", wallType, wallType + num2str(i) + num2str(k));
            end
        end
    end
    % Similar process for SCW (secondary core walls) if needed
end

function createFrameElements(FrameObj, i, typ_story_H, CB_str, CB_end, Col_str, Col_end, MB_str, MB_end, OB_str, OB_end)
    % Create coupling beams
    createBeams(FrameObj, CB_str, CB_end, i, typ_story_H, 'CB');
    
    % Create columns
    createColumns(FrameObj, Col_str, Col_end, i, typ_story_H);
    
    % Create main beams
    createBeams(FrameObj, MB_str, MB_end, i, typ_story_H, 'MB');
    
    % Create outer beams
    createBeams(FrameObj, OB_str, OB_end, i, typ_story_H, 'OB');
end

function createBeams(FrameObj, start_coords, end_coords, i, typ_story_H, beam_type)
    for l = 1:size(start_coords, 1)
        coords = [start_coords(l,:), i*typ_story_H, end_coords(l,:), i*typ_story_H];
        if i/N_story < 1/2
            beam_prop = "L_" + beam_type;
        else
            beam_prop = "H_" + beam_type;
        end
        ret = FrameObj.AddByCoord(coords(1), coords(2), coords(3), coords(4), coords(5), coords(6), '1', beam_prop);
    end
end

function createColumns(FrameObj, start_coords, end_coords, i, typ_story_H)
    for l = 1:size(start_coords, 1)
        coords = [start_coords(l,:), (i-1)*typ_story_H, end_coords(l,:), i*typ_story_H];
        if i/N_story < 1/2
            col_prop = "L_COL";
        else
            col_prop = "H_COL";
        end
        ret = FrameObj.AddByCoord(coords(1), coords(2), coords(3), coords(4), coords(5), coords(6), "1", col_prop, "Col" + num2str((i-1)*size(start_coords,1) + l));
    end
end

function setupGroupsAndAssignments(SapModel, N_story, D_bld, B_bld, typ_story_H)
    GroupDef = NET.explicitCast(SapModel.GroupDef, 'SAP2000v1.cGroup');
    PointObj = NET.explicitCast(SapModel.PointObj, 'SAP2000v1.cPointObj');
    
    ret = GroupDef.SetGroup("Story");
    for i = 1:N_story
        groupName = "Story" + num2str(i);
        ret = GroupDef.SetGroup(groupName);
        ret = PointObj.AddCartesian(D_bld/2, B_bld/2, i*typ_story_H, groupName, groupName);
        ret = PointObj.SetGroupAssign(groupName, groupName, false, SAP2000v1.eItemType.Objects);
        ret = PointObj.SetGroupAssign(groupName, "Story", false, SAP2000v1.eItemType.Objects);
    end
    
    % Set up other groups (SLAB, WALL, MW)
    otherGroups = {"SLAB", "WALL", "MW"};
    for i = 1:length(otherGroups)
        ret = GroupDef.SetGroup(otherGroups{i});
    end
    
    % Assign objects to groups
    assignObjectsToGroups(SapModel);
end

function assignObjectsToGroups(SapModel)
    SelectObj = NET.explicitCast(SapModel.SelectObj, 'SAP2000v1.cSelect');
    AreaObj = NET.explicitCast(SapModel.AreaObj, 'SAP2000v1.cAreaObj');
    
    groupAssignments = struct('SLAB', 'SLAB', 'WALL', {'LMW', 'LSW', 'HMW', 'HSW'}, 'MW', {'LMW', 'HMW'});
    fields = fieldnames(groupAssignments);
    
    for i = 1:length(fields)
        ret = SelectObj.ClearSelection;
        if iscell(groupAssignments.(fields{i}))
            for j = 1:length(groupAssignments.(fields{i}))
                ret = SelectObj.PropertyArea(groupAssignments.(fields{i}){j});
            end
        else
            ret = SelectObj.PropertyArea(groupAssignments.(fields{i}));
        end
        ret = AreaObj.SetGroupAssign(fields{i}, fields{i}, false, SAP2000v1.eItemType.SelectedObjects);
    end
end

function setupDiaphragms(SapModel, typ_story_H)
    PointObj = NET.explicitCast(SapModel.PointObj, 'SAP2000v1.cPointObj');
    SelectObj = NET.explicitCast(SapModel.SelectObj, 'SAP2000v1.cSelect');
    
    ret = PointObj.SetSelected("ALL", true, SAP2000v1.eItemType.Group);
    [~, ~, ~, ptNamelist] = SelectObj.GetSelected;
    ret = SelectObj.ClearSelection;
    
    for i = 1:length(ptNamelist)
        [~, ~, ~, z] = PointObj.GetCoordCartesian(ptNamelist(i));
        if z ~= 0
            slabind = round(z/typ_story_H);
            slName = "Slab" + num2str(slabind);
            ret = PointObj.SetConstraint(ptNamelist(i), slName, SAP2000v1.eItemType.Objects);
            ret = PointObj.SetGroupAssign(ptNamelist(i), "Story_pt" + num2str(slabind), false, SAP2000v1.eItemType.Objects);
        else
            ret = PointObj.SetRestraint(ptNamelist(i), true(1,6));
            ret = PointObj.SetGroupAssign(ptNamelist(i), "Base", false, SAP2000v1.eItemType.Objects);
        end
    end
    
    % Set edge constraints
    AreaObj = NET.explicitCast(SapModel.AreaObj, 'SAP2000v1.cAreaObj');
    ret = AreaObj.SetEdgeConstraint("SLAB", true, SAP2000v1.eItemType.Group);
    ret = AreaObj.SetEdgeConstraint("WALL", true, SAP2000v1.eItemType.Group);
end
```

### MDOF_sap_running.m
```matlab
run('tpu_test.m');
%%

ret = SapModel.SetModelIsLocked(false());
%메쉬 변경
AreaObj = NET.explicitCast(SapModel.AreaObj,'SAP2000v1.cAreaObj');
mesh_size=1.2;
ret=AreaObj.SetAutoMesh("SLAB", 4, 1, 1,mesh_size,mesh_size, false(), false(), true(), 0, 2, false(), false(),false() ,false() , "ALL",true() ,mesh_size ,SAP2000v1.eItemType.Group);
ret=AreaObj.SetAutoMesh("WALL", 4, 1, 1,mesh_size,mesh_size, false(), false(), true(), 0, 2, false(), false(),false() ,false() , "ALL",true() ,mesh_size ,SAP2000v1.eItemType.Group);
%%
LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
ModalEigen = NET.explicitCast(LoadCases.ModalEigen,'SAP2000v1.cCaseModalEigen');

%모달 최대 개수 정의
ret=ModalEigen.SetNumberModes("MODAL",Modal_num,Modal_num-10);

%% 물성치 수정
PropMaterial = NET.explicitCast(SapModel.PropMaterial,'SAP2000v1.cPropMaterial');
ret=PropMaterial.SetWeightAndMass("4000Psi",1,Weight_concrete);
ret=PropMaterial.SetMPIsotropic("4000Psi",E_concrete,0.2,9.9*10^(-6));
%% 단면 수정
along_bool=1;
across_bool=0;
torsion_bool=0;
%%
AreaObj = NET.explicitCast(SapModel.AreaObj,'SAP2000v1.cAreaObj');
ret = AreaObj.SetEdgeConstraint("SLAB", true(), SAP2000v1.eItemType.Group);
ret = AreaObj.SetEdgeConstraint("WALL", true(), SAP2000v1.eItemType.Group);
%% Time history data 생성
cd 'E:\result_data\TPU_file'
for i=1:N_story
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');

    Windloadname ="U1" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,1,0,0,0,0);
end


for i=1:N_story
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');

    Windloadname ="U2" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,0,1,0,0,0);
end

for i=1:N_story
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');

    Windloadname ="R3" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,0,0,1,0,0);
end

cd 'C:\Users\chanho\Desktop\research\MATLAB'
%% 풍하중 load pattern 생성
for i=1:N_story
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');

    Windloadname ="U1" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,1,0,0,0,0);
end


for i=1:N_story
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');

    Windloadname ="U2" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,0,1,0,0,0);
end

for i=1:N_story
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');

    Windloadname ="R3" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,0,0,1,0,0);
end

%% 풍하중 load case 생성
relative_height=PH_N_story/N_story;
applied_WL_X=zeros(length(BX_ensenbled(:,1)),N_story,t_Numb);
applied_WL_Y=zeros(length(BX_ensenbled(:,1)),N_story,t_Numb);
applied_WL_Z=zeros(length(BX_ensenbled(:,1)),N_story,t_Numb);
for l=1:t_Numb
    lcount=0;
    %%
    if along_bool==1
        for i=1:N_story
            lcount=lcount+1;
            if B_bld/D_bld==2 && H_bld/B_bld==2
                relative_height=(PH_N_story+1)/N_story;
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end

                if assigned_PHN>=3
                    assigned_PHN=assigned_PHN-1;
                else
                    assigned_PHN=1;
                end

            elseif B_bld/D_bld==2 && H_bld/B_bld==5
                relative_height=(PH_N_story+3)/N_story;
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end

                if assigned_PHN>=7
                    assigned_PHN=assigned_PHN-3;
                elseif assigned_PHN<=6 && assigned_PHN>=5
                    assigned_PHN=3;
                elseif assigned_PHN<=4 && assigned_PHN>=3
                    assigned_PHN=2;
                elseif assigned_PHN<=2 && assigned_PHN>=1
                    assigned_PHN=1;
                end
            else
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end
            end
            funcname=num2str(l)+"th_case_"+num2str(assigned_PHN)+"th_V"+num2str(V_0)+"_B_"+num2str(B_bld)+"_A_"+anglestr+"_along";
            applied_WL_X(:,i,l)=BX_ensenbled(:,assigned_PHN,l);

            Windloadname ="U1" + strtrim(num2str(i));
            MyLoadType(lcount)="Load";
            MyLoadName(lcount)=Windloadname;
            MyFunc(lcount)=funcname;
            MySF(lcount)=1;
            MyTF(lcount)=1;
            MyAT(lcount)=0;
            MyCSys(lcount)="Global";
            MyAng(lcount)=0;
        end
    end
    %%
    if across_bool==1
        for i=1:N_story
            lcount=lcount+1;
            if B_bld/D_bld==2 && H_bld/B_bld==2
                relative_height=(PH_N_story+1)/N_story;
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end

                if assigned_PHN>=3
                    assigned_PHN=assigned_PHN-1;
                else
                    assigned_PHN=1;
                end

            elseif B_bld/D_bld==2 && H_bld/B_bld==5
                relative_height=(PH_N_story+3)/N_story;
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end

                if assigned_PHN>=7
                    assigned_PHN=assigned_PHN-3;
                elseif assigned_PHN<=6 && assigned_PHN>=5
                    assigned_PHN=3;
                elseif assigned_PHN<=4 && assigned_PHN>=3
                    assigned_PHN=2;
                elseif assigned_PHN<=2 && assigned_PHN>=1
                    assigned_PHN=1;
                end
            else
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end
            end
            funcname=num2str(l)+"th_case_"+num2str(assigned_PHN)+"th_V"+num2str(V_0)+"_B_"+num2str(B_bld)+"_A_"+anglestr+"_across";

            applied_WL_Y(:,i,l)=BY_ensenbled(:,assigned_PHN,l);
            Windloadname ="U2" + strtrim(num2str(i));
            MyLoadType(lcount)="Load";
            MyLoadName(lcount)=Windloadname;
            MyFunc(lcount)=funcname;
            MySF(lcount)=1;
            MyTF(lcount)=1;
            MyAT(lcount)=0;
            MyCSys(lcount)="Global";
            MyAng(lcount)=0;
        end
    end


    %%
    if torsion_bool==1
        for i=1:N_story
            lcount=lcount+1;
            if B_bld/D_bld==2 && H_bld/B_bld==2
                relative_height=(PH_N_story+1)/N_story;
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end

                if assigned_PHN>=3
                    assigned_PHN=assigned_PHN-1;
                else
                    assigned_PHN=1;
                end

            elseif B_bld/D_bld==2 && H_bld/B_bld==5
                relative_height=(PH_N_story+3)/N_story;
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end

                if assigned_PHN>=7
                    assigned_PHN=assigned_PHN-3;
                elseif assigned_PHN<=6 && assigned_PHN>=5
                    assigned_PHN=3;
                elseif assigned_PHN<=4 && assigned_PHN>=3
                    assigned_PHN=2;
                elseif assigned_PHN<=2 && assigned_PHN>=1
                    assigned_PHN=1;
                end
            else
                assigned_PHN=round(i*relative_height);
                if assigned_PHN==0
                    assigned_PHN=1;
                end
            end
            funcname=num2str(l)+"th_case_"+num2str(assigned_PHN)+"th_V"+num2str(V_0)+"_B_"+num2str(B_bld)+"_A_"+anglestr+"_torsion";

            applied_WL_Z(:,i,l)=BZ_ensenbled(:,assigned_PHN,l);
            Windloadname ="R3" + strtrim(num2str(i));
            MyLoadType(lcount)="Load";
            MyLoadName(lcount)=Windloadname;
            MyFunc(lcount)=funcname;
            MySF(lcount)=1;
            MyTF(lcount)=1;
            MyAT(lcount)=0;
            MyCSys(lcount)="Global";
            MyAng(lcount)=0;
        end
    end
    LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
    ModHistLinear = NET.explicitCast(LoadCases.ModHistLinear,'SAP2000v1.cCaseModalHistoryLinear');
    case_name=num2str(l)+"th_case_"+"Wind_TH";
    ret = ModHistLinear.SetCase(case_name);
    ret = ModHistLinear.SetLoads(case_name, N_story*(along_bool+across_bool+torsion_bool), MyLoadType, MyLoadName, MyFunc, MySF, MyTF, MyAT, MyCSys, MyAng);
    ret= ModHistLinear.SetTimeStep(case_name,count,dt);
    ret=ModHistLinear.SetDampConstant(case_name,damp);
    ret=ModHistLinear.SetMotionType(case_name,2);


    for i=1:N_story
        Windloadname ="U1" + strtrim(num2str(i));
        LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
        ret=LoadCases.Delete(Windloadname);

        Windloadname ="U2" + strtrim(num2str(i));
        LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
        ret=LoadCases.Delete(Windloadname);

        Windloadname ="R3" + strtrim(num2str(i));
        LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
        ret=LoadCases.Delete(Windloadname);
    end
end


%% 코드 하중 생성


GD_val_along=GD_along(B_bld,H_bld,V_H_bld,alpha,T_along,damp,N_story,asp_r,angle,betax); %GD(B,H,V,alpha,T,zd,M,N_story)
WD_val_along=WD_along(alpha,N_story,B_bld,D_bld,H_bld,V_H_bld,GD_val_along,asp_r,angle);
indices = 1:length(WD_val_along);
multiplied_values = indices .* typ_story_H .* WD_val_along;
WD_val_along_SDOF = sum(multiplied_values);
WD_val_along_SDOF=WD_val_along_SDOF;

WD_val_along_mean=WD_along(alpha,N_story,B_bld,D_bld,H_bld,V_H_bld,1,asp_r,angle);

WD_val_across=WD_across(alpha,T_across,N_story,B_bld,D_bld,H_bld,V_H_bld,damp,betay);
multiplied_values = indices .* typ_story_H .* WD_val_across;
WD_val_across_SDOF = sum(multiplied_values);
WD_val_across_SDOF=WD_val_across_SDOF;

WD_val_torsion=WD_torsion(alpha,T_torsion,N_story,B_bld,D_bld,H_bld,V_H_bld,damp,betaz);
WD_val_torsion_SDOF = sum(WD_val_torsion);

%%
Code_load(res_cyc_count,1,1:N_story)=WD_val_along(1:N_story);
Code_load(res_cyc_count,2,1:N_story)=WD_val_across(1:N_story);
Code_load(res_cyc_count,3,1:N_story)=WD_val_torsion(1:N_story);
%%
lcount=0;
for i=1:N_story
    lcount=lcount+1;
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    Windloadname = "Wind_along_code" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,WD_val_along(i),0,0,0,0);
    MyLoadType_code(lcount)="Load";
    MyLoadName_code(lcount)=Windloadname;
    MySF_code(lcount)=1;
end



%%
for i=1:N_story
    lcount=lcount+1;
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    Windloadname = "Wind_across_code" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,0,WD_val_across(i),0,0,0);
    MyLoadType_code(lcount)="Load";
    MyLoadName_code(lcount)=Windloadname;
    MySF_code(lcount)=1;
end


%%
for i=1:N_story
    lcount=lcount+1;
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    Windloadname = "Wind_torsion_code" + strtrim(num2str(i));
    slName = "Slab" + strtrim(num2str(i));
    LoadPatterns = NET.explicitCast(SapModel.LoadPatterns,'SAP2000v1.cLoadPatterns');
    ret=LoadPatterns.Add(Windloadname, SAP2000v1.eLoadPatternType.Wind);
    AutoWind = NET.explicitCast(LoadPatterns.AutoWind,'SAP2000v1.cAutoWind');
    ret = AutoWind.SetUserLoad(Windloadname,slName,0,0,WD_val_torsion(i),0,0);
    MyLoadType_code(lcount)="Load";
    MyLoadName_code(lcount)=Windloadname;
    MySF_code(lcount)=1;
end

LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
StaticLinear = NET.explicitCast(LoadCases.StaticLinear,'SAP2000v1.cCaseStaticLinear');
ret=StaticLinear.SetCase("Wind_code");
ret=StaticLinear.SetLoads("Wind_code",N_story*3,MyLoadType_code,MyLoadName_code,MySF_code);

for i=1:N_story
    Windloadname = "Wind_torsion_code" + strtrim(num2str(i));
    LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
    ret=LoadCases.Delete(Windloadname);
end

for i=1:N_story
    Windloadname = "Wind_along_code" + strtrim(num2str(i));
    LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
    ret=LoadCases.Delete(Windloadname);
end

for i=1:N_story
    Windloadname = "Wind_across_code" + strtrim(num2str(i));
    LoadCases = NET.explicitCast(SapModel.LoadCases,'SAP2000v1.cLoadCases');
    ret=LoadCases.Delete(Windloadname);
end


%%
Analyze = NET.explicitCast(SapModel.Analyze,'SAP2000v1.cAnalyze');
ret = Analyze.RunAnalysis();

```
