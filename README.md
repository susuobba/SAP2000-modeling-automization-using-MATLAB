# SAP2000-modeling-automization-using-MATLAB

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

--------------------------------------------------------------------------------------------------------------
%set the parameter or get the plan informations
----------------------------------------------------------------------------------------
if model_mode=="MDOF"
    run('MDOF_modeling.m')
else
    run('modeling.m')
end

Kzr=1.87566687625052*exp(-6.4856276*alpha)*H_bld^alpha;
V_H_bld=V_0*Kzr;
Modal_num=N_story*3;
T_along=Period_arr(res_cyc_count,1);
T_across=Period_arr(res_cyc_count,2);
T_torsion=Period_arr(res_cyc_count,3);
run('MDOF_pre_sap_running.m');
run('find_beta_MDOF.m')
run('get_period.m')

if model_mode=="MDOF"
    run('MDOF_sap_running.m');
else
    run('sap_running.m')
end
run('MDOF_correction_run.m');
run('get_total_mass.m')


```


