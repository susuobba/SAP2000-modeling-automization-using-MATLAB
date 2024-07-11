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
```

get the response
- KDS code response
```matlab
case_name="Wind_code";
for jj=1:N_story
    run('Wind_Code_data.m')

    Code_disp(res_cyc_count,1,jj)=joint_disp_U1;
    Code_disp(res_cyc_count,2,jj)=joint_disp_U2;
    Code_disp(res_cyc_count,3,jj)=joint_disp_R3;
end
Code_stress_base(res_cyc_count,1)=Vx_TH;
Code_stress_base(res_cyc_count,2)=Vy_TH;
Code_stress_base(res_cyc_count,3)=Vz_TH;

Code_moment_base(res_cyc_count,1)=abs(Mx_TH);
Code_moment_base(res_cyc_count,2)=abs(My_TH);
Code_moment_base(res_cyc_count,3)=abs(Mz_TH);

run('section_analyze_Code_design.m')
for jj=1:N_story-1

    Code_story_M2(res_cyc_count,jj)=Code_shear_M2(res_cyc_count,jj+1)-Code_shear_M2(res_cyc_count,jj);
    Code_story_M3(res_cyc_count,jj)=Code_shear_M3(res_cyc_count,jj+1)-Code_shear_M3(res_cyc_count,jj);

    Code_story_V2(res_cyc_count,jj)=Code_shear_V2(res_cyc_count,jj+1)-Code_shear_V2(res_cyc_count,jj);
    Code_story_V3(res_cyc_count,jj)=Code_shear_V3(res_cyc_count,jj+1)-Code_shear_V3(res_cyc_count,jj);
    Code_story_T(res_cyc_count,jj)=Code_shear_T(res_cyc_count,jj+1)-Code_shear_T(res_cyc_count,jj);

end
Code_story_M2(res_cyc_count,N_story)=Code_story_M2(res_cyc_count,N_story-1);
Code_story_M3(res_cyc_count,N_story)=Code_story_M3(res_cyc_count,N_story-1);
Code_story_V2(res_cyc_count,N_story)=Code_story_V2(res_cyc_count,N_story-1);
Code_story_V3(res_cyc_count,N_story)=Code_story_V3(res_cyc_count,N_story-1);
Code_story_T(res_cyc_count,N_story)=Code_story_T(res_cyc_count,N_story-1);

Code_story_M2=abs(Code_story_M2);
Code_story_M3=abs(Code_story_M3);
Code_story_V2=abs(Code_story_V2);
Code_story_V3=abs(Code_story_V3);
Code_story_T=abs(Code_story_T);
```

- time history response
```matlab
run('section_analyze_Code_design.m')


%%
clear Mx_TH My_TH Mz_TH Vx_TH Vy_TH Vz_TH


ensenbeled_base_Mx=[];
ensenbeled_base_My=[];
ensenbeled_base_Mz=[];

ensenbeled_joint_disp_along=[];
ensenbeled_joint_disp_across=[];
ensenbeled_joint_disp_torsion=[];

ensenbeled_joint_disp=[];
temp_dx=[];
temp_dy=[];
temp_dz=[];
temp_My=[];
Mx_TH_data=[];
My_TH_data=[];
Mz_TH_data=[];

Mx_peak_MDOF(res_cyc_count,:)=zeros(1,2);
My_peak_MDOF(res_cyc_count,:)=zeros(1,2);
Mz_peak_MDOF(res_cyc_count,:)=zeros(1,2);

Vx_peak_MDOF(res_cyc_count,:)=zeros(1,2);
Vy_peak_MDOF(res_cyc_count,:)=zeros(1,2);
Base_moment_std(res_cyc_count,:)=[0 0 0];
for jj=1:N_story
    joint_disp=[];
    clear joint_disp_along joint_disp_across joint_disp_torsion
    for l=1:t_Numb
        run('Wind_TH_data.m');

        ensenbeled_joint_disp_along(:,l)=joint_disp_along;
        ensenbeled_joint_disp_across(:,l)=joint_disp_across;
        ensenbeled_joint_disp_torsion(:,l)=joint_disp_torsion;
        if jj==1
            run('Wind_TH_base.m');
            ensenbeled_base_My(:,l)=My_TH;
            ensenbeled_base_Mx(:,l)=Mx_TH;
            ensenbeled_base_Mz(:,l)=Mz_TH;
            ensenbeled_base_Vy(:,l)=Vy_TH;
            ensenbeled_base_Vx(:,l)=Vx_TH;

            % Mx_TH_res(res_cyc_count,1:length(My_TH))=Mx_TH;
            % My_TH_res(res_cyc_count,1:length(My_TH))=My_TH;
            % Mz_TH_res(res_cyc_count,1:length(My_TH))=Mz_TH;
            % Vy_TH_res(res_cyc_count,1:length(My_TH))=Vy_TH;
            % Vx_TH_res(res_cyc_count,1:length(My_TH))=Vx_TH;


            joint_disp=Mx_TH;
            Mx_TH_data(:,l)=joint_disp;

            gng=Across_pf(alpha,T_across,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
            Mx_peak_MDOF(res_cyc_count,1)=Mx_peak_MDOF(res_cyc_count,1)+(std(joint_disp)*gng)/t_Numb;
            Base_moment_std(res_cyc_count,2)=Base_moment_std(res_cyc_count,2)+std(joint_disp)/t_Numb;

            joint_disp=My_TH;
            My_TH_data(:,l)=joint_disp;
            temp_My(:,l)=My_TH;



            gng=Along_pf(B_bld,H_bld,V_H_bld,alpha,T_along,damp,N_story,asp_r,angle);
            My_peak_MDOF(res_cyc_count,1)=My_peak_MDOF(res_cyc_count,1)+(mean(joint_disp)+std(joint_disp)*gng)/t_Numb;
            My_peak_MDOF(res_cyc_count,2)=My_peak_MDOF(res_cyc_count,2)+(mean(joint_disp)-std(joint_disp)*gng)/t_Numb;
            Base_moment_std(res_cyc_count,1)=Base_moment_std(res_cyc_count,1)+std(joint_disp)/t_Numb;
            Base_moment_mean(res_cyc_count,1)=Base_moment_mean(res_cyc_count,1)+mean(joint_disp)/t_Numb;


            joint_disp=Mz_TH;
            Mz_TH_data(:,l)=joint_disp;

            gng=Torsion_pf(alpha,T_torsion,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
            Mz_peak_MDOF(res_cyc_count,1)=Mz_peak_MDOF(res_cyc_count,1)+(std(joint_disp)*gng)/t_Numb;
            Base_moment_std(res_cyc_count,3)=Base_moment_std(res_cyc_count,3)+std(joint_disp)/t_Numb;

            joint_disp=Vx_TH;

            gng=Along_pf(B_bld,H_bld,V_H_bld,alpha,T_along,damp,N_story,asp_r,angle);
            Vx_peak_MDOF(res_cyc_count,1)=Vx_peak_MDOF(res_cyc_count,1)+(mean(joint_disp)+std(joint_disp)*gng)/t_Numb;
            Vx_peak_MDOF(res_cyc_count,2)=Vx_peak_MDOF(res_cyc_count,2)+(mean(joint_disp)-std(joint_disp)*gng)/t_Numb;

            joint_disp=Vy_TH;

            gng=Across_pf(alpha,T_across,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
            Vy_peak_MDOF(res_cyc_count,1)=Vy_peak_MDOF(res_cyc_count,1)+(std(joint_disp)*gng)/t_Numb;


        end
        % Along

        % joint_disp_along_res(res_cyc_count,jj,1:length(joint_disp_along))=joint_disp_along;
        % joint_disp_across_res(res_cyc_count,jj,1:length(joint_disp_along))=joint_disp_across;
        % joint_disp_torsion_res(res_cyc_count,jj,1:length(joint_disp_along))=joint_disp_torsion;


        joint_disp=joint_disp_along;
        TPU_mean_disp_arr_along_MDOF(res_cyc_count,jj)=mean(joint_disp);


        if jj==1 && l==1
            TPU_disp_arr_along_MDOF=zeros(length(joint_disp),l);
            TPU_disp_arr_across_MDOF=zeros(length(joint_disp),l);
            TPU_disp_arr_torsion_MDOF=zeros(length(joint_disp),l);

            TPU_peak_disp_arr_along_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
            TPU_peak_disp_arr_across_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
            TPU_peak_disp_arr_torsion_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
        end
        TPU_disp_arr_along_MDOF(:,l)=joint_disp;

        gng=Along_pf(B_bld,H_bld,V_H_bld,alpha,T_along,damp,N_story,asp_r,angle);
        peak_factor(res_cyc_count,1)=gng;
        pf=mean(joint_disp)+std(joint_disp)*gng;
        TPU_peak_disp_arr_along_MDOF(res_cyc_count,jj,1)=TPU_peak_disp_arr_along_MDOF(res_cyc_count,jj,1)+pf/t_Numb;

        pf=mean(joint_disp)-std(joint_disp)*gng;
        TPU_peak_disp_arr_along_MDOF(res_cyc_count,jj,2)=TPU_peak_disp_arr_along_MDOF(res_cyc_count,jj,2)+pf/t_Numb;
        if jj==N_story
            gng=Along_pf(B_bld,H_bld,V_H_bld,alpha,T_along,damp,N_story,asp_r,angle);
            pf=mean(joint_disp)+std(joint_disp)*gng;
            maxdp(l)=pf;


            temp_dx(:,l)=joint_disp;
        end

        % Across
        joint_disp=joint_disp_across;
        TPU_mean_disp_arr_across_MDOF(res_cyc_count,jj)=mean(joint_disp);
        TH_data_across(res_cyc_count,jj,1)=TH_data_across(res_cyc_count,jj,1)+mean(joint_disp)/t_Numb;
        TH_data_across(res_cyc_count,jj,2)=TH_data_across(res_cyc_count,jj,2)+std(joint_disp)/t_Numb;

        TPU_disp_arr_across_MDOF(:,l)=joint_disp;
        gng=Across_pf(alpha,T_across,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
        peak_factor(res_cyc_count,2)=gng;
        pf=std(joint_disp)*gng;
        TPU_peak_disp_arr_across_MDOF(res_cyc_count,jj,1)=TPU_peak_disp_arr_across_MDOF(res_cyc_count,jj,1)+pf/t_Numb;
        if jj==N_story
            gng=Across_pf(alpha,T_across,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
            pf=mean(joint_disp)+std(joint_disp)*gng;
            temp_dy(:,l)=joint_disp;
        end

        %Torsion
        joint_disp=joint_disp_torsion;
        TPU_mean_disp_arr_torsion_MDOF(res_cyc_count,jj)=mean(joint_disp);
        TH_data_torsion(res_cyc_count,jj,1)=TH_data_torsion(res_cyc_count,jj,1)+mean(joint_disp)/t_Numb;
        TH_data_torsion(res_cyc_count,jj,2)=TH_data_torsion(res_cyc_count,jj,2)+std(joint_disp)/t_Numb;

        TPU_disp_arr_torsion_MDOF(:,l)=joint_disp;
        gng=Torsion_pf(alpha,T_torsion,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
        peak_factor(res_cyc_count,3)=gng;
        pf=std(joint_disp)*gng;
        TPU_peak_disp_arr_torsion_MDOF(res_cyc_count,jj,1)=TPU_peak_disp_arr_torsion_MDOF(res_cyc_count,jj,1)+pf/t_Numb;


        if jj==N_story
            gng=Torsion_pf(alpha,T_across,N_story,B_bld,D_bld,H_bld,V_H_bld,damp);
            pf=mean(joint_disp)+std(joint_disp)*gng;
            temp_dz(:,l)=joint_disp;
        end

    end
end

clearvars peak_disp_Mx mean_joint_disp_Mx peak_factor_Mx std_joint_disp_Mx
clearvars peak_disp_My mean_joint_disp_My peak_factor_My std_joint_disp_My
clearvars peak_disp_Mz mean_joint_disp_Mz peak_factor_Mz std_joint_disp_Mz
clearvars peak_disp_Vx mean_joint_disp_Vx peak_factor_Vx std_joint_disp_Vx
clearvars peak_disp_Vy mean_joint_disp_Vy peak_factor_Vy std_joint_disp_Vy




joint_disp=[];
clear joint_shear_along joint_shear_across TPU_shear_V2 TPU_shear_V3 TPU_shear_M2 TPU_shear_M3 TPU_shear_T



for l=1:t_Numb
    run('get_max_overturnM2_time.m')
    run('section_analyze_TPU_design.m')
    for jj=1:N_story

        joint_shear_along(1:floor(time_end_pt/dt)-ceil(time_str_pt/dt)+1)=(TPU_shear_V2(jj,ceil(time_str_pt/dt):floor(time_end_pt/dt)));
        joint_shear_across(1:floor(time_end_pt/dt)-ceil(time_str_pt/dt)+1)=(TPU_shear_V3(jj,ceil(time_str_pt/dt):floor(time_end_pt/dt)));
        joint_moment_x=[];
        joint_moment_y=[];
        joint_moment_z=[];

        joint_moment_x(1:floor(time_end_pt/dt)-ceil(time_str_pt/dt)+1)=(TPU_shear_M2(jj,ceil(time_str_pt/dt):floor(time_end_pt/dt)));
        joint_moment_y(1:floor(time_end_pt/dt)-ceil(time_str_pt/dt)+1)=(TPU_shear_M3(jj,ceil(time_str_pt/dt):floor(time_end_pt/dt)));
        joint_moment_z(1:floor(time_end_pt/dt)-ceil(time_str_pt/dt)+1)=(TPU_shear_T(jj,ceil(time_str_pt/dt):floor(time_end_pt/dt)));


        if jj==1 && l==1
            TPU_shear_arr_along_MDOF=zeros(N_story,length(joint_disp));
            TPU_shear_arr_across_MDOF=zeros(N_story,length(joint_disp));

            TPU_peak_shear_arr_along_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
            TPU_peak_shear_arr_across_MDOF(res_cyc_count,:,:)=zeros(1,100,2);

            TPU_moment_arr_x_MDOF=zeros(N_story,length(joint_disp));
            TPU_moment_arr_y_MDOF=zeros(N_story,length(joint_disp));
            TPU_moment_arr_z_MDOF=zeros(N_story,length(joint_disp));

            TPU_peak_moment_arr_x_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
            TPU_peak_moment_arr_y_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
            TPU_peak_moment_arr_z_MDOF(res_cyc_count,:,:)=zeros(1,100,2);
        end
        % Along
        joint_disp=joint_shear_along;
        gng=peak_factor(res_cyc_count,1);
        pf=mean(joint_disp)+std(joint_disp)*gng;
        TPU_peak_shear_arr_along_MDOF(res_cyc_count,jj,1)=TPU_peak_shear_arr_along_MDOF(res_cyc_count,jj,1)+pf/t_Numb;

        pf=mean(joint_disp)-std(joint_shear_along)*gng;
        TPU_peak_shear_arr_along_MDOF(res_cyc_count,jj,2)=TPU_peak_shear_arr_along_MDOF(res_cyc_count,jj,2)+pf/t_Numb;


        % Across
        joint_disp=joint_shear_across;
        gng=peak_factor(res_cyc_count,2);
        pf=mean(joint_disp)+std(joint_disp)*gng;
        TPU_peak_shear_arr_across_MDOF(res_cyc_count,jj,1)=TPU_peak_shear_arr_across_MDOF(res_cyc_count,jj,1)+pf/t_Numb;

        pf=mean(joint_disp)-std(joint_disp)*gng;
        TPU_peak_shear_arr_across_MDOF(res_cyc_count,jj,2)=TPU_peak_shear_arr_across_MDOF(res_cyc_count,jj,2)+pf/t_Numb;


        % Mx
        joint_disp=joint_moment_x;
        gng=peak_factor(res_cyc_count,2);
        pf=mean(joint_disp)+std(joint_disp)*gng;
        TPU_peak_moment_arr_x_MDOF(res_cyc_count,jj,1)=TPU_peak_moment_arr_x_MDOF(res_cyc_count,jj,1)+pf/t_Numb;

        pf=mean(joint_disp)-std(joint_disp)*gng;
        TPU_peak_moment_arr_x_MDOF(res_cyc_count,jj,2)=TPU_peak_moment_arr_x_MDOF(res_cyc_count,jj,2)+pf/t_Numb;

        % My
        joint_disp=joint_moment_y;
        gng=peak_factor(res_cyc_count,1);
        pf=mean(joint_disp)+std(joint_disp)*gng;
        TPU_peak_moment_arr_y_MDOF(res_cyc_count,jj,1)=TPU_peak_moment_arr_y_MDOF(res_cyc_count,jj,1)+pf/t_Numb;

        pf=mean(joint_disp)-std(joint_disp)*gng;
        TPU_peak_moment_arr_y_MDOF(res_cyc_count,jj,2)=TPU_peak_moment_arr_y_MDOF(res_cyc_count,jj,2)+pf/t_Numb;

        % Mz
        joint_disp=joint_moment_z;
        gng=peak_factor(res_cyc_count,3);
        pf=mean(joint_disp)+std(joint_disp)*gng;
        TPU_peak_moment_arr_z_MDOF(res_cyc_count,jj,1)=TPU_peak_moment_arr_z_MDOF(res_cyc_count,jj,1)+pf/t_Numb;

        pf=mean(joint_disp)-std(joint_disp)*gng;
        TPU_peak_moment_arr_z_MDOF(res_cyc_count,jj,2)=TPU_peak_moment_arr_z_MDOF(res_cyc_count,jj,2)+pf/t_Numb;
    end
end

shear_story_x=zeros(1,N_story);
shear_story_y=zeros(1,N_story);
shear_story_z=zeros(1,N_story);
shear_story_z_temp=max(squeeze(abs(TPU_peak_moment_arr_z_MDOF(res_cyc_count,1:N_story,:))),[],2);

for i=2:N_story

    shear_story_x(i-1)=abs(TPU_peak_shear_arr_along_MDOF(res_cyc_count,i-1,2))-abs(TPU_peak_shear_arr_along_MDOF(res_cyc_count,i,2));
    shear_story_y(i-1)=abs(TPU_peak_shear_arr_across_MDOF(res_cyc_count,i-1,1))-abs(TPU_peak_shear_arr_across_MDOF(res_cyc_count,i,1));
    shear_story_z(i-1)=abs(shear_story_z_temp(i-1))-abs(shear_story_z_temp(i));

end

shear_story_x(N_story)=shear_story_x(N_story-1);
shear_story_y(N_story)=shear_story_y(N_story-1);
shear_story_z(N_story)=shear_story_z(N_story-1);
shear_story_x=abs(shear_story_x);
shear_story_y=abs(shear_story_y);
```

