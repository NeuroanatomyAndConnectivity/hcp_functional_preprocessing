#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO DO REGRESS OUT NUISANCE COVARIATES FROM RESTING_STATE SCAN
## nuisance covariates are: global signal, white matter (WM, WHITE), CSF, and
## 6 motion parameters obtained during motion correction step (see 2_funcpreproc.sh)
##
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################

## subject
subject=$1
## analysisdirectory/subject
dir=$2
## resting-state filename (no extension)
rest=$3
## TR
TR=$4
## number of timepoints in the resting-state scan
n_vols=$5
## full path to template nuisance feat .fsf file; e.g. /full/path/to/template.fsf
nuisance_template=$6
scan=$7
hp=0.005
lp=0.1

## directory setup
anat_dir=${dir}/mprage_1
func_dir=${dir}
reg_dir=${dir}/reg
segment_dir=${dir}/segment
nuisance_dir=${func_dir}/nuisance
rest=${subject}_rfMRI_${scan}_RAI


##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################



echo --------------------------------------------
echo !!!! RUNNING NUISANCE SIGNAL REGRESSION !!!!
echo --------------------------------------------

cd ${func_dir}

echo ${func_dir}
rm -rf ${nuisance_dir}
## 1. make nuisance directory
mkdir -p ${nuisance_dir}

##11.Create Mask
echo "Generating mask of preprocessed data for ${subject}"
fslmaths ${rest}.nii.gz -Tmin -bin ${rest}_bin.nii.gz -odt char


# 2. Seperate motion parameters into seperate files
echo "Splitting up ${subject} motion parameters"
awk '{print $1}' Movement_Regressors.txt > ${nuisance_dir}/mc1.1D
awk '{print $2}' Movement_Regressors.txt > ${nuisance_dir}/mc2.1D
awk '{print $3}' Movement_Regressors.txt > ${nuisance_dir}/mc3.1D
awk '{print $4}' Movement_Regressors.txt > ${nuisance_dir}/mc4.1D
awk '{print $5}' Movement_Regressors.txt > ${nuisance_dir}/mc5.1D
awk '{print $6}' Movement_Regressors.txt > ${nuisance_dir}/mc6.1D

# Extract signal for global, csf, and wm
## 3. Global
echo "Extracting global signal for ${subject}"
#3dmaskave -mask ${segment_dir}/global_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/global.1D

## 4. csf
echo "Extracting signal from csf for ${subject}"
3dmaskave -mask /scr/lessing3/cb_method/segment/${subject}_T1w_acpc_dc_restore1_clone_transform_mask_cp_strip_norm_mems_csf_binmask_def_binmask.nii.gz -quiet ${func_dir}/${rest}.nii.gz > ${nuisance_dir}/csf.1D

## 5. wm
echo "Extracting signal from white matter for ${subject}"
3dmaskave -mask /scr/lessing3/cb_method/segment/ero_${subject}_T1w_acpc_dc_restore1_clone_transform_mask_cp_strip_norm_mems_wm_binmask_def_binmask.nii.gz -quiet ${func_dir}/${rest}.nii.gz > ${nuisance_dir}/wm.1D

## 6. Generate mat file (for use later)
## create fsf file
echo "Modifying model file"
sed -e s:nuisance_dir:"${nuisance_dir}":g <${nuisance_template} >${nuisance_dir}/temp1
sed -e s:nuisance_model_outputdir:"${nuisance_dir}/residuals.feat":g <${nuisance_dir}/temp1 >${nuisance_dir}/temp2
sed -e s:nuisance_model_TR:"${TR}":g <${nuisance_dir}/temp2 >${nuisance_dir}/temp3
sed -e s:nuisance_model_numTRs:"${n_vols}":g <${nuisance_dir}/temp3 >${nuisance_dir}/temp4
sed -e s:nuisance_model_input_data:"${func_dir}/${rest}_pp.nii.gz":g <${nuisance_dir}/temp4 >${nuisance_dir}/nuisance.fsf 

rm ${nuisance_dir}/temp*

echo "Running feat model"
feat_model ${nuisance_dir}/nuisance

minVal=`3dBrickStat -min -mask ${func_dir}/${rest}_bin.nii.gz ${func_dir}/${rest}.nii.gz`

## 7. Get residuals
echo "Running film to get residuals"
film_gls -rn ${nuisance_dir}/stats -noest -sa -ms 5 ${func_dir}/${rest}.nii.gz ${nuisance_dir}/nuisance.mat ${minVal}

## 8. Demeaning residuals and ADDING 100
3dTstat -mean -prefix ${nuisance_dir}/stats/res4d_mean.nii.gz ${nuisance_dir}/stats/res4d.nii.gz
3dcalc -a ${nuisance_dir}/stats/res4d.nii.gz -b ${nuisance_dir}/stats/res4d_mean.nii.gz -expr '(a-b)+100' -prefix ${func_dir}/${rest}_res.nii.gz

##10.Detrending
echo "Removing linear trends for ${subject}"
3dTstat -mean -prefix ${rest}_res_mean.nii.gz ${rest}_res.nii.gz
3dDetrend -polort 1 -prefix ${rest}_dt.nii.gz ${rest}_res.nii.gz
3dcalc -a ${rest}_res_mean.nii.gz -b ${rest}_dt.nii.gz -expr 'a+b' -prefix ${rest}_pp.nii.gz

##9. Temporal filtering
echo "Band-pass filtering ${subject}"
3dFourier -lowpass ${lp} -highpass ${hp} -retrend -prefix ${rest}_proc.nii.gz ${rest}_pp.nii.gz

rm -f ${rest}_res_mean.nii.gz ${rest}_dt.nii.gz ${rest}_res.nii.gz ${rest}_pp.nii.gz
