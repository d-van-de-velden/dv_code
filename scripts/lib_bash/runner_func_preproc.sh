#! /bin/bash
which bash
echo "~~~##########    fMRI preprocessing pipeline (offline)    ##########~~~"
echo "                              [START]                                  "
while getopts ":s:e:f:u:l:d:t:a:o:" option; do
  case $option in
    s) subject_id="$OPTARG" ;;
    e) session="$OPTARG" ;;
    f) file_name="$OPTARG" ;;
    u) func_dir="$OPTARG" ;;
    l) fmap_dir="$OPTARG" ;;
    d) anat_dir="$OPTARG" ;;
    t) totReTime="$OPTARG" ;;
    a) anat_img="$OPTARG" ;;  # Here it recommended to use the nu.mgz from FREESURFEr recon-all
    o) output_dir="$OPTARG" ;;  
    *)
      echo "Usage: $0 [-s subject_id] [-e session] [-f file_name] [-u func_dir] [-l fmap_dir] [-d anat_dir] [-t totReTime] [-a struc_img] [-o output_dir]"
      exit 1
      ;;
  esac
done
#######################################################################
file_name0=$(basename "$file_name")
f0="${file_name0%.*}"
echo "$f0"

echo "Cding in functional data directory: $func_dir"
cd $func_dir
echo "Copying functional data to derivatives directory: $func_dir"
cp $file_name "${output_dir}"

echo "Original functional file:"
echo "$file_name"
echo "Temporary functional file:"
file_name="${output_dir}/${f0}"
echo "$file_name"


echo ""
echo "#########################################################"
echo ""
echo "Starting to apply FSL - SliceTime Correction"

filename_=$(basename "$file_name")
f="${filename_%.*}"
echo "$f"
file_name1="${output_dir}/${f}_st.nii" 

slicetimer -i $file_name -r 2 -o $file_name1

echo ""
echo "#########################################################"
echo ""

echo "Starting to apply FSL - MCFLIRT (Movement)"

echo "Applying MCFLIRT on functional data: $file_name1"
mcflirt -in $file_name1 -plots

echo ""
echo "#########################################################"
echo ""

echo "Starting to apply FSL - topup (Distortion correction)"
echo "Cding in phase encoding data directory: $fmap_dir"
cd $fmap_dir

echo "Applying TopUp on phase encoding map: ${subject_id}_${session}_dir-PA_epi.nii.gz"
fslroi "${subject_id}_${session}_dir-PA_epi.nii.gz" PA 0 1

echo "Applying TopUp on phase encoding map: ${subject_id}_${session}_dir-AP_epi.nii.gz"
fslroi "${subject_id}_${session}_dir-AP_epi.nii.gz" AP 0 1

echo "Merge phase encoding maps: AP & PA"
fslmerge -t AP-PA AP.nii.gz PA.nii.gz


# the readout time in seconds for the parameter file will be: For singleband images :
# readout = (echospacing * (acquisitionmatrix[0] * (percentsampling/100))) / 1e6 
# For mutli-band images the readout calculation is more complex:
# readout = ( ( ceil ((1/Round_factor) * AcquisitionMatrixPE / Asset_R_factor ) * Round_factor) - 1 ) * EchoSpacing * 0.000001
echo "Setting up aquisition parameter map: $file_name1"
rm acq_param.txt -f
cat << EOF >> acq_param.txt
0 1 0 $totReTime
0 -1 0 $totReTime
EOF

# Get distrortion maps
echo "Make distortion maps from phase encoding map:"
topup --imain=AP-PA.nii.gz --datain=acq_param.txt --config=b02b0.cnf --out=AP_PA_topup

echo ""
echo "#########################################################"
echo ""

# Apply distortion maps to functional data
# Might take a while
file_name2=$(basename "$file_name1")
f2="${file_name2%.*}"
echo "$f2"
file_name2b="${output_dir}/${f2}_mcf.nii.gz"
echo "${output_dir}/${file_name2b}"
file_name3="${output_dir}/${f2}_mcf_topUP.nii.gz"
echo "${output_dir}/${file_name3}"

echo "Applying top on functional data: $file_name2b"
applytopup --imain=$file_name2b --inindex=1 --datain=acq_param.txt --topup=AP_PA_topup --method=jac --out=$file_name3

echo ""
echo "#########################################################"
echo ""

file_name4=$(basename "$anat_img")
f3="${file_name4%.*}"
file_name_anat=$(basename "$f3")
f_anat="${file_name_anat%.*}"
anat_in_brain="${anat_dir}/${f_anat}_brain.nii"
echo ""
echo "Apply Brain Extraction: functional --> anatomical "
echo "IN  -> $anat_img"
echo "OUT -> $anat_in_brain"
bet "${anat_img}" "${anat_in_brain}"

#########################################################################################
echo ""
echo "Applying coregistration of functional data to FreeSurfer anatomical data: $file_name3"
#export SUBJECTS_DIR=/data/p_02825/cocoa/data/FREESURFER/
#bbregister --s ${subject_id} --mov {$file_name3} --init-fsl --bold --fslmat flirt.mtx --reg register.dat --12
# https://surfer.nmr.mgh.harvard.edu/fswiki/bbregister
file_name_MNI="${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz"
file_name_MNI_func="${FSLDIR}/data/standard/MNI152_T1_2mm.nii.gz"
file_name_MNI_brain="${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz"
file_name_MNI_brainmask="${FSLDIR}/data/standard/MNI152_T1_1mm_brain_mask_dil.nii"
file_name_MNI_func_brainmask="${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask_dil.nii"

echo ""
echo "Spatially normalizing to MNI152 template"
echo $file_name_MNI
echo ""

file_name_anat=$(basename "$f3")
f_anat="${file_name_anat%.*}"
anat_out_brain="${anat_dir}/${f_anat}_brain_w.nii"

anat_in="${anat_dir}/${f_anat}.nii"
anat_out="${anat_dir}/${f_anat}_w.nii"

anat_warp="${anat_dir}/${f_anat}_warp_anat2mni.nii.gz"
aff_guess="${anat_dir}/${f_anat}_aff_guess.mat"

file_name2=$(basename "$file_name3")
f_func="${file_name2%.*}"
file_name_func=$(basename "$f_func")
f_func="${file_name_func%.*}"

file_func_out1="${output_dir}/${f_func}_r.nii.gz"
file_func_out_matrix="${output_dir}/${f_func}_func2struct.mat"


echo ""
echo "Apply Coregistration: functional --> anatomical "
echo "IN  -> $file_name3"
echo "OUT -> $file_func_out1"
echo "OUT -> $file_func_out_matrix"
flirt -ref $anat_in_brain -in $file_name3 -out $file_func_out1 -omat $file_func_out_matrix -dof 6


file_name_func=$(basename "$file_func_out1")
f_func="${file_name_func%.*}"
echo "$f_func"
file_func_out="${output_dir}/${f_func}"
file_name_func=$(basename "$f_func")
f_func="${file_name_func%.*}"
echo "$f_func"
file_func_out="${output_dir}/${f_func}_w.nii.gz"

# Apply spatial normalization
# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide#Deformation_model
echo ""
echo "Applying non-linear spatial normalization to anatomical data and use same transformation matrix for functional data:"


echo "$anat_warp"

if test ! -f $anat_warp; then
  echo ""
  echo "Using FLIRT to get an initial linear transformation matrix for further advanced spatial normalization..."
  echo "Using betted data..."
  echo "IN  -> $anat_in_brain"
  echo "OUT -> $anat_out_brain"
  flirt -ref $file_name_MNI_brain -in $anat_in_brain -out $anat_out_brain  -omat $aff_guess
  echo ""
  echo "Using FNIRT with the initial linear transformation matrix to find and apply a non-linear transformation vector field for spatial normalization..."
  echo "Using betted data..."
  echo "IN  -> $anat_in"
  echo "OUT -> $anat_out"
  fnirt --ref=$file_name_MNI --in=$anat_in --iout=$anat_out --aff=$aff_guess --cout=$anat_warp --refmask=$file_name_MNI_brainmask --verbose
  echo ""
else
  echo "Already did FLIRT + FNIRT to get a non-linear transformation vector field for spatial normalization..."
fi

# Apply warp field to functional data
echo ""
echo "Apply warp field to functional data"
echo "IN  -> $file_name3"
echo "OUT -> $file_func_out"
applywarp -i $file_name3 -o $file_func_out -r $file_name_MNI_func -w $anat_warp -m $file_name_MNI_func_brainmask --premat=func2struct.mat

# Apply warp field to c1 data
echo ""
echo "Apply warp field to c1 data"
in_c1="${anat_dir}/${f_anat}_c1.nii"
out_c1="${anat_dir}/${f_anat}_c1_w.nii"
echo "IN  -> $in_c1"
echo "OUT -> $out_c1"
applywarp -r $file_name_MNI -i $in_c1 -o $out_c1 -w $anat_warp -m $file_name_MNI_brainmask

# Apply warp field to c2 data
echo ""
echo "\nApply warp field to c2 data"
in_c2="${anat_dir}/${f_anat}_c2.nii"
out_c2="${anat_dir}/${f_anat}_c2_w.nii"
echo "IN  -> $in_c2"
echo "OUT -> $in_c2"
applywarp -r $file_name_MNI -i $in_c2 -o $out_c2 -w $anat_warp -m $file_name_MNI_brainmask

# Apply warp field to brain data
echo ""
echo "Apply warp field to brain data"
in_brain="${anat_dir}/${f_anat}_brain.nii"
out_brain="${anat_dir}/${f_anat}_brain_w.nii"
echo "IN  -> $in_brain"
echo "OUT -> $out_brain"
applywarp -r $file_name_MNI -i $in_brain -o $out_brain -w $anat_warp -m $file_name_MNI_brainmask

# Apply warp field to aseg+aparc data
echo ""
echo "Apply warp field to aseg+aparc data"
red_brain1="${anat_dir}/${f_anat}"
red_brain="${red_brain1::-3}"
in_aseg="${red_brain}_aparc+aseg.nii"
out_aseg="${red_brain}_aparc+aseg_w.nii"
echo "IN  -> $in_aseg"
echo "OUT -> $out_aseg"
applywarp -r $file_name_MNI -i $in_aseg -o $out_aseg -w $anat_warp -m $file_name_MNI_brainmask


echo ""
echo "                                [END]                                  "
echo "~~~##########    fMRI preprocessing pipeline (offline)    ##########~~~"
kill $BASHPID
exit 0
q
quit
