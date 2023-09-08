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
echo "$f3"
bet "${anat_img}" "${f3}_brain.nii"


echo "Applying coregistration of functional data to FreeSurfer anatomical data: $file_name2b"
bbregister --s ${subject_id} --mov {$file_name3} --init-fsl --bold --fslmat flirt.mtx --reg register.dat --12
# https://surfer.nmr.mgh.harvard.edu/fswiki/bbregister


# Apply spatial normalization
# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FNIRT/UserGuide#Deformation_model
echo "Applying non-linear spatial normalization to anatomical data and use same transformation matrix for functional data:"
file_name4=$(basename "$anat_img")
f3="${file_name4%.*}"
file_name_MNI="/afs/.cbs.mpg.de/software/fsl/6.0.3/debian-bullseye-amd64/data/standard/MNI152_T1_1mm.nii.gz"
file_name_MNI_brain="/afs/.cbs.mpg.de/software/fsl/6.0.3/debian-bullseye-amd64/data/standard/MNI152_T1_1mm_brain.nii.gz"
file_name_MNI_brainmask="/afs/.cbs.mpg.de/software/fsl/6.0.3/debian-bullseye-amd64/data/standard/MNI152_T1_1mm_brain_mask_dil.nii"

anat_in_brain="${anat_dir}${f3}_brain.nii"
anat_out_brain="${anat_dir}${f3}_brain_w.nii"

anat_in="${anat_dir}${f3}.nii"
anat_out="${anat_dir}${f3}_w.nii"

anat_warp="${anat_dir}${f3}_warp_anat2mni.nii"
aff_guess="${anat_dir}${f3}_aff_guess.mat"

echo $file_name3
echo "$file_name_MNI"
echo "$anat_in"
echo "$anat_out"
echo "$anat_warp"
echo "$aff_guess"

if [ ! -f $anat_warp ]
then
  flirt -ref $file_name_MNI_brain -in $anat_in_brain -out $anat_out_brain  -omat $aff_guess
  fnirt --ref=$file_name_MNI --in=$anat_in --iout=$anat_out --aff=$aff_guess --cout=$anat_warp --refmask=$file_name_MNI_brainmask --verbose
else
  # Apply warp field to functional data 
  applywarp -i $file_name3 -o $file_name5 -r $file_name_MNI  -w $anat_warp -m $file_name_MNI_brainmask

  # Apply warp field to c1 data 
  in_c1="${anat_dir}${f3}_c1.nii"
  out_c1="${anat_dir}${f3}_c1_w.nii"
  applywarp -r $file_name_MNI -i $in_c1 -o $out_c1 -w $anat_warp -m $file_name_MNI_brainmask

  # Apply warp field to c2 data 
  in_c2="${anat_dir}${f3}_c2.nii"
  out_c2="${anat_dir}${f3}_c2_w.nii"
  applywarp -r $file_name_MNI -i $in_c2 -o $out_c2 -w $anat_warp -m $file_name_MNI_brainmask

  # Apply warp field to brain data 
  in_brain="${anat_dir}${f3}_brain.nii"
  out_brain="${anat_dir}${f3}_brain_w.nii"
  applywarp -r $file_name_MNI -i $in_brain -o $out_brain -w $anat_warp -m $file_name_MNI_brainmask

  # Apply warp field to aseg+aparc data 
  red_brain1="${anat_dir}${f3}"
  red_brain="${red_brain1::-3}"
  in_aseg="${red_brain}_aparc+aseg.nii"
  out_aseg="${red_brain}_aparc+aseg_w.nii"
  applywarp -r $file_name_MNI -i $in_aseg -o $out_aseg -w $anat_warp -m $file_name_MNI_brainmask

fi

echo "                                [END]                                  "
echo "~~~##########    fMRI preprocessing pipeline (offline)    ##########~~~"
kill $BASHPID
exit 0
q
quit