from os import listdir

def fetch_all_subjects(data_dir):


	subject_ids = listdir(data_dir)

	subject_list = []
	for iSubj in range(len(subject_ids)):
		if "sub-" in subject_ids[iSubj]:
			subject_list.append(subject_ids[iSubj])

	subject_list = sorted(subject_list)
	print(subject_list)
	
	return subject_list
