import os

eos_prefix = "root://cmsxrootd-site.fnal.gov/"
eos_suffix = "\");"
shared_final = "/eos/uscms"#/store/user/dokstolp/darkpho/Gen/df_1/"
use_path = {'1en0':["/store/group/lpcgg/dokstolp/RECO/df1en0/161222_233804/0000/"],'1en1':["/store/group/lpcgg/dokstolp/RECO/df1en1/161222_222822/0000/"],'1en2':["/store/group/lpcgg/dokstolp/RECO/df1en2/161222_222955/0000/"],'1en3':["/store/group/lpcgg/dokstolp/RECO/df1en3/161222_223119/0000/"]}
for dfs in use_path:
	number_files = -1
	for path in use_path[dfs]:
		directory_path = shared_final+path
		files = os.listdir(directory_path)
		input = 0
		for file in files:
			if '.root' not in file:
				continue
			else:
				number_files+=1
				if number_files is 0:
					print "about to make file"
					input+=1
					output_file = open('RecoOuts/reco_file_df'+dfs+'.txt', 'w')
				code_line = eos_prefix+path+file
				output_file.write(code_line+"\n")
