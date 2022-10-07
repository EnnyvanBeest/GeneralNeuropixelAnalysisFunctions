from pathlib import Path
#import matlab
import numpy as np
from mtscomp import compress

def MTSCOMPRESSION_From_Matlab(datapath,JsonPath,savepath):
	
	# Use decompress from mtscomp
	compress(datapath,outmeta=JsonPath,out=savepath, sample_rate=30000., n_channels=385,dtype=np.int16, check_after_compress = True)
	
	success=1
	#print(chunk)
	return success

#datapath = '//128.40.198.18/Subjects/EB014/2022-05-05/ephys/2022-05-05_EB014_1_g0/2022-05-05_EB014_1_g0_imec0/2022-05-05_EB014_1_g0_t0.imec0.ap.cbin'
#start_time = 0
#end_time = 5000
#channel = 384
success = MTSCOMPRESSION_From_Matlab(datapath,JsonPath,savepath)
#print(chunk+2)
