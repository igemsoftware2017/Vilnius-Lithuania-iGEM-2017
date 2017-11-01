import os
import cairosvg
from IPython import embed

INPUT_DIR = '/home/aurimas/test/pdf'

for (dirpath, dirnames, filenames) in os.walk(INPUT_DIR):
    for filename in filenames:
    	if '.svg' in filename:
	    	if not os.path.exists(dirpath + '/pdf'):
	        	os.mkdir(dirpath + '/pdf')

	    	cairosvg.svg2pdf(url=dirpath+'/'+filename, write_to=dirpath+'/pdf/'+filename[:-3] + '.pdf')