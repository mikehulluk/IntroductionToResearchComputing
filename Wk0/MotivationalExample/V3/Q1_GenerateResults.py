
import os
import util


def Create_HOC_For_Diam( diam, hocfilename ):
	tmpl = util.ReadFile( "Model3.hoc.tmpl")
	substFile = tmpl.replace("__DIAM__", str(diam) )
	util.WriteFile( hocfilename, substFile )
	
util.ClearDirectory("tmp")

for diam in [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6]:
	hoc_filename="tmp/Model3_%2.2f.hoc"%diam
	
	# Create the hoc-file
	Create_HOC_For_Diam( diam, hoc_filename )
	
	# Run the hoc-file
	os.system("nrngui " + hoc_filename)

	
