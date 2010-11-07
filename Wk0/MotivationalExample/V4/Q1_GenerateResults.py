
import os
import util


def Create_HOC_For_Diam( diam, hocfilename ):
	tmpl = util.ReadFile( "Model4_Q1.hoc.tmpl")
	substFile = tmpl.replace("__DIAM__", str(diam) )
	util.WriteFile( hocfilename, substFile )
	
util.ClearDirectory("tmp")

for diam in [0.5,1.0,1.5]:
	hoc_filename="tmp/Model4_Q1_%2.2f.hoc"%diam
	
	# Create the hoc-file
	Create_HOC_For_Diam( diam, hoc_filename )
	
	# Run the hoc-file
	os.system("nrngui " + hoc_filename)

	
