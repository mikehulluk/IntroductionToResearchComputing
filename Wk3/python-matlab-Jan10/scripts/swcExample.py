#
#swcFile = open( "data/091113 dlc 293 LHS1-short.swc")
#
#for line in swcFile.readlines():
#	print line,
#
#
#swcFile = open( "data/091113 dlc 293 LHS1.swc")
#for line in swcFile.readlines():
#	lineSplit = line.split()
#	print lineSplit
#	
#	
#
#
#swcFile = open( "data/091113 dlc 293 LHS1-short.swc")
#for line in swcFile.readlines():
#	lineSplit = line.split()
#	
#	id = int( lineSplit[0] )
#	type = int( lineSplit[1] )
#	x = float( lineSplit[2] ) 
#	y = float( lineSplit[3] ) 
#	z = float( lineSplit[4] ) 
#	rad = float( lineSplit[5] ) 
#	parent_id = int( lineSplit[6] ) 
#
#	print """ID: %d (%2.2f,%2.2f,%2.2f) Rad:%2.2f Parent:%d""" % (id,x,y,z,rad,parent_id)
#
#
#
#
#
#
#
#
#
#
#import pprint
#idDictionary = {}
#
#swcFile = open( "data/091113 dlc 293 LHS1-short.swc")
#for line in swcFile.readlines():
#	lineSplit = line.split()
#	
#	id = int( lineSplit[0] )
#	type = int( lineSplit[1] )
#	x = float( lineSplit[2] ) 
#	y = float( lineSplit[3] ) 
#	z = float( lineSplit[4] ) 
#	rad = float( lineSplit[5] ) 
#	parent_id = int( lineSplit[6] ) 
#
#	idDictionary[id] = [ x,y,z,parent_id]
#
#
#pprint.pprint( idDictionary )
#
#
#
#
#
#
#
#
#idDictionary = {}
#
#swcFile = open( "data/091113 dlc 293 LHS1-short.swc")
#for line in swcFile.readlines():
#	lineSplit = line.split()
#	
#	id = int( lineSplit[0] )
#	type = int( lineSplit[1] )
#	x = float( lineSplit[2] ) 
#	y = float( lineSplit[3] ) 
#	z = float( lineSplit[4] ) 
#	rad = float( lineSplit[5] ) 
#	parent_id = int( lineSplit[6] ) 
#
#	idDictionary[id] = [ x,y,z,rad,parent_id]
#
#
#for id in idDictionary.keys():
#	swcData = idDictionary[id]
#	print id, " -> ",swcData
#	x,y,z,rad,parent_id = swcData
#
#
#
#
#
#
#
#
#
#from math import pi,sqrt
#
#def calcLength( x1,y1,z1, x2,y2,z2):
#	x = x1-x2
#	y = y1-y2
#	z = z1-z2
#	return sqrt( x*x + y*y + z*z )
#
#def calcSurfaceAreaOfChoppedCone(r1,r2,l):
#	return (r1+r2)*pi*l
#
#
#idDictionary = {}
#
#swcFile = open( "data/091113 dlc 293 LHS1-short.swc")
#for line in swcFile.readlines():
#	lineSplit = line.split()
#	
#	id = int( lineSplit[0] )
#	type = int( lineSplit[1] )
#	x = float( lineSplit[2] ) 
#	y = float( lineSplit[3] ) 
#	z = float( lineSplit[4] ) 
#	rad = float( lineSplit[5] ) 
#	parent_id = int( lineSplit[6] ) 
#
#	idDictionary[id] = [ x,y,z,rad,parent_id]
#
#
#for id, swcData in idDictionary.iteritems():
#	x,y,z,rad,parent_id = swcData
#	
#	if parent_id == -1:
#		print "Root node - nothing to do"
#	
#	else:
#		parent_swc_data = idDictionary[ parent_id ]
#		p_x,p_y,p_z,p_rad,p_parent_id = parent_swc_data
#	
#		l = calcLength(	x,y,z, p_x,p_y,p_z)
#		sa = calcSurfaceAreaOfChoppedCone(rad,p_rad,l)
#		print "ID: %d - Length %2.2f, SA: %2.2f"%(id, l,sa)	
#
#
#
#
#
#
#
#
#
#
#
#
#
#




from math import pi,sqrt
from numpy import array
from scipy.linalg import norm

swcFile = open( "data/091113 dlc 293 LHS1-short.swc")

linesSplit = [ line.split() for line in swcFile.readlines() ] 
linesSplitForDict = [ (int(l[0]), (array( ( float(l[2]), float(l[3]), float(l[4]) ) ), float(l[5]), int(l[6]) ) )  for l in linesSplit ] 
idDictionary = dict(linesSplitForDict )  

for id, (xyz,rad,parentid) in idDictionary.iteritems():
	if parentid == -1: continue

	p_xyz,p_rad,p_parentid = idDictionary[ parentid ]
	l = norm( p_xyz - xyz )
	sa = (rad + p_rad) * pi * l
	print "ID: %d - Length %2.2f, SA: %2.2f"%(id, l,sa)	

