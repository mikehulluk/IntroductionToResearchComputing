
import shutil,os
import numpy as np

def ReadFile(fname):
	f = open(fname)
	d = f.read()
	f.close()
	return d

def WriteFile(fname, data):
	f = open(fname,"w")
	f.write(data)
	f.close()

def ClearDirectory(dName):
	if os.path.exists(dName):		
		shutil.rmtree(dName)
	os.makedirs(dName)

def DeleteIfExists(fName):
	if os.path.exists(fName):
		os.unlink(fName)


def EnsureDirectoryExists(dName):
	if not os.path.exists(dName):		
		os.makedirs(dName)

def BungleData(d):
    d[:,1] = np.log( d[:,1] )
    return d
