import os,hashlib
from os.path import join, getsize

def fileMD5(filename):
  with open(filename) as f:
    m = hashlib.md5(f.read()).hexdigest()
    return m
    
#Calculate MD5 and size of every file:
md5sums = {}
for root, dirs, files in os.walk('.'):
    if '.git' in dirs:
        dirs.remove('.git')  
    for f in files:
      fName = join(root, f)
      k = ( getsize(fName), fileMD5(fName) )
      if not k in md5sums: 
        md5sums[k] = [fName]
      else: 
        md5sums[k].append( fName )
      
# Find Duplicates:
duplicates = []
for k,v in md5sums.iteritems():
  if len(v) > 1: duplicates.append( (k,v) ) 
duplicates.sort()

#Print out duplicates
print "Duplicates Found:"
for k,v in duplicates:
  print k[0],"bytes:", v

