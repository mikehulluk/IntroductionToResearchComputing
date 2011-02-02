import os
from os.path import join, getsize
for root, dirs, files in os.walk('.'):
    print root, "consumes",
    print sum(getsize(join(root, name)) for name in files),
    print "bytes in", len(files), "files"
    if '.git' in dirs:
        dirs.remove('.git')  # don't visit CVS directories
