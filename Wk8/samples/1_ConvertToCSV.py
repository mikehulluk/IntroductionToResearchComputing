import datetime

fdate = None
for line in open("1_testdata.txt"):
  if not line: continue
  date,value = line.strip().rsplit(" ",1)
  d = datetime.datetime.strptime(date.strip(),"%Y-%m-%d %H:%M:%S")
  if fdate  is None: fdate = d
  t_elapsed = (d - fdate).total_seconds()
  print "%d, %s"%(t_elapsed, value)
