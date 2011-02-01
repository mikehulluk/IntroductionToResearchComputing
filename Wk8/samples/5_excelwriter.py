from xlwt import Workbook
book = Workbook()
sheet1 = book.add_sheet('Sheet 1')
book.add_sheet('Sheet 2')
sheet1.write(0,0,'Parameter')
sheet1.write(0,1,'Value')

sheet1.row(1).write(0,'Input Resistance')
sheet1.row(1).write(1,'2000')
sheet1.row(2).write(0,'Capacitance')
sheet1.row(2).write(1,'1.0')
book.save('simple.xls')

