
f = open('ARID1B_filled.csv')
fR = open('ARID1B_reduced.csv','w')

zeroCt = 0
for line in f:
	# print(line.split(','))
	if(line.split(',')[1].strip() == '0'):
		zeroCt+=1
		if(zeroCt > 9):
			zeroCt = 0
			fR.write(line)
	else:
		fR.write(line)

