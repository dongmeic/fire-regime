. /Users/dongmeichen/GitHub/fire-regime/CMG_2.sh 2001
. /Users/dongmeichen/GitHub/fire-regime/CMG_3.sh 2002
. /Users/dongmeichen/GitHub/fire-regime/CMG_4.sh 2002
for year in `seq 2003 2016`
do
	bash /Users/dongmeichen/GitHub/fire-regime/CMG_1.sh $year
done

