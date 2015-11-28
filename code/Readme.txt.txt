
##################
Risk分级实验
ScoreSystem(form,FileRead,FileOut,bnum,model.num)
ScoreSystem总共五个参数，各参数含义如下：
	form：输入公式，如：form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
	FileRead：数据文件的路径，如本例中是FileReader=’原始数据.csv’
	FileOut：输出的几个文件（bootstrap.csv,final result.csv,model.csv,thraval.csv）所在的文件目录
	bnum：是bootstrap抽样的次数，
	model.num是评分级别所选用的方法，
	1代表处理混合数据的第一种按cutoff点划分的方法。
	2代表处理混合数据的第二种按比例划分的方法
	3代表处理原始数据的第一种按cutoff点划分的方法
	4代表处理原始数据的第二种按比例划分的方法
	
	

###########################

全模型自动二分类
AllModelAllBinary.R文件  这个文件里面的是代码是用来实现全模型的全部变量自动二分类的代码。
其中最重要的函数是AllModelAllBinary (form,FileRead,FileOut,bnum)
各参数含义如下：
	form：输入公式，如：form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
	FileRead：数据文件的路径，如本例中是FileReader=’原始数据.csv’
	FileOut：输出的几个文件（bootstrap.csv,final result.csv,model.csv,thraval.csv）所在的文件目录
	bnum：是bootstrap抽样的次数，
	

##########################
全模型混合数据不自动二分类
相应代码是在AllModelNotBinary.R文件里面。
其中重要的函数是：hayida.bootstrap(form,FileRead,FileOut,bnum)
各参数含义如下：
	form：输入公式，如：form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
	FileRead：数据文件的路径，如本例中是FileReader=’混合数据.csv’
	FileOut：输出的几个文件（bootstrap.csv,final result.csv,model.csv,thraval.csv）所在的文件目录
	bnum：是bootstrap抽样的次数，


