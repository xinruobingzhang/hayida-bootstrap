
##################
Risk�ּ�ʵ��
ScoreSystem(form,FileRead,FileOut,bnum,model.num)
ScoreSystem�ܹ�����������������������£�
	form�����빫ʽ���磺form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
	FileRead�������ļ���·�����籾������FileReader=��ԭʼ����.csv��
	FileOut������ļ����ļ���bootstrap.csv,final result.csv,model.csv,thraval.csv�����ڵ��ļ�Ŀ¼
	bnum����bootstrap�����Ĵ�����
	model.num�����ּ�����ѡ�õķ�����
	1�����������ݵĵ�һ�ְ�cutoff�㻮�ֵķ�����
	2�����������ݵĵڶ��ְ��������ֵķ���
	3������ԭʼ���ݵĵ�һ�ְ�cutoff�㻮�ֵķ���
	4������ԭʼ���ݵĵڶ��ְ��������ֵķ���
	
	

###########################

ȫģ���Զ�������
AllModelAllBinary.R�ļ�  ����ļ�������Ǵ���������ʵ��ȫģ�͵�ȫ�������Զ�������Ĵ��롣
��������Ҫ�ĺ�����AllModelAllBinary (form,FileRead,FileOut,bnum)
�������������£�
	form�����빫ʽ���磺form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
	FileRead�������ļ���·�����籾������FileReader=��ԭʼ����.csv��
	FileOut������ļ����ļ���bootstrap.csv,final result.csv,model.csv,thraval.csv�����ڵ��ļ�Ŀ¼
	bnum����bootstrap�����Ĵ�����
	

##########################
ȫģ�ͻ�����ݲ��Զ�������
��Ӧ��������AllModelNotBinary.R�ļ����档
������Ҫ�ĺ����ǣ�hayida.bootstrap(form,FileRead,FileOut,bnum)
�������������£�
	form�����빫ʽ���磺form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
	FileRead�������ļ���·�����籾������FileReader=���������.csv��
	FileOut������ļ����ļ���bootstrap.csv,final result.csv,model.csv,thraval.csv�����ڵ��ļ�Ŀ¼
	bnum����bootstrap�����Ĵ�����


