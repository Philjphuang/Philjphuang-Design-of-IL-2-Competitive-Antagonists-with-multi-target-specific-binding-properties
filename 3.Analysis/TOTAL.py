# 需要的文件是所有的resfile和一个所有任务输出的文件
#只要修改下面两个地方
#1. 修改数字的范围
name_list = range(10,29)
i = 0
fw = open('total.out','w')
#2. 修改整个输出文件的名字
fr = open('con.out')
data = fr.readlines()
dict = 0
for line in data:
    if 'TOTAL' in line:
        if dict !=1:
            fw.write(str(name_list[i]))
            with open('resfile.'+str(name_list[i])+'.in') as fa:
                for index in fa.readlines():
                    if 'PIKAA' in index:
                        fw.write(' '+index.strip().split()[-1])
        fw.write(' '+line.strip().split()[-1])
        dict += 1
        if dict == 2:
            i += 1
            fw.write('\n')
            dict = 0



