from collections import defaultdict
import numpy as np
import pandas as pd
import time

class Reader:  # 读文件，返回一个字典{序列名:序列}
    def __init__(self, infasta=None):
        self.infasta = infasta

        self.dom = {}
        self.data = []

    def _read_data(self):

        """Sets data to stripped lines from the fasta file
        """

        with open(self.infasta) as infasta:
            self.data = [l.strip() for l in infasta]

    def _upper_seq_per_line(self):  # self.data中字符串统统变大写
        """Sets data to upper case, single line sequences for each header
        """
        new_data = []
        dom1 = {}
        seq = ""
        for i, line in enumerate(self.data):
            if line[0] == ">":
                if seq:
                    new_data.append(seq.upper())  # 小写字母转换为大写字母
                    seq = ""
                else:
                    assert i == 0, "There may be a header without a sequence at line {}.".format(i)
                new_data.append(line)
            else:
                seq += line
        new_data.append(seq.upper())
        k = 0
        for index, item in enumerate(new_data):
            if index % 2 == 0:
                dom1[item] = new_data[index + 1]
                k += 1
        self.dom = dom1

    def get_lines(self):
        self._read_data()
        self._upper_seq_per_line()
        return self.dom


class cut:  # 第二个字典
    def __init__(self, **dom):
        self.dom = dom

    def kmer(self):

        dom2 = defaultdict(int)

        for key, values in self.dom.items():
            seq = values
            length = len(values)
            for k in range(1,7):#*********************************KMER
                for c in range(length - k + 1):
                    kmer = seq[c:c + k]
                    dom2[kmer] += 1
            # print('标准化之前', dom2)
            # f = open(r"C:\Users\shilinyue\Desktop\x.txt", 'a')
            # f.write(str(dom2))
            # f.write("\n")
            # l = 0
            # for keys, values in dom2.items():
            #     l += dom2[keys]

            for keys, values in dom2.items():
                dom2[keys] = (values / length)  # / l
            # print('标准化之后', dom2)
            # f.write(str(dom2))
            # f.write("\n")
            self.dom[key] = dom2
            # dom2 = {}
            dom2 = defaultdict(int)
        for key, values in self.dom.items():
            v = self.dom[key]
            jichu = ['A', 'G', 'C', 'T','N']
            list = sorted(v.keys(), key=lambda i: len(i), reverse=True)
            # print('字符从长到短排序：', list)
            # f.write('字符从长到短排序：'+str(list))
            # f.write("\n")

            for i in list:
                if i in jichu:
                    break
                else:
                    # print('字符排序后从最长依次出来', i)
                    # print('父字符和子字符', v[i], v[i[:-1]])
                    # f.write('字符排序后从最长依次出来' + i)
                    # f.write('父字符和子字符' +str(v[i])+ ' '+str(v[i[:-1]]))
                    # f.write("\n")
                    # f.close()
                    v[i] = v[i] / v[i[:-1]]
                    self.dom[key][i] = v[i]

        # print('更新后的字典', self.dom)

        return self.dom


def jihe():
    s = list(hum.keys())
    h=list(hum.keys())
    ls = len(s)
    lh=len(h)
    test = pd.DataFrame(np.zeros((lh, ls)), index=h, columns=s)
    for i, va in enumerate(s):

        for j, vb in enumerate(h):
            if vb != va:
                value1 = defaultdict(int, hum[va])
                value2 = defaultdict(int, hum[vb])
                a = list(value1.keys())
                b = list(value2.keys())
                # print(type(value2.keys()))

                # print('当前谁是擂主', a)

                # print('打擂者', b)
                e = list(set(a).union(set(b)))  # 并
                f = list(set(a).intersection(set(b)))  # 交
                # print('先并后交', e, f)
                # 求交集数据集合，并集数据集合
                # 交----求最小值
                inter = []
                for i in f:
                    if value1[i] < value2[i]:
                        inter.append(value1[i])
                    else:
                        inter.append(value2[i])
                # print('交对应的数值集合', inter)
                # 并---求最大值
                u = []
                for i in e:

                    # ----------------------------------------------------------------------
                    if value1[i] > value2[i]:
                        # -----------------------------------------------------------------
                        u.append(value1[i])
                    else:
                        u.append(value2[i])
                # print('并对应的数值集合', u)
                # 求相似----向量模之比
                x = np.linalg.norm(inter)
                y = np.linalg.norm(u)
                # print('交的模，并的', x, y)
                sim = x / y
                test.loc[vb, va] = sim

            else:
                test.loc[va, vb] = 1

        # print('相似性矩阵', test)
    test.to_csv(r"v21.csv")  # C:\Users\shilinyue\Desktop\1.csv#C:\Users\DELL\Desktop***********

a_sta=time.time()
def file(infasta):
    ss = Reader(infasta)  # /home/zxteng/shily/list/1000.fa
    c = ss.get_lines()
    dd = cut(**c)
    c = dd.kmer()
    return c







# print("这是c")
# print(c)
hum=file(infasta=r"v21-01.fa")
#mos=file(infasta=r"1000.fa")
s=time.time()
jihe()
e=time.time()
print('Tooki %f second' % (e - s))
a_end=time.time()
print('Tooko %f second' % (a_end - a_sta))
