USAGE
=====

msort, sort file rows by multiple field, written by  Ruan Jue <ruanjue@genomics.org.cn>
Usage: msort [-hrfn] [-t <filed_separators>]  [-k [rfmne]<col>[start-end]{enum1,...},...] [<file>]
Options:
   -h  display this document
   -l  specify line brokers, you can define more than one characters
        defaultly, they are <newline>
   -L  treat N line as a block, defaultly N = 1
   -t  specify field separators, you can define more than one characters
        defaultly, they are <space>, <tab>
   -r  reverse sort, it will be overwritten if fileds are specified
   -f  ignore character`s case, ...
   -n  treat fileds as number, ...
   -m  treat fileds as number or string, ...
   -k  specify fileds to be sorted, eg.
         rn10[2-6]: reverse sort 2-6 (include 6) of column 10, treat it as number
         n11: sort column 11, treat it as number
         f2[6-]: sort 6-end of column 2, treat it as string, ignore case
         rfm3: reverse sort column 3, treat it as string or number, ignore string case
         f3{red green blue}: sort column 3 , treat it as enum, ignore string case
If given file, read text from file, otherwise from STDIN

INSTALL
=======

> tar -xzf msort.tar.gz
> cd msort
> ./autogen.sh
>./configure
>make

---------------------------------------------------------------------------------------------
中文文档

msort 是一个类似于sort的软件，它的特点是：
(1) msort可以对数字和字母混合的列进行排序, 修饰符m
(2) msort可以对枚举类型进行排序，2{red green blue}
(3) msort可以对一个列中指定的范围进行排序，例如 2[4-9]表示2列的4到9一共不大于6个字符排序
(4) msort可以指定多个行分割符，多个列分割符，这样可以方便地处理多种格式的文本文件
(5) msort可以对每个用来排序的字段进行修饰，修饰符包括 rfnme
	  i，r 表示从大到小，默认的顺序是小到大
	 ii，f 表示忽略字符的大小写
	iii，n 表示该字段为数字型
	 iv，m 表示该字段为数字或字符，数字排在字符前面，如果存在修饰符r，结果是反向的
	  v，e 表示该字段是枚举类型，一般情况下，程序发现用户给出了枚举后，就自动标识
	        该字段为枚举类型，所以e不推荐主动使用

语法：
(1) msort -l '>\n\t' -t '_:|'
    行分割符使用 '>', '\n', '\t'
	列分割符使用 '_', ':', '|'
(2) msort -r -f 
    对输入文件的所有列进行排序，把所有列都看作字符，忽略大小写，结果是反向的
(3) msort -k 'rf1[1-3]{Jan Feb Mar Apr May Jun Jul Aug Seq Oct Nov Dec}'
    对第一列的前3个字符按枚举排序，忽略大小写，排序结果是月份的顺序
(4) msort -k '3,2,10'
    多列排序，首先根据3列，3列相同的在按2列，...10列
(5) msort test.table
    对文件test.table进行排序，结果输出到标准控制台
(6) msort -
    从标准输入读取数据，进行排序，结果输出到标准控制台
(7) msort
    从标准输入读取数据，进行排序，结果输出到标准控制台
例子：
> cat test.table
chrX:+:loci100/4
chrX:+:loci110/2
chrX:+:loci110/3
chr1:+:loci110/3
chr1:-:loci110/3
> msort -t ':/' -k 'm1[4-],2{+ -},n3[5-],n4' test.table
chr1:+:loci110/3
chr1:-:loci110/3
chrX:+:loci100/4
chrX:+:loci110/2
chrX:+:loci110/3

