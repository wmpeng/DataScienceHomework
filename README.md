# 1 Implement Smith-Waterman Algorithm

This is a amplement of Smith-Waterman algorithm, maintaned by:
  Mingpeng Wang(1511212) <wmpeng@outlook.com>

For more information about this homework please visit github repository <https://github.com/wmpeng/DataScienceHomework> (that will be visible at 24:00, July 5th.)

# 2 Contents

```
homework
 ├── DataScienceHomework                vs2015 porject
 │    ├── DataScienceHomework
 │    │    └ smithWaterman.cpp          source code file
 │    └── DataScienceHomework.sln
 ├── DataScienceHomework.exe            executable file of release version
 ├── A.txt                              sample file : the first DNA sequence file
 ├── B.txt                              sample file : the second DNA sequence file
 ├── S.txt                              sample file : score matrix and gap penalty file
 └── README.md
```

My homework include a Vistual Studio 2015 project folder, a An executable file, three parameter files named "A.txt", "B.txt", "S.txt", and a readme file that you are reading.

The source code's path is ".\DataScienceHomework\DataScienceHomework\smithWaterman.cpp".

The executable file is "DataScienceHomework.exe".

# 3 Run Process

The process can run without installation.

## 3.1 Input and parameters

You can run the executable file in Powershell by command line like this : DataScienceHomework <1stFile> <2ndFile> <3rdFile>

<1stFile> and <2ndFile> are the first and the second DNA sequences file that only contain upper case letter {A,T,G,C} (line breaks are allowed).

<3rdFile> talls the scoring matrix and gap penalty ,whose format likes that：

```
scoreMatrix:
	A	T	G	C
A	1	-2	-2	-2
T	-2	1	-2	-2
G	-2	-2	1	-2
C	-2	-2	-2	1

GapPenalty:
-2
```

↑ When the SingleGapPenalty isn't define ,which is default ↑

```
scoreMatrix:
...

GapPenalty:
-1 -2 -3 -4
```
↑ When the SingleGapPenalty isn't define ,which is default ↑

## 3.2 Output

The result will be output at screen. When the the parameter files are the sample files (different and gap -2, the same +1), the result is :

```
highest value is 18.000000
sequences 1 from  1758 to  1784 : GAAGTTGAGCACGTAAGGCCTACGTTC
sequences 2 from  2004 to  2030 : GAAGTTGAGCACGTGATGCCTAAGTTC
```

# Contact

Contact me at anytime :

  github <https://github.com/wmpeng/DataScienceHomework>

  email <wmpeng@outlook.com>
