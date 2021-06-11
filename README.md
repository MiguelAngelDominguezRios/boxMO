# Compiling the Sources
There is a makefile file to compile the sources. This makefile mut be edited prior to the compilation to set the correct location of the CPLEX path and the name of the CPLEX libraries, which depend on the version installed. The current makefile is prepared for CPLEX 12.6.2 running on a Linux distribution.

Once you edited the makefile, run it with the `make` command to generate an executable linked with the static library of CPLEX. 

# Executing TPA and holzmann

To execute algorithm TPA, type:

`./BOXES (arg1) (arg2) (arg3) 0 p-partition holzmann reduced_scaled alternate RE 1 1`

To execute algorithm holzmann, type:

`./BOXES (arg1) (arg2) (arg3) 0 full holzmann volume 1 RE 1 1`

- (arg1) is the objective costs file
- (arg2) is the `.lp` file
- (arg3) is the maximum total execution time in seconds (0 for unlimited time)

# Other options

The software provide more options to execute, such as limitation of the number of solutions. To see all of them, type `./BOXES`

`./BOXES (arg1) (arg2) (arg3) (arg4) (arg5) (arg6) (arg7) (arg8) (arg9) (arg10) (arg11)`

- (arg1) is the objective costs file
- (arg2) is the `.lp` file
- (arg3) is the maximum total execution time in seconds (0 for unlimited time)
- (arg4) is the maximum size of the Pareto front (0 for unlimited size)
- (arg5) is the type of partition. Type 1(`full`) or 2(`p-partition`)
- (arg6) is the parameterization model. Type 1(`chalmet`), 2(`tchebycheff`) or 4(`holzmann`)
- (arg7) is the box value. Type 11(`volume`) or 14(`reduced_scaled`)
- (arg8) is the number of set of boxes used. Type `1` or `alternate`.
- (arg9) is the filtering proccess. Type 1(`RE`).
- (arg10) is the value of CPX_PARAM_PARALLEL. Type `-1`, `0` or `1`
- (arg11) is the value of CPX_PARAM_THREADS. Type `0` or (a positive integer)

# Objective file and lp file

All the objective cost files and .lp files are provided in the folder `Instances`

## Objective file archive must have the following format:
```
p
n  m
c_11 c_12 ... c_1n
c_21 c_22 ... c_2n
......
c_p1 c_p2 ... c_pn
```
Instance example: (AP_p-3_n-5_ins-1)	Assignment problem with dimension 3; 25 variables and 10 constraints. The costs of the objective functions are in the three last lines.

```
3
25 10
6 1 20 2 3 2 6 9 10 18 1 6 20 5 9 6 8 6 9 6 7 10 10 6 2 
17 20 8 8 20 10 13 1 10 15 4 11 1 13 1 19 13 7 18 17 15 3 5 1 11 
10 7 1 19 12 2 15 12 10 3 11 20 16 12 9 10 15 20 11 7 1 9 20 7 6
```

## LP file must be the following format (provided by CPLEX):

```\ENCODING=ISO-8859-1
\Problem name: AP_p-3_n-5_ins-1
Minimize
 obj:
Subject To
 c1:  x1 + x2 + x3 + x4 + x5  = 1
 c2:  x6 + x7 + x8 + x9 + x10  = 1
 c3:  x11 + x12 + x13 + x14 + x15  = 1
 c4:  x16 + x17 + x18 + x19 + x20  = 1
 c5:  x21 + x22 + x23 + x24 + x25  = 1
 c6:  x1 + x6 + x11 + x16 + x21  = 1
 c7:  x2 + x7 + x12 + x17 + x22  = 1
 c8:  x3 + x8 + x13 + x18 + x23  = 1
 c9:  x4 + x9 + x14 + x19 + x24  = 1
 c10: x5 + x10 + x15 + x20 + x25  = 1
Bounds
 0 <= x1 <= 1	 
 0 <= x2 <= 1	 
 0 <= x3 <= 1	 
 0 <= x4 <= 1	 
 0 <= x5 <= 1
 0 <= x6 <= 1	 
 0 <= x7 <= 1	 
 0 <= x8 <= 1	 
 0 <= x9 <= 1	 
 0 <= x10 <= 1
 0 <= x11 <= 1	 
 0 <= x12 <= 1	 
 0 <= x13 <= 1	 
 0 <= x14 <= 1	 
 0 <= x15 <= 1
 0 <= x16 <= 1	 
 0 <= x17 <= 1	
 0 <= x18 <= 1	
 0 <= x19 <= 1	
 0 <= x20 <= 1
 0 <= x21 <= 1	 
 0 <= x22 <= 1	
 0 <= x23 <= 1	 
 0 <= x24 <= 1	 
 0 <= x25 <= 1
Binaries
 x1  x2  x3  x4  x5  x6  x7  x8  x9  x10  x11  x12  x13  x14  x15  x16  x17  x18  x19  x20  x21  x22  x23  x24  x25 
End
```

