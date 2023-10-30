# GACR - a genetic algorithm based channel router

## Description
This project provides a genetic algorithm based channel router implemented by C++  
Paper: https://www.cs.ou.edu/~thulasi/misc/genetic_channel_routing.pdf

## Usage

* ### Build
```
make
```
(The project can not be compiled under c++98 or earlier versions)
 
* ### Run
```
./GACR --i [input_file] --o [output_file] --p 500 --g 100 --d 500 --h 5 --r 0.5 --opt 1
```
 `--i`  path of input file (necessary)  
 `--o`  path of output file (necessary, will generate layout guide file, gdt file, gds file)  
 `--p`  population of the genetic algorithm (default to 500)    
 `--g`  number of generation to perform in the genetic algorithm (default to 150)  
 `--d`  number of descendant generated in each generation (default to 500)  
 `--h`  initial number of rows in the channel (default to 5)  
 `--r`  return r'th best solution in the population (default to 1, which is the best)  
 `--opt`  will not perform post optimization if set to zero (default to 1, which will do post optimization)  
 
* ### Input Format  
   * Please refer to: https://github.com/yoyojs200602/PDA/blob/main/Lab4/Lab4.pdf  

* ### Output Format  
   * Output will contain a layout guide file, a gdt file and a gds file

* ### Example
   * case2 (.gds file)  
     * ![image](https://github.com/yoyojs200602/GACR/blob/baaf57f3da20f9d4aae0552e24dc6530fd0a0158/pic/case2.png)

   
 ## 
