# Usage 
## Preparing histograms

For joint fit,  
you need to prepare histograms of **time since last muon(dtlSH)**,  
their are **3 catogories** of muon which classified by muon.    
Each muon catogory should have a dtlSH histogram,  
but we further slice the histogram by **N_slices = 5** for example.  
This slicing is based on certain IBD quantity,  
and we have **4 IBD quantities** Ep, Ed, dT, dist.

As a summary, 
+ Joint fit a site should have **3(muon energy) x 4(IBD quantities) x 5(N_slices)** histograms.
- **Example - format of histogram name:** 
   + **dtlSH_EH1_Nslice5_0_0_0**  
   _*EH = 1  
   N_slices = 5  
   Muon energy range = 0(low)  
   IBD quantity = 0(Ep)  
   Index of slice = 0*_  
   (Numbers are following the order in format.)
+ Program should be able to access all the histograms in a single rootfile.
+ Change the parameter of file name in main.cc `string rootFileName = "example.root"`, which contains all histograms.
   
For people who want to do a joint fit in configuration of N_slices = 5,  
please do remenber to prepare those 3 x 4 x 5 histograms and follow the format above.

## Modify some important parameters in main.cc

+ First is making sure that what configuration of N_slices,  
if you prepared histograms for N_slices = 5,  
please check the `int N_slices = 5` in main.cc.

+ Second is to check the `int livetime[8]`,  
which will be used to transfer livetime to daily rates,  
*ps.The livetime here should consider the effect of muon veto efficiency.*

## Modify the Makefile

You can just type the comment below to generate execution file: 
```bash
make
```

This is defined in Makefile, so what important is to make sure the name of executable file is what you want,
for example that I will define the executable file name for joint fit of N_slices = 5 as **unified5** 
```makefile
all: main.cc
    g++ -g -o unified5 main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer
```

## Run the joint fit
+ Format of command:
    `./unified5 directory name 5`
    The directory name will be used in ./FitPlot/_*directory name*_/, then the final option _*5*_ will further specify the sub-directory, for saving plots, _*s_5*_ which corresponding to the number of slices.
    As summary, plots will be saved in ./FitPlot/_*directory name*_/_*s_5*_/, and ./Slices/_*directory name*_/_*s_5*_/.
