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
   
For people who want to do a joint fit in configuration of N_slices = 5,  
please do remenber to prepare those 3 x 4 x 5 histograms and follow the format above.
