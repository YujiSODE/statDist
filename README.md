# statDist
Web Worker interface for statistical test (F-test, t-test, and U-test).  
https://github.com/YujiSODE/statDist

>Copyright (c) 2016 Yuji SODE \<yuji.sode@gmail.com\>  
>This software is released under the Artistic License  
>See LICENSE or http://dev.perl.org/licenses/ and http://dev.perl.org/licenses/artistic.html
______

## Script
* statDist.js

## How to use
1. To create Worker object with "statDist.js"  
   
2. The following formatted data is available for the created Worker object:  
   two csv formatted numerical arrays separated with '@',  
  e.g., '3,4,5@0,1,2'  
   
3. The Worker object returns an object with upper probabilities as follows:  
  {f: p-value for F, t: p-value for t, u: p-value for U}.

## Library list
* statistics-distributions-js_Ymdf.js (Ben Tilly,2008; modified by Yuji SODE, 2016):  
  the Artistic License or the GPL v1 or later;  
  https://github.com/YujiSODE/statistics-distributions-js_Ymdf
  
* calcStat_FtU.js (Yuji SODE, 2016):  
  the MIT License;  
  https://github.com/YujiSODE/calcStat_FtU
