/* calcStat_FtU.js
*
*    Copyright (c) 2016 Yuji SODE <yuji.sode@gmail.com>
*
*    This software is released under the MIT License.
*    See LICENSE or http://opensource.org/licenses/mit-license.php
*/
//============================================================================
//it returns calculated statistical values of F, t, and U, with given numerical arrays: X and Y.
function _FtU(X,Y){
  var Xav,Yav,vrX,vrY,F,t,U,nX=X.length,nY=Y.length,slf=this;
  //it returns average of a given numerical array A.
  var av=function(A){
    var s=0,i=0,n=A.length;
    while(i<n){s+=A[i],i+=1;}
    return s/n;};
  //it returns [variance,difference] with a given numerical array A, and average of A: av, while difference=total of (x-average)*(x-average), and variance=difference/(n-1).
  var vr=function(A,av){
    var b=0,s=0,i=0,n=A.length;
    while(i<n){b=A[i]-av,s+=b*b,i+=1;}
    return [s/(n-1),s];};
  //it returns a rank object available for the statistical rank test, by a given numerical array: A.
  var rk=function(A){
    var i=0,j=0,s=0,Oe,N=A.length,O={},B=slf.JSON.parse(slf.JSON.stringify(A)).sort(function(a,b){return a-b;});
    while(i<N){
      //O is rank object; O[B[I]]=[rank,count,true rank].
      if(!O[B[i]]){O[B[i]]=[i+1,1,i+1];}else{O[B[i]][1]+=1;}
      i+=1;}
    for(var el in O){
      Oe=O[el];
      if(Oe[1]>1){
        j=0,s=0;while(j<Oe[1]){s+=Oe[0]+j,j+=1;}
        Oe[2]=s/Oe[1];}
    }
    return O;};
  //------------ <functions estimating values of statistic> ------------
  //Zav=av(Z), and vZ=vr(Z,average of Z), while Z is a numerical array (X and Y), and nZ is Z.length.
  //it returns statistic F as [value of F,degree of freedom in numerator,degree of freedom in denominator].
  var fCal=function(vX,nX,vY,nY){
    if(vX[0]/vY[0]<1){
      return [vY[0]/vX[0],nY-1,nX-1];}
    else{
      return [vX[0]/vY[0],nX-1,nY-1];}
  };
  //it returns statistic t as [value of t,degree of freedom].
  var tCal=function(Xav,vX,nX,Yav,vY,nY){
    return [Math.abs(Xav-Yav)/(Math.sqrt((vX[1]+vY[1])/(nX+nY-2))*Math.sqrt((1/nX)+(1/nY))),nX+nY-2];};
  //it returns statistic U as [value of U, U approximated with standard normal distribution, 1st degree of freedom, 2nd degree of freedom].
  var uCal=function(A,nA,B,nB){
    var C=[],R2=0,I=0,u=0,u1=0,u2=0,Obj;
    I=0;while(I<nA){C.push(A[I]),I+=1;}
    I=0;while(I<nB){C.push(B[I]),I+=1;}
    Obj=rk(C);
    //rank total by B.
    I=0;while(I<nB){R2+=Obj[B[I]][2],I+=1;}
    u1=nA*nB+(nB*(nB+1))/2-R2,u2=nA*nB-u1,u=Math.min(u1,u2);
    return [u,Math.abs(u-(nA*nB/2))/Math.sqrt(nA*nB*(nA+nB+1)/12),nA,nB];};
  //------------ </functions estimating values of statistic> ------------
  Xav=av(X),Yav=av(Y),vrX=vr(X,Xav),vrY=vr(Y,Yav);
  //=== estimation of F value; F=[F,n1,n2]. ===
  F=fCal(vrX,nX,vrY,nY);
  //=== estimation of t value; t=[t,n]. ===
  t=tCal(Xav,vrX,nX,Yav,vrY,nY);
  //=== estimation of U value; U=[U,(U-m)/s,n1,n2]. ===
  U=uCal(X,nX,Y,nY);
  return {f:F,t:t,u:U};}
