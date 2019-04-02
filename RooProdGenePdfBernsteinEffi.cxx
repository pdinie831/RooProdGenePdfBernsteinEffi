/***************************************************************************** 
 * Project: RooProdGenePdfBernsteinEffi                                      * 
 *                                                                           * 
 * P.Dini fecit, Anno Domini MMXVIII                                         *
 * "Est modus in rebus."                                                     * 
 *                                                                           * 
 * Class to compute 3D angular efficiency x P-wave in B0->K*MuMU Analysis    * 
 *                                                                           * 
 *****************************************************************************/ 


#include "Riostream.h" 

#include "RooProdGenePdfBernsteinEffi.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h"
#include "RooRealProxy.h" 
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include <math.h> 
#include "TMath.h"  
#include <algorithm>  

ClassImp(RooProdGenePdfBernsteinEffi); 
   
  
   RooProdGenePdfBernsteinEffi::RooProdGenePdfBernsteinEffi(const char *name, const char *title, 
                        RooAbsReal& x, RooAbsReal& y, RooAbsReal& z, 
			const RooArgList& coefList1,
			const RooArgList& coefList2,
			int maxDegree1, int maxDegree2, int maxDegree3
			) :
   RooAbsPdf(name,title), 
   _x("_x","_x",this,x),
   _y("_y","_y",this,y),
   _z("_z","_z",this,z), 
   _coefList1("_coefList1","_coefList1",this),
   _coefList2("_coefList2","_coefList2",this),
   _maxDegree1(maxDegree1),
   _maxDegree2(maxDegree2),
   _maxDegree3(maxDegree3)
{ 
    TIterator* coefIter1 = coefList1.createIterator() ;
     RooAbsArg* coef1 ;
     while((coef1 = (RooAbsArg*)coefIter1->Next())) {
       if (!dynamic_cast<RooAbsReal*>(coef1)) {
    	 std::cout << "RooProdGenePdfBernsteinEffi::ctor(" << GetName() << ") ERROR: coefficient " << coef1->GetName()
    	      << " is not of type RooAbsReal" << std::endl ;
    	 assert(0) ;
       }
       _coefList1.add(*coef1) ;
//   	 std::cout << "RooProdGenePdfBernsteinEffi::ctor(" << GetName() << ") coefficient " << coef->GetName()<<std::endl ;
      }
    TIterator* coefIter2 = coefList2.createIterator() ;
     RooAbsArg* coef2 ;
     while((coef2 = (RooAbsArg*)coefIter2->Next())) {
       if (!dynamic_cast<RooAbsReal*>(coef2)) {
    	 std::cout << "RooProdGenePdfBernsteinEffi::ctor(" << GetName() << ") ERROR: coefficient " << coef2->GetName()
    	      << " is not of type RooAbsReal" << std::endl ;
    	 assert(0) ;
       }
       _coefList2.add(*coef2) ;
//   	 std::cout << "RooProdGenePdfBernsteinEffi::ctor(" << GetName() << ") coefficient " << coef->GetName()<<std::endl ;
      }
     delete coefIter1;
     delete coefIter2;
     _maxDegreeV = std::max(maxDegree3,std::max(maxDegree1,maxDegree2));
     printf("RooBernstein: Max(Numbers of degree) = %d\n",_maxDegreeV);
     _numParamsE = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);
     printf("RooBernstein: Numbers of parameters = %d\n",_numParamsE);
     for(int ii = 0; ii <NormInteg  ; ++ii) {
      int i = ii+_numParamsE;;
      Norm[ii] = ((RooAbsReal&) _coefList2[i]).getVal();
       printf("RooBernstein: Norm[%d] = %f\n",ii,Norm[ii]);
     }
 } 
//==============================================================================================
RooProdGenePdfBernsteinEffi::RooProdGenePdfBernsteinEffi(const RooProdGenePdfBernsteinEffi& other, const char* name) :  
   RooAbsPdf(other,name), 
   _x("_x",this,other._x),
   _y("_y",this,other._y),
   _z("_z",this,other._z), 
   _coefList1("_coefList1",this,other._coefList1),
   _coefList2("_coefList2",this,other._coefList2),
   _maxDegree1(other._maxDegree1),
   _maxDegree2(other._maxDegree2),
   _maxDegree3(other._maxDegree3),
   _maxDegreeV(other._maxDegreeV),
   _numParamsE(other._numParamsE)
   
{ 
     for(int ii = 0; ii <NormInteg  ; ++ii) {
      int i = ii+_numParamsE;;
      Norm[ii] = ((RooAbsReal&) _coefList2[i]).getVal();
       printf("RooBernstein: Norm[%d] = %f\n",ii,Norm[ii]);
     }
} 
//==============================================================================================
Double_t RooProdGenePdfBernsteinEffi::evaluate() const 
{ 
   fptype x    = _x; // x, y, z...
   fptype y    = _y; // x, y, z...
   fptype z    = _z; // x, y, z...
   double F_L	= ((RooAbsPdf&) _coefList1[0]).getVal();
   double P_1	= ((RooAbsPdf&) _coefList1[1]).getVal();
   double P_2	= ((RooAbsPdf&) _coefList1[2]).getVal();
   double P_3	= ((RooAbsPdf&) _coefList1[3]).getVal();
   double P4p	= ((RooAbsPdf&) _coefList1[4]).getVal();
   double P5p	= ((RooAbsPdf&) _coefList1[5]).getVal();
   double P6p	= ((RooAbsPdf&) _coefList1[6]).getVal();
   double P8p	= ((RooAbsPdf&) _coefList1[7]).getVal();
   double F_T	= (1.-F_L);
   double cosxq =     x*x;
   double cosyq =     y*y;
   double sinxq = (1.-cosxq);
   double sinyq = (1.-cosyq);
   double cos2x = (2.*cosxq-1.);
   double sin2x =  2.*x*sqrt(sinxq);
   double sin2y =  2.*y*sqrt(sinyq);
   double sinx  =  sqrt(sinxq);
  
   double pdf = 9./(32.*TMath::Pi())*(3./4.*F_T*sinyq+F_L*cosyq+\
   (1./4.*F_T*sinyq-F_L*cosyq)*cos2x+0.5*P_1*F_T*sinyq*sinxq*TMath::Cos(2.*z)+\
   sqrt(F_L*F_T)*TMath::Cos(z)*(P4p*sin2y*sin2x+P5p*sin2y*sinx)+\
   sqrt(F_L*F_T)*TMath::Sin(z)*(P6p*sin2y*sinx+P8p*sin2y*sin2x)+\
   2.*P_2*F_T*sinyq*x-P_3*F_T*sinyq*sinxq*TMath::Sin(2.*z));
 
 
       double xmin = _x.min();
       double xdif = _x.max()-xmin;
              x    =(_x-xmin)/xdif;
       double ymin = _y.min();
       double ydif = _y.max()-ymin;
              y    =(_y-ymin)/ydif;
       double zmin = _z.min();
       double zdif = _z.max()-zmin;
              z    =(_z-zmin)/zdif;
       
       double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
       double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
       double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       
       
       sx[0]=1.0;
       sy[0]=1.0;
       sz[0]=1.0;
       for( int i = 1; i <= _maxDegreeV ; ++i){
        sx[i]= sx[i-1]*(1.-x);
        sy[i]= sy[i-1]*(1.-y);
        sz[i]= sz[i-1]*(1.-z);
       }
       
   
// 
       double bernknvalx = 0.;
       double bernknvaly = 0.;
       double bernknvalz = 0.;
       
       
       int ipar =0;
       double func =0.0;
       
       double tx = 1.;
//       long double sx = sx_ini;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
       
         double ty = 1.;
//         long double sy =sy_ini;
         for(int j = 0; j <= _maxDegree2 ; ++j) {
	  bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
	 
          double tz = 1.;
//          long double sz = sz_ini;
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
	   func += ((RooAbsReal&) _coefList2[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) _coefList[ipar]).getVal()<<std::endl;
	   ipar++;
	   tz *= z;
	  }
	   ty *= y;
         }
	   tx *= x;
       }
//        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//      if(func<1.E-30)  std::cout<<"evaluate = "<<func<<" x,y,z ="<<_x<<" "<<_y<<_z<<" "<<std::endl;
 
       if(func<1.E-30) func=1.E-30; 
//       func = func/(xdif*ydif*zdif); 
       func=pdf*func;
//        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//       exit(0);
//      func=func/(1.+maxDegree1)/(1.+maxDegree2)/(1.+maxDegree3);
//      std::cout<<"eval = "<<func/(xdif*ydif*zdif)<<std::endl;
      return  func/(xdif*ydif*zdif);
} 

fptype RooProdGenePdfBernsteinEffi::device_coeffbinomial  (fptype enne, fptype kappa) const {
        fptype factor=1.;
        for(fptype i = 1; i <=kappa; ++i) {
          factor *= (enne+1-i)/i; 
        }	 
        if (factor<=0 ){
	 printf("Error in RooProdGenePdfBernsteinEffi coeffbinomial=> factor = %f enne=%f kappa=%f",factor,enne,kappa);
         return 0;
	} 
       return factor;
}
//==============================================================================================
fptype  RooProdGenePdfBernsteinEffi::device_bernsteinkn_func  (fptype x, fptype enne, fptype kappa) const{
   return RooProdGenePdfBernsteinEffi::device_coeffbinomial(enne,kappa)*pow(x,kappa)*pow(1.0-x,enne-kappa);
}
//==============================================================================================
Int_t RooProdGenePdfBernsteinEffi::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if       (matchArgs(allVars, analVars, _x,_y,_z)){
   std::cout<<"getAnalyticalIntegral==1: Analytic integral over x y z => Yes, available"<<std::endl;
   return 1;
  }else if (matchArgs(allVars, analVars, _y,_z)){
   std::cout<<"getAnalyticalIntegral==2: Analytic integral over y z => Not available"<<std::endl;
   return 0;
  }else if (matchArgs(allVars, analVars, _x,_z)){
   std::cout<<"getAnalyticalIntegral==3: Analytic integral over x z => Not available"<<std::endl;
   return 0;
  }else if (matchArgs(allVars, analVars, _x,_y)){
   std::cout<<"getAnalyticalIntegral==4: Analytic integral over x y => Not available"<<std::endl;
   return 0;
  }else{  
   std::cout << "Error in RooBernsteinEffi::analyticalIntegral" << std::endl;
   return 0;
 } 
}
// //==============================================================================================
Double_t RooProdGenePdfBernsteinEffi::analyticalIntegral(Int_t code, const char* rangeName) const
{
    if(code==1){
     double xdif  = _x.max()-_x.min();
     double ydif  = _y.max()-_y.min();
     double zdif  = _z.max()-_z.min();
     double F_L   = ((RooAbsPdf&) _coefList1[0]).getVal();
     double P_1   = ((RooAbsPdf&) _coefList1[1]).getVal();
     double P_2   = ((RooAbsPdf&) _coefList1[2]).getVal();
     double P_3   = ((RooAbsPdf&) _coefList1[3]).getVal();
     double P4p   = 2.*((RooAbsPdf&) _coefList1[4]).getVal();
     double P5p   = ((RooAbsPdf&) _coefList1[5]).getVal();
     double P6p   = -((RooAbsPdf&) _coefList1[6]).getVal();
     double P8p   = 2.*((RooAbsPdf&) _coefList1[7]).getVal();
//      double F_L = 6.04202e-01;
//      double P4p = 0.00000e+00;
//      double P5p = 0.00000e+00;
//      double P6p = 0.00000e+00;
//      double P8p = 0.00000e+00;
//      double P_1 = -1.88576e-01;
//      double P_2 = 4.16427e-01;
//      double P_3 = 0.00000e+00;
//      double F_T   = (1.-F_L);
//      double F_L = 1.;
//      double P4p = 0.00000e+00;
//      double P5p = 0.00000e+00;
//      double P6p = 0.00000e+00;
//      double P8p = 0.00000e+00;
//      double P_1 = 0.;
//      double P_2 = 0.;
//      double P_3 = 0.00000e+00;
//      double F_T   = 0.;
     double F_T   = (1.0-F_L);

// 
     double norm = (F_T*Norm[0] + F_L*Norm[1] +\
     F_T*(-Norm[2]) + F_L*Norm[3] + P_1*F_T*Norm[4] +\
     sqrt(F_L*F_T)*(P4p*(-Norm[5]) + P5p*(-Norm[6]))+\
     sqrt(F_L*F_T)*(P6p*(-Norm[7]) + P8p*(-Norm[8]))+\
     P_2*F_T*Norm[9] - P_3*F_T*(-Norm[10]));
// 
 
     norm=norm/(xdif*ydif*zdif);

//     double     norm1=0.398921;
//     double     norm1=0.398921/(xdif*ydif*zdif);
     
//     
//     std::cout<<std::scientific << std::setprecision(5)<<"norm = "<< norm <<"   norm1 = "<< norm1 <<std::endl;
     return norm;
//     return norm/(xdif*ydif*zdif);
    }else{
    
    return 0;
    }

    return 0;
}    
