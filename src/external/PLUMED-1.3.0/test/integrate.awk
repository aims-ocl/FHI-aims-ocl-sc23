#!/usr/bin/awk -f
#
# usage:
# integrate.awk time=1 cv=2 constr=3 kappa=100.0 COLVAR 
#
BEGIN{tot=0;i=0;
   if(ARGC<6){
    print "usage: "
    print "integrate.awk time=1 cv=2 constr=3 kappa=100.0 COLVAR "
    exit
   }
}
{

   if($2!="FIELDS"){
    if(i==0){valoldconstr=$constr;valoldtime=$time}
    else{
      valtime=$time  
      valconstr=$constr
      valcv=$cv  
      vel=(valconstr-valoldconstr)/(valtime -valoldtime)
      if (vel!=0.){
        tot+=-kappa*vel*(valcv-valconstr)*(valtime -valoldtime)
        valoldtime=valtime
        valoldconstr=valconstr   
        printf("%12.6f %12.6f %12.6f\n",valtime,valconstr,tot) 
      }
    }
    i++
   }
}

