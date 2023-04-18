#!/usr/bin/python

def leps(rab, x):
    import sys
    from math import exp
    
    a=0.05
    b=0.30
    c=0.05
    dab=4.746
    dbc=4.746
    dac=3.445
    r0=0.742
    alpha=1.942
    rac=3.742
    kc=0.2025
    rbc=rac-rab
    
    
    Qab= dab/2*(3/2*exp(-2*alpha*(rab-r0))-exp(-alpha*(rab-r0)))
    
    
    Qbc= dbc/2*(3/2*exp(-2*alpha*(rbc-r0))-exp(-alpha*(rbc-r0)))
    
    Qac= dac/2*(3/2*exp(-2*alpha*(rac-r0))-exp(-alpha*(rac-r0)))
    
    
    Jab= dab/4*(exp(-2*alpha*(rab-r0))-6*exp(-alpha*(rab-r0)))
    
    
    Jbc= dbc/4*(exp(-2*alpha*(rbc-r0))-6*exp(-alpha*(rbc-r0)))
    
    Jac= dac/4*(exp(-2*alpha*(rac-r0))-6*exp(-alpha*(rac-r0)))
    

    result = Qab/(1+a)+Qbc/(1+b)+Qac/(1+c) - ((Jab/(1+a))**2 + (Jbc/(1+b))**2 + (Jac/(1+c))**2-Jab*Jbc/((1+a)*(1+b))- Jbc*Jac/((1+b)*(1+c)) - Jab*Jac/((1+a)*(1+c)))**0.5 + 2*kc*(rab-(rac/2-x/1.154))**2

    return result
 

def lepsdr(rab, x):
    import sys
    from math import exp
    
    a=0.05
    b=0.30
    c=0.05
    dab=4.746
    dbc=4.746
    dac=3.445
    r0=0.742
    alpha=1.942
    rac=3.742
    kc=0.2025
    rbc=rac-rab
    
    
    Qab= dab/2*(3/2*exp(-2*alpha*(rab-r0))-exp(-alpha*(rab-r0)))
    dQab= dab/2*(3/2*(-2*alpha)*exp(-2*alpha*(rab-r0))-(-alpha)*exp(-alpha*(rab-r0)))
    
    
    Qbc= dbc/2*(3/2*exp(-2*alpha*(rbc-r0))-exp(-alpha*(rbc-r0)))
    dQbc= dbc/2*(3/2*2*alpha*exp(-2*alpha*(rbc-r0))-alpha*exp(-alpha*(rbc-r0)))
    
    Qac= dac/2*(3/2*exp(-2*alpha*(rac-r0))-exp(-alpha*(rac-r0)))
    dQac=0
    
    
    Jab= dab/4*(exp(-2*alpha*(rab-r0))-6*exp(-alpha*(rab-r0)))
    dJab= dab/4*((-2*alpha)*exp(-2*alpha*(rab-r0))-6*(-alpha)*exp(-alpha*(rab-r0)))
        
    Jbc= dbc/4*(exp(-2*alpha*(rbc-r0))-6*exp(-alpha*(rbc-r0)))
    dJbc= dbc/4*((2*alpha)*exp(-2*alpha*(rbc-r0))-6*alpha*exp(-alpha*(rbc-r0)))
    
    Jac= dac/4*(exp(-2*alpha*(rac-r0))-6*exp(-alpha*(rac-r0)))
    dJac= 0

    result1 = dQab/(1+a)+dQbc/(1+b) - 0.5*((Jab/(1+a))**2 + (Jbc/(1+b))**2 + (Jac/(1+c))**2-Jab*Jbc/((1+a)*(1+b))- Jbc*Jac/((1+b)*(1+c)) - Jab*Jac/((1+a)*(1+c)))**(-0.5)*(2*Jab*dJab/(1+a)**2+2*Jbc*dJbc/(1+b)**2-dJab*Jbc/((1+a)*(1+b))-Jab*dJbc/((1+a)*(1+b))-dJbc*Jac/((1+b)*(1+c))-dJab*Jac/((1+a)*(1+c))) + 4*kc*(rab-(rac/2-x/1.154))
    result2 = 4*kc/1.154*(rab-rac/2+x/1.154)

    return result1, result2
 


geo = open('geometry.in')
lines = geo.readlines()[3]
geo.close()
lines = lines.split()
d1,d2 = lepsdr(float(lines[1]), float(lines[2]))
d1 = -d1
d2 = -d2
ener = leps(float(lines[1]), float(lines[2]))

print "  | Number of atoms                   :        1 "
print "  | Total energy corrected        :         %10.10f \n" % ener
print "Total atomic forces (unitary forces cleaned) "
print " |  1   %10.10f   %10.10f   0.0000000" % (d1, d2)
