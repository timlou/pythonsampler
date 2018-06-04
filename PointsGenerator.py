import numpy as np
import sympy as sp
import vegas
import math
from sympy.parsing.sympy_parser import parse_expr
from sympy.solvers.inequalities import reduce_rational_inequalities
from sympy import Interval, S
from functools import reduce
import csv

#define some useful functions

#replace string by a dictionary
def str_replace(st,d):
    outst=st
    for k,v in d.items():
        outst = outst.replace(k, v)
    return outst

#convert strings such as -xi or xi into i
#useful for sympy expression parser
def var_ind(st):
    return int(st.strip('-')[1:])




# a class to contain all the expressions for generating points
class PointsGenerator():
    dim = 0
    constraints = []
    # list of variable names
    var_list = None
    # list of sympy variables
    spvar_list = None
    npt=1
    
    def __init__(self, filename):
        with open(filename, "r") as infile:
            iline = 0
            var_sub = None

            for line in infile:
                # ignore white space 
                li=line.strip()
                if li.startswith("#"):
                    continue

                # first line should be dimensions
                if iline==0:
                    self.dim=int(li)
                    
                    if self.dim > 40:
                        ValueError('Currently sobol sequence is used in the algorithm and only d<=40 is implemented'+ \
                                   ' an invalid input dimension {0:d} is found'.format(self.dim))
                    
                    self.dim_low=[0]*self.dim
                    self.dim_high=[1]*self.dim
                    
                # second line for testing
                elif iline==1:
                    iline+=1
                    # implement later
                    continue
                    
                # all other lines should be conditions
                else:

                    # construct list of variables for symbolic computation
                    if self.var_list is None or var_sub is None:
                        # get names of all variables
                        self.var_list = ['x'+str(i) for i in range(self.dim)]
                        
                        # substitution rules for convert x[i] to xi
                        var_sub = {'x[{0:d}]'.format(i):'x{0:d}'.format(i) for i in range(self.dim)}
                        
                        # declare sympy variable
                        self.spvar_list = sp.symbols(reduce(lambda x,y: ""+x+' '+y, self.var_list))
                        
                    if '>=' in li:
                        pieces=li.split('>=')
                    elif '<=' in li:
                        pieces=reverse(li.split('<='))                       
                    
                    
                    expr_st='('+pieces[0].strip()+')-('+pieces[1].strip()+')'
                    expr_st=str_replace(expr_st,var_sub)
                    expr=parse_expr(expr_st)
                    
                    
                    # now see if the constraint can be converted to simple boundaries on individual variables
                    try:
                        if len(expr.free_symbols) == 1:
                            # find which variable this expr refers to
                            n=var_ind(str(list(expr.free_symbols)[0]))
                            xn=self.spvar_list[n]

                            # error may occur if inconsistent boundaries exist!
                            sol=sp.solve([str_replace(li,var_sub), self.dim_low[n]<=xn, self.dim_high[n]>=xn],xn)
                                                        
                            #get the final interval
                            interval_list=[]
                            for l in sol.args:
                                interval_list.append(sp.solveset(l,xn,domain=S.Reals))
                            final_interval=reduce(lambda x,y:x.intersection(y),interval_list)
                            
                            self.dim_low[n]=final_interval.args[0]
                            self.dim_high[n]=final_interval.args[1]
                            iline+=1
                            continue

                        
                        # try to simplify the constraints
                        #else:
                        #    print("trying to simplify")
                        #    solexpr = [str_replace(li,var_sub)]
                        #    print(solexpr)
                        #    sol=sp.solve(solexpr,list(expr.free_symbols)[0])
                        #    print(sol)
                        # TODO: sympy cannot solve multivariate inequalities currently
                        # so multivariate inequalities do not work
                            
                    except AttributeError:
                        # if inconsistent inequalities exist, an error is thrown
                        raise ValueError('the condition '+li+' is an inconsistent constraint! no solution available')
                    except NotImplementedError:
                        print("cannot simplify "+li+", continue.")
                    
                    self.constraints.append(expr)
                    
                iline+=1
            print("{0:d} non-trivial constraint(s), all > 0:".format(len(self.constraints)))
            for i in self.constraints:
                print(i)
            
            #convert to numpy arrays
            self.dim_low=np.asarray(self.dim_low, dtype=np.double)
            self.dim_high=np.asarray(self.dim_high, dtype=np.double)
        
    # test a sample pt on non-trivial constraint
    def valid(self,pt):
        for expr in self.constraints:
            val=expr.evalf(subs={ s:v for s,v in zip(self.spvar_list,pt) })
            if val < 0:
                return False
        return True
    
    def gen(self, n=1, outfile=None, nstop=1000000, reset=False):
        # maximum number of steps to try
        nmax=n*1000+nstop
        
        import sobol_seq as so

        #output file specified
        if outfile is not None:
            fout = open(outfile, 'w')
            writer = csv.writer(fout, lineterminator='\n')
            
        # generate the same points again
        if reset:
            self.npt=1
            
        ngen=0
        #cutoff to prevent infinite loop
        nstep=0
        result=[]

        #adaptive strategy flag
        adapt=False
        
        while ngen<n and nstep < nmax and not adapt:
            pt_test,self.npt = so.i4_sobol(self.dim, self.npt)
            
            # scale to boundary
            pt_test=self.dim_low + np.asarray(pt_test,dtype=np.double)*(self.dim_high-self.dim_low)
            
            #test if working
            if self.valid(pt_test):
                result.append(pt_test)
                # write to output
                if outfile is not None:
                    writer.writerow(pt_test)
                    
                ngen+=1
            nstep+=1

            # check and see if generation is going too slowly...
            if nstep >= 1000 and ngen <100:
                print("basic algorithm is too slow, switching to an adaptive algorithm.")
                adapt=True

        if not adapt:
            print("basic algorithm generated {0:d} proper points.".format(ngen))
        
        else:
            print("Beginning adaptive strategy...")
                
            # define an integrand that will aid sampling
            # utilize the vegas package for adaptive sampling
            def integrand(x):
                pt=np.asarray(x)
                fx=1.0
                for expr in self.constraints:
                    val=expr.evalf(subs={ s:v for s,v in zip(self.spvar_list,pt) })
                    if val < 0:
                        return 0
                    # cap function off at 1
                    fx*= (val if val < 1.0 else 1.0)
                return fx

            # declare vegas integration, will be used for random number generation
            integrator = vegas.Integrator([ [l,h] for l,h in zip(self.dim_low,self.dim_high) ])

            # try performing the integral and give points generated
            # max number of iterations
            max_it = 3
            nval=5*n
            it=0
            while it < max_it:
                # the integrator will generate points for us
                int_reult=integrator(integrand, nitn=5, neval=nval)
                nval*=5

                integrator_list = [ np.asarray(x) for x,wt in integrator.random() if integrand(x) > 0]

                print("adaptive algorithm generated {0:d} proper points.".format(len(integrator_list)))
                
                # if not enough good sampled points, try higher neval                
                if len(integrator_list) < n and it < max_it-1:
                    it+=1
                    continue
                else:
                    if len(integrator_list) < n:
                        ngen=len(integrator_list)
                    else:
                        ngen=n

                    result= integrator_list[:ngen]


                    for i in result:
                        if not self.valid(i):
                            print("NOT valid!")
                        
                    if outfile is not None:
                        for pt_test in result:
                            writer.writerow(pt_test)
                    break
                    
                

        if outfile is not None:
            fout.close()

        if ngen < n:
            print('only {0:d} points are successfully generated.'.format(ngen))
            
        return np.asarray(result, dtype=np.double)


if __name__ == "__main__":
    import sys

    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    n_event = int(sys.argv[3])
    #load the file
    mygen = PointsGenerator(input_file)
    mygen.gen(n=n_event, outfile=output_file)
    print("exiting...")
    sys.exit()
