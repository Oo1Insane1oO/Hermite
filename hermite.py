########################################################################
## Script for generating hermite polynomials upto order n (with H_0=1) #
########################################################################

import math
import sys
import sympy as sp

def fac(n):
    """ calculate factorial recursively """
    if n == 0:
        """ default power to the zeroth """
        return 1
    else:
        return n*fac(n-1)
    # and ifelse
# end function fac

def hermiteCoefficients(n):
    """ use recursive relation for coefficients and return a table with
    coefficients each from 0 to n """
    table = [[1], [0,1]]
    while (n>=len(table)):
        i = len(table)
        preva = table[i-1]
        ppreva = table[i-2]
        tmp = [0 for i in range(i+1)]
        m = 1 - i
        tmp[0] =  m*ppreva[0]
        tmp[i-1] = preva[i-2]
        tmp[i] = preva[i-1]
        for k in range(1,i-1):
            tmp[k] = preva[k-1] + m*ppreva[k]
        # end fork
        table.append(tmp)
    # end while

    # convert to physicists'
    for i,t in enumerate(table):
        tmp = sp.sympify(2**(sp.Rational(i,2)))
        for k,tt in enumerate(t):
            table[i][k] = tmp*sp.sympify(2**(sp.Rational(k,2)))*tt
    return table
# end function hermiteCoefficients

def hermite(x,n):
    """ Use recursive relation """
    if n<0:
        """ ignore negative idices """
        return 0
    elif n==0:
        """ first value """
        return 1
    else:
      return 2*(x*hermite(x,n-1)-(n-1)*hermite(x,n-2))
# end function hermite

def choose(nom, dom):
    """ calculate ratio of integer factorials (n!/(n-m)!)"""
    res = 1
    for i in range(dom):
        res *= nom - i
    # end fori
    return res
# end cFacRatio

def hermiteDerivative(x,n,m=1):
    """ return the derivative """
#     if m > n:
#         return 0
#     else:
#         return 2**m*math.factorial(m)*math.factorial(n)/math.factorial(n-m) * hermite(x,n-m)
    return sp.diff(hermite(x,n),x,m)
# end function hermiteDerivative

def usage():
    """ print usage """
    print "USAGE: python hermite.py 'order' 'filename'\n"
# end function usage

def turnToCPP(n,H,funcName):
    """ turn hermite polynomial H of degree n into C++ template function """
    return ''' template<typename T> static T ''' + funcName + '''%(n)i(const T& x) {return %(expr)s;}''' % (
            {'n':n, 'expr':sp.printing.ccode(H)})
# end functtion turnToCPP

def turnCoeffsToCPP(coeffs, funcName):
    """ turn coefficients to lists in template function """
    codes = []
    for n,c in enumerate(coeffs):
#         s = "   std::vector<long int> coeffs = std::vector<long int>{"
        s = "std::vector<long int>{"
        for i in range(len(c)):
            """ create string with vector of coefficients """
            if (i==len(c)-1):
                s += "%i};\n" % c[i]
            else:
                s += "%i," % c[i]
            # end ifelse
        # end fori
#         codes.append('''static std::vector<long int> ''' + funcName + '''%i() {\n ''' % n + s + '''\n    return coeffs;\n}''')
        codes.append(s);
    return codes
# end function turnCoeffsToCPP

def appendToFile(codes, fname, funcName, o="w+"):
    """ append polynomial templates to file """
    with open(fname, o) as ofile:
        ofile.write("#pragma GCC diagnostic push\n")
        ofile.write('#pragma GCC diagnostic ignored "-Wunused-parameter"\n')
        for c in codes:
            ofile.write(c+"\n")
        # end for c
        ofile.write("#pragma GCC diagnostic pop\n")
        ofile.write("template<typename T> static T " + funcName + "(const T& x, const int& n) {\n")
        ofile.write("   if (n > %i) {\n" % (len(codes)-1))
        ofile.write("       return -1;\n")
        ofile.write("   }\n")
        ofile.write("   switch(n) {\n")
        for i in range(len(codes)):
            ofile.write(("       case %i: return " + funcName + "%i(x);\n") % (i,i))
        ofile.write("       default: return 0;\n")
        ofile.write("   }\n}\n")
    # end ofile
# end function appendToFile

def appendCoeffsToFile(coeffs, fname, funcName, o="w+"):
    with open(fname, o) as ofile:
        sizes = [0]
        totalsize = sum(len(c) for c in coeffs)
        for i in range(0,len(coeffs)-1):
            sizes.append(sizes[i] + len(coeffs[i]))
        # end fori
        ofile.write("struct HC {\n"
                    "   private:\n"
                    "      static constexpr std::array<long int, %i> coeffs = "
                    "{\n" % totalsize)
        for c in range(len(coeffs)):
            for s in range(len(coeffs[c])):
                if ((s == len(coeffs[c])-1) and (c == len(coeffs)-1)):
                    ofile.write("%i" % coeffs[c][s])
                else:
                    ofile.write("%i," % coeffs[c][s])
            # fors
            ofile.write("           \n")
        # end forc
        ofile.write("       };\n")
        ofile.write("       static constexpr std::array<unsigned int, %i> "
                "displ = {" % len(coeffs))
        for c in range(len(coeffs)):
            if c != (len(coeffs)-1):
                ofile.write(str(sizes[c]) + ",")
            else:
                ofile.write(str(sizes[c]) + "};\n")
            # end ifelse
        # end forc
        ofile.write("\n     public:\n")
        ofile.write("           static constexpr long int coeff(unsigned int n, unsigned "
                    "int i) {\n"
                    "               return coeffs[displ[n] + i];\n"
                    "           }\n};\n\n")
    # end ofile
# end function appendCoeffsToFile

if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        filename = sys.argv[2]
    except IndexError:
        usage()
        raise IndexError
    # end try-except

    x = sp.symbols('x')
    Hlist = [hermite(x,0), hermite(x,1)]
    for i in range(2, n+1):
        Hlist.append(2*(x*Hlist[i-1] - (i-1)*Hlist[i-2]));
    # endif

    Hcodes = [turnToCPP(i,sp.powsimp(sp.simplify(Hlist[i])), "H") for i in range(n+1)]
#     Hdercodes = [turnToCPP(i,sp.powsimp(sp.simplify(sp.diff(Hlist[i],x,1))), "dH") for i in range(n+1)]
#     Hderdercodes = [turnToCPP(i,sp.powsimp(sp.simplify(sp.diff(Hlist[i],x,2))), "ddH") for i in range(n+1)]
    coeffCodes = turnCoeffsToCPP(hermiteCoefficients(n+1), "HC")

    with open(filename, "w+") as ofile:
        hname = filename.replace(".h", "_H").upper()
        ofile.write("#ifndef " + hname + "\n")
        ofile.write("#define " + hname + "\n\n")
        ofile.write("#include <vector>\n")
        ofile.write("#include <cmath>\n")
        ofile.write("#include <array>\n\n")
    # end with open
    appendCoeffsToFile(hermiteCoefficients(n+1), filename, "HC", "a")
    appendToFile(Hcodes,filename,"H", "a")
#     appendToFile(Hdercodes,filename,"dH","a")
#     appendToFile(Hderdercodes,filename,"ddH","a")
    with open(filename, "a") as ofile:
        ofile.write("#endif /* " + hname + " */\n")
    # end withopen
# end ifmain
