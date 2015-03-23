import sys
import re

def spliter(x):
    if len(x)>0:
        return x.split(",")
    else:
        return []

def getHeader(string):
    headerRE = re.compile(r"""contents=(.*)\n\[BAYES-NETWORK:BAYES-PROPERTIES\]""", re.MULTILINE)
    header   = headerRE.findall(string)[0]
    headerRE = re.compile(r"P\(([^\s\|]*)\|?(\S*)\)")
    header   = headerRE.findall(header)
    spliter1 = re.compile(r"[\[,]")
    header   = [[spliter1.split(i[0][0:-1]),spliter(i[1])] for i in header]
    return header

def getBody(string):
    body = string.split("[BAYES-NETWORK:BAYES-TABLE]")
    if len(body) < 2:
        return []
    body = body[1]
    bodyRE = re.compile(r"(?:P\((.*)\)=([\d.]*))")
    body = bodyRE.findall(body)
    spliter = re.compile(r"[,|]")
    body = map(lambda x: spliter.split(x[0])+[x[1]], body)
    return body

def stractValue(string, values):
    if string.find("+") == 0:
        return string[1:], 0
    if string.find("-") == 0:
        return string[1:], 1
    out = string.split("=")
    try:
        out[1] = values[out[0]].index(out[1])
    except ValueError:
        out[1] = int(out[1])
    return out


def stractInfo(string):
    header = getHeader(string)
    values = {}
    cpds   = {}
    #find out the values for each random variable
    for line in header:
        cpds[tuple([line[0][0]]+line[1])] = []
        values[line[0][0]] = line[0][1:]

    #now, find the cpds, please
    #ok, first find the size of each cpd
    for line in header:
        curr_cpd =[line[0][0]]+line[1]
        sizes = map(lambda x: len(values[x]), curr_cpd)
        cpds[tuple(curr_cpd)] = [0.0]*reduce(lambda x, y: x*y, sizes)

    body = getBody(string)

    #finally, fill the cpds with the correct prob. distribution
    for line in body:
        #print "::: ", line
        vars  = map(lambda x: stractValue(x, values), line[:-1])
        sizes = map(lambda x: len(values[x[0]]), vars)
        #print "--->", vars
        pos   = int(vars[0][1])
        for i in range(1,len(vars)):
            pos = pos*(sizes[i]) + vars[i][1]
        cpds[tuple(map(lambda x: x[0], vars))][pos] = float(line[-1])

    #thanks, return it
    return values, cpds

def convert2BBNF(string_in, file_out = open("/Users/brenowcarvalho/Documents/BBNet", "w")):
    #given the variables and factors, it is time to write the file in my format
    vars, cpds = stractInfo(string_in)
    vars_keys = vars.keys()
    cpds_keys = cpds.keys()

    str_vars = {}

    for var_count in range(len(vars_keys)):
        name = vars_keys[var_count]
        values = " ".join(vars[vars_keys[var_count]])
        str_vars[name] = "( %s | %s )" %(name, values)
        file_out.write("%02d Variable %s\n" %(var_count, str_vars[name]))

    for cpd_count in range(len(cpds_keys)):
        #f_vars = " ".join(cpds_keys[cpd_count])
        f_vars = map(str_vars.get, cpds_keys[cpd_count])
        file_out.write("%02d Factor [< %s > %s]\n" %(cpd_count, " ".join(f_vars), " ".join(map(str,cpds[cpds_keys[cpd_count]]))))

    len_vars = len(vars)
    for sc_count in range(len(cpds_keys)):
        fact_id = len_vars + sc_count
        for var in cpds_keys[sc_count]:
            var_id  = vars_keys.index(var)
            file_out.write("%02d -> %02d\n" %(fact_id, var_id))
    file_out.flush()

def main(args):
    if(len(args) < 3):
        print """
    USAGE: python %s <destination_folder> <file_to_convert>+
    This command converts a Bayesian Network in a EG format to my temporary format
        """ %(args[0])
        return

    for name in args[2:]:
        print name
        prefix = "".join(name.split(".")[:-1])
        f_out  = open(args[1]+"BBNet-"+prefix+".txt", "w")
        f_in   = open(name, "r")
        str_in = "".join(f_in.readlines())
        convert2BBNF(str_in, f_out)
        f_out.close()
        f_in.close()

if(__name__ == "__main__"):
    main(sys.argv)
