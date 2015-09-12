from math import log, ceil
from random import choice
import re
import sys
from os import makedirs

#REGEX list
double_str = r"(?:\d+(?:\.\d*)(?:e[-]?\d+)?)";
double_list_str = r"(?:"+double_str + r"(?:\s*,\s*"+double_str+r")*)"

word_str = r"(?:[\w_][\w_\d]*)"
word_list_str = r"(?:"+word_str + r"(?:\s*,\s*"+word_str+r")*)"

body_str = r"\s*network\s+"+word_str+r"\s*{\s*}\s*(.*)" #gets the body of the file (variables and probabilities)

variable_str = r"\s*variable\s*("+word_str+r")\s*{\s*type\s+discrete\s*\[\s*\d*\s*\]\s*{\s*("+word_list_str+r")\s*}\s*;\s*}"

values_str = r"(?:\(\s*("+word_str+"+(?:\s*,\s*"+word_str+r")*)\s*\)\s*(" + double_list_str + r")*\s*;\s*)"
prob_value_str = r"(?:\(\s*(?:"+word_str+r"+(?:\s*,\s*"+word_str+r")*)\s*\)\s*" + double_list_str + r"*\s*;\s*)"

probability_simple_str = r"probability\s*\(\s*("+word_str+r")()\s*\)\s*{"+ \
                                    r"\s*table\s+(" + double_list_str + r")\s*;\s*}" # gets the variable name and its distribution
probability_str = r"probability\s*\(\s*(?:("+word_str+r")\s*\|\s*(" + word_list_str + r"))\s*\)\s*{\s*(" + \
                    prob_value_str + r"+)\s*}" # gets the variable name and its distribution
#End of REGEX list

def break_list(l, foo = str.strip):
    return map(lambda x: foo(x), l.split(","))

def varsFromBIFBody(body):
    vars_str_list = re.compile(variable_str).findall(body)
    vars = map(lambda x: [x[0], break_list(x[1])], vars_str_list)
    return vars

def distributionFromBIFBody(body):
    factors_str_list = re.compile(probability_simple_str).findall(body)
    factors_str_list = map (lambda x: [x[0],[], break_list(x[2], float)], factors_str_list)
    def extract_cpd(vals_seq): #improve it to make sure it is taking the values in the right order
        valuesRE = re.compile(values_str)
        vals = valuesRE.findall(vals_seq)
        cpd = []
        for line in vals:
            cpd += break_list(line[1], float)
        return cpd
    factors_str_list1 = re.compile(probability_str).findall(body)
    factors_str_list += map(lambda x: [x[0], break_list(x[1]), extract_cpd(x[2])], factors_str_list1)

    arrows = []

    for factor in factors_str_list:
        for parent in factor[1]:
            arrows.append((parent, factor[0]))
    return factors_str_list, arrows

def structureFromBIF(in_str):
    body = re.compile(body_str).findall(in_str);
    if len(body) == 0:
        raise Exception("Invalid file format!\n")
    vars = varsFromBIFBody(body[0])
    factors, arrows = distributionFromBIFBody(body[0])
    return {"vars": vars, "factors": factors, "arrows" : arrows}
    #return {"vars": {}, "factors": {}, "arrows" : {}}

#----------------------------------------------------------

def createFakeData(variables, name, step = 100, max = 10000, suffix = ".csv"):
    header = ", ".join(map(lambda x: x[0], variables))+"\n"

    bag = []
    for line in range(step):
        vals = []
        for var in variables:
            vals.append(choice(var[1])) #get one of the variable values
        bag.append(", ".join(vals))

    for file_number in range(step, max+1, step):
        fout = open(name+"-"+str(file_number)+suffix, "w")
        fout.write(header)
        for _ in range(step, file_number+1, step):
            fout.write("\n".join(bag)+"\n")
        fout.close()

def convert2BBNF(string_in, file_out = open("/Users/brenowcarvalho/Documents/BBNet", "w")):
    #given the variables and factors, it is time to write the file in my format

    string_in = str.replace(string_in, "\r", " ")
    net = structureFromBIF(string_in)

    vars_str = {}
    vars = net["vars"]#[:10]
    vars_ids = {}
    #['var', [var1, ..., varN]]  to "Variable ( var | var1, ..., varN )"
    tuple2var_str = lambda x: "( %s | %s )" %(x[0], " ".join(x[1]))

    id_size= str(int(ceil(log(2*len(vars),10))))
    for count in range(len(vars)):
        file_out.write(("%0"+str(id_size)+"d Variable %s\n") %(count+1, tuple2var_str(vars[count])))
        vars_str[vars[count][0]] = tuple2var_str(vars[count])
        vars_ids[vars[count][0]] = count+1

    factors = net["factors"]
    factors_str = {}
    factors_ids = {}

    tuple2factor_str = lambda fact, vars: "[< %s > %s ]" %(" ".join([vars_str[fact[0]]] + map(lambda x: vars_str[x], fact[1])),\
                                                                    " ".join(map(str, fact[2])))
    for count in range(len(factors)):
        file_out.write(("%0"+str(id_size)+"d Factor %s\n") %(count+len(vars), tuple2factor_str(factors[count], vars)))
        factors_str[factors[count][0]] = tuple2factor_str(factors[count], vars)
        factors_ids[factors[count][0]] = count+len(vars)

    arrows = net["arrows"]
    arrows_str = map(lambda arrow: "%s -> %s\n" %(factors_ids[arrow[0]], vars_ids[arrow[1]]), arrows)

    for arrow in arrows_str:
        #arrows_str.append("%s -> %s\n" %(factors_ids[var[0]], vars_ids[var[0]]))
        file_out.write(arrow)

    for var in vars:
        #arrows_str.append("%s -> %s\n" %(factors_ids[var[0]], vars_ids[var[0]]))
        file_out.write(("%"+id_size+"d -> %"+id_size+"s\n") \
                         %(factors_ids[var[0]], vars_ids[var[0]]))
    file_out.flush()
    return net

def main(args):
    if(len(args) < 2):
        print """
    USAGE: python %s <input_files>+
    This command converts a Bayesian Network in a EG format to my temporary format
        """ %(args[0])
        return

    for name in args[1:]:
        print name
        prefix = "".join(name.split(".")[:-1])
        f_out  = open("BBNet-"+prefix+".txt", "w")
        f_in   = open(name, "r")
        str_in = "".join(f_in.readlines()).replace("\n", " ")
        net = convert2BBNF(str_in, f_out)
        #makedirs(prefix)
        #createFakeData(net["vars"], prefix+"/"+prefix) #TODO make it safe using os.[something] xD
        f_out.close()
        f_in.close()



if(__name__ == "__main__"):
    main(sys.argv)
