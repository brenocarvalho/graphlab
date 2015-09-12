import sys, itertools
from random import choice, sample, seed
from os import path, makedirs
from reading_bnets import createBBNF


'''
This method returns a bayesian network, a pair (vars,cpds) of maps that encompass the structure of the network and the cpd's
It is useful as a generic method of generation that permit others, more particular, methods to be easily implemented
nodes_qtd is the total of nodes in the network
choose_cardinality is a method that returns a integer, the quantity of values of a variable
choose children is a method that returns a list with all the node's children, it receives as parameter the id of the node
'''
def simpleGenerator1(nodes_qtd, choose_cardinality, choose_children):
    vars = {}
    #setting variables:
    for id in range(1,nodes_qtd +1):
        values = map(lambda x: "val%d" % x, range(1, choose_cardinality(id)+1))
        vars["Var%d" %id] = values

    #print "Vars:", vars
    scopes = [[]]*len(vars)
    cpds = {}

    for id in range(1,nodes_qtd+1):
        scopes[id-1] = ["Var%d" %id]

    for id in range(1,nodes_qtd+1):
        children = choose_children(id)
        for child in children:
            #print "Child:", child, "Parent:", id, len(scopes)
            scopes[child-1].append("Var%d" %id)

    #print "Scopes", scopes
    for var_id in range(1, len(scopes)+1):
        scope = tuple(scopes[var_id-1])
        #print "scope:", scope
        cpds[scope] = [0.0]*(reduce(lambda x,y: x*y, map(lambda z: len(vars[z]),scope)))
    #print "->", vars, cpds
    #print "GENERATOR1 DONE"
    return vars, cpds

def deterministic(nodes_level, cardinality, out_degree, levels):
    def choose_children_determ(id):
        n_children = min(nodes_level, out_degree)
        #get the closest child to the father that can afford having n_children/2 siblings on its both sides
        mean = id
        if(id / nodes_level >= levels -1):
            children = []
        else:
            if( mean - nodes_level*((id-1)/nodes_level) <= n_children/2 ):
                mean = nodes_level*((id-1)/nodes_level) + n_children/2+1
            if( nodes_level*(1+(id-1)/nodes_level) - mean  < (n_children-1)/2+ n_children%2):
                mean = nodes_level*(1+(id-1)/nodes_level) - n_children/2 - n_children%2
            mean = mean + nodes_level
            #print "Mean:", mean, "id:", id
            children = range(mean - n_children/2, mean + n_children/2 + n_children%2)
            #print "Children", children
        return children

    choose_cardinality = lambda x: cardinality
    return simpleGenerator1(nodes_level*levels, choose_cardinality, choose_children_determ)


'''
This method returns a bayesian network, a pair (vars,cpds) of maps that encompass the structure of the network and the cpd's
nodes_qtd is the total of nodes in the network
cardinality is the quantity of values each variable has
out_degree is the maximum quantity of out edges a node can have, this maximum will be reached whenever possible
r_seed is the random seed, setting it permits us to replicate experiments
'''
def random1(nodes_qtd, cardinality, out_degree, r_seed = 10):
    seed(r_seed)
    def choose_children_random(id):
        n_children = min(nodes_qtd-id, out_degree)
        if(n_children > 0):
            return sample(range(id+1, nodes_qtd+1), n_children)
        else:
            return []

    choose_cardinality = lambda x: cardinality
    return simpleGenerator1(nodes_qtd, choose_cardinality, choose_children_random)


'''
This method returns a bayesian network, a pair (vars,cpds) of maps that encompass the structure of the network and the cpd's
nodes_qtd is the total of nodes in the network
max_cardinality is the maximum quantity of values each variable has, the minimum is one
out_degree is the maximum quantity of out edges a node can have, this maximum will be reached whenever possible
r_seed is the random seed, setting it permits us to replicate experiments
'''
def random2(nodes_qtd, max_cardinality, out_degree, r_seed = 10):
    seed(r_seed)
    def choose_children_random(id):
        n_children = min(nodes_qtd-id, out_degree)
        if(n_children > 0):
            return sample(range(id+1, nodes_qtd+1), n_children)
        else:
            return []

    def choose_cardinality(id):
        if(max_cardinality > 0):
            return choice(range(1,max_cardinality+1))
        return 0
    return simpleGenerator1(nodes_qtd, choose_cardinality, choose_children_random)

def createDummyData(vars, n_lines, f_out, steps = 100):
    vars_names = vars.keys()
    #vars_names.sort()
    if(len(vars_names) > 0):
        f_out.write("%s" % vars_names[0])
    for name in vars_names[1:]:
        f_out.write(", %s" % name)
    f_out.write("\n")
    f_out.flush

    for i in xrange(n_lines):
        if(i%steps == 0):
            seed(100) #repeats the values for every "steps" lines
        if(len(vars_names) > 0):
            f_out.write("%s" % choice(vars[vars_names[0]]))
        for name in vars_names[1:]:
            f_out.write(", %s" % choice(vars[name]))
        f_out.write("\n")
        f_out.flush

def createTestNets(dest_folder, file_name = "test", card_list = (2, 3, 4, 5),
                   nodes_per_level_list = (100,), levels_list = (100,),
                   out_degree_list = (1, 2, 5), max_lines = 100000):
    if(dest_folder[-1] != path.sep):
        dest_folder = dest_folder + path.sep
    for parameters in itertools.product(nodes_per_level_list,card_list,out_degree_list,levels_list): #itertools.product generates the combination of all those lists
        params_names = (dest_folder, file_name) + parameters
        b_net_out_name = "%ssimple-%s-n%d-c%d-d%d-l%d" %params_names
        data_out_name = "%sfake_data_simple-%s-n%d-c%d-d%d-l%d" %params_names
        vars,cpds = deterministic(*parameters)

        f_out_net  = open(b_net_out_name+".txt", "w")
        createBBNF(vars, cpds, f_out_net)
        #print "FILE OUTPUT:",
        print b_net_out_name
        f_out_net.close()

        f_out_data  = open(data_out_name+".csv", "w")
        createDummyData(vars, max_lines, f_out_data)
        #print "DATA FILE OUTPUT:", data_out_name
        f_out_data.close()
        del(vars, cpds)
    return

def main(args):
    man = """
USAGE: python %s <destination_folder> <file_name> <nodes_per_level> <cardinality> <out_degree> <levels> <lines>+

This command creates a Bayesian Network in my temporary format.
<destination_folder>: Destination folder of all the generated files
<file_name>         : Base for creating the name of the files
<nodes_per_level>   : Quantity of nodes in each "level"
<cardinality>       : Quantity of values of each variable
<out_degree>        : Out degree of any non-leaf node
<levels>            : Depth of the DAG
<lines>?            : Number of lines in the data file created
""" %(args[0])

    #Check if we have the right number of parameters
    if(len(args) < 7):
        print man
        return

    #Check if all the parameters that shoould be numeric are numeric
    try:
        folder = args[1]
        if(folder[-1] != path.sep):
            folder = folder + path.sep
        file_name = args[2]
        nodes_per_level = int(args[3])
        card = int(args[4])
        o_degree = int(args[5])
        qtd_levels = int(args[6])
        if len(args) > 7:
            n_lines = int(args[7])
    except ValueError:
        print "Invalid Arguments!"
        print man
        return

    #Set output file name
    out_name = "%ssimple-%s-n%d-c%d-d%d-l%d.txt" %(folder, file_name, nodes_per_level,card,o_degree,qtd_levels)
    vars, cpds = deterministic(nodes_per_level,card,o_degree,qtd_levels)

    #Create path if needed
    if not path.exists(folder):
        makedirs(folder)

    #Open file, create the bnet, save it and close
    f_out_net  = open(out_name, "w")
    createBBNF(vars, cpds, f_out_net)
    print "FILE OUTPUT:", out_name
    f_out_net.close()

    #For each of the
    if(n_lines in dir()):
        out_name = "%sfake_data_simple-%s-n%d-c%d-d%d-l%d.csv" %(folder, file_name, nodes_per_level,card,o_degree,qtd_levels)
        f_out_data  = open(out_name, "w")
        createDummyData(vars, n_lines, f_out_data)
        print "DATA FILE OUTPUT:", out_name
        f_out_data.close()
    print "DONE"

if(__name__ == "__main__"):
    main(sys.argv)
