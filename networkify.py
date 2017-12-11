#!/usr/bin/env python
# https://www.biostars.org/p/81795/
import sys
import string

cluster_filename=sys.argv[1]
network_filename=sys.argv[2]


interactions=[]

# read interactions
num_interactions=0
with open(network_filename) as in_handle:
    for line in in_handle:
        parts=line.strip().split("\t")
        ac1=parts[0]
        ac2=parts[1]
        interactions.append(ac1+"#"+ac2)
        interactions.append(ac2+"#"+ac1)
        num_interactions=num_interactions+1

print "%s interactions" % (num_interactions)

cluster=0

with open(cluster_filename) as in_handle:
    for line in in_handle:
        cluster = cluster+1
        acs=line.strip().split("\t")
        print "number of proteins: %s" % (len(acs))
        if len(acs) > 2:
            output=open(cluster_filename+ "." + str(cluster) + ".network.txt" ,  "w")

            i = 0
            while i < len(acs) :
                j=i
                while j < len(acs):
                    if acs[i] + "#" + acs[j] in interactions:
                        output.write(acs[i] + "\t" + acs[j] + "\n")
                    j=j+1
                i=i+1

            output.close()
