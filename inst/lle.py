#!/usr/local/bin/python

import sys
import sklearn
import sklearn.manifold
import numpy
import getopt

def main(argv):
    global inputfile
    global outputfile
    global neighbors
    global components
    global method
    inputfile = ''
    outputfile = ''
    neighbors = 30
    components = 2
    method = ''
    try:
        opts, args = getopt.getopt(argv,"h:i:o:n:c:m:",["help","in=","out=","neighbors","components","method"])
    except getopt.GetoptError:
        print('lle.py -i <inputfile> -o <outputfile> -n <neighbors> -c <components> -m <method>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h","--help"):
            print('lle.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--in"):
            inputfile = arg
        elif opt in ("-o", "--out"):
            outputfile = arg
        elif opt in ("-n", "--neighbors"):
            neighbors = int(arg)
        elif opt in ("-c", "--components"):
            components = int(arg)
        elif opt in ("-m", "--method"):
            method = arg
            
def lle(path,out,n_neighbors, n_components, s_method):
    Y = numpy.loadtxt(path,delimiter=",")
    if s_method == '':
        m = sklearn.manifold.LocallyLinearEmbedding(n_neighbors=n_neighbors, n_components=n_components)
    else:
        m = sklearn.manifold.LocallyLinearEmbedding(n_neighbors=n_neighbors, n_components=n_components,method=s_method)
    X = m.fit_transform(Y)
    numpy.savetxt(out,X,fmt='%.10f',delimiter=',',newline='\n')      
    return(X)

 
if __name__ == "__main__":
    main(sys.argv[1:])
    X=lle(inputfile,outputfile,neighbors,components,method)
