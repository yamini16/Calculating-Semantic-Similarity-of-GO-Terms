import copy;
import re;
from collections import defaultdict
import untangle #use the command <sudo pip install untangle> to install the library
import time;
import pickle;
import math;
import itertools;
import sys;
 
class Node:
    def __init__(self, GO_ID, f):
        self.go_id = GO_ID
        self.parents = [] #List of immediate parents
        self.f = f #Freq of the node
        self.cf = 0 #c=Cumulative freq of the node
        self.cfUpdated = False
        self.ancestors = {} #Mapping of ancestors to distance 

    def getID(self):
        return self.go_id

    def namespace(self, ns):
        self.namespace = ns

    def addParent(self, node):
        self.parents.append(node)

class Cache:
    def __init__(self, term_f, gene_annotation, term_list, term_graph, roots):
        self.term_f = term_f
        self.gene_annotation = gene_annotation
        self.term_list = term_list
        self.term_graph = term_graph
        self.roots = roots

class GeneSimiliarity:
    def __init__(self):
        try:
            cache = pickle.load(open("cache", "rb")) #cache=pickled graph file
            self.term_f = cache.term_f
            self.gene_annotation = cache.gene_annotation
            self.term_list = cache.term_list
            self.term_graph = cache.term_graph
            self.roots = cache.roots
        except:
            print "Parsing Gene Association..."
            self.term_f, self.gene_annotation = self.parseGeneAssociation()
            print "Parsing GO Terms..."
            self.term_list, self.term_graph, self.roots = self.parseAnnotations()
            print "Computing Frequency and Ancestors for GO Terms..."
            for root in self.roots.values():
                self.compute_cf_and_ancestors(root)
            print "Saving Generated Cache.."
            pickle.dump(Cache(self.term_f, self.gene_annotation, self.term_list, self.term_graph, self.roots), open("cache", "wb"))

    def parseGeneAssociation(self, fileName = "gene_association.goa_ref_human"):
        try:
            f = open(fileName)
            data = f.read()
            #term_f = {}
            term_f = defaultdict(int)
            gene_annotation = defaultdict(list) #GO_REF -> [go_id1,...go_idn]
            for match in re.finditer('(?:.*?\t){2}(?P<gene_id>\S*).*?(?P<go_id>GO:\d+).*', data, re.MULTILINE):
            #for match in re.finditer('(?:.*?\t)+(?P<gene_id>\S*).*?(?P<go_id>GO:\d+).*', data, re.MULTILINE):
                go_id, gene_id = match.group('go_id'), match.group('gene_id')
                #print go_id + " => " + go_id + ", " + gene_id
                '''match group GeneID and GoId
                '''
                '''
                if term_f.has_key(go_id):
                    term_f[go_id] += 1
                else:
                    term_f[go_id] = 1
                '''
                '''make a dictionary of GoId and frequency
                compute its frequency
                '''
                term_f[go_id] += 1
                '''in each annotation dictionary of gene add Go_Id
                '''
                if (go_id=='GO:0043241'):
                    print "***********"
                    print gene_id
                if (gene_id=='Q9UBN7'):
                    print "***********"
                if (gene_id=='Q9Y616'):
                    print "***********"
                gene_annotation[gene_id].append(go_id)
                #print gene_id
        except Exception as e:
            print e
            return
        return term_f, gene_annotation

    def parseAnnotations(self, file_name = 'go_daily-termdb.obo-xml'):
            term_list = {} #Str GO_ID-> Node object
            term_graph = defaultdict(list) #Graph str-> str list
            roots = {} #Namespace -> root
            obj = untangle.parse(file_name)
            '''parse the file with Antangle'''
            f = open(file_name)
            data = f.read()

            #Need regex to capture multiple matches for a group (is_a in our case)...
            #for m in re.finditer('<term>.*?<id>(?P<go_id>GO:\d+)</id>.*?(?:(?:<is_root>(?P<is_root>1)</is_root>)|(?:<is_a>(?P<parents>GO:\d+)</is_a>))+?.*?</term>', data, re.DOTALL):
            #    print m.group("go_id")
            #   print m.group("is_root")
            #   print m.group("parents")

            for term in obj.obo.term:
                child = term.id.cdata
                '''specify the root  each of them is the root of different graph
                '''
                if len(term.get_elements("is_root"))>0:
                    roots[term.namespace.cdata] = child
                if self.term_f.has_key(child):
                    f = self.term_f[child]
                else:
                    f = 0
                '''specify each node in each name space'''
                node = Node(child, f)
                node.namespace(term.namespace.cdata)
                try:
                    for parent in term.is_a:
                        '''in each node the parent is specified by is_a
                        '''
                        term_graph[parent.cdata].append(child)
                        node.addParent(parent.cdata)
                except:
                    pass #No parent
                term_list[child] = node
            return term_list, term_graph, roots #go_id -> node, graph, roots

    def compute_cf_and_ancestors(self,a):
        self.term_list[a].cf = self.term_list[a].f
        '''the accumulative frequency of the last child is equal to its frequency'''

        try:
            childList = self.term_graph[a]
        except Exception as ex:
            print(ex)
            return #added

        #print str(childList)
        '''the accumulative frequency of each node is the summation of its children cf  and its frequency (but this number divided by the number of its parents)        '''
        for child in childList:
            childObj = self.term_list[child]
            temp = copy.deepcopy(self.term_list[a].ancestors)
            for k,v in temp.iteritems():
                temp[k] = v + 1
                if k in childObj.ancestors.keys():
                    if v+1 < childObj.ancestors[k]:
                        childObj.ancestors[k] = temp[k]
                else:
                        childObj.ancestors[k] = temp[k]
            #childObj.ancestors.update(temp)
            childObj.ancestors[a] = 1
            self.compute_cf_and_ancestors(child)
            self.term_list[a].cf += float(childObj.cf) / len(childObj.parents)
        return True

    def semantic_similiarity(self,a, b):

        if self.term_list[a].namespace != self.term_list[b].namespace:
            #print "InputError:Genes from diferent Namespace..."
            return 0
        if a==b:
            #print "InputError:Genes from diferent Namespace..."
            return 1
        '''both of node should be from the same namespace'''
        closestAncestors = []
        maxIC = -1 * float("inf")

        aA = self.term_list[a].ancestors
        bA = self.term_list[b].ancestors

        root_f = self.term_list[self.roots[self.term_list[a].namespace]].cf

        for id in set(aA).intersection(set(bA)):
            p = float(self.term_list[id].cf) / root_f
            ic = -1*math.log(p)
            if ic == maxIC:
                closestAncestors.append(id)
            elif ic > maxIC:
                maxIC = ic
                closestAncestors = [id]

        '''Compute similarity
        '''
        #print "ClosestAncestors: " + str(closestAncestors)
        sim = -1 * float('inf')

        #print "Root CF" + str(root_f)
        a_ic = -1*math.log(self.term_list[a].cf/root_f)
        b_ic = -1*math.log(self.term_list[b].cf/root_f)

        for ancestor in closestAncestors:
            sim_t = (2*maxIC)/(a_ic+b_ic)
            #print sim_t
            if sim_t > sim:
                sim = sim_t
        return sim

if __name__=='__main__':
    start_time = time.time()
    gs = GeneSimiliarity()
    total_time = time.time() - start_time
    print "total time taken = ",total_time

    #total_f = 0
    #roots = [u'GO:0003674', u'GO:0005575', u'GO:0008150']
    #for root in roots:
    #    print gs.term_list[root].cf
    #    total_f += gs.term_list[root].cf
    #print total_f
    #for go_id, node in gs.term_list.iteritems():
    #x = gs.term_list['GO:0043241']
    #print str(x.f)
    #print go_id + " => " + str(node.f) + ", " + str(node.cf)

    gene1 = sys.argv[1]

    gene2 = sys.argv[2]
    #print gene2


    go_list1 = gs.gene_annotation[gene1]
    print 'No of go_id in gene1:', len(go_list1)
    go_list2 = gs.gene_annotation[gene2]
    print 'No of go_id in gene2:', len(go_list2)
    go_id_pairs = list(itertools.product(go_list1, go_list2))
    print 'No of go_id pairs in gene product:', len(go_id_pairs)
    i = 0
    maxsim=0

    for go_pair in go_id_pairs:
        gssim=gs.semantic_similiarity(go_pair[0], go_pair[1])
        #print go_pair,  gssim
        if (gssim>maxsim ):
            print go_pair,  gssim
            maxsim =gssim;


        #print go_pair,  gssim
    print '\n\nSemantic similarity:', maxsim

    #print gs.semantic_similiarity("GO:0043241", "GO:0043242")

