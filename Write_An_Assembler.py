# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Write A Genome Assembler with ``blasr`` and ``(I)Python``
# ==========================================================

# <markdowncell>

# Introduction
# ============
# 
# I am not sure it is right thing to do some coding in a vacation. (See this tweet
# https://twitter.com/pm/status/284102222416056320, not sure if I would agree.)
# Anyway, before the vacation, I decided to organize whole bunch of random papers
# I collected in the last few years in my laptop. I eventually felt that I should
# read some of some theoretical papers about genome assembly again. I grabbed Gene
# Meyer's paper (
# http://bioinformatics.oxfordjournals.org/content/21/suppl_2/ii79.abstract ),
# and started to think about the problem about constructing unitigs
# (high-confident contigs). Before I tried to read the paper in detail, I just
# felt maybe that it was useful to write some quick code to check out what real
# data looked like. I started writing some code visualizing the local overlapping
# and generating the global overlapping graph. It was pretty straight forward and quite
# inspiring from the visualization. The visualization motived me to write
# more ad-hoc code in IPython to go beyond generating simple visualization during the
# vacation. I started to implement a very simple greedy algorithm to connect the
# input DNA fragments.  Eventually, I found I can atucally get the whole genome 
# assembly right!! This Notebook shows the work step by step toward a simple genome 
# assembler for long read data using ``IPython`` and ``Python``.
# 
# The input data was the E. coli K12 data from
# https://github.com/PacificBiosciences/DevNet/wiki/Hierarchical-Genome-Assembly-Process-%28HGAP%29.
# I generated a new set of pre-assembled reads as the input DNA fragment using
# HBAR-DTK ( https://github.com/PacificBiosciences/HBAR-DTK ). The statistics of
# the input reads is shown below:
# 
#     Pre-assembled Read Stats
#     #Seqs   17497
#     Min     448
#     1st Qu. 5365
#     Median  6149
#     Mean    6193
#     3rd Qu. 7189
#     Max     14762
#     Total   108363627
#     n50     6521
#     n90     4694
#     n95     3864
# 
# This is about 23x of average length of 6.2kb reads of E. coli K12 genome. The
# mean sequence identity when aligning the pre-assembled reads to the canonical
# E. coli K12 reference is 99.787% (median identity = 99.921%). The pre-assembled
# read data is available for download from https://www.dropbox.com/s/ybsxsfe5ppb04hx/pre_assembled_reads.fa. According the
# Lander-Waterman statistics, this set of reads is enough to get single contig
# for the single chromosome. The only thing one needs to worry about is that
# whether the repeats in genome can be uniquely. Since it is shown one can use
# Celera Assembler to assemble this data set into one single contig, this implies
# that all of the repeats of the rRNA operon can be resolved from this data set.
# Typically, one might need multiple libraries, especially those long jumping
# libraries, to resolve long repeats from short read (~100bp) data. For example,
# ALLPATHS-LG can use multiple short-read libraries and long read data to get
# perfect assemblies for bacteria
# ( http://genome.cshlp.org/content/early/2012/07/24/gr.141515.112.abstract and
# http://www.broadinstitute.org/software/allpaths-lg/blog/?page_id=14 ). However,
# when one can get long range information from one single library, the assembly
# process can be probably greatly simplified.  As demonstrated here, one can 
# actually get pretty good assembly using Python code with blasr for overlapping and
# Quiver for polishing the residual errors.
# 
# For long read data, the de Bruijn graph with some sophisticated graph
# operations to get around the limit of short read data to assemble genomes may
# not be the right approach. The Overlap-Layout-Consensus can be much more
# efficient for assembling long reads. Also, when the reads are getting longer,
# the overlapping between the reads becomes more and more specific. This will
# significantly reduce the complexity of the overlapping graph between the reads.
# (Here we define an overlapping graph as a graph where the nodes are reads and
# the edges are constructed from the pairs of proper overlapped reads.) When the
# overlapping graph is simple, it becomes possible to use some naive algorithm to
# get a reasonably good assembly. To make this clear, one can image if the read
# lengths are as long as the genome, the number of overlapping one needs to worry
# about will be in single digit and there will be no ambiguity at all. The
# overlapping and layout will be trivial in such extreme case.
# 
# When I use Celera Assembler (CA) to assemble genomes, I always try to peek the
# intermediate data in order to understand the process better and maybe fix few
# things manually. Michael Schatz (http://schatzlab.cshl.edu/ , twitter:
# @mike_schatz ) and Adam Phillippy ( twitter: @aphillippy ) both pointed to me
# some very useful information about CA.  (Thanks, Michael and Adam.)
# ``tigStore`` is a very useful command to see the relationship between the
# generated unitigs and the initial DNA fragments.  The file
# ``4-unitigger/best.edges`` contains the information best overlapped reads which
# one can use to construct the overlapping graph and one can learn a lot subtle
# points about genome assembly by visualizing those.  This notebook is definitely
# inspired by the great work done by the CA's developers.
# 
# This IPython notebook is just a "Notebook". I write these code for fun and for
# (self-)educational purposes. It shows that given the good long read data and
# the proper tools ``blasr`` ( https://github.com/PacificBiosciences/blasr ),
# ``pbdagcon/HBAR-DTK`` ( https://github.com/PacificBiosciences/pbdagcon ,
#  https://github.com/PacificBiosciences/HBAR-DTK ), and ``quiver``
# ( https://github.com/PacificBiosciences/GenomicConsensus)  one can learn how to
# do genome assembly with a "not-so-high-performance" language ``Python``.  The
# great ``IPython Notebook`` ( http://ipython.org/ ) allows me to try out different
# code snippet and different algorithm ideas without a lot of overhead. One
# should notice this code in this ``IPython Notebook`` is mainly for fun, not for
# "production". I chose the least resistant path for coding.  Coding style and
# code performance were considered but not taking seriously as I was doing a lot
# of experiments for different algorithm approaches when I coded this notebook.
# 
# I consider this Notebook is as a personal project for fun. The data and the
# used software can be all downloaded publicly.  Although I am an employee of
# PacBio(R), this notebook is developed in my own time. PacBio will not be 
# responsible for any information associated with this Notebook.  
# 
# --Jason Chin, Apr 2, 2013
# 

# <markdowncell>

# Get Data
# --------
# I put the pre-assembled read file in DropBox. The URL is https://www.dropbox.com/s/ybsxsfe5ppb04hx/pre_assembled_reads.fa.
# For most of part of the Notebook, one does not need the fasta file. You can also visit https://www.dropbox.com/sh/m2bjpwdua6t3iuk/3xgt1Zi9Go to see other data files in the same DropBox folder.

# <markdowncell>

# Overlapping
# -----------
# The first step to assemble the genome is to use ``blasr`` to find the overlapping between the reads. This can be done by executing this commend:
# 
#     blasr pre_assembled_reads.fa pre_assembled_reads.fa -m 4 -nproc 32 -bestn 32 -maxLCPLength 15 -nCandidates 36 -out pr_pr_bn32.m4 -noSplitSubreads
# 
# We ask ``blasr`` to output the coordinates of the alignment between the reads by specified the ``-m 4`` option.
# The ``-nproc 32`` will let ``blasr`` to use 32 threads.  For each reads, we allow ``blasr`` to find the other 32 best overlaps by using ``-bestn 32`` and scaning 36 potential hit candidateswith ``-nCandidates 36``. The ``-noSplitSubreads`` has only cosmatic effect here. It makes the ouput query identifier cleaner. 
# 
# ``blasr`` is a general purpose aligner. It is not really tuned to do overlapping work. One will see sometimes it does not handle the ends of the alignment as one would like it to do. The ``containment_tolerance`` variables used below was introduced for being less picky at the ends of the alignment. However, one would need to worry about potential mis-assemblies if the tolerance is too big.
# 
# The aligned results can be found in the DropBox folder mentioned above. Or, one can download it from this link https://www.dropbox.com/s/64rq78l9t74afkb/pr_pr_bn32.m4.

# <markdowncell>

# "``blasr -m 4``" output format
# ------------------------------
# The cells below show the content of the ``m4``  file. It containes the alignment information between the reads. The fields are ``query identifier``, ``target identifier``, ``alignment score (more negative -> better alignment)``, ``alignment identity``, ``query strand``, ``query start``, ``query end``, ``query length``, ``target strand``, ``target start``, ``target end``, ``target length``.   The rest of the fields are not used.  One trap here is that the "``target start``" and "``target end``" are not relative to the original sequence if "``target strand``"  is ``1`` .  When "``target strand= 1``", it means that the reverse-complimentary strand is aligned and the "``target start``" and "``target end``" are the coordinates in the reverse complimentary strand.  In the stand ``blasr`` output, "``query strand``" is always ``1``.

# <codecell>

#!head -10 pr_pr_bn32.m4

# <markdowncell>

# Convert ``blasr``'s output to python internal data
# --------------------------------------------------
# 
# We deinfe ``read_node`` class for the objects that holds the local overlapping information of each read. It simply holds a ``python`` dictionary (hash_map) where the key is the ``target identifier`` and the elements from ``m4`` line are converted to the right types and stored. Only "non-contaminated"" alignments are stored. If one read is full aligned within the other read, such alignment is not stored in a ``read_node``. The ``get_overlap_data`` function processes the ``m4`` file and determine whether the alignments are ``contaminated`` with some tolerence of the ends of the alignments. 

# <codecell>

import sys
import os
class read_node_0:  #experimenting on keeping the data in a list rather than a dictionary, not used

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.hits = [] 
        self.hit_map = None
        
    def add_hit(self, aln_info):
        self.hits.append ( aln_info )

    def __getitem__(self, target_name):
        if self.hit_map == None:
            self.hit_map = dict(zip( [h[0] for h in self.hits], self.hits ) )
        return self.hit_map[target_name]

class read_node:

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.hit_map = {}
        
    def add_hit(self, aln_info):
        self.hit_map[aln_info[0]] = aln_info

    def __getitem__(self, target_name):
        return self.hit_map[target_name]
    
    @property
    def hits(self):  #another convenient representation of the same data
        return self.hit_map.values()


def get_overlap_data(m4_filename):
    containment_tolerance = 50 
    permitted_error_pct = 2
    overlap_data = {}
    contained_reads = set()
    with open(m4_filename) as m4_f:
        for l in m4_f:
            l = l.strip().split()
            q_name, t_name =l[0:2]
            if q_name == t_name:
                continue
            aln_score = int(l[2])
            aln_idt = float(l[3])
            
            if aln_idt < 100-permitted_error_pct:
                continue
    
            q_strand, q_start, q_end, q_len = ( int(x) for x in l[4:8])
            t_strand, t_start, t_end, t_len = ( int(x) for x in l[8:12])
    
            if q_len - (q_end - q_start) < containment_tolerance:
                contained_reads.add(q_name)
    
            if t_len - (t_end - t_start) < containment_tolerance:
                contained_reads.add(t_name)
    
            if q_name not in overlap_data:
                overlap_data[ q_name ] = read_node(q_name, q_len)
    
            assert q_strand == 0
            if t_name not in [x[0] for x in overlap_data[ q_name ].hits]:
                overlap_data[ q_name ].add_hit( (t_name, 
                                                 aln_score, 
                                                 aln_idt, 
                                                 (q_strand, q_start, q_end, q_len),
                                                 (t_strand, t_start, t_end, t_len) ) )
    
            #symmetrized the alignment record
            if t_name not in overlap_data:
                overlap_data[ t_name ] = read_node(t_name, t_len)
    
            if q_name not in [x[0] for x in overlap_data[ t_name ].hits]:
                if t_strand == 1: 
                    t_start, t_end = t_len - t_end, t_len - t_start
                    q_start, q_end = q_len - q_end, q_len - t_start
                    overlap_data[ t_name ].add_hit( (q_name, 
                                                 aln_score, 
                                                 aln_idt, 
                                                 (0, t_start, t_end, t_len),
                                                 (1, q_start, q_end, q_len) ) )
                else:
                    overlap_data[ t_name ].add_hit( (q_name, 
                                                 aln_score, 
                                                 aln_idt, 
                                                 (0, t_start, t_end, t_len),
                                                 (0, q_start, q_end, q_len) ) )
    return overlap_data, contained_reads

# <codecell>

overlap_data, contained_reads = get_overlap_data("pr_pr_bn32.m4")

# <markdowncell>

# A simple function to convert the overlapping data to a GML file
# ---------------------------------------------------------------
# 
# The ``generate_overlap_gml`` takes the ``overlap_data`` and the ``contained_reads`` to generate a overlap graph in ``GML`` format. (You will need to install ``networkx`` in your ``Python`` environment to use this code.)  

# <codecell>

def generate_overlap_gml(overlap_data, contained_reads, gml_filename):
    containment_tolerance = 50 
    permitted_error_pct = 2
    import networkx as nx
    G=nx.DiGraph()
    node_in_graph = set()
    for q_name in overlap_data:
        if q_name in contained_reads:
            continue
        if q_name not in node_in_graph:
            G.add_node(q_name)
        targets = overlap_data[ q_name ].hits
        targets_3prime = [ h for h in targets if h[4][1] < containment_tolerance and h[0] not in contained_reads]
        targets_5prime = [ h for h in targets if h[3][1] < containment_tolerance and h[0] not in contained_reads]
        targets_3prime.sort(key = lambda k:k[1])
        targets_5prime.sort(key = lambda k:k[1])

        if len(targets_3prime) > 0:
            t = targets_3prime[0]
            t_name = t[0]
            if t_name not in node_in_graph:
                G.add_node(t_name)
            G.add_edge(q_name, t_name)

        if len(targets_5prime) > 0:
            t = targets_5prime[0]
            t_name = t[0]
            if t_name not in node_in_graph:
                G.add_node(t_name)
            G.add_edge(q_name, t_name)
    nx.write_gml(G, gml_filename)

# <codecell>

generate_overlap_gml(overlap_data, contained_reads, "overlap_graph.gml")

# <markdowncell>

# Overlapping graph visualized by Gephi
# -------------------------------------
# 
# Ther generated overlap_graph.gml can be visualized by Gephi ( https://gephi.org/
# ). One can use "YifanHu' Multilevel" graph layout algorithm to do the rough
# layout followed by "ForceAtlas 2" algorithm to smooth the graph. The cell below shows the
# global overlapping structure. One can see that it is topologically a circle plus
# a number of singletons. Namely, the read data does capture the whole genome in
# one single circular layout.

# <codecell>

from IPython.display import Image
Image(filename = "overlap_graph01.png")

# <markdowncell>

# We can also see some local "bubbles" in the graph.

# <codecell>

Image(filename = "overlap_graph02.png")

# <codecell>

Image(filename = "overlap_graph03.png")

# <markdowncell>

# Some utility constants, functions and a class for viewing the local overlaping with ``SVG`` files
# --------------------------------------------------------------------------------------------------
# 
# The cell below defines some contants and utility funtions so we can generate ``SVG`` to show the local overlapping alignments.

# <codecell>

arrow_defs = """<defs>
    <marker id="Triangle_green" stroke="green" fill="green"
      viewBox="0 0 10 10" refX="0" refY="5"
      markerUnits="strokeWidth"
      markerWidth="4" markerHeight="3"
      orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" />
    </marker>
    <marker id="Triangle_blue" stroke="blue" fill="blue"
      viewBox="0 0 10 10" refX="0" refY="5"
      markerUnits="strokeWidth"
      markerWidth="4" markerHeight="3"
      orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" />
    </marker>
    <marker id="Triangle_black" stroke="black" fill="black"
      viewBox="0 0 10 10" refX="0" refY="5"
      markerUnits="strokeWidth"
      markerWidth="4" markerHeight="3"
      orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" />
    </marker>
    <marker id="Triangle_red" stroke="red" fill="red"
      viewBox="0 0 10 10" refX="0" refY="5"
      markerUnits="strokeWidth"
      markerWidth="4" markerHeight="3"
      orient="auto">
      <path d="M 0 0 L 10 5 L 0 10 z" />
    </marker>
  </defs>"""

def svg_arrow( x1, y1, x2, y2, col, w):
    return """<g stroke = "%s" stroke-width = "%f" fill = "none">\
<path d = "M %f %f L %f %f" marker-end="url(#Triangle_%s)"/> </g>""" % (col, w, x1, y1, x2, y2, col)

def svg_line( x1, y1, x2, y2, col, w):
    return """<g stroke = "%s" stroke-width = "%f" fill = "none">\
<path d = "M %f %f L %f %f"/> </g>""" % (col, w, x1, y1, x2, y2)

class Overlap_SVG_view:
    def __init__(self):
        self.header = """<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" version="1.1" width="900px" height="300px">
"""
        self.footer = """</svg>"""
        self.elements = []
        self._svg_data = None
        
    def add_element(self, elm):
        self.elements.append(elm)
        
    def _construct_SVG(self):
        SVG_str = []
        SVG_str.append( self.header )
        for elm in self.elements:
            SVG_str.append( elm )
        SVG_str.append( self.footer )
        self._svg_data = "\n".join(SVG_str)
        self._svg_data = self._svg_data.decode('utf-8')

    def output(self, filename):
        if self._svg_data is None:
            self._construct_SVG()
        with open(filename, "w") as out_f:
            print >> out_f, self._svg_data
            
    def _repr_svg_(self):
        if self._svg_data is None:
            self._construct_SVG()
        return self._svg_data
        
    @property
    def svg(self):
        return SVG(self._repr_svg_())
            

# <markdowncell>

# ``generate_overlap_view``: a function taking the overlapping data to generate a ``SVG`` file
# ---------------------------------------------------------------------------------------------
# 
# The ``generate_overlap_view`` simple does several coordinate transformations to create a ``SVG`` file so we can see the local overlap around a given DNA fragment. The blank arrow indicates the fragment in interesting. Green and blue arrows indicate the forward and reversed proper overlapped alignments. The red arrows indicate containment alignments or improper alignments (locally inconsistent alignments due to repeats, mis-mapping or sequencing artifacts). 

# <codecell>

from IPython.display import SVG

import os

def generate_overlap_view(overlap_data, q_name):
    containment_tolerance = 10
    
    hits = overlap_data[q_name].hits
    SV = Overlap_SVG_view()
    SV.add_element(arrow_defs)
    x_offset = 10000
    x_scaling = 0.02
    x1 = 0
    x2 = overlap_data[q_name].length
    x1 += x_offset 
    x2 += x_offset
    x1 *= x_scaling
    x2 *= x_scaling
    y = 200
    SV.add_element(svg_arrow(x1, y, x2, y, "black", 4))
    hits.sort(key = lambda k: (k[3][1],k[3][3],k[3][2]))
    for h in hits:
        t_name, t_score, t_idt = h[0], h[1], h[2]
        if t_name in contained_reads:
            continue
        q_strand, q_start, q_end, q_len = h[3]
        t_strand, t_start, t_end, t_len = h[4]
        y -=10
        x1 = q_start
        x2 = q_end
        col = None
        
        if q_start < containment_tolerance and t_len - t_end < containment_tolerance:
            x1 = -t_start
        elif q_len - q_end < containment_tolerance and t_start < containment_tolerance:
            x2 = q_end + (t_len - t_end)
        else:
            col = "red"
            if t_strand == 1:
                x1, x2 = x2, x1
            print h
        if t_strand == 0 and col == None:
            col = "green"
        elif col == None:
            col = "blue"
            x1, x2 = x2, x1
        #print h
        #print min(x1, x2), max(x1, x2), abs(x1-x2)
        x1 += x_offset
        x2 += x_offset
        x1 *= x_scaling
        x2 *= x_scaling
        SV.add_element(svg_arrow(x1, y, x2, y, col, 3))
    return SV
    
q_name = overlap_data.keys()[0]
OSV = generate_overlap_view(overlap_data, q_name)
OSV.output("test.svg")
#display(OSV.svg)

# <markdowncell>

# Generate the "best overlap graph"
# ---------------------------------
# 
# Given that we see a beautiful circle overlapping graph, we like to see if can choose a seeding frgament and extend the overlapping to form a contigs.
# The code below scans through the ``overlap_data`` and assign the best score alignments for both 3' and 5' of a read as the best partner. The return of the ``get_best_overlap_graph`` is a dictionary where the key is a fragment identifier and the value is the alignment information of the 3' and 5' best overlapping alignments.

# <codecell>

def filter_overlap_hit(overlap_hit):
    containment_tolerance = 50
    q_strand, q_start, q_end, q_len = overlap_hit[3]
    t_strand, t_start, t_end, t_len = overlap_hit[4]
    if q_len - (q_end - q_start) < containment_tolerance and  t_len - (t_end - t_start) / t_len < containment_tolerance:
        return False
    elif overlap_hit[0] in contained_reads:
        return False
    else:
        return True

def get_best_overlap_graph(overlap_data):
    best_overlap_graph = {}
    containment_tolerance = 50
    
    for q_name in overlap_data:
        if q_name in contained_reads:
            continue
        
        if q_name not in best_overlap_graph:
            best_overlap_graph[q_name] = {"5p":None, "3p":None}
        
        
        targets = overlap_data[ q_name ].hits
        
        targets_5prime = [ h for h in targets if h[3][1] < containment_tolerance and filter_overlap_hit(h)]
        targets_3prime = [ h for h in targets if h[4][1] < containment_tolerance and filter_overlap_hit(h)]
    
        targets_5prime.sort(key = lambda k:k[1])
        targets_3prime.sort(key = lambda k:k[1])
    
        if len(targets_5prime) > 0:
            best_overlap_graph[q_name]["5p"] = targets_5prime[0]
            t_name = targets_5prime[0][0]
    
                
        if len(targets_3prime) > 0:
            best_overlap_graph[q_name]["3p"] = targets_3prime[0]
            t_name = targets_3prime[0][0]
    
    return best_overlap_graph

# <markdowncell>

# Check the output of ``get_best_overlap_graph``
# ----------------------------------------------
# The cell below find the best alignment from the read of the first key in the ``best_overlap_graph``. It shows the best 3' and 5' partners.

# <codecell>

best_overlap_graph = get_best_overlap_graph(overlap_data)
q_name = best_overlap_graph.keys()[0]
print q_name, "has best partners:", best_overlap_graph[q_name]

# <markdowncell>

# Find a best overlapping path from a seeding fragment
# ----------------------------------------------------
# 
# The two functions defined below allows us to extend the overlapping through the chain of overlapped fragments. The extension is terminated when (1) some inconsistent overlapping is detected (2) the overlapped fragment is already used by other contigs (3) no overlapped fragment is found and (4) a circle is detected in the path. Such rules may not be the best way to generate unitigs but they are very simple to implement. At least, these rules do work well for this E. coli data set. One can image some more complicated scenarios that it will be necessary to refine these rules or add more rules to get the best results and avoid mis-assemblies. 

# <codecell>

def proper_overlap_hit(overlap_hit):
    containment_tolerance = 50
    q_strand, q_start, q_end, q_len = overlap_hit[3]
    t_strand, t_start, t_end, t_len = overlap_hit[4]
    if q_start < containment_tolerance:
        if t_len - t_end > containment_tolerance:
            return False
        else:
            return True
    if t_start < containment_tolerance:
        if q_len - q_end > containment_tolerance:
            return False
        else:
            return True
        
def find_path( q_name_end, frag_used = set() ):
    reverse_end = {"3p":"5p", "5p":"3p"}
    path = []
    path_q_name = set()
    q_name, end = q_name_end 
    if end == "5p":
        reversing_path = True
    else:
        reversing_path = False
        
    while 1:
        if q_name in frag_used:
            if reversing_path:
                path.reverse()
            return path, "frag_used"
        
        path.append( (q_name, end) )
        path_q_name.add(q_name)
        if q_name not in best_overlap_graph:
            if reversing_path:
                path.reverse()
            return path, "terminate_1" 
        next_hit = best_overlap_graph[q_name][end]
        #print next_hit
        if next_hit == None:
            if reversing_path:
                path.reverse()
            return path, "terminate_2"
        
        if next_hit[0] in best_overlap_graph: #if not mutual good hits, break the path
            
            # Using mutual best hit might be to strigent, 
            # bh = []
            #if best_overlap_graph[next_hit[0]]["5p"]:
            #    bh.append( best_overlap_graph[next_hit[0]]["5p"][0] )
            #if best_overlap_graph[next_hit[0]]["3p"]:
            #    bh.append( best_overlap_graph[next_hit[0]]["3p"][0] )
            
            bh = [h[0] for h in overlap_data[next_hit[0]].hits if proper_overlap_hit(h)]
       
            if q_name not in bh:
                if reversing_path:
                    path.reverse()
                return path, "branch"
        
        q_name = next_hit[0]
        
        if q_name in path_q_name:
            if reversing_path:
                path.reverse()
            return path, "circle"
        
        if next_hit[4][0] == 1: #reversed strand
            end = reverse_end[end]


# <markdowncell>

# Examples of the best overlapping paths
# --------------------------------------
# 
# Here we show how to call the ``find_path`` function to get the overlapping path for both 5' and 3' from a seed fragment. We pick the longest fragment in the data set as the seed.  One will notice that the total fragments in from both 5'-parh and 3'-path are greater than the total number of fragments in the ``best_overlap_graph``. It mean that there are significant overlapping between the 5'-path and 3'-path. This is consistent with the genome is circular.

# <codecell>

len_qname = [ (overlap_data[x].length, x) for x in best_overlap_graph.keys() ]
len_qname.sort()
q_name = len_qname[-1][1]
#q_name = "0xa4919d66_0061381"
print "The longest fragment is", q_name

path_5p, s_5p = find_path( (q_name, "5p") )
path_3p, s_3p = find_path( (q_name, "3p") )
print "The number of the fragment of the 5' path is", len(path_5p)
print "The number of the fragment of the 3' path is", len(path_3p)
print "The total number of framgent of both 3' and 5' path is %d." % (len(path_3p)+len(path_5p)-1)
print "%d is greater than the total number of fragments %d in the best_overlap_graph." % (len(path_3p)+len(path_5p)-1, len(best_overlap_graph))
print "The begin of the 5'-path is", path_5p[0]
print "The end of the 3'-path is", path_3p[-1]

# <markdowncell>

# Build a squence database
# ------------------------
# Since we get the paths that seems covering the whole genome, we can use the paths to construct draft contigs. In order to do that, we need to get the fragment sequences.

# <codecell>

seq_db = {}
with open("pre_assembled_reads.fa") as seq_file:
    for l in seq_file:
        l = l.strip()
        if l[0] == ">":
            name = l[1:]
            continue
        else:
            seq_db[name] = l

# <markdowncell>

# Build a "draft" contig using the best overlapping paths
# -------------------------------------------------------
# 
# The ``layout_path`` function goes through a layout path and fill in the base from the fragment sequences. Nothing is really special here. However, it is alwaye tricky to track the right orientation of the fragements. The function also tracks the fragments used in the layout. Each fragement can only be assigned to one layout/contig.

# <codecell>

def rc_seq(seq):
    rev_map = dict(zip("acgtACGTNn-","tgcaTGCANn-"))
    return "".join([rev_map[c] for c in seq[::-1]])

def layout_path(full_path, frag_used, out_fn, tig_name):
    
    if len(full_path) == 0 or full_path[0][0] in frag_used:
        return None
    
    if len(full_path) == 1:
        with open("singleton_"+out_fn, "w") as out_seq_file:
            seq = seq_db[full_path[0][0]]
            print >>out_seq_file, ">%s" % ( "singleton_" + tig_name  )
            print >>out_seq_file, seq
            frag_used.add(full_path[0][0])
        return None
    
    first = full_path[0]
    offset = 0
    current_orientation = first[1]
    revserse_orientation = {"5p":"3p", "3p":"5p"}
    tig_seq = []
    frag_in_layout = set()
    frag_in_layout.add(first[0])
    for second in full_path[1:]:
        overlap_hit = overlap_data[first[0]][second[0]]
        
        q_strand, q_start, q_end, q_len = overlap_hit[3]
        t_strand, t_start, t_end, t_len = overlap_hit[4]
        #print overlap_hit, offset, current_orientation
        seq = seq_db[first[0]]
        if current_orientation == "3p":
            seq = rc_seq(seq)
        del tig_seq[offset:-1]
        tig_seq.extend(seq)
        
        if current_orientation == "5p":
            offset += q_start
        elif current_orientation == "3p":
            offset += q_len - q_end
        if t_strand == 1:
           current_orientation  = revserse_orientation[current_orientation] 
        frag_in_layout.add(first[0])
        first = second
        if second[0] in frag_in_layout:
            break
        if second[0] in frag_used:
            break

    seq = seq_db[first[0]]
    if current_orientation == "3p":
        seq = rc_seq(seq)
    del tig_seq[offset:-1]
    tig_seq.extend(seq)
    frag_in_layout.add(first[0])
    
    frag_used.update( frag_in_layout )
    
    tig_seq = tig_seq[0:offset+len(seq)]
    tig_seq = "".join(tig_seq)
    with open(out_fn, "w") as out_seq_file:
        print >>out_seq_file, ">%s" % tig_name
        print >>out_seq_file, tig_seq

# <markdowncell>

# Loop through all fragments to generate all contigs and singletons
# -----------------------------------------------------------------
# 
# Ok, ready for the primer time!!
# 
# The code below will pick the longest fragment as the first seed and generate a contig from its overlapping path. Then, we will pick another longest one from those reads that have not been used in a contig. This process repeats until no more fragement is avaiable. 
# 
# In the data set, we generate one big contig (4650011 bp) and 49 singleton reads.

# <codecell>

frag_used = set()

all_best_overlap_frags = set(best_overlap_graph.keys())

unused_frag = all_best_overlap_frags - frag_used
i = 0
while len(unused_frag) != 0:

    len_qname = [ (overlap_data[x].length, x) for x in unused_frag  ]
    len_qname.sort()
    q_name = len_qname[-1][1]
    print "iteration: ",i
    print "frag is used?", q_name in frag_used
    if q_name in frag_used:
        continue

    path_5p, s_5p = find_path( (q_name, "5p"), frag_used  )
    path_3p, s_3p = find_path( (q_name, "3p"), frag_used  )
    print "seeding frag:", q_name, "Length:", overlap_data[q_name].length
    print "number of unused frag:", len(unused_frag), "total overlapped frag:", len(path_3p)+len(path_5p)
    print len(path_5p), s_5p, len(path_3p), s_3p
    print "--------------------------"
    #print path_5p[0], path_3p[-1]
    if len(path_5p) + len(path_3p) == 0:
        out_fn = "tig_%05d.fa" % i
        tig_name ="tig%05d" % i
        with open("singleton_"+out_fn, "w") as out_seq_file:
            seq = seq_db[q_name]
            print >>out_seq_file, ">%s" % ( "singleton_" + tig_name  )
            print >>out_seq_file, seq
        frag_used.add(q_name)
        unused_frag = all_best_overlap_frags - frag_used
        i += 1
        continue
        
    if len(path_5p) > 0:
        #assert path_5p[-1] == path_3p[0]

        full_path = path_5p + path_3p[1:]
    else:
        full_path = path_3p
    layout_path(full_path, frag_used, "tig_%05d.fa" % i, "tig%05d" % i)
    unused_frag = all_best_overlap_frags - frag_used
    i += 1


# <markdowncell>

# What are the singletons?
# ================================
# 
# One can align the contigs to the canoincal E. coli K12 reference to evaluate the quality of the assembly. One can see from the alignment data below, it seems that some fragments have lower accuracy that does not pass the 98% identity threshold to be called "overlapped" with other fragments. Three of them does not have full alginments. In this case, one will need to examine the raw reads that corresponds to these three fragments to see whether they are minor variants or some lower quality reads are not filtered.

# <codecell>

#!cat singleton_tig_000* > singletons.fa
#!blasr singletons.fa ecoli_k12_MG1655.fasta -m 4 -noSplitSubreads -nCandidates 20 -bestn 1 | awk '{print $1" "$2" "$4" "$6" "$7" "$8" "$8-($7-$6)}'

# <markdowncell>

# How good is the main contig?
# ============================
# 
# The main contig is 4650011bp. It seems that it should cover the whole genome.  We can see the whole genome alignment using ``gepard`` ( http://www.helmholtz-muenchen.de/en/mips/services/analysis-tools/gepard/index.html ). We can see there is no larger scale mis-assembly. Promising!!
# (The strand of the K12 for this data is actuall slight different from the canonical one. We do expect to see some small structure variations.)

# <codecell>

Image(filename="tig_align.jpeg")

# <markdowncell>

# Accuracy?
# ---------
# We don't expect to get high accuracy from the "draft" assembly. The way we construct the contig is just to copy-and-paste from the pre-assembled reads which have mean identity 99.789% to the reference genome. For a 4.7Mb genome, we will expect 4700000 * (1-0.99787) ~ 10k differences. One can do another round of consensus using the pre-assembled reads. Or, we can align the raw reads on top of the draft contig and apply the "Quiver" algorithm to get the best accuracy from all raw data. Let's see how "Quiver" helps to remove the errors in the draft contig.

# <markdowncell>

# First, let's evalute the differece between the draft contigs to the canonical reference with the ``dnadiff`` from Mummer3 package ( http://mummer.sourceforge.net/ )

# <codecell>

#!dnadiff ecoli_k12_MG1655.fasta tig_00000.fa -p diff_tig0
#!echo the number of SNPs is `cat diff_tig0.snps  | wc | awk '{print $1}'`

# <markdowncell>

# The draft contig has about 10k SNPs as what we expect from the estimation earlier. We apply ``Quvier`` consensus algorithm on the the draft contig using the following script::
# 
#     #!/bin/bash
#     export SEYMOUR_HOME=/PathTo/PacBio/SMRTAnalysis
#     . $SEYMOUR_HOME/etc/setup.sh 
#     cd /Path/To/WorkingDirectory
#     cp tig_00000.fa asm.ctg.fasta
#     referenceUploader -c -p $PWD -n assembly -f asm.ctg.fasta --skipIndexUpdate
# 
#     compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=24 --noXML --h5mode=w --h5fn=out.cmp.h5 --minAccuracy=0.70 --minLength=200 \
#     -x -nCandidates 50 -x -minMatch 12 -x -bestn 1 -x -minPctIdentity 70.0 /Path/To/WorkingDirectory/input.fofn assembly/
#     
#     loadPulses /Path/To/WorkingDirectory/input.fofn out.cmp.h5\
#     -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread \
#     
#     cmph5tools.py sort out.cmp.h5 --tmp /scratch
#     
#     variantCaller.py --algorithm quiver -j 16 --referenceFilename assembly/sequence/assembly.fasta \
#     --parameters best -o output.gff  -o output.fasta -o output.fastq -q 0  -X 80 -x 5 --mapQvThreshold 0 out.cmp.h5
# 
# This will generate the consensus sequence using Quiver. The consensus output is in ``output.fasta``. We find 26 SNPs in the final consensus results. This gives a lower bound of the final assembly phred QV about -10*log10(26/4639675) = 52.5.  One can see that the ``Quiver`` algorithm which utilizes all important information from the raw trace signal is able to correct the ~10k errors in the
# draft assembly. This allows us more freedom while constructing the initial draft assembly.  We do not need to construct perfect draft assembly from perfect reads. The goal for constructing the draft is to get the contiguity right. Given the long reads, one can expect that most mapping is unambiguous and that will lead to almost perfect consensus using ``Quiver``.

# <codecell>

#dnadiff ecoli_k12_MG1655.fasta output.fasta -p diff_quiver_tig0
#!echo the number of SNPs is `cat diff_quiver_tig0.snps | wc | awk '{print $1}'`
from math import log
print "The phred scale QV of the assembly is about %0.1f" % (-10*log(26.0/4639675)/log(10))
print "The concordnace is about %0.5f%%" % (100*(1-26.0/4639675))

# <markdowncell>

# Final Thought
# =============
# 
# I had some fun with this exercise. I think this approach is general but it will
# need more tests for sure. We only have n=1 successful story for now. However,
# the reality is when we can generate longer and longer reads, the assembly
# problem will become easier and easier. DNA sequencing is about connecting
# different bases across distance. Some problems can be only solved by enough
# read length spans. I do like to see more fundamental theories published about
# genome assembly and haplotyping (e.g. the classical Lander-Waterman paper and
# the recent careful analysis done by Davie Tse's group
# http://arxiv.org/abs/1301.0068 ).  Getting good consensus from long reads is
# probably a relative simpler problem than dealing with combinatorial explosion
# when one tries to infer long range information using short reads.  Taking
# bacteria assembly as an example, it is likely impossible to assemble the E.
# coli K12 genome to single contig with single-end short read library.  Without
# any long pair-end or jumping libraries, the single-end short read data simply
# does not the power to resolve various scales of repeat structure in the genome.
# Such constraints can not be solved by simply sequencing more (see Tse's argument).
# Actually, regardless how cheap DNA sequencing is, without proper understanding
# the fundamental theories and mathematical constraints, one could just waste
# money collecting tons of useless data.  Cheap data still costs money. Cheap data
# that can not solve problems is actually expansive. One can sequence E coli. K12 genome
# to 10000x coverage with short read technologies cheaply. However, without
# proper libraries to get long range information, 10000x coverage from short fragments 
# will still lead to fragmented assemblies. One the other hand, one can see from
# the example in this ``IPython`` Notebook. Once we can generate good quality long reads,
# even an assembly newbie like me can put together a prototype of an assembler with python
# in a few days. By the way, it only takes 20 seconds on my 2012 MacBook Air (Intel i7) to 
# get the draft assembly from the ``m4`` file.

# <codecell>

#!time python Write_An_Assembler.py > /dev/null

# <codecell>


