from pymol import cmd,util
import inspect
from glob import glob
from spectrumany import spectrumany

# Pymol commands used by the Das Lab
# (C) R. Das 2010-2013.
#
# Some documentation and sample images available at:
#
# https://docs.google.com/document/d/1uWeEEGPjAceaw07ESf9bec-FrxW4Bx6jGaBqoHbSXuo/edit
#


def sa(intra=False,rainbow=True):
  """
  Superimpose all open models onto the first one.
  This may not work well with selections.
  Option intra can be set to True to enable intra_fit first, for working with multi-state (nmr) pdbs.
  [Thanks to Kyle Beauchamp for this one]
  """
  AllObj=cmd.get_names("all")
  for x in AllObj:
    print(AllObj[0],x)
    if intra==True:
      cmd.intra_fit(x)
    if rainbow==True:
      cmd.util.chainbow(x)
    cmd.align(x,AllObj[0])
    cmd.zoom()

def superimpose_all(intra=False,rainbow=True):
  sa( intra, rainbow );

def chainbow():
  """
  run chainbow on all molecules, one by one.
  """
  AllObj=cmd.get_names("all")
  for x in AllObj:
    print(AllObj[0],x)
    cmd.util.chainbow(x)

def color_by_data( filename, offset = 0, min_val=-1.0, max_val = 0.0, palette = "rainbow" ):
  """
  Read in a text file with rows like:

  125 0.12
  126 1.50

  and color specified residue numbers by scalar values.
  Takes advantage of B-factor column, and color by temperature
  function in pymol. Note that coloring is scaled/offset based
  on lowest/highest scalar value.
  """
  lines = open( filename ).readlines()
  data = {}
  data_backbone = {}

  avg_data = 0.0
  min_data = 0.0
  max_data = 0.0
  for line in lines:
    cols = string.split( line )
    dataval = float( cols[1] )
    if min_val >= 0 and dataval < min_val: dataval = min_val
    if max_val > 0 and dataval > max_val: dataval = max_val
    data[ int( cols[0] )  ] = dataval
    avg_data = avg_data + dataval
    if ( dataval < min_data ): min_data = dataval
    if ( dataval > max_data ): max_data = dataval

    if len( cols ) > 2:
      dataval2 = float( cols[2] )
      if min_val >= 0 and dataval2 < min_val: dataval2 = min_val
      if max_val > 0 and dataval2 > max_val: dataval2 = max_val
      data_backbone[ int( cols[0] ) ] = dataval2


  avg_data /= len( data.keys() )

  cmd.alter( 'all', 'b=%6.3f' % avg_data )

  for i in data.keys():
    cmd.alter( 'resi  \\%d' % (i+int(offset)),  'b=%6.3f' % data[i] )

  backbone_tag = " and (name o1p+o2p+o3p+p+op1+op2+'c1*'+'c2*'+'c3*'+'c5*'+'o2*'+'o3*'+'o4*'+'o5*'+'c1*'+'c2*'+'c3*'+'c4*'+'o2*'+'o4*'+c1'+c2'+c3'+c5'+o2'+o3'+o4'+o5'+c1'+c2'+c3'+c4'+o2'+o4') and (not name c1+c2+c3+c4+c5+o2+o3+o4+o5)"
  for i in data_backbone.keys():
    cmd.alter( 'resi  \\%d %s' % (i+int(offset),backbone_tag),  'b=%6.3f' % data_backbone[i] )

  if ( min_val < 0 ): min_val = min_data
  if ( max_val < 0 ): max_val = max_data

  if palette == 'rainbow':
    cmd.spectrum( "b", palette,"all",min_val,max_val )
  else:
    spectrumany( "b", palette,"all",min_val,max_val )
  #cmd.ramp_new("ramp_obj", "1gid_RNAA", range=[0, 0, max_val], color="[blue, white, red ]")


def align_all( subset = [] ):
  """
  Superimpose all open models onto the first one.
  This may not work well with selections.
  """
  AllObj=cmd.get_names("all")
  for x in AllObj[1:]:
    print(AllObj[0],x)

    subset_tag = ''
    if isinstance( subset, int ):
      subset_tag = ' and resi %d' % subset
    elif isinstance( subset, list ) and len( subset ) > 0:
      subset_tag = ' and resi %d' % (subset[0])
      for m in range( 1,len(subset)): subset_tag += '+%d' % subset[m]
    elif isinstance( subset, str ) and len( subset ) > 0:
      subset_tag = ' and %s' % subset
    cmd.align(x+subset_tag,AllObj[0]+subset_tag)
    cmd.zoom()

def render_molecules():
  rd()

def rd():
  """
  rhiju's favorite coloring of proteins and generic molecules
  side chains are all-heavy-atom and colored CPK, backbone is
  rainbow cartoon from N to C terminus.
  """
  cmd.bg_color( "white" )
  AllObj=cmd.get_names("all")

  for x in AllObj:
    #print(AllObj[0],x)
    print x
    cmd.show( "cartoon", x )
    cmd.hide( "line", x )
    cmd.color( "white", x+" and elem C" )
    cmd.color( "blue", x+" and elem N" )
    cmd.color( "red", x+" and elem O" )
    cmd.color( "yellow", x+" and elem S" )
    cmd.spectrum( "count", "rainbow", x+" and name CA+C" )
    cmd.show( "sticks", x +" and not elem H and not name C+N+O" )
    cmd.show( "sticks", x +" and resn PRO and name N" )
    cmd.show( "sticks", x + " and name NR+CR+CS+CP+CQ" )
    cmd.show( "sticks", x + " and not elem H and neighbor name NR+CQ+CR+CS+CP" )
    cmd.show( "sticks", x + " and not elem H and neighbor neighbor name NR+CQ+CR+CS+CP" )
    cmd.set( "cartoon_oval_width", 0.1 )
    cmd.set( "cartoon_oval_length", 0.5 )

def rx():
  """
  rhiju's favorite coloring of proteins, more details --
  no cartoon; heavy backbone
  """
  cmd.bg_color( "white" )
  AllObj=cmd.get_names("all")

  for x in AllObj:
    #print(AllObj[0],x)
    print x
    cmd.hide( "line", x )
    cmd.color( "white", x+" and elem C" )
    cmd.color( "blue", x+" and elem N" )
    cmd.color( "red", x+" and elem O" )
    cmd.color( "yellow", x+" and elem S" )
    cmd.spectrum( "count", "rainbow", x+" and name CA+C" )
    #cmd.show( "sticks", x +" and not elem H and not name C+N+O" )

    cmd.select('backbone','name o+c+ca+n')
    cmd.show('sticks','not elem H')

    if not x.count( 'BACKBONE' ):
      cmd.create( x+"_BACKBONE", x+" and not element H and backbone" )


    cmd.set('stick_radius', '0.5', "*BACKBONE" )

def render_x():
  rx()

def rj():
  """
  rhiju's residue-level favorite coloring of proteins
  """
  cmd.bg_color( "white" )
  AllObj=cmd.get_names("all")

  for x in AllObj:
    #print(AllObj[0],x)
    print x
    cmd.show( "cartoon", x )
    #cmd.hide( "line", x )
    cmd.show( "line", x )
    cmd.color( "gray", x+" and resn trp+phe+ala+val+leu+ile+pro+met" )
    cmd.color( "orange", x+" and resn gly" )
    cmd.color( "red", x+" and resn asp+glu" )
    cmd.color( "blue", x+" and resn lys+arg+his" )
    cmd.color( "purple", x+" and resn cys" )
    cmd.color( "forest", x+" and resn tyr+thr+ser+gln+asn" )
    #cmd.spectrum( "count", "rainbow", x+" and name CA" )
    cmd.show( "sticks", x +" and not elem H and not name C+N+O" )
    cmd.show( "sticks", x +" and resn PRO and name N" )
    cmd.hide( "sticks", x + " and name NR+CR+CS+CP+CQ" )
    cmd.show( "sticks", x + " and not elem H and neighbor name NR+CQ+CR+CS+CP" )
  cmd.set( "cartoon_rect_length", 0.75 )
  cmd.set( "cartoon_rect_width", 0.1 )
  cmd.set( "cartoon_oval_length", 0.6 )
  cmd.set( "cartoon_oval_width", 0.2 )

def render_rhiju():
  rj()


def rr():
  """
  rhiju's favorite coloring of RNA
  with 2' OH as spheres,
  bases as filled rings, and backbone as cartoon
  ribbons, rainbow colored from 5' to 3'. No hydrogens,
  white background.
  """
  cmd.bg_color( "white" )

  cmd.hide( 'everything' )
  cmd.show('sticks','not elem H')

  cmd.color( 'red','resn rG+G+DG')
  cmd.color( 'forest','resn rC+C+DC')
  cmd.color( 'orange','resn rA+A+DA')
  cmd.color( 'blue','resn rU+U+DT+BRU')

  #cmd.set( 'cartoon_ring_color',  'red','resn rG+G+DG')
  #cmd.set( 'cartoon_ring_color',  'forest','resn rC+C+DC')
  #cmd.set( 'cartoon_ring_color',  'orange','resn rA+A+DA')
  #cmd.set( 'cartoon_ring_color',  'blue','resn rU+U+DT+BRU')

  #cmd.select('bases','name c2+c4+c5+c6+c8+n1+n2+n3+n4+n6+n7+n9+o2+o4+o6+n1p')
  #cmd.select('backbone', 'name o1p+o2p+o3p+p+c1*+c2*+c3*+c5*+o2*+o3*+o4*+o5*')
  #cmd.select('sugar', 'name c1*+c2*+c3*+c4*+o2*+o4*')
  AllObj=cmd.get_names("all")

  cmd.color( 'red','resn rG+G and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2')
  cmd.color( 'forest','resn rC+C and name n1+c2+o2+n3+c4+n4+c5+c6')
  cmd.color( 'orange','resn rA+A and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2')
  cmd.color( 'blue','resn rU+U and name n3+c4+o4+c5+c6+n1+c2+o2')


  cmd.select( 'backbone', " (name o1p+o2p+o3p+p+op1+op2+'c1*'+'c2*'+'c3*'+'c5*'+'o2*'+'o3*'+'o4*'+'o5*'+'c1*'+'c2*'+'c3*'+'c4*'+'o2*'+'o4*'+c1'+c2'+c3'+c5'+o2'+o3'+o4'+o5'+c1'+c2'+c3'+c4'+o2'+o4') and (not name c1+c2+c3+c4+c5+o2+o3+o4+o5) ")

  for x in AllObj:
    cmd.show( "cartoon", x )
    cmd.spectrum( "count", "rainbow", x+" and backbone" )
    #cmd.color( 'white', 'backbone' )

  cmd.cartoon( "tube", "backbone" )

  cmd.set( "cartoon_ring_mode", 3 )
  cmd.set( "cartoon_ring_transparency", 0.0 )
  cmd.set( "cartoon_tube_radius", 0.2 )

  cmd.hide( "sticks", "backbone" )

  cmd.delete('backbone')


  cmd.alter( "name o2*","vdw=0.5" )
  cmd.show( "spheres", "name o2'+'o2*' and not name o2" )
  cmd.show( "sticks", "name o2'+c2'+'o2*'+'c2*' " )
  cmd.show( "sticks", "resn hoh" )

  cmd.alter( "resn mg", "vdw=1.0")
  cmd.alter( "resn hoh", "vdw=0.5")
  cmd.show( "spheres", "resn mg+sr+co+zn+hoh and not elem H")
  cmd.hide( "ev","name RP*")

def render_rna():
  rr()

def rrs():
  """)
  rhiju's favorite coloring of RNA, showing
  all heavy atoms as sticks -- more detail than rr().
  """
  rr()
  cmd.show( 'sticks', 'not elem H' )
  cmd.hide( "ev","name RP*")

def render_rna_sticks():
  rr()

def rr2():
  """
  rhiju's favorite coloring of RNA, showing
  all heavy atoms as sticks -- more detail than rr().
  """
  rr()
  cmd.hide( 'spheres' )
  cmd.hide( 'sticks' )
  cmd.set( "cartoon_ring_mode", 0 )

def render_rna2():
  rr2()


def get_residue_colors( sele ):
  """
  Get RGB color values associated with a selection.
  Useful if you want to exactly match coloring of 3D models
  with coloring in, say, a MATLAB script.
  """
  pymol.stored.colors = []
  cmd.iterate( sele, "stored.colors.append( (chain, resi, name, color))")
  res_colors = {}
  for chain, resi, name, color in pymol.stored.colors:
    if name == 'CA': # c-alpha atom
      res_colors[(chain, resi)] = cmd.get_color_tuple(color)
  print res_colors
  return res_colors

def spr():
  """
  Load up these commands again after, say, an edit.
  """

  cmd.do( 'run '+inspect.getfile(inspect.currentframe()) )

def source_pymol_rhiju():
  """
  Load up these commands again after, say, an edit.
  """
  spr()


def loop_color( start, end, native=None, zoom=False ):
  """
  Used for rendering protein loop modeling puzzles.
  White in background, colored red/blue over loop.
  """

  rd()

  cmd.color( "white", "not resi %d-%d" % (start,end) )
  #cmd.hide( "cartoon", "resi %d-%d" % (start,end) )
  #cmd.show( "sticks", "not elem H and resi %d-%d" % (start,end) )

  #before_start = start - 1
  #cmd.show( "sticks", "name C and resi %d" % (before_start) )
  #after_end = end + 1
  #cmd.show( "sticks", "name N and resi %d" % (after_end) )

  cmd.color( "salmon",  "elem C and resi %d-%d" % (start,end) )

  #cmd.show( "lines", "not elem H" )
  #cmd.hide( "cartoon",  "resi %d-%d" % (start,end) )
  #cmd.show( "sticks",  "name C+N+CA+O and resi %d-%d" % (start,end) )
  cmd.hide( "sticks", "resi %d-%d and name C+N+O" % (start,end) )
  cmd.show( "sticks", "resn PRO and name N")
  cmd.show( "sticks", x +" and ( not elem H and neighbor name NR+CR+CS+CP+CQ )" )


  if native:

    # reassign colors based on native -- spectrum colors by atom count and
    # messes up loop coloring on small loop subsegments.
    #colors = get_residue_colors( "%s and resi %d-%d" % (native,start,end) )
    #for x in AllObj:
      #for m in range( start, end+1):
        #cmd.set_color( 'color%d' % m, colors[ ('','%d' % m) ] )
        #cmd.color( 'color%d' % m, 'elem C and resi %d' % m )


    cmd.color( "white", native + " and not resi %d-%d" % (start,end) )
    #cmd.color( "palecyan", native+" and not name C+N+CA+O")
    cmd.color( "skyblue", native+" and elem C and resi %d-%d" % (start,end) )

  if zoom: cmd.zoom( "resi %d-%d" % (start,end) )


def rb():
  """
  basic cartoon coloring
  """

  AllObj=cmd.get_names("all")
  cmd.bg_color( "white" )
  cmd.hide( "ev" )
  cmd.show( "cartoon" )
  cmd.cartoon( "rectangle" )
  cmd.set( "cartoon_ring_mode", 1 )
  cmd.set( "cartoon_rect_length", 0.7 )
  cmd.set( "cartoon_rect_width", 0.2 )
  for x in AllObj:
    print(AllObj[0],x)
    cmd.spectrum( "count", "rainbow", x )

def atomcolor():
  """
  atom coloring
  """

  cmd.bg_color( "white" )
  cmd.hide( "ev" )
  cmd.show( "sticks", "not elem H" )
  cmd.show( "lines", "elem H" )
  util.cbag()
  cmd.color( "white", "elem C" )

def rc( selection = "all" ):
  """
  tube coloring for large RNA comparisons
  """
  cmd.bg_color( "white" )

  cmd.hide( 'everything',selection )

  cmd.color( 'red','resn rG+G+DG and '+selection )
  cmd.color( 'forest','resn rC+C+DC and '+selection)
  cmd.color( 'orange','resn rA+A+DA and '+selection)
  cmd.color( 'blue','resn rU+U+DT+BRU and '+selection)

  cmd.select( 'backbone', " (name o1p+o2p+o3p+p+op1+op2+'c1*'+'c2*'+'c3*'+'c5*'+'o2*'+'o3*'+'o4*'+'o5*'+'c1*'+'c2*'+'c3*'+'c4*'+'o2*'+'o4*'+c1'+c2'+c3'+c5'+o2'+o3'+o4'+o5'+c1'+c2'+c3'+c4'+o2'+o4') and (not name c1+c2+c3+c4+c5+o2+o3+o4+o5) ")

  #for x in AllObj:
  #print x
  cmd.show( "cartoon", selection )
  cmd.spectrum( "resi", "rainbow", selection+" and backbone" )

  cmd.cartoon( "tube", "backbone and "+selection )

  cmd.set( "cartoon_ring_mode", 0 )
  cmd.set( "cartoon_ring_transparency", 0.0 )
  cmd.set( "cartoon_tube_radius", 1.0 )

  cmd.color( 'red','resn rG+G and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2 and '+selection)
  cmd.color( 'forest','resn rC+C and name n1+c2+o2+n3+c4+n4+c5+c6 and '+selection)
  cmd.color( 'orange','resn rA+A and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2 and '+selection)
  cmd.color( 'blue','resn rU+U and name n3+c4+o4+c5+c6+n1+c2+o2 and '+selection)

  cmd.delete('backbone')

def rcd( selection = "all" ):
  """
  fancy ribbon coloring for large RNA comparisons
  """
  rc( selection )
  cmd.cartoon( 'dumbbell')
  cmd.set( 'cartoon_dumbbell_radius', 0.5 )

def render_cartoon( selection = "all" ):
  rc( selection )

def load_movie( filename_string, movie_name = "mov" ):
  lst = glob( filename_string )
  lst.sort()
  for fil in lst: cmd.load(fil, movie_name )
