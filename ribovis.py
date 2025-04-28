from __future__ import print_function
from pymol import cmd,util
import inspect
import math
from glob import glob
#from spectrumany import spectrumany

# Pymol commands used by the Das Lab
# (C) R. Das 2010-2013, 2017-2018
#
# Some documentation and sample images available at:
#
# https://docs.google.com/document/d/1uWeEEGPjAceaw07ESf9bec-FrxW4Bx6jGaBqoHbSXuo/edit
#

@cmd.extend
def all_vs_all(filename='rms_matrix.txt'):
  AllObj=cmd.get_names("public_objects") # not temps created by alignment or selections
  matrix = []
  for x in AllObj:
    x_to_all = []
    for y in AllObj:
      rmsd = cmd.rms( x, y, cutoff=10)#[0] # super-inclusive
      x_to_all.append(rmsd)
    matrix.append(x_to_all)
  print(matrix)
  outfile = open(filename,'w')
  # No names 'on top' -- we can figure these out from the vertical ordering
  # and they're way too long...
  #for name in AllObj:
  #  outfile.write("\t%s" % name)
  #outfile.write("\n")
  maxlen = 0
  for obj in AllObj:
    if maxlen < len(obj):
      maxlen = len(obj)

  print_allobj = []
  for obj in AllObj:
    while len(obj) < maxlen:
      obj += ' '
    print_allobj.append(obj)
  ni = 0
  for x in matrix:
    outfile.write(print_allobj[ni])
    ni += 1
    for y in x:
      outfile.write("\t%0.2f" % y)
    outfile.write("\n")

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

  If you specify a third column, that will set color of the backbone,
  separately from the base.

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
  nan_res = []
  for line in lines:
    cols = line.split()
    dataval = float( cols[1] )
    if math.isnan(dataval) or not cols[0].isnumeric() or int(cols[0])==-1e18:
      if cols[0].isnumeric(): nan_res.append( int(cols[0]) )
      continue

    if min_val >= 0 and dataval < min_val: dataval = min_val
    if max_val > 0 and dataval > max_val: dataval = max_val
    try:
      data[ int( cols[0] )  ] = dataval
    except:
      pass
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


  #if palette == 'rainbow':
  cmd.spectrum( "b", palette,"all",min_val,max_val )
  #else:
  #  spectrumany( "b", palette,"all",min_val,max_val )
  #cmd.ramp_new("ramp_obj", "1gid_RNAA", range=[0, 0, max_val], color="[blue, white, red ]")

  for i in nan_res: cmd.color( 'white', 'resi  \\%d' % (i+int(offset)) )


def color_by_rgb( filename, selection = "all" ):
  """
  Read in a text file with rows like:

  125 0.1 0.1 0.1
  126 0.5 0.2 0.2

  and color specified residue numbers by RGB values.

  """
  lines = open( filename ).readlines()
  print("Applying RGB values from: " , filename)
  for line in lines:
    cols = string.split( line )
    j = cols[0]
    colorname = string.split(filename)[0] + str(j)
    cmd.set_color( colorname, (float(cols[1]),float(cols[2]),float(cols[3])) )
    #print(colorname)
    cmd.color( colorname, 'resi %s and %s' % (cols[0],selection) )

def blue_red( which_property = "resi", selection="all", minval = 0.0, maxval = 2.0 ):
  """
  mainly a mnemonic for a command I use all the time.
  Note that which_property can be "p.bisulfite" or something,
   if the user-defined property is defined!
  """
  cmd.spectrum( which_property, "blue_white_red", selection, minval, maxval )


def align_all( subset = [], cealign = False ):
  """
  Superimpose all open models onto the first one.
  This may not work well with selections.
  """
  AllObj=cmd.get_names("public_objects")
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

    if cealign:
      cmd.cealign(AllObj[0]+subset_tag, x+subset_tag)
    else:
      cmd.align(x+subset_tag,AllObj[0]+subset_tag)
    cmd.zoom()

def cealign_all( subset = [] ):
  align_all( subset, cealign = True )

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
    print(x)
    cmd.show( "cartoon", x )
    cmd.hide( "line", x )
    cmd.color( "white", x+" and elem C" )
    cmd.color( "blue", x+" and elem N" )
    cmd.color( "red", x+" and elem O" )
    cmd.color( "yellow", x+" and elem S" )
    cmd.spectrum( "resi", "rainbow", x+" and name CA+C" )
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
    print(x)
    cmd.hide( "line", x )
    cmd.color( "white", x+" and elem C" )
    cmd.color( "blue", x+" and elem N" )
    cmd.color( "red", x+" and elem O" )
    cmd.color( "yellow", x+" and elem S" )
    cmd.spectrum( "resi", "rainbow", x+" and name CA+C" )
    #cmd.show( "sticks", x +" and not elem H and not name C+N+O" )

    cmd.select('backbone_','name o+c+ca+n')
    cmd.show('sticks','not elem H')

    if not x.count( 'BACKBONE' ):
      cmd.create( x+"_BACKBONE", x+" and not element H and backbone_" )


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
    print(x)
    cmd.show( "cartoon", x )
    #cmd.hide( "line", x )
    cmd.show( "line", x )
    cmd.color( "gray", x+" and resn trp+phe+ala+val+leu+ile+pro+met" )
    cmd.color( "orange", x+" and resn gly" )
    cmd.color( "red", x+" and resn asp+glu" )
    cmd.color( "blue", x+" and resn lys+arg+his" )
    cmd.color( "purple", x+" and resn cys" )
    cmd.color( "forest", x+" and resn tyr+thr+ser+gln+asn" )
    #cmd.spectrum( "resi", "rainbow", x+" and name CA" )
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


def rr( selection = "all" ):
  """
  rhiju's favorite coloring of RNA
  with 2' OH as spheres,
  bases as filled rings, and backbone as cartoon
  ribbons, rainbow colored from 5' to 3'. No hydrogens,
  white background.
  """
  cmd.bg_color( "white" )

  cmd.hide( "everything",selection )
  cmd.show('sticks','not elem H and ' + selection )
  cmd.hide( "everything","resn HOH" )

  cmd.color('white','elem C and ' + selection )
  cmd.color( 'red','resn rG+G+GTP+DG+GUA and ' + selection)
  cmd.color( 'forest','resn rC+C+DC+CYT and ' + selection)
  cmd.color( 'orange','resn rA+A+DA+ADE and ' + selection)
  cmd.color( 'blue','resn rU+U+DT+BRU+URA+THY and ' + selection)

  cmd.color( 'red','resn rG+G+GTP+DG+GUA and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2 and ' + selection)
  cmd.color( 'forest','resn rC+C+DC+CYT and name n1+c2+o2+n3+c4+n4+c5+c6 and ' + selection)
  cmd.color( 'orange','resn rA+A+DA+ADE and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2 and ' + selection)
  cmd.color( 'blue','resn rU+U+URA+THY and name n3+c4+o4+c5+c6+n1+c2+o2 and ' + selection)

  cmd.select( 'backbone_', " (name o1p+o2p+o3p+p+op1+op2+'c1*'+'c2*'+'c3*'+'c5*'+'o2*'+'o3*'+'o4*'+'o5*'+'c1*'+'c2*'+'c3*'+'c4*'+'o2*'+'o4*'+c1'+c2'+c3'+c5'+o2'+o3'+o4'+o5'+c1'+c2'+c3'+c4'+o2'+o4') and (not name c1+c2+c3+c4+c5+o2+o3+o4+o5) and " + selection)
  cmd.spectrum( "resi", "rainbow", "backbone_" )
  cmd.cartoon( "tube", "backbone_" )
  cmd.hide( "sticks", "backbone_" )
  cmd.delete('backbone_')

  cmd.show( "cartoon", selection )
  cmd.set( "cartoon_ring_mode", 3 )
  cmd.set( "cartoon_ring_transparency", 0.0 )
  cmd.set( "cartoon_tube_radius", 0.2 )

  cmd.alter( "name o2* and "+selection,"vdw=0.5" )
  cmd.show( "spheres", "name o2'+'o2*' and not name o2+o2p and "+selection)
  cmd.show( "sticks", "name o2'+c2'+'o2*'+'c2*' and "+selection )
  cmd.show( "sticks", "resn hoh and "+selection )

  cmd.alter( "name MG and "+selection,"vdw=0.5")
  cmd.show( "spheres", "name MG and "+selection )

  cmd.alter( "resn mg and "+selection, "vdw=1.0")
  cmd.alter( "resn hoh and "+selection, "vdw=0.5")
  cmd.show( "spheres", "resn mg+sr+co+zn and not elem H and "+selection)
  cmd.hide( "ev","name RP* and "+selection)

def render_rna( selection = "all" ):
  rr( selection )

def rrs( selection = "all" ):
  """)
  rhiju's favorite coloring of RNA, showing
  all heavy atoms as sticks -- more detail than rr().
  """
  rr( selection )
  cmd.show( 'sticks', 'not elem H and '+selection )
  cmd.hide( "ev","name RP* and "+selection)

def render_rna_sticks( selection = "all" ):
  rr( selection )

def rf( selection = "all" ):
  """
  rhiju's favorite coloring of RNA, with fewer details
  """
  rr( selection )
  cmd.hide( 'spheres', selection )
  cmd.hide( 'sticks', selection )
  cmd.hide( 'lines', selection )
  cmd.set( 'cartoon_ring_mode', 1 )
  cmd.show( 'lines', "not elem H and not name *' and not name *P* and "+selection )
  cmd.show( 'lines', "name O2'+C2' and "+selection )

def render_rna_fine( selection = "all" ):
  rf( selection )

def get_residue_colors( sele = "all", outfile = "colors.txt" ):
  """
  Get RGB color values associated with a selection.
  Useful if you want to exactly match coloring of 3D models
  with coloring in, say, a MATLAB script.
  """
  pymol.stored.colors = []
  cmd.iterate( sele, "stored.colors.append( (chain, resi, name, color))")
  res_colors = {}
  stored_colors = pymol.stored.colors;
  fid = open( outfile, 'w' )
  for chain, resi, name, color in stored_colors:
    if name == 'CA' or name == 'P': # c-alpha atom
      color_tuple = cmd.get_color_tuple(color)
      res_colors[chain+resi] = color_tuple
      fid.write( '%s%s %f %f %f\n' % (chain,resi,color_tuple[0],color_tuple[1],color_tuple[2]) )
  print("Outputted RGB colors to: ", outfile)
  return res_colors

def fade_color( sele = "all", fade_extent = 0.7, by_chain = True ):
  """
  Fade out RGB color values associated with a selection.
  The by_chain variable assumes that chains are uniform color, to
    allow speedy color setting based on chain, rather than based on
    residue (which can take a long time!)
  """
  pymol.stored.colors = []
  cmd.iterate( sele, "stored.colors.append( (chain, resi, name, color))")
  res_colors = {}
  stored_colors = pymol.stored.colors;

  cols = []
  colorname = ''
  chain = ''
  prev_chain = ''
  for chain, resi, name, color in stored_colors:
    if name == 'CA' or name == 'P': # c-alpha atom
      if by_chain and chain != prev_chain and len( colorname ) > 0:
        cmd.color( colorname, 'chain %s and %s' % (prev_chain,sele) )
      color_tuple = cmd.get_color_tuple(color)
      cols = [];
      for i in range(3):
        cols.append( fade_extent  + ( 1.0  - fade_extent ) * float(color_tuple[ i ]) )
      colorname = 'temp_chain%s' % chain
      cmd.set_color( colorname, cols )
      if not by_chain: cmd.color( colorname, 'resi %s and chain %s and %s' % (resi,chain,sele) )
      prev_chain = chain

  if by_chain: cmd.color( colorname, 'chain %s and %s' % (chain,sele) )

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
    cmd.spectrum( "resi", "rainbow", x )

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

  cmd.color( 'red','resn rG+G+GTP+DG and '+selection )
  cmd.color( 'forest','resn rC+C+DC and '+selection)
  cmd.color( 'orange','resn rA+A+DA and '+selection)
  cmd.color( 'blue','resn rU+U+DT+BRU and '+selection)

  cmd.select( 'backbone_', " (name o1p+o2p+o3p+p+op1+op2+'c1*'+'c2*'+'c3*'+'c5*'+'o2*'+'o3*'+'o4*'+'o5*'+'c1*'+'c2*'+'c3*'+'c4*'+'o2*'+'o4*'+c1'+c2'+c3'+c5'+o2'+o3'+o4'+o5'+c1'+c2'+c3'+c4'+o2'+o4') and (not name c1+c2+c3+c4+c5+o2+o3+o4+o5) ")

  #for x in AllObj:
  #print(x)
  cmd.show( "cartoon", selection )
  cmd.spectrum( "resi", "rainbow", selection+" and backbone_" )

  cmd.cartoon( "tube", "backbone_ and "+selection )

  cmd.set( "cartoon_ring_mode", 0 )
  cmd.set( "cartoon_ring_transparency", 0.0 )
  cmd.set( "cartoon_tube_radius", 1.0 )

  cmd.color( 'red','resn rG+G+GTP and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2 and '+selection)
  cmd.color( 'forest','resn rC+C and name n1+c2+o2+n3+c4+n4+c5+c6 and '+selection)
  cmd.color( 'orange','resn rA+A and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2 and '+selection)
  cmd.color( 'blue','resn rU+U and name n3+c4+o4+c5+c6+n1+c2+o2 and '+selection)

  cmd.delete('backbone_')

def rcd( selection = "all" ):
  """
  fancy ribbon coloring for large RNA comparisons
  """
  rc( selection )
  cmd.cartoon( 'dumbbell')
  cmd.set( 'cartoon_dumbbell_radius', 0.5 )

def render_cartoon( selection = "all" ):
  rc( selection )

def re( selection = "all" ):
  """
  eterna-like coloring
  """
  cmd.bg_color( "white" )

  cmd.hide( "everything",selection )
  cmd.show('sticks','not elem H and ' + selection )
  cmd.hide( "everything","resn HOH" )

  cmd.color('white','elem C and ' + selection )
  cmd.color( '0xA02C28','resn rG+G+GTP+DG+GUA and ' + selection)
  cmd.color( '0x458147','resn rC+C+DC+CYT and ' + selection)
  cmd.color( '0xF4C25C','resn rA+A+DA+ADE and ' + selection)
  cmd.color( '0x3577AF','resn rU+U+DT+BRU+URA+THY and ' + selection)

  #cmd.color( 'red','resn rG+G+GTP+DG+GUA and name n1+c6+o6+c5+c4+n7+c8+n9+n3+c2+n1+n2 and ' + selection)
  #cmd.color( 'forest','resn rC+C+DC+CYT and name n1+c2+o2+n3+c4+n4+c5+c6 and ' + selection)
  #cmd.color( 'orange','resn rA+A+DA+ADE and name n1+c6+n6+c5+n7+c8+n9+c4+n3+c2 and ' + selection)
  #cmd.color( 'blue','resn rU+U+URA+THY and name n3+c4+o4+c5+c6+n1+c2+o2 and ' + selection)

  cmd.select( 'backbone_', " (name o1p+o2p+o3p+p+op1+op2+'c1*'+'c2*'+'c3*'+'c5*'+'o2*'+'o3*'+'o4*'+'o5*'+'c1*'+'c2*'+'c3*'+'c4*'+'o2*'+'o4*'+c1'+c2'+c3'+c5'+o2'+o3'+o4'+o5'+c1'+c2'+c3'+c4'+o2'+o4') and (not name c1+c2+c3+c4+c5+o2+o3+o4+o5) and " + selection)
  # cmd.spectrum( "resi", "rainbow", "backbone_" )
  cmd.cartoon( "tube", "backbone_" )
  cmd.hide( "sticks", "backbone_" )
  cmd.delete('backbone_')

  cmd.show( "cartoon", selection )
  cmd.set( "cartoon_ring_mode", 3 )
  cmd.set( "cartoon_ring_transparency", 0.0 )
  cmd.set( "cartoon_tube_radius", 0.2 )

  cmd.bg_color("0x10213B" );

def render_eterna( selection = "all" ):
  re( selection )


def load_movie( filename_string, movie_name = "mov" ):
  lst = glob( filename_string )
  lst.sort()
  for fil in lst: cmd.load(fil, movie_name )

def raytracemode():
  """
  Compliments of Possu Huang, updated with ray_tracr_gain by Rhiju.
  Use with rcd() for good visualization of ribosome.
  """
  cmd.set( 'cartoon_fancy_helices', 1 )
  cmd.set( 'cartoon_dumbbell_width', 0.6 )
  #cmd.set( 'cartoon_dumbbell_length', 1.80000 )
  cmd.set( 'cartoon_dumbbell_length', 1.00 ) # better for ribosome
  cmd.set( 'cartoon_dumbbell_radius', 0.42 )
  cmd.set( 'cartoon_loop_radius', 0.3 )
  cmd.set( 'ray_trace_mode', 1 )
  cmd.set( 'ray_opaque_background', 0 )
  cmd.set( 'ray_trace_gain', 0.5 )  # better for ribosome
  print("Now type: ray 1200,1200")
